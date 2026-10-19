reference_maxmin <- function(counts, coordinates, gene_rows) {
    distances <- lapply(gene_rows, function(rows) {
        sources <- which(colSums(counts[rows, , drop = FALSE]) > 0L)
        if (!length(sources)) return(rep(Inf, ncol(counts)))
        vapply(seq_len(ncol(counts)), function(bin) {
            displacement <- sweep(
                coordinates[sources, , drop = FALSE], 2L,
                coordinates[bin, ], `-`
            )
            ceiling(min(sqrt(rowSums(displacement^2))))
        }, numeric(1L))
    })
    do.call(pmax, distances)
}

reference_umi_null <- function(counts, coordinates, gene_rows, seed) {
    category_totals <- vapply(gene_rows, function(rows) {
        sum(counts[rows, , drop = FALSE])
    }, numeric(1L))
    total <- sum(counts)
    labels <- withr::with_rng_version("4.0.0", withr::with_seed(seed, {
        selected <- sample.int(total, sum(category_totals), replace = FALSE)
        shuffled_labels <- sample(rep(seq_along(gene_rows), category_totals))
        list(selected = selected, labels = shuffled_labels)
    }))
    # Explicitly expand each UMI for this small reference case rather than
    # using the package's cumulative-slot mapping or its KD-tree query.
    umi_bins <- rep(seq_len(ncol(counts)), colSums(counts))
    assigned_bins <- split(umi_bins[labels$selected], labels$labels)
    null_presence <- matrix(0L, nrow = length(gene_rows), ncol = ncol(counts))
    for (k in seq_along(assigned_bins)) {
        null_presence[k, unique(assigned_bins[[k]])] <- 1L
    }
    reference_maxmin(null_presence, coordinates, as.list(seq_along(gene_rows)))
}

test_that("MaxMin and within-sample null agree with an expanded-UMI reference", {
    spe <- make_toy_spe()
    counts <- SummarizedExperiment::assay(spe, "counts")
    coordinates <- as.matrix(SpatialExperiment::spatialCoords(spe))
    gene_rows <- lapply(toy_gene_sets, function(genes) match(genes, rownames(spe)))
    observed <- reference_maxmin(counts, coordinates, gene_rows)
    expect_equal(unname(maxMinScore(spe, toy_gene_sets)$maxmin), observed)

    prep <- SEAC:::.seac_prepare_sample(spe, toy_gene_sets)
    null_scores <- lapply(seq_len(7L), function(id) {
        seed <- 100L + id
        reference <- reference_umi_null(counts, coordinates, gene_rows, seed)
        expect_equal(unname(SEAC:::.seac_within_null(prep, seed)), reference)
        reference
    })
    null_matrix <- do.call(rbind, null_scores)
    expected_less <- colSums(sweep(null_matrix, 2L, observed, `<`))
    expected_leq <- colSums(sweep(null_matrix, 2L, observed, `<=`))

    output <- tempfile("seac-reference-")
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 7L, chunk_size = 3L,
        seed = 100L, output_dir = output
    )
    results <- loadSEACResults(output)
    expect_equal(results$n_null_less, as.integer(expected_less))
    expect_equal(results$n_null_leq, as.integer(expected_leq))
    expect_equal(results$p_value, (1 + expected_leq) / 8)
    expect_equal(results$p_adjusted, p.adjust(results$p_value, "BH"))
})

test_that("pooled mode agrees with a cohort-wide expanded-UMI reference", {
    first <- make_toy_spe()
    second <- make_toy_spe()
    SummarizedExperiment::assay(second, "counts")[1L, 1L] <- 2L
    samples <- list(first = first, second = second)
    counts <- lapply(samples, SummarizedExperiment::assay, "counts")
    coordinates <- lapply(samples, function(spe) {
        as.matrix(SpatialExperiment::spatialCoords(spe))
    })
    gene_rows <- lapply(toy_gene_sets, function(genes) {
        match(genes, rownames(first))
    })
    observed <- Map(reference_maxmin, counts, coordinates,
                    MoreArgs = list(gene_rows = gene_rows))
    category_totals <- vapply(gene_rows, function(rows) {
        sum(vapply(counts, function(x) sum(x[rows, , drop = FALSE]), numeric(1L)))
    }, numeric(1L))
    all_bin_totals <- unlist(lapply(counts, colSums), use.names = FALSE)
    umi_bins <- rep(seq_along(all_bin_totals), all_bin_totals)
    sample_sizes <- vapply(counts, ncol, integer(1L))
    sample_ends <- cumsum(sample_sizes)
    sample_starts <- c(0L, head(sample_ends, -1L))

    reference <- lapply(seq_len(7L), function(id) {
        draws <- withr::with_rng_version("4.0.0", withr::with_seed(200L + id, {
            list(
                selected = sample.int(length(umi_bins), sum(category_totals)),
                labels = sample(rep(seq_along(gene_rows), category_totals))
            )
        }))
        assigned <- split(umi_bins[draws$selected], draws$labels)
        Map(function(sample_index, sample_counts, sample_coordinates) {
            presence <- matrix(0L, nrow = length(gene_rows),
                               ncol = ncol(sample_counts))
            for (k in seq_along(assigned)) {
                local <- assigned[[k]] - sample_starts[[sample_index]]
                local <- local[local >= 1L & local <= sample_sizes[[sample_index]]]
                presence[k, unique(local)] <- 1L
            }
            reference_maxmin(
                presence, sample_coordinates, as.list(seq_along(gene_rows))
            )
        }, seq_along(samples), counts, coordinates)
    })

    output <- tempfile("seac-pooled-reference-")
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    maxMinPermutationTest(
        samples, toy_gene_sets, mode = "pooled", n_permutations = 7L,
        chunk_size = 3L, seed = 200L, output_dir = output
    )
    for (s in seq_along(samples)) {
        null_matrix <- do.call(rbind, lapply(reference, `[[`, s))
        expected_less <- colSums(sweep(null_matrix, 2L, observed[[s]], `<`))
        expected_leq <- colSums(sweep(null_matrix, 2L, observed[[s]], `<=`))
        result <- loadSEACResults(output, sample_id = names(samples)[[s]])
        expect_equal(result$n_null_less, as.integer(expected_less))
        expect_equal(result$n_null_leq, as.integer(expected_leq))
    }
})
