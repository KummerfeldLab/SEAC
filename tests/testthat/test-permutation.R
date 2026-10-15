test_that("permutation output is deterministic and has compact result columns", {
    spe <- make_toy_spe()
    d1 <- tempfile("seac_perm_")
    scored1 <- maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 5L, chunk_size = 2L,
        seed = 42L, output_dir = d1
    )
    expect_true(file.exists(file.path(d1, "results.h5")))
    expect_true(all(c("p_value", "p_adjusted", "legacy_hotspot") %in%
        names(SummarizedExperiment::colData(scored1))))
    result <- loadSEACResults(d1)
    expect_equal(nrow(result), ncol(spe))
    expect_true(all(result$p_value > 0 & result$p_value <= 1))
})

test_that("a completed run can be resumed without changing results", {
    spe <- make_toy_spe()
    d1 <- tempfile("seac_resume_")
    first <- maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 4L, chunk_size = 2L,
        seed = 9L, output_dir = d1
    )
    second <- maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 4L, chunk_size = 1L,
        seed = 9L, output_dir = d1, resume = TRUE
    )
    expect_equal(first$p_value, second$p_value)
})

test_that("a zero-progress checkpoint and seed validation protect restarts", {
    checkpoint_dir <- tempfile("seac_checkpoint_")
    dir.create(checkpoint_dir)
    SEAC:::.seac_write_checkpoint(
        checkpoint_dir, integer(), list(integer(2L)), list(integer(2L)),
        list(version = 1L)
    )
    expect_true(file.exists(file.path(checkpoint_dir, "checkpoint_000000.rds")))

    spe <- make_toy_spe()
    expect_error(
        maxMinPermutationTest(
            spe, toy_gene_sets, n_permutations = 1L, seed = NA_integer_,
            output_dir = tempfile("seac_bad_seed_")
        ),
        "seed"
    )
})

test_that("checkpoint reading falls back from an interrupted atomic write", {
    checkpoint_dir <- tempfile("seac_checkpoint_fallback_")
    dir.create(checkpoint_dir)
    specification <- list(version = 2L)
    SEAC:::.seac_write_checkpoint(
        checkpoint_dir, integer(), list(integer(2L)), list(integer(2L)),
        specification
    )
    writeLines("incomplete", file.path(checkpoint_dir, "checkpoint_000002.rds"))

    recovered <- SEAC:::.seac_read_checkpoint(checkpoint_dir, specification)
    expect_identical(recovered$completed, integer())
})

test_that("resume signatures detect rearranged per-bin input", {
    prep <- SEAC:::.seac_prepare_sample(make_toy_spe(), toy_gene_sets)
    changed <- prep
    changed$bin_totals <- rev(changed$bin_totals)

    expect_equal(sum(prep$bin_totals), sum(changed$bin_totals))
    expect_false(identical(
        SEAC:::.seac_input_signature(prep)$hash,
        SEAC:::.seac_input_signature(changed)$hash
    ))
})

test_that("permutation testing rejects non-integer expression totals", {
    spe <- make_toy_spe()
    SummarizedExperiment::assay(spe, "counts")[1L, 1L] <- 0.5
    expect_error(
        maxMinPermutationTest(
            spe, toy_gene_sets, n_permutations = 1L,
            output_dir = tempfile("seac_non_integer_")
        ),
        "integer UMI totals"
    )
})

test_that("UMI-slot draws preserve category totals without reusing a UMI", {
    prep <- SEAC:::.seac_prepare_sample(make_toy_spe(), toy_gene_sets)
    slots <- SEAC:::.seac_draw_category_slots(
        sum(prep$bin_totals), prep$category_totals, seed = 33L
    )

    expect_equal(unname(lengths(slots)), unname(prep$category_totals))
    expect_false(anyDuplicated(unlist(slots, use.names = FALSE)) > 0L)
    expect_true(all(unlist(slots, use.names = FALSE) >= 1L))
    expect_true(all(unlist(slots, use.names = FALSE) <= sum(prep$bin_totals)))
})

test_that("worker inputs omit the full SpatialExperiment", {
    prep <- SEAC:::.seac_prepare_sample(make_toy_spe(), toy_gene_sets)
    worker_input <- SEAC:::.seac_permutation_input(prep)

    expect_named(worker_input, c(
        "coords", "observed", "bin_totals", "cumulative",
        "category_totals", "n_bins"
    ))
    expect_false("spe" %in% names(worker_input))
})

test_that("chunking is deterministic and pooled mode returns every sample", {
    spe <- make_toy_spe()
    d1 <- tempfile("seac_chunk1_")
    d2 <- tempfile("seac_chunk3_")
    one <- maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 3L, chunk_size = 1L,
        seed = 17L, output_dir = d1
    )
    three <- maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 3L, chunk_size = 3L,
        seed = 17L, output_dir = d2
    )
    expect_equal(one$p_value, three$p_value)

    pooled <- maxMinPermutationTest(
        list(a = spe, b = spe), toy_gene_sets, mode = "pooled",
        n_permutations = 2L, output_dir = tempfile("seac_pooled_")
    )
    expect_named(pooled, c("a", "b"))
    expect_true(all(vapply(pooled, inherits, logical(1L), "SpatialExperiment")))
})

test_that("permutation testing preserves the caller's random state", {
    spe <- make_toy_spe()
    stats::runif(1L)
    state_before <- .Random.seed
    kind_before <- RNGkind()

    maxMinPermutationTest(
        spe, toy_gene_sets, n_permutations = 2L, seed = 27L,
        output_dir = tempfile("seac_rng_")
    )

    expect_identical(.Random.seed, state_before)
    expect_identical(RNGkind(), kind_before)
})
