test_that("permutation reads saved post-QC data without running BLADE", {
    qc_dir <- tempfile("seac-post-qc-")
    permutation_dir <- tempfile("seac-permutation-")
    on.exit(unlink(c(qc_dir, permutation_dir), recursive = TRUE), add = TRUE)
    dir.create(qc_dir)

    spe_dir <- file.path(qc_dir, "spe")
    HDF5Array::saveHDF5SummarizedExperiment(
        make_toy_spe(), dir = spe_dir, as.sparse = TRUE
    )
    saved_object <- file.path(spe_dir, "se.rds")
    original_hash <- unname(tools::md5sum(saved_object))

    testthat::local_mocked_bindings(
        runBladeQC = function(...) stop("BLADE was called unexpectedly"),
        .package = "SEAC"
    )
    tested <- maxMinPermutationTest(
        x = qc_dir,
        gene_sets = toy_gene_sets,
        n_permutations = 3L,
        chunk_size = 2L,
        output_dir = permutation_dir
    )

    result <- loadSEACResults(permutation_dir)
    expect_s4_class(tested, "SpatialExperiment")
    expect_equal(nrow(result), ncol(tested))
    expect_equal(result$barcode, colnames(tested))
    expect_equal(result$maxmin, tested$maxmin)
    expect_identical(unname(tools::md5sum(saved_object)), original_hash)
    expect_true(file.exists(file.path(permutation_dir, "results.h5")))
})

test_that("permutation rejects raw directories rather than running BLADE", {
    raw_dir <- tempfile("seac-raw-")
    output_dir <- tempfile("seac-permutation-")
    dir.create(raw_dir)
    on.exit(unlink(c(raw_dir, output_dir), recursive = TRUE), add = TRUE)

    expect_error(
        maxMinPermutationTest(
            x = raw_dir,
            gene_sets = toy_gene_sets,
            n_permutations = 1L,
            output_dir = output_dir
        ),
        "Expected an HDF5-backed SpatialExperiment"
    )
    expect_false(file.exists(output_dir))
})

test_that("upstream lower-tail rank rule includes ties", {
    null_scores <- c(1L, 1L, 2L, 3L, 4L)
    threshold <- sort(null_scores)[[2L]]
    observed <- 0:4
    rank_rule <- vapply(observed, function(value) {
        sum(null_scores < value) < 2L
    }, logical(1L))
    prep <- list(spe = make_toy_spe()[, seq_along(observed)], observed = observed)
    result <- SEAC:::.seac_result_table(
        prep = prep,
        n_less = vapply(observed, function(value) sum(null_scores < value), integer(1L)),
        n_leq = vapply(observed, function(value) sum(null_scores <= value), integer(1L)),
        n_permutations = length(null_scores),
        legacy_alpha = 0.4,
        fdr_cutoff = 0.05
    )

    expect_identical(rank_rule, observed <= threshold)
    expect_identical(result$legacy_hotspot, observed <= threshold)
})
