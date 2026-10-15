test_that("runSEAC validates workflow switches before writing output", {
    output_dir <- tempfile("seac-workflow-")

    expect_error(
        runSEAC(
            tempfile("missing-input-"), output_dir = output_dir,
            permutation = TRUE
        ),
        "Supply `gene_sets`"
    )
    expect_false(file.exists(output_dir))

    expect_error(
        runSEAC(tempfile("missing-input-"), permutation = NA),
        "must be TRUE or FALSE"
    )
})

test_that("scored post-QC metadata remains available after reloading", {
    output_dir <- tempfile("seac-persist-")
    spe_dir <- file.path(output_dir, "spe")
    on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
    dir.create(output_dir)

    HDF5Array::saveHDF5SummarizedExperiment(
        make_toy_spe(), dir = spe_dir, as.sparse = TRUE
    )
    qc <- HDF5Array::loadHDF5SummarizedExperiment(spe_dir)
    scored <- maxMinScore(qc, toy_gene_sets)
    reloaded <- SEAC:::.seac_persist_post_qc_metadata(scored, output_dir)

    expect_equal(reloaded$maxmin, scored$maxmin)
    expect_equal(
        reloaded$maxmin_distance_first,
        scored$maxmin_distance_first
    )
    expect_equal(
        reloaded$maxmin_distance_second,
        scored$maxmin_distance_second
    )
    expect_equal(
        S4Vectors::metadata(reloaded)$SEAC$maxmin,
        S4Vectors::metadata(scored)$SEAC$maxmin
    )
    expect_equal(
        as.matrix(SummarizedExperiment::assay(reloaded)),
        as.matrix(SummarizedExperiment::assay(scored))
    )
    expect_gt(file.info(file.path(spe_dir, "assays.h5"))$size, 0)
})

test_that("runSEAC persists scores in its post-QC output", {
    input_dir <- tempfile("seac-input-")
    output_dir <- tempfile("seac-workflow-")
    dir.create(input_dir)
    on.exit(unlink(c(input_dir, output_dir), recursive = TRUE), add = TRUE)

    mock_run_blade <- function(input_dir, output_dir, bin_size, sample_id) {
        spe_dir <- file.path(output_dir, "spe")
        dir.create(output_dir, recursive = TRUE)
        HDF5Array::saveHDF5SummarizedExperiment(
            make_toy_spe(), dir = spe_dir, as.sparse = TRUE
        )
        HDF5Array::loadHDF5SummarizedExperiment(spe_dir)
    }
    testthat::local_mocked_bindings(
        runBladeQC = mock_run_blade,
        .package = "SEAC"
    )

    result <- runSEAC(
        input_dir,
        output_dir = output_dir,
        gene_sets = toy_gene_sets
    )
    reloaded <- HDF5Array::loadHDF5SummarizedExperiment(
        file.path(output_dir, "post_qc", "spe")
    )

    expect_equal(reloaded$maxmin, result$post_qc$maxmin)
    expect_true(all(c(
        "maxmin_distance_first", "maxmin_distance_second", "maxmin"
    ) %in% colnames(SummarizedExperiment::colData(reloaded))))
})
