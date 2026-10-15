test_that("runBladeQC imports and writes the bundled Visium HD example", {
    skip_if_not_installed("VisiumIO")
    input_dir <- system.file("extdata", package = "VisiumIO", mustWork = TRUE)
    output_dir <- tempfile("seac-post-qc-")

    qc <- runBladeQC(
        input_dir = input_dir,
        output_dir = output_dir,
        bin_size = "002",
        sample_id = "visiumio-example"
    )

    expect_s4_class(qc, "SpatialExperiment")
    expect_gt(ncol(qc), 0L)
    expect_true(file.exists(file.path(output_dir, "blade_mask.h5")))
    expect_true(file.exists(file.path(output_dir, "spe", "se.rds")))
    expect_true(file.exists(file.path(output_dir, "manifest.json")))
})

test_that("the Visium HD importer remains disk-backed without Arrow", {
    skip_if_not_installed("VisiumIO")
    input_dir <- system.file("extdata", package = "VisiumIO", mustWork = TRUE)
    paths <- SEAC:::.seac_resolve_hd_paths(input_dir, "002")
    positions <- SEAC:::.seac_import_positions(paths$positions)
    spe <- SEAC:::.seac_import_visium_hd(paths, "example", positions)

    expect_s4_class(SummarizedExperiment::assay(spe), "TENxMatrix")
    expect_setequal(colnames(spe), positions$barcode)
    expect_equal(as.character(spe$barcode), colnames(spe))
    expect_equal(
        colnames(SpatialExperiment::spatialCoords(spe)),
        c("array_col", "array_row")
    )
    identifiers <- SEAC:::.seac_feature_identifiers(spe)
    expect_equal(
        identifiers$aliases,
        as.character(SummarizedExperiment::rowData(spe)$Symbol)
    )
    expect_true(any(identifiers$primary != identifiers$aliases))

    synthetic <- matrix(0L, nrow(spe), ncol(spe), dimnames = dimnames(spe))
    synthetic[1L, 1L] <- 1L
    synthetic[2L, 2L] <- 1L
    SummarizedExperiment::assay(spe, "counts") <- synthetic
    scored <- maxMinScore(spe, list(first = "Xkr4", second = "Rp1"))
    expect_length(scored$maxmin, ncol(spe))
})
