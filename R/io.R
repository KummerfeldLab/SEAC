.seac_resolve_hd_paths <- function(input_dir, bin_size = "002") {
    .seac_assert_scalar_character(input_dir, "input_dir")
    .seac_assert_scalar_character(bin_size, "bin_size")
    valid_integer <- grepl("^[0-9]+$", bin_size) && nchar(bin_size) <= 9L
    if (!valid_integer) {
        .seac_abort("`bin_size` must be a positive integer encoded as text, e.g. `002`.")
    }
    bin_number <- as.integer(bin_size)
    if (bin_number <= 0L) {
        .seac_abort("`bin_size` must be a positive integer encoded as text, e.g. `002`.")
    }
    if (!dir.exists(input_dir)) {
        .seac_abort("Input directory does not exist: ", input_dir)
    }
    bin_name <- paste0("square_", sprintf("%03d", bin_number), "um")
    direct <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
    candidates <- c(
        direct,
        file.path(direct, "binned_outputs", bin_name),
        file.path(direct, "outs", "binned_outputs", bin_name)
    )
    is_bin_dir <- vapply(candidates, function(path) {
        file.exists(file.path(path, "filtered_feature_bc_matrix.h5")) &&
            file.exists(file.path(path, "spatial", "tissue_positions.parquet"))
    }, logical(1L))
    if (!any(is_bin_dir)) {
        .seac_abort(
            "Could not find filtered_feature_bc_matrix.h5 and ",
            "spatial/tissue_positions.parquet for ", bin_name, "."
        )
    }
    bin_dir <- candidates[which(is_bin_dir)[1L]]
    resolved_bin_size <- sub("^square_([0-9]+)um$", "\\1", basename(bin_dir))
    if (!grepl("^[0-9]+$", resolved_bin_size)) {
        resolved_bin_size <- sprintf("%03d", bin_number)
    }
    list(
        bin_dir = bin_dir,
        h5 = file.path(bin_dir, "filtered_feature_bc_matrix.h5"),
        positions = file.path(bin_dir, "spatial", "tissue_positions.parquet"),
        bin_size = resolved_bin_size
    )
}

.seac_import_positions <- function(path) {
    .seac_as_data_frame(nanoparquet::read_parquet(path))
}

.seac_import_visium_hd <- function(paths, sample_id, positions) {
    # Import the HDF5 count matrix directly so reading the one flat positions
    # table does not require the much larger Arrow C++ runtime. TENxIO keeps
    # the assay as a disk-backed TENxMatrix.
    sce <- BiocIO::import(TENxIO::TENxH5(paths$h5, ranges = NA_character_))
    idx <- match(colnames(sce), positions$barcode)
    if (anyNA(idx)) {
        .seac_abort("Some count-matrix barcodes are absent from tissue positions.")
    }
    spatial_data <- positions[idx, , drop = FALSE]
    rownames(spatial_data) <- colnames(sce)

    SpatialExperiment::SpatialExperiment(
        assays = SummarizedExperiment::assays(sce),
        rowData = SummarizedExperiment::rowData(sce),
        mainExpName = SingleCellExperiment::mainExpName(sce),
        altExps = SingleCellExperiment::altExps(sce),
        sample_id = sample_id,
        colData = S4Vectors::DataFrame(spatial_data),
        spatialCoordsNames = c("array_col", "array_row")
    )
}

.seac_sample_id_from_input <- function(input_dir) {
    absolute <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
    base <- basename(absolute)
    if (grepl("^square_[0-9]+um$", base)) {
        return(basename(dirname(dirname(dirname(absolute)))))
    }
    if (identical(base, "binned_outputs")) {
        return(basename(dirname(dirname(absolute))))
    }
    if (identical(base, "outs")) {
        return(basename(dirname(absolute)))
    }
    base
}

.seac_write_mask_h5 <- function(mask, path) {
    if (file.exists(path)) {
        .seac_abort("Refusing to overwrite existing mask file: ", path)
    }
    rhdf5::h5createFile(path)
    for (nm in names(mask)) {
        value <- mask[[nm]]
        if (is.factor(value)) value <- as.character(value)
        rhdf5::h5write(value, path, nm)
    }
    invisible(path)
}

#' Import, mask, and save Visium HD data
#'
#' @param input_dir A Space Ranger root, `outs`, or `square_XXXum` directory.
#' @param output_dir Destination directory. By default a new directory is made
#'   below the current working directory.
#' @param bin_size Visium HD square bin size, such as `"002"`.
#' @param sample_id Optional sample identifier.
#'
#' @return A HDF5-backed `SpatialExperiment` containing non-BLADE bins.
#' @export
#' @examples
#' if (requireNamespace("VisiumIO", quietly = TRUE)) {
#'     input <- system.file("extdata", package = "VisiumIO")
#'     output <- tempfile("seac-qc-")
#'     qc <- runBladeQC(input, output, bin_size = "002",
#'                      sample_id = "example")
#'     dim(qc)
#' }
runBladeQC <- function(input_dir, output_dir = NULL, bin_size = "002",
                       sample_id = NULL) {
    paths <- .seac_resolve_hd_paths(input_dir, bin_size)
    if (is.null(sample_id)) {
        sample_id <- .seac_sample_id_from_input(input_dir)
    }
    .seac_assert_scalar_character(sample_id, "sample_id")
    if (is.null(output_dir)) {
        output_dir <- file.path(getwd(), paste0(sample_id, "_SEAC_post_qc"))
    }
    .seac_assert_scalar_character(output_dir, "output_dir")
    if (dir.exists(output_dir) || file.exists(output_dir)) {
        .seac_abort("Output already exists: ", output_dir)
    }

    positions <- .seac_import_positions(paths$positions)
    mask <- bladeMaskHD(positions)
    spe <- .seac_import_visium_hd(paths, sample_id, positions)
    idx <- match(colnames(spe), mask$barcode)
    if (anyNA(idx)) {
        .seac_abort("Some count-matrix barcodes are absent from tissue positions.")
    }
    keep <- mask$kept[idx]
    qc_spe <- spe[, keep, drop = FALSE]
    cd <- SummarizedExperiment::colData(qc_spe)
    cd$blade <- FALSE
    qc_spe <- SummarizedExperiment::`colData<-`(qc_spe, value = cd)
    md <- S4Vectors::metadata(qc_spe)
    md$SEAC <- list(
        sample_id = sample_id,
        bin_size = paths$bin_size,
        blade = list(algorithm = "kdtree-radius", radius = 1,
                     n_positions = nrow(mask), n_removed = sum(!mask$kept))
    )
    qc_spe <- S4Vectors::`metadata<-`(qc_spe, value = md)

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    mask_path <- file.path(output_dir, "blade_mask.h5")
    .seac_write_mask_h5(mask, mask_path)
    spe_dir <- file.path(output_dir, "spe")
    HDF5Array::saveHDF5SummarizedExperiment(
        qc_spe, dir = spe_dir, replace = FALSE, as.sparse = TRUE
    )
    .seac_write_manifest(file.path(output_dir, "manifest.json"), list(
        package = "SEAC",
        stage = "post_qc",
        sample_id = sample_id,
        input = normalizePath(paths$bin_dir, winslash = "/", mustWork = TRUE),
        bin_size = paths$bin_size,
        blade = list(algorithm = "kdtree-radius", radius = 1),
        n_positions = nrow(mask),
        n_retained = sum(mask$kept),
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ))
    HDF5Array::loadHDF5SummarizedExperiment(spe_dir)
}
.seac_load_post_qc <- function(path) {
    if (inherits(path, "SpatialExperiment")) {
        return(path)
    }
    .seac_assert_scalar_character(path, "x")
    spe_dir <- if (dir.exists(file.path(path, "spe"))) file.path(path, "spe") else path
    if (!file.exists(file.path(spe_dir, "se.rds"))) {
        .seac_abort("Expected an HDF5-backed SpatialExperiment at: ", spe_dir)
    }
    HDF5Array::loadHDF5SummarizedExperiment(spe_dir)
}
