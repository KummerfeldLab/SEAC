.seac_abort <- function(...) {
    stop(..., call. = FALSE)
}

.seac_assert_scalar_character <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        .seac_abort("`", name, "` must be one non-empty character string.")
    }
    invisible(x)
}

.seac_is_scalar_number <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.seac_assert_whole_number <- function(x, name, minimum = 0,
                                      maximum = .Machine$integer.max) {
    if (!.seac_is_scalar_number(x) || x < minimum || x > maximum ||
        x != floor(x)) {
        .seac_abort(
            "`", name, "` must be one whole number between ", minimum,
            " and ", maximum, "."
        )
    }
    invisible(as.integer(x))
}

.seac_as_data_frame <- function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
}

.seac_coordinates <- function(x) {
    cd <- .seac_as_data_frame(SummarizedExperiment::colData(x))
    required <- c("array_row", "array_col")
    missing <- setdiff(required, names(cd))
    if (!length(missing)) {
        coords <- as.matrix(cd[, required, drop = FALSE])
    } else {
        # SpatialExperiment commonly moves named coordinate columns out of
        # colData and into spatialCoords() during construction/import.
        spatial <- SpatialExperiment::spatialCoords(x)
        spatial_missing <- setdiff(required, colnames(spatial))
        if (length(spatial_missing)) {
            .seac_abort(
                "`SpatialExperiment` must provide array_row and array_col in ",
                "colData or spatialCoords()."
            )
        }
        coords <- spatial[, required, drop = FALSE]
    }
    coords <- as.matrix(coords)
    storage.mode(coords) <- "double"
    if (anyNA(coords) || any(!is.finite(coords))) {
        .seac_abort("Array coordinates must be finite and non-missing.")
    }
    coords
}

.seac_counts <- function(x) {
    nms <- SummarizedExperiment::assayNames(x)
    if (!length(nms)) {
        .seac_abort("The object has no assay.")
    }
    selected <- if ("counts" %in% nms) "counts" else nms[[1L]]
    SummarizedExperiment::assay(x, selected)
}

.seac_col_sums <- function(x) {
    as.numeric(DelayedMatrixStats::colSums2(x, useNames = FALSE))
}

.seac_row_sums <- function(x) {
    as.numeric(DelayedMatrixStats::rowSums2(x, useNames = FALSE))
}

.seac_with_seed <- function(seed, code) {
    withr::with_rng_version(
        "4.0.0",
        withr::with_seed(as.integer(seed), code)
    )
}

.seac_hash <- function(x) {
    path <- tempfile("seac-hash-", fileext = ".rds")
    on.exit(if (file.exists(path)) unlink(path), add = TRUE)
    saveRDS(x, path, version = 2L)
    unname(tools::md5sum(path))
}

.seac_write_manifest <- function(path, value) {
    jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE)
    invisible(path)
}
