#' Create the Visium HD KD-tree BLADE mask
#'
#' @param positions A data frame with `barcode`, `in_tissue`, `array_row`, and
#'   `array_col` columns from a Visium HD tissue-positions file.
#'
#' @return A data frame in the original order with seed indicators,
#'   `nearest_seed_distance`, `blade`, and `kept` columns.
#' @export
#'
#' @examples
#' positions <- data.frame(
#'     barcode = paste0("b", 1:9),
#'     in_tissue = c(0L, rep(1L, 8L)),
#'     array_row = rep(1:3, each = 3),
#'     array_col = rep(1:3, 3)
#' )
#' bladeMaskHD(positions)
bladeMaskHD <- function(positions) {
    positions <- .seac_as_data_frame(positions)
    required <- c("barcode", "in_tissue", "array_row", "array_col")
    missing <- setdiff(required, names(positions))
    if (length(missing)) {
        .seac_abort(
            "`positions` is missing required columns: ",
            paste(missing, collapse = ", "), "."
        )
    }
    if (!nrow(positions)) {
        .seac_abort("`positions` must contain at least one row.")
    }
    if (!is.character(positions$barcode) || anyNA(positions$barcode) ||
        any(!nzchar(positions$barcode))) {
        .seac_abort("`positions$barcode` must contain non-empty character values.")
    }
    if (anyDuplicated(positions$barcode)) {
        .seac_abort("`positions$barcode` must be unique.")
    }
    if (anyNA(positions$in_tissue) ||
        any(!positions$in_tissue %in% c(0, 1))) {
        .seac_abort("`in_tissue` must contain only zero and one, without missing values.")
    }

    coords <- as.matrix(positions[, c("array_row", "array_col"), drop = FALSE])
    storage.mode(coords) <- "double"
    if (anyNA(coords) || any(!is.finite(coords))) {
        .seac_abort("Array coordinates must be finite and non-missing.")
    }

    capture_border_seed <-
        positions$array_row %in% range(positions$array_row) |
        positions$array_col %in% range(positions$array_col)
    outside_tissue_seed <- positions$in_tissue == 0
    seed <- capture_border_seed | outside_tissue_seed
    seed_coords <- coords[seed, , drop = FALSE]

    # This is one exact nearest-seed KD-tree query, equivalent to the Python
    # query_ball_tree(..., r = 1) union used by SEAC.py::add_blade_column().
    nearest <- RANN::nn2(
        data = seed_coords,
        query = coords,
        k = 1L,
        eps = 0
    )$nn.dists[, 1L]
    blade <- nearest <= 1

    cbind(
        positions,
        capture_border_seed = capture_border_seed,
        outside_tissue_seed = outside_tissue_seed,
        nearest_seed_distance = as.numeric(nearest),
        blade = blade,
        kept = !blade
    )
}
