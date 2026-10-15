.seac_feature_identifiers <- function(spe) {
    primary <- rownames(spe)
    if (is.null(primary) || anyNA(primary) || any(!nzchar(primary)) ||
        anyDuplicated(primary)) {
        .seac_abort("The input must have unique, non-empty feature row names.")
    }
    rd <- .seac_as_data_frame(SummarizedExperiment::rowData(spe))
    alias_column <- intersect(c("Symbol", "symbol", "gene_symbol"), names(rd))
    aliases <- rep.int(NA_character_, length(primary))
    if (length(alias_column)) {
        aliases <- as.character(rd[[alias_column[[1L]]]])
        aliases[!nzchar(aliases)] <- NA_character_
    }
    list(primary = primary, aliases = aliases)
}

.seac_validate_gene_sets <- function(gene_sets, identifiers) {
    if (!is.list(gene_sets) || !length(gene_sets) || is.null(names(gene_sets)) ||
        anyNA(names(gene_sets)) || any(!nzchar(names(gene_sets))) ||
        anyDuplicated(names(gene_sets))) {
        .seac_abort("`gene_sets` must be a named list with unique, non-empty names.")
    }
    if (length(gene_sets) < 2L) {
        .seac_abort("`gene_sets` must contain at least two gene sets.")
    }
    valid_sets <- vapply(gene_sets, function(x) {
        is.character(x) && length(x) && !anyNA(x) && all(nzchar(x))
    }, logical(1L))
    if (any(!valid_sets)) {
        .seac_abort(
            "Every gene set must contain non-empty, non-missing character values."
        )
    }
    gene_sets <- lapply(gene_sets, unique)
    flat <- unlist(gene_sets, use.names = FALSE)
    duplicate_genes <- unique(flat[duplicated(flat)])
    if (length(duplicate_genes)) {
        .seac_abort(
            "Genes may not belong to multiple sets: ",
            paste(duplicate_genes, collapse = ", "), "."
        )
    }
    available <- unique(c(identifiers$primary, identifiers$aliases))
    available <- available[!is.na(available)]
    matched <- lapply(gene_sets, function(gs) intersect(gs, available))
    missing <- lapply(seq_along(gene_sets), function(i) setdiff(gene_sets[[i]], matched[[i]]))
    names(missing) <- names(gene_sets)
    missing_text <- unlist(lapply(names(missing), function(nm) {
        if (!length(missing[[nm]])) return(NULL)
        paste0(nm, ": ", paste(missing[[nm]], collapse = ", "))
    }), use.names = FALSE)
    if (length(missing_text)) {
        warning(
            "Genes absent from the input were removed: ",
            paste(missing_text, collapse = "; "), call. = FALSE
        )
    }
    empty <- names(matched)[!lengths(matched)]
    if (length(empty)) {
        .seac_abort(
            "No genes from these sets occur in the input: ",
            paste(empty, collapse = ", "), "."
        )
    }
    rows <- lapply(matched, function(gs) {
        which(identifiers$primary %in% gs | identifiers$aliases %in% gs)
    })
    resolved_rows <- unlist(rows, use.names = FALSE)
    if (anyDuplicated(resolved_rows)) {
        .seac_abort(
            "Gene sets resolve to overlapping features through row names ",
            "and gene-symbol aliases."
        )
    }
    list(gene_sets = matched, rows = rows)
}

.seac_nearest_distance <- function(source_coords, query_coords) {
    if (!nrow(source_coords)) {
        .seac_abort("Cannot query a KD-tree with no source coordinates.")
    }
    out <- RANN::nn2(
        data = source_coords,
        query = query_coords,
        k = 1L,
        eps = 0
    )$nn.dists[, 1L]
    ceiling(as.numeric(out))
}

.seac_score_components <- function(spe, gene_sets) {
    if (!inherits(spe, "SpatialExperiment")) {
        .seac_abort("`spe` must be a SpatialExperiment.")
    }
    counts <- .seac_counts(spe)
    identifiers <- .seac_feature_identifiers(spe)
    resolved <- .seac_validate_gene_sets(gene_sets, identifiers)
    gene_sets <- resolved$gene_sets
    coords <- .seac_coordinates(spe)
    distances <- vector("list", length(gene_sets))
    positives <- vector("list", length(gene_sets))
    names(distances) <- names(positives) <- names(gene_sets)

    for (nm in names(gene_sets)) {
        rows <- resolved$rows[[nm]]
        present <- .seac_col_sums(counts[rows, , drop = FALSE]) > 0
        if (!any(present)) {
            .seac_abort("Gene set `", nm, "` has no positive post-QC bins.")
        }
        positives[[nm]] <- present
        distances[[nm]] <- .seac_nearest_distance(coords[present, , drop = FALSE], coords)
    }
    distance_matrix <- do.call(cbind, distances)
    colnames(distance_matrix) <- names(gene_sets)
    list(
        gene_sets = gene_sets,
        gene_rows = resolved$rows,
        coords = coords,
        distances = distance_matrix,
        positives = positives,
        maxmin = apply(distance_matrix, 1L, max)
    )
}

#' Calculate the MaxMin gene-set spatial co-localization statistic
#'
#' @param spe A post-QC `SpatialExperiment`.
#' @param gene_sets Named list of non-overlapping gene vectors.
#'
#' @return `spe` with one distance column per gene set and a `maxmin` column in
#'   `colData`.
#' @export
#'
#' @examples
#' # See the package vignette for a SpatialExperiment example.
maxMinScore <- function(spe, gene_sets) {
    components <- .seac_score_components(spe, gene_sets)
    cd <- SummarizedExperiment::colData(spe)
    distance_names <- paste0(
        "maxmin_distance_", make.names(colnames(components$distances), unique = TRUE)
    )
    for (i in seq_along(distance_names)) {
        cd[[distance_names[[i]]]] <- components$distances[, i]
    }
    cd$maxmin <- as.integer(components$maxmin)
    spe <- SummarizedExperiment::`colData<-`(spe, value = cd)
    md <- S4Vectors::metadata(spe)
    md$SEAC <- utils::modifyList(md$SEAC %||% list(), list(
        maxmin = list(
            gene_sets = components$gene_sets,
            distance_columns = stats::setNames(distance_names, colnames(components$distances)),
            distance = "exact-kdtree-euclidean-ceiling"
        )
    ))
    spe <- S4Vectors::`metadata<-`(spe, value = md)
    spe
}

`%||%` <- function(x, y) if (is.null(x)) y else x
