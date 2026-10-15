make_toy_spe <- function() {
    counts <- rbind(
        A = c(1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L),
        B = c(0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L),
        C = c(0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L, 0L),
        D = c(0L, 0L, 0L, 1L, 0L, 0L, 0L, 1L, 0L)
    )
    colnames(counts) <- paste0("b", 1:9)
    cd <- S4Vectors::DataFrame(
        array_row = rep(1:3, each = 3),
        array_col = rep(1:3, 3),
        row.names = colnames(counts)
    )
    SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = cd,
        spatialCoordsNames = c("array_col", "array_row")
    )
}

toy_gene_sets <- list(first = c("A", "B"), second = c("C", "D"))
