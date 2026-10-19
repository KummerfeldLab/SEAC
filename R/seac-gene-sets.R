#' Example gene sets used by the original SEAC analysis
#'
#' This object is example data only. Users must supply `gene_sets` explicitly
#' to MaxMin functions so the statistic is not tied to this list.
#'
#' @format A named list of character vectors.
#' @export
#' @examples
#' names(seacGeneSets)
#' seacGeneSets[["SASP"]]
seacGeneSets <- list(
    SASP = c("IL1B", "IL6", "CXCL12", "CCL2", "CXCL8", "TGFB1", "VEGFA"),
    EMC = c("PLAUR", "PLAU", "PAPPA", "MMP14", "ADAMTS4", "COL6A1", "THBS1"),
    Anti_Apoptotic = c("SERPINE1", "SERPINB2", "SERPINE2", "BCL2L1", "MCL1", "BCL2"),
    Cell_Cycle_Related = c("CDKN1A", "TP53", "GADD45A", "CDKN2A", "GADD45B")
)
