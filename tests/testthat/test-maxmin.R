test_that("MaxMin adds exact KD-tree distance columns", {
    spe <- make_toy_spe()
    expect_false(all(c("array_row", "array_col") %in%
        names(SummarizedExperiment::colData(spe))))
    scored <- maxMinScore(spe, toy_gene_sets)
    expect_true("maxmin" %in% names(SummarizedExperiment::colData(scored)))
    expect_equal(length(scored$maxmin), ncol(scored))
    expect_equal(scored$maxmin, c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L))
})

test_that("gene-set overlap and empty matches are rejected", {
    spe <- make_toy_spe()
    expect_error(maxMinScore(spe, list(one = "A", two = "A")), "multiple")
    expect_warning(
        expect_error(
            maxMinScore(spe, list(one = "A", two = "NO_GENE")), "No genes"
        ),
        "Genes absent from the input were removed"
    )
    expect_error(
        maxMinScore(spe, list(one = "A", two = NA_character_)),
        "non-empty, non-missing"
    )
    expect_error(
        maxMinScore(spe, list(one = "A", two = numeric())),
        "non-empty, non-missing"
    )
})

test_that("gene sets match TENx gene-symbol aliases and all duplicate symbols", {
    spe <- make_toy_spe()
    SummarizedExperiment::rowData(spe)$Symbol <- c("X", "X", "Y", "Y")
    components <- SEAC:::.seac_score_components(
        spe, list(first = "X", second = "Y")
    )
    counts <- SummarizedExperiment::assay(spe, "counts")

    expect_equal(
        components$positives$first,
        unname(colSums(counts[c("A", "B"), , drop = FALSE]) > 0)
    )
    expect_equal(
        components$positives$second,
        unname(colSums(counts[c("C", "D"), , drop = FALSE]) > 0)
    )
})

test_that("aliases cannot place one feature in multiple gene sets", {
    spe <- make_toy_spe()
    SummarizedExperiment::rowData(spe)$Symbol <- c("X", "B", "C", "D")
    expect_error(
        maxMinScore(spe, list(primary_id = "A", symbol = "X")),
        "overlapping features"
    )
})
