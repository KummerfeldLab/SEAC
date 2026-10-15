test_that("KD-tree BLADE radius one agrees with the grid contract", {
    pos <- expand.grid(array_row = 1:7, array_col = 1:7)
    pos$barcode <- paste0("b", seq_len(nrow(pos)))
    pos$in_tissue <- 1L
    pos$in_tissue[pos$array_row == 4 & pos$array_col == 4] <- 0L
    out <- bladeMaskHD(pos[, c("barcode", "in_tissue", "array_row", "array_col")])
    center <- out$array_row == 4 & out$array_col == 4
    axial <- out$array_row == 3 & out$array_col == 4
    diagonal <- out$array_row == 3 & out$array_col == 3
    expect_true(out$blade[center])
    expect_true(out$blade[axial])
    expect_false(out$blade[diagonal])
    expect_equal(out$nearest_seed_distance[center], 0)
    expect_equal(out$nearest_seed_distance[axial], 1)
    expect_equal(out$nearest_seed_distance[diagonal], sqrt(2))
})

test_that("BLADE validates positions", {
    expect_error(bladeMaskHD(data.frame(barcode = "x")), "missing required")
    bad <- data.frame(
        barcode = c("x", "y"), in_tissue = c(0L, 2L),
        array_row = c(1L, 2L), array_col = c(1L, 2L)
    )
    expect_error(bladeMaskHD(bad), "only zero and one")

    invalid_barcode <- transform(bad, in_tissue = c(0L, 1L))
    invalid_barcode$barcode[1L] <- NA_character_
    expect_error(bladeMaskHD(invalid_barcode), "non-empty character")

    invalid_barcode$barcode <- c("x", "x")
    expect_error(bladeMaskHD(invalid_barcode), "must be unique")
})
