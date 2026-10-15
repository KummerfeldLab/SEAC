#' Run the SEAC Visium HD workflow
#'
#' @param input_dir A Space Ranger root, `outs`, or `square_XXXum` directory.
#' @param output_dir A new root output directory.
#' @param bin_size Visium HD bin size.
#' @param permutation Whether to run the optional MaxMin permutation test.
#' @param gene_sets Optional named gene sets. When supplied, MaxMin scores are
#'   added to `post_qc`; they are required when `permutation = TRUE`.
#' @param ... Arguments forwarded to [maxMinPermutationTest()].
#'
#' @return A list containing `post_qc` and, when requested, `permutation`.
#'   When gene sets are supplied, the scored post-QC metadata are also
#'   persisted in `output_dir/post_qc/spe`.
#' @export
runSEAC <- function(input_dir, output_dir = NULL, bin_size = "002",
                    permutation = FALSE, gene_sets = NULL, ...) {
    .seac_assert_scalar_character(input_dir, "input_dir")
    if (!is.logical(permutation) || length(permutation) != 1L ||
        is.na(permutation)) {
        .seac_abort("`permutation` must be TRUE or FALSE.")
    }
    if (permutation && is.null(gene_sets)) {
        .seac_abort("Supply `gene_sets` when `permutation = TRUE`.")
    }
    sample_id <- .seac_sample_id_from_input(input_dir)
    if (is.null(output_dir)) output_dir <- file.path(getwd(), paste0(sample_id, "_SEAC"))
    .seac_assert_scalar_character(output_dir, "output_dir")
    if (dir.exists(output_dir) || file.exists(output_dir)) {
        .seac_abort("Output already exists: ", output_dir)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    qc <- runBladeQC(
        input_dir = input_dir,
        output_dir = file.path(output_dir, "post_qc"),
        bin_size = bin_size,
        sample_id = sample_id
    )
    scored <- if (is.null(gene_sets)) {
        qc
    } else {
        .seac_persist_post_qc_metadata(
            maxMinScore(qc, gene_sets),
            file.path(output_dir, "post_qc")
        )
    }
    result <- list(post_qc = scored)
    if (permutation) {
        result$permutation <- maxMinPermutationTest(
            x = scored,
            gene_sets = gene_sets,
            output_dir = file.path(output_dir, "permutation"),
            ...
        )
    }
    .seac_write_manifest(file.path(output_dir, "manifest.json"), list(
        package = "SEAC",
        stage = "workflow",
        input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
        bin_size = as.character(bin_size),
        maxmin_scored = !is.null(gene_sets),
        permutation = permutation,
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ))
    result
}
