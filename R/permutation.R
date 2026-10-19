.seac_prepare_sample <- function(spe, gene_sets) {
    components <- .seac_score_components(spe, gene_sets)
    counts <- .seac_counts(spe)
    category_totals <- vapply(components$gene_rows, function(rows) {
        sum(.seac_row_sums(counts[rows, , drop = FALSE]))
    }, numeric(1L))
    if (any(category_totals <= 0)) {
        .seac_abort("Every gene set must have at least one UMI after QC.")
    }
    bin_totals <- .seac_col_sums(counts)
    if (any(bin_totals < 0) || any(!is.finite(bin_totals)) ||
        any(bin_totals != floor(bin_totals)) ||
        any(category_totals != floor(category_totals))) {
        .seac_abort("Counts must produce finite, non-negative integer UMI totals.")
    }
    list(
        spe = spe,
        coords = components$coords,
        gene_sets = components$gene_sets,
        observed = as.integer(components$maxmin),
        bin_totals = bin_totals,
        cumulative = cumsum(bin_totals),
        category_totals = category_totals,
        n_bins = ncol(spe)
    )
}

.seac_permutation_input <- function(prep) {
    prep[c(
        "coords", "observed", "bin_totals", "cumulative",
        "category_totals", "n_bins"
    )]
}

.seac_input_signature <- function(prep) {
    coords <- prep$coords
    list(
        hash = .seac_hash(list(
            barcodes = colnames(prep$spe),
            coordinates = coords,
            bin_totals = prep$bin_totals,
            category_totals = prep$category_totals,
            observed = prep$observed
        )),
        n_features = nrow(prep$spe),
        n_bins = prep$n_bins,
        total_umis = sum(prep$bin_totals),
        category_totals = unname(prep$category_totals),
        coordinate_sum = unname(colSums(coords)),
        coordinate_range = unname(c(range(coords[, 1L]), range(coords[, 2L])))
    )
}

.seac_slots_to_bins <- function(slots, cumulative) {
    findInterval(slots - 1, cumulative) + 1L
}

.seac_draw_category_slots <- function(total_umis, category_totals, seed) {
    category_total <- sum(category_totals)
    if (category_total > total_umis) {
        .seac_abort("Gene-set UMI totals exceed the total UMI count.")
    }
    .seac_with_seed(seed, {
        selected_slots <- sample.int(total_umis, category_total, replace = FALSE)
        labels <- sample(rep.int(seq_along(category_totals), as.integer(category_totals)))
        split(selected_slots, labels)
    })
}

.seac_null_maxmin <- function(prep, category_bins) {
    n_sets <- length(category_bins)
    distances <- vector("list", n_sets)
    for (k in seq_len(n_sets)) {
        bins <- unique(category_bins[[k]])
        if (!length(bins)) {
            # In pooled mode a rare category can be assigned to no bin of one
            # sample. Its nearest distance is infinite, so the lower-tail
            # event cannot be significant for that sample in this iteration.
            distances[[k]] <- rep.int(Inf, prep$n_bins)
            next
        }
        distances[[k]] <- .seac_nearest_distance(
            prep$coords[bins, , drop = FALSE], prep$coords
        )
    }
    do.call(pmax, distances)
}

.seac_within_null <- function(prep, seed) {
    slots <- .seac_draw_category_slots(
        sum(prep$bin_totals), prep$category_totals, seed
    )
    category_bins <- lapply(slots, .seac_slots_to_bins, cumulative = prep$cumulative)
    .seac_null_maxmin(prep, category_bins)
}

.seac_pooled_null <- function(preps, seed) {
    bin_totals <- unlist(lapply(preps, `[[`, "bin_totals"), use.names = FALSE)
    category_totals <- Reduce(`+`, lapply(preps, `[[`, "category_totals"))
    slots <- .seac_draw_category_slots(sum(bin_totals), category_totals, seed)
    global_bins <- lapply(slots, .seac_slots_to_bins, cumulative = cumsum(bin_totals))

    sample_ends <- cumsum(vapply(preps, `[[`, numeric(1L), "n_bins"))
    sample_starts <- c(0L, head(sample_ends, -1L))
    output <- vector("list", length(preps))
    for (s in seq_along(preps)) {
        category_bins <- lapply(global_bins, function(bins) {
            owner <- findInterval(bins - 1L, sample_ends) + 1L
            bins[owner == s] - sample_starts[[s]]
        })
        output[[s]] <- .seac_null_maxmin(preps[[s]], category_bins)
    }
    output
}

.seac_permutation_block <- function(preps, ids, mode, seed) {
    n_samples <- length(preps)
    n_less <- lapply(preps, function(x) integer(x$n_bins))
    n_leq <- lapply(preps, function(x) integer(x$n_bins))
    for (id in ids) {
        # Per-permutation seed makes results independent of block boundaries.
        # Use double arithmetic before reducing into set.seed()'s integer range
        # so a valid large seed cannot overflow to NA.
        this_seed <- as.integer((as.double(seed) + id) %% .Machine$integer.max)
        null_scores <- if (identical(mode, "pooled")) {
            .seac_pooled_null(preps, this_seed)
        } else {
            # Match upstream per-sample permutations: permutation `id` has the
            # same seed for each sample, irrespective of list order or whether
            # a sample is run alone. Each sample draws only from its own UMIs.
            lapply(preps, .seac_within_null, seed = this_seed)
        }
        for (s in seq_len(n_samples)) {
            n_less[[s]] <- n_less[[s]] + as.integer(null_scores[[s]] < preps[[s]]$observed)
            n_leq[[s]] <- n_leq[[s]] + as.integer(null_scores[[s]] <= preps[[s]]$observed)
        }
    }
    list(n_less = n_less, n_leq = n_leq, completed = ids)
}

.seac_checkpoint_files <- function(checkpoint_dir) {
    list.files(checkpoint_dir, pattern = "^checkpoint_[0-9]+\\.rds$", full.names = TRUE)
}

.seac_read_checkpoint <- function(checkpoint_dir, specification) {
    files <- .seac_checkpoint_files(checkpoint_dir)
    if (!length(files)) return(NULL)
    ids <- as.integer(sub(".*checkpoint_([0-9]+)\\.rds$", "\\1", files))
    for (file in files[order(ids, decreasing = TRUE)]) {
        checkpoint <- tryCatch(readRDS(file), error = function(e) NULL)
        if (is.null(checkpoint)) next
        if (!identical(checkpoint$specification, specification)) {
            .seac_abort("Checkpoint parameters do not match the requested permutation run.")
        }
        return(checkpoint)
    }
    NULL
}

.seac_write_checkpoint <- function(checkpoint_dir, completed, n_less, n_leq,
                                   specification) {
    checkpoint <- list(
        completed = sort(unique(as.integer(completed))),
        n_less = n_less,
        n_leq = n_leq,
        specification = specification
    )
    checkpoint_id <- if (length(checkpoint$completed)) {
        max(checkpoint$completed)
    } else {
        0L
    }
    destination <- file.path(
        checkpoint_dir, sprintf("checkpoint_%06d.rds", checkpoint_id)
    )
    temporary <- tempfile("checkpoint-", tmpdir = checkpoint_dir, fileext = ".rds")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    saveRDS(checkpoint, temporary)
    if (!file.rename(temporary, destination)) {
        .seac_abort("Could not finalize checkpoint file: ", destination)
    }
    invisible(checkpoint)
}

.seac_write_results_h5 <- function(results, path) {
    if (file.exists(path)) {
        .seac_abort("Refusing to overwrite existing result file: ", path)
    }
    temporary <- tempfile("seac-results-", tmpdir = dirname(path), fileext = ".h5")
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
    rhdf5::h5createFile(temporary)
    for (sample_name in names(results)) {
        rhdf5::h5createGroup(temporary, sample_name)
        for (nm in names(results[[sample_name]])) {
            rhdf5::h5write(
                results[[sample_name]][[nm]], temporary,
                paste0(sample_name, "/", nm)
            )
        }
    }
    if (!file.rename(temporary, path)) {
        .seac_abort("Could not finalize result file: ", path)
    }
    invisible(path)
}

.seac_result_table <- function(prep, n_less, n_leq, n_permutations,
                               legacy_alpha, fdr_cutoff) {
    p_value <- (1 + n_leq) / (1 + n_permutations)
    adjusted <- stats::p.adjust(p_value, method = "BH")
    data.frame(
        barcode = colnames(prep$spe),
        maxmin = prep$observed,
        n_null_less = as.integer(n_less),
        n_null_leq = as.integer(n_leq),
        p_value = p_value,
        p_adjusted = adjusted,
        fdr_significant = adjusted <= fdr_cutoff,
        legacy_hotspot = n_less < ceiling(legacy_alpha * n_permutations),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

.seac_normalize_inputs <- function(x) {
    if (inherits(x, "SpatialExperiment") || (is.character(x) && length(x) == 1L)) {
        return(list(sample1 = .seac_load_post_qc(x)))
    }
    if (!is.list(x) || !length(x)) {
        .seac_abort("`x` must be a SpatialExperiment, post-QC directory, or named list.")
    }
    if (is.null(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
        .seac_abort("A list input must have unique, non-empty sample names.")
    }
    if (any(grepl("/", names(x), fixed = TRUE))) {
        .seac_abort("Sample names may not contain `/` because they name HDF5 groups.")
    }
    lapply(x, .seac_load_post_qc)
}

#' Run the MaxMin UMI-level permutation test
#'
#' @param x A post-QC `SpatialExperiment`, saved output directory from
#'   [runBladeQC()], or named list of either for multi-sample testing. This
#'   function does not run BLADE; perform QC and save the output first.
#' @param gene_sets Named list of non-overlapping gene vectors.
#' @param mode Whether to permute category labels within each sample or pooled
#'   across multiple samples. Within-sample results for one sample do not
#'   depend on which other samples are included or their order.
#' @param n_permutations Number of UMI-level shuffles.
#' @param legacy_alpha Lower-tail threshold used by the original pipeline.
#' @param fdr_cutoff BH adjusted p-value cutoff.
#' @param seed Integer seed.
#' @param chunk_size Number of permutations per checkpoint.
#' @param output_dir New destination directory.
#' @param resume Resume a compatible incomplete output directory.
#' @param BPPARAM A BiocParallel parameter object.
#'
#' @return A `SpatialExperiment`, or named list of them, with result columns in
#'   `colData`.
#' @export
maxMinPermutationTest <- function(x, gene_sets,
                                  mode = c("within", "pooled"),
                                  n_permutations = 5000L,
                                  legacy_alpha = 0.01,
                                  fdr_cutoff = 0.05,
                                  seed = 1L,
                                  chunk_size = 100L,
                                  output_dir = NULL,
                                  resume = TRUE,
                                  BPPARAM = BiocParallel::SerialParam()) {
    mode <- match.arg(mode)
    .seac_assert_whole_number(n_permutations, "n_permutations", minimum = 1)
    .seac_assert_whole_number(chunk_size, "chunk_size", minimum = 1)
    .seac_assert_whole_number(seed, "seed")
    if (!.seac_is_scalar_number(legacy_alpha) || legacy_alpha <= 0 ||
        legacy_alpha >= 1 || !.seac_is_scalar_number(fdr_cutoff) ||
        fdr_cutoff <= 0 || fdr_cutoff >= 1) {
        .seac_abort("`legacy_alpha` and `fdr_cutoff` must be between zero and one.")
    }

    samples <- .seac_normalize_inputs(x)
    if (identical(mode, "pooled") && length(samples) < 2L) {
        .seac_abort("`mode = 'pooled'` requires at least two samples.")
    }
    preps <- lapply(samples, .seac_prepare_sample, gene_sets = gene_sets)
    permutation_inputs <- lapply(preps, .seac_permutation_input)
    if (is.null(output_dir)) {
        output_dir <- file.path(getwd(), "SEAC_permutation")
    }
    .seac_assert_scalar_character(output_dir, "output_dir")
    checkpoint_dir <- file.path(output_dir, "checkpoints")
    specification <- list(
        version = if (identical(mode, "within")) 3L else 2L,
        sample_names = names(samples),
        input_signature = lapply(preps, .seac_input_signature),
        n_bins = vapply(preps, `[[`, numeric(1L), "n_bins"),
        gene_sets = lapply(preps, `[[`, "gene_sets"),
        mode = mode,
        n_permutations = as.integer(n_permutations),
        legacy_alpha = legacy_alpha,
        fdr_cutoff = fdr_cutoff,
        seed = as.integer(seed),
        rng_version = "4.0.0",
        kdtree = list(k = 1L, eps = 0, distance = "euclidean-ceiling")
    )
    if (identical(mode, "within")) {
        specification$within_seed_scheme <-
            "permutation-id-independent-of-sample-order"
    }

    if (dir.exists(output_dir)) {
        if (!isTRUE(resume)) {
            .seac_abort("Output directory exists; use `resume = TRUE` for a matching run.")
        }
        checkpoint <- .seac_read_checkpoint(checkpoint_dir, specification)
        if (is.null(checkpoint)) {
            .seac_abort("Output directory exists but has no compatible checkpoint.")
        }
    } else {
        dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
        checkpoint <- list(
            completed = integer(),
            n_less = lapply(preps, function(z) integer(z$n_bins)),
            n_leq = lapply(preps, function(z) integer(z$n_bins)),
            specification = specification
        )
        .seac_write_manifest(file.path(output_dir, "manifest.json"), list(
            package = "SEAC", stage = "permutation", specification = specification,
            created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
        ))
        # A zero-progress checkpoint makes even an interruption during the
        # first block resumable.
        .seac_write_checkpoint(
            checkpoint_dir, checkpoint$completed, checkpoint$n_less,
            checkpoint$n_leq, specification
        )
    }

    pending <- setdiff(seq_len(as.integer(n_permutations)), checkpoint$completed)
    if (length(pending)) {
        blocks <- split(pending, ceiling(seq_along(pending) / as.integer(chunk_size)))
        # Process one checkpoint block at a time so completed work is durable.
        # Independent per-permutation seeds make worker scheduling irrelevant.
        for (ids in blocks) {
            workers <- BiocParallel::bpworkers(BPPARAM)
            if (length(workers) != 1L || is.na(workers) || workers < 1L) {
                workers <- 1L
            }
            workers <- min(length(ids), as.integer(workers))
            task_ids <- if (workers == 1L) {
                list(ids)
            } else {
                split(ids, cut(seq_along(ids), workers, labels = FALSE))
            }
            task_results <- BiocParallel::bplapply(task_ids, function(block_ids) {
                .seac_permutation_block(
                    permutation_inputs, block_ids, mode, as.integer(seed)
                )
            }, BPPARAM = BPPARAM)
            for (result in task_results) {
                checkpoint$n_less <- Map(`+`, checkpoint$n_less, result$n_less)
                checkpoint$n_leq <- Map(`+`, checkpoint$n_leq, result$n_leq)
                checkpoint$completed <- sort(unique(c(checkpoint$completed, result$completed)))
            }
            .seac_write_checkpoint(
                checkpoint_dir, checkpoint$completed, checkpoint$n_less,
                checkpoint$n_leq, specification
            )
        }
    }

    results <- Map(.seac_result_table, preps, checkpoint$n_less, checkpoint$n_leq,
                   MoreArgs = list(
                       n_permutations = as.integer(n_permutations),
                       legacy_alpha = legacy_alpha,
                       fdr_cutoff = fdr_cutoff
                   ))
    result_path <- file.path(output_dir, "results.h5")
    if (!file.exists(result_path)) {
        .seac_write_results_h5(results, result_path)
    }

    output <- Map(function(prep, result) {
        spe <- prep$spe
        cd <- SummarizedExperiment::colData(spe)
        for (nm in setdiff(names(result), "barcode")) cd[[nm]] <- result[[nm]]
        spe <- SummarizedExperiment::`colData<-`(spe, value = cd)
        spe
    }, preps, results)
    if (length(output) == 1L) output[[1L]] else output
}

#' Load compact SEAC permutation results
#'
#' @param output_dir Permutation output directory.
#' @param sample_id Sample name; defaults to the only sample in the result file.
#' @param barcodes Optional barcode vector for subsetting.
#'
#' @return A data frame of per-bin test results.
#' @export
loadSEACResults <- function(output_dir, sample_id = NULL, barcodes = NULL) {
    path <- file.path(output_dir, "results.h5")
    if (!file.exists(path)) .seac_abort("No results.h5 found in: ", output_dir)
    listing <- rhdf5::h5ls(path, recursive = FALSE)
    samples <- listing$name[listing$otype == "H5I_GROUP"]
    if (is.null(sample_id)) {
        if (length(samples) != 1L) {
            .seac_abort("Specify `sample_id`; the result contains multiple samples.")
        }
        sample_id <- samples[[1L]]
    }
    if (!sample_id %in% samples) .seac_abort("Unknown sample_id: ", sample_id)
    columns <- rhdf5::h5ls(path, recursive = TRUE)
    columns <- columns[columns$group == paste0("/", sample_id) &
                           columns$otype == "H5I_DATASET", "name"]
    result <- lapply(columns, function(nm) rhdf5::h5read(path, paste0(sample_id, "/", nm)))
    names(result) <- columns
    result <- as.data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
    if (!is.null(barcodes)) result <- result[result$barcode %in% barcodes, , drop = FALSE]
    result
}
