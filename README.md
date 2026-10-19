# SEAC

SEAC is an R package for 10x Genomics Visium HD data. It provides:

- a KD-tree BLADE mask that removes capture-border, non-tissue, and adjacent bins;
- a MaxMin score for spatial co-localization of two or more gene sets; and
- an optional UMI-level permutation test with checkpoints and HDF5 results.

The R package and its GitHub repository are both named `SEAC`.

## Install from GitHub

SEAC is not currently distributed through Bioconductor. Install it directly
from GitHub. It requires R 4.4 or later and several Bioconductor dependencies;
the following setup lets R resolve those dependencies before installing SEAC:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
options(repos = BiocManager::repositories())

if (!requireNamespace("remotes", quietly = TRUE)) {
    BiocManager::install("remotes")
}
remotes::install_github(
    "KummerfeldLab/SEAC",
    build_vignettes = FALSE
)
```

Then load the package:

```r
library(SEAC)
```

## Install a local copy

After cloning or downloading this repository, use the same setup above and run:

```r
remotes::install_local(
    "/path/to/SEAC",
    build_vignettes = FALSE
)
```

## Input

Pass a Space Ranger run directory, its `outs` directory, or one Visium HD bin
directory. The selected bin directory must contain:

```text
square_016um/
├── filtered_feature_bc_matrix.h5
└── spatial/
    └── tissue_positions.parquet
```

Use the matching bin size, for example `"002"`, `"008"`, or `"016"`.

## Minimal workflow

Run BLADE explicitly. `runBladeQC()` saves the full mask and a disk-backed,
filtered `SpatialExperiment` to `sample_SEAC_post_qc`. You can inspect this
saved result before starting permutations.

```r
library(SEAC)

qc_dir <- "sample_SEAC_post_qc"
qc <- runBladeQC(
    input_dir = "/path/to/spaceranger/outs",
    output_dir = qc_dir,
    bin_size = "016"
)

dim(qc)
```

The saved object is at `sample_SEAC_post_qc/spe`; the complete BLADE mask is
at `sample_SEAC_post_qc/blade_mask.h5`. In a later R session, pass `qc_dir`
directly to the permutation function. It reads the saved object and does not
run BLADE again.

Define at least two non-overlapping gene sets. Replace the mouse small-intestine
examples below with genes present in your data. Gene names may be Ensembl
feature IDs or 10x gene symbols.

```r
gene_sets <- list(
    Paneth = "Lyz1",
    Goblet = c("Muc2", "Clca1"),
    Enterocyte = c("Fabp1", "Fabp2", "Apoa4")
)

tested <- maxMinPermutationTest(
    x = qc_dir,
    gene_sets = gene_sets,
    output_dir = "sample_SEAC_permutation"
)

results <- loadSEACResults("sample_SEAC_permutation")
head(results)
```

The default is 5,000 permutations. `legacy_hotspot` applies the upstream
VisHD per-bin lower-tail threshold; `p_adjusted` and `fdr_significant` provide
the additional BH-adjusted result. For a quick workflow smoke test, reduce
`n_permutations` and `chunk_size`, using a new output directory:

```r
maxMinPermutationTest(
    x = qc_dir,
    gene_sets = gene_sets,
    n_permutations = 20L,
    chunk_size = 5L,
    seed = 1L,
    output_dir = "sample_SEAC_smoke"
)
```

Twenty permutations only check that the workflow runs; they cannot support
fine-grained p-values or biological claims. For multiple saved samples, pass
a named list of QC directories:

```r
qc_dirs <- list(sample_A = "sample_A_post_qc", sample_B = "sample_B_post_qc")
maxMinPermutationTest(
    x = qc_dirs,
    gene_sets = gene_sets,
    mode = "within",  # use "pooled" for a shared cross-sample null
    output_dir = "two_sample_permutation"
)
```

In `within` mode, each sample's category totals are preserved. Running one
sample alone or alongside other samples with the same seed gives the same
result for that sample. In `pooled` mode, category totals are preserved across
the full cohort, so results depend on the cohort supplied.

Interrupted permutation runs are checkpointed. Repeat the same call with the
same output directory and `resume = TRUE` to continue it. The updated
within-sample seed scheme cannot resume within-sample checkpoints made with
SEAC 0.99.1; start those analyses in a new output directory. Pooled-mode
checkpoints are unchanged.

## Real-data validation

The staged workflow was tested with the official 10x Genomics Visium HD
Mouse Small Intestine binned outputs. The 16 µm test input contained 19,059
features, 91,033 filtered bins, and 138,666,717 UMIs. SEAC retained 83,408 bins
and successfully wrote and reloaded its post-QC object, checkpoints, and
permutation results. A full run with SEAC 0.99.3 then completed 5,000
within-sample permutations on the retained bins using the three example gene
sets above (seed `20261019`, chunk size 100, eight workers). Validation checked
that all 5,000 permutations were recorded, all 83,408 result rows were
reloadable, and the score, permutation-count, and p-value fields met their
expected bounds. This run produced 19,319 `legacy_hotspot` calls and zero
BH-FDR-significant bins. These are workflow-validation results for the selected
gene sets, not evidence of a biological discovery or validation of every
possible gene-set analysis. No 10x data are included in this repository.

The source data are described in the [10x Genomics Visium HD analysis
guide](https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd).

For details, see the package vignette:

```r
browseVignettes("SEAC")
```
