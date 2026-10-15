# SEAC

SEAC is an R package for 10x Genomics Visium HD data. It provides:

- a KD-tree BLADE mask that removes capture-border, non-tissue, and adjacent bins;
- a MaxMin score for spatial co-localization of two or more gene sets; and
- an optional UMI-level permutation test with checkpoints and HDF5 results.

The R package is named `SEAC`; its GitHub repository is `VisHD_R`.

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
    "KummerfeldLab/VisHD_R",
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
    "/path/to/VisHD_R",
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

Define at least two non-overlapping gene sets. Replace the mouse small-intestine
examples below with genes present in your data. Gene names may be Ensembl
feature IDs or 10x gene symbols.

```r
library(SEAC)

gene_sets <- list(
    Paneth = "Lyz1",
    Goblet = c("Muc2", "Clca1"),
    Enterocyte = c("Fabp1", "Fabp2", "Apoa4")
)

fit <- runSEAC(
    input_dir = "/path/to/spaceranger/outs",
    output_dir = "sample_SEAC",
    bin_size = "016",
    gene_sets = gene_sets,
    permutation = FALSE
)

head(as.data.frame(SummarizedExperiment::colData(fit$post_qc)))
```

The scored post-QC object is saved at `sample_SEAC/post_qc/spe` and can be
reopened in a new R session:

```r
post_qc <- HDF5Array::loadHDF5SummarizedExperiment(
    "sample_SEAC/post_qc/spe"
)
```

## Permutation test

For a quick smoke test, start with a small number of permutations:

```r
fit <- runSEAC(
    input_dir = "/path/to/spaceranger/outs",
    output_dir = "sample_SEAC_smoke",
    bin_size = "016",
    gene_sets = gene_sets,
    permutation = TRUE,
    n_permutations = 20L,
    chunk_size = 5L,
    seed = 1L
)

results <- loadSEACResults("sample_SEAC_smoke/permutation")
head(results)
```

Use substantially more permutations, such as the default 5,000, for a formal
analysis. A small smoke test only checks that the workflow runs; it cannot
support fine-grained p-values or biological claims.

Interrupted permutation runs are checkpointed. Repeat the same call with the
same output directory and `resume = TRUE` to continue it.

## Real-data smoke test

The end-to-end workflow was tested with the official 10x Genomics Visium HD
Mouse Small Intestine binned outputs. The 16 µm test input contained 19,059
features, 91,033 filtered bins, and 138,666,717 UMIs. SEAC retained 83,408 bins
and successfully wrote and reloaded its post-QC object, checkpoints, and
permutation results. No 10x data are included in this repository.

The source data are described in the [10x Genomics Visium HD analysis
guide](https://www.10xgenomics.com/analysis-guides/segmentation-visium-hd).

For details, see the package vignette:

```r
browseVignettes("SEAC")
```
