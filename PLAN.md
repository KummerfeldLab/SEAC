# SEAC implementation plan

## Purpose

SEAC is a pure-R Bioconductor-oriented package for Visium HD data.  It
implements the KD-tree BLADE mask used by the supplied `VisHD-main` code and a
general, user-supplied gene-set MaxMin spatial co-localisation statistic with a
UMI-level Monte Carlo permutation test.

## Confirmed method contract

### KD-tree BLADE mask

For every row in the complete `tissue_positions.parquet`, a seed is any
position on the rectangular capture boundary or with `in_tissue == 0`.
Coordinates of all seeds form one KD-tree.  Every position is queried against
that tree using exact Euclidean nearest-neighbour distance.  It is masked when
that distance is at most one array-coordinate unit.  This is the behaviour of
`SEAC.py::add_blade_column()` and deliberately is **not** a one-spot iterative
neighbour operation or a k-nearest-neighbour (`k = 8`) mode.

### MaxMin

For each named gene set, find post-QC bins with at least one UMI from that
set, build a KD-tree from their array coordinates, and query every retained
bin.  Distances are rounded up.  The MaxMin statistic at a bin is the maximum
of those per-set distances; smaller values indicate tighter spatial
co-occurrence of all sets.

### Permutation null

Each UMI is labelled as one of the supplied, non-overlapping gene-set
categories or `Other`.  Category labels are shuffled while bin-level UMI
totals remain fixed.  `within` permutations preserve category totals per
sample; `pooled` permutations preserve category totals across all supplied
samples, matching the newer Python pipeline.  The default is 5,000
permutations in blocks of 100.  Results are accumulated rather than retaining
the full null matrix.

## User workflow

1. `runBladeQC()` imports a Space Ranger root, `outs`, or `square_XXXum`
   directory, calculates the mask, filters counts, and writes a HDF5-backed
   `SpatialExperiment` under `post_qc/`.
2. `maxMinScore()` calculates a statistic on a post-QC object with any valid
   named gene-set list.
3. `maxMinPermutationTest()` optionally runs the 5,000 UMI-level
   permutations, supports content-validated atomic restart checkpoints, and
   writes compact HDF5 results under `permutation/`.
4. `runSEAC()` combines steps 1--3, with permutation disabled by default.

## Package and review principles

- No Python runtime, reticulate, or custom C/C++ source is shipped.
- Counts remain sparse/HDF5-backed where supported by the upstream importer.
- `RANN` supplies exact KD-tree queries; `nanoparquet` reads the flat Visium HD
  position table; `TENxIO` imports the disk-backed HDF5 count matrix; and
  `HDF5Array` stores the post-QC object. This avoids a runtime dependency on
  Apache Arrow for a single small metadata table.
- The built-in `seacGeneSets` dataset is an example, never an implicit
  analysis default.
- Cross-set duplicate genes are errors; missing genes warn and are removed;
  an empty or unobserved set is an error.
- Gene sets resolve against both feature row names and 10x `Symbol` aliases;
  every row sharing a requested symbol contributes, while alias-based overlap
  between sets is an error.
- The package describes its BLADE component as the local KD-tree HD mask, not
  as a replacement for the full published BLADE artefact-detection suite.
- Permutations use the R 4.0.0 RNG contract with independent per-permutation
  seeds, so results remain stable across chunking, worker scheduling, and
  future supported R versions.

## Deliverables

- Exported functions: `bladeMaskHD`, `runBladeQC`, `maxMinScore`,
  `maxMinPermutationTest`, `runSEAC`, and `loadSEACResults`.
- Unit tests for KD-tree masking, MaxMin, permutation invariants, resume,
  validation, and HDF5 results.
- A small bundled SEAC gene-set example and an evaluated Bioconductor-style
  vignette using an in-memory toy `SpatialExperiment`.
- `R CMD check` and BiocCheck configuration using the fixed MSI R module.
