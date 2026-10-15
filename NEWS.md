# SEAC 0.99.1 (2026-10-15)

- Validates the complete workflow on non-empty official 10x Genomics Visium HD
  Mouse Small Intestine 16 micrometer binned output.
- Persists MaxMin score columns and metadata in the HDF5-backed post-QC object,
  so they remain available after a new R session reloads the output.
- Adds direct GitHub and local installation instructions and a minimal workflow
  tutorial without bundling validation data.

# SEAC 0.99.0 (2026-10-01)

- Initial Bioconductor submission candidate.
- Adds KD-tree radius-one BLADE QC for Visium HD.
- Adds general MaxMin scoring and UMI-level within/pooled permutation tests.
- Imports Visium HD positions with lightweight `nanoparquet` and counts with
  `TENxIO`, avoiding the Apache Arrow C++ runtime while retaining disk-backed
  assays.
- Resolves gene sets against both feature row names and 10x `Symbol` aliases,
  including duplicate-symbol feature rows and rejecting cross-set alias
  overlap.
- Uses content-hashed, atomically written permutation checkpoints with fallback
  from an interrupted latest write.
- Reduces parallel-worker transfer size by excluding the full expression object
  from permutation tasks, and tightens barcode and gene-set validation.
