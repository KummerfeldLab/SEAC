\name{SEAC}
\alias{SEAC}
\alias{bladeMaskHD}
\alias{runBladeQC}
\alias{maxMinScore}
\alias{maxMinPermutationTest}
\alias{loadSEACResults}
\alias{seacGeneSets}
\title{KD-tree quality control and MaxMin testing for Visium HD}
\description{SEAC imports Visium HD Space Ranger outputs, applies a KD-tree radius-one BLADE mask, and tests general gene-set spatial co-localization with a UMI-level permutation null model.}
\usage{
bladeMaskHD(positions)
runBladeQC(input_dir, output_dir = NULL, bin_size = "002", sample_id = NULL)
maxMinScore(spe, gene_sets)
maxMinPermutationTest(x, gene_sets, mode = c("within", "pooled"), n_permutations = 5000L, legacy_alpha = 0.01, fdr_cutoff = 0.05, seed = 1L, chunk_size = 100L, output_dir = NULL, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
loadSEACResults(output_dir, sample_id = NULL, barcodes = NULL)
seacGeneSets
}
\arguments{
\item{positions}{A tissue positions data frame.}
\item{input_dir}{A Space Ranger root, \code{outs}, or Visium HD bin directory.}
\item{output_dir}{A new output directory.}
\item{bin_size}{A square Visium HD bin size.}
\item{sample_id}{Optional sample name.}
\item{spe}{A post-QC \code{SpatialExperiment}.}
\item{gene_sets}{Named non-overlapping character vectors of feature row
names or 10x gene symbols from \code{rowData(x)$Symbol}.}
\item{x}{A post-QC object, a saved \code{runBladeQC()} output directory,
or a named list of either. The permutation function never runs BLADE.}
\item{mode}{Within-sample or pooled-sample label permutations.}
\item{n_permutations}{Number of Monte Carlo permutations.}
\item{legacy_alpha}{Original lower-tail hotspot threshold.}
\item{fdr_cutoff}{BH adjusted p-value cutoff.}
\item{seed}{Integer random seed.}
\item{chunk_size}{Permutations per checkpoint.}
\item{resume}{Whether a matching interrupted run may resume.}
\item{BPPARAM}{A BiocParallel backend.}
\item{barcodes}{Optional result barcode subset.}
}
\value{\code{bladeMaskHD} returns a mask data frame. \code{runBladeQC} saves and returns a post-QC \code{SpatialExperiment}; \code{maxMinScore} returns a scored object. The permutation function returns a scored object or named list and saves results separately, while \code{loadSEACResults} returns a data frame.}
\details{The BLADE mask constructs a KD-tree from capture-border and non-tissue seeds, then marks every bin whose nearest seed is at Euclidean array-coordinate distance at most one. Call \code{runBladeQC()} explicitly and pass its saved output directory to \code{maxMinPermutationTest()}. MaxMin is the maximum across input gene sets of the ceil-rounded nearest positive-bin distance. Lower values are more spatially co-localized. Gene sets resolve against both feature row names and 10x \code{Symbol} aliases; all rows sharing a requested symbol contribute. In within-sample mode, each sample preserves its own category totals and its result is independent of the other samples supplied. Within-sample checkpoints made with version 0.99.1 cannot be resumed under the updated seed scheme; pooled checkpoints are unchanged.}
\examples{
positions <- data.frame(
    barcode = paste0("b", 1:4), in_tissue = c(0L, 1L, 1L, 1L),
    array_row = c(1L, 1L, 2L, 2L), array_col = c(1L, 2L, 1L, 2L)
)
bladeMaskHD(positions)

# The small VisiumIO example exercises import and BLADE QC.
if (requireNamespace("VisiumIO", quietly = TRUE)) {
    input <- system.file("extdata", package = "VisiumIO")
    qc_output <- tempfile("seac-qc-")
    qc <- runBladeQC(input, qc_output, bin_size = "002",
                     sample_id = "example")
    dim(qc)
}

# A post-QC toy object with positive counts exercises MaxMin and permutations.
counts <- rbind(A = c(1L, 0L, 0L, 1L),
                B = c(0L, 1L, 1L, 0L))
colnames(counts) <- paste0("bin", 1:4)
coords <- S4Vectors::DataFrame(
    array_col = c(1L, 2L, 1L, 2L),
    array_row = c(1L, 1L, 2L, 2L),
    row.names = colnames(counts)
)
spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts), colData = coords,
    spatialCoordsNames = c("array_col", "array_row")
)
sets <- list(first = "A", second = "B")
maxMinScore(spe, sets)$maxmin
permutation_output <- tempfile("seac-permutations-")
maxMinPermutationTest(spe, sets, n_permutations = 3L,
                      output_dir = permutation_output)
loadSEACResults(permutation_output, barcodes = "bin1")
names(seacGeneSets)
}
\references{Kummerfeld E, et al. (2025). Artifacts in spatial transcriptomics data: their detection, importance, prevalence, and prevention. \emph{Briefings in Bioinformatics}, 26(4), bbaf306. \doi{10.1093/bib/bbaf306}.}
