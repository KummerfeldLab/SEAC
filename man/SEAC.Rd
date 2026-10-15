\name{SEAC}
\alias{SEAC}
\alias{bladeMaskHD}
\alias{runBladeQC}
\alias{maxMinScore}
\alias{maxMinPermutationTest}
\alias{runSEAC}
\alias{loadSEACResults}
\alias{seacGeneSets}
\title{KD-tree quality control and MaxMin testing for Visium HD}
\description{SEAC imports Visium HD Space Ranger outputs, applies a KD-tree radius-one BLADE mask, and tests general gene-set spatial co-localization with a UMI-level permutation null model.}
\usage{
bladeMaskHD(positions)
runBladeQC(input_dir, output_dir = NULL, bin_size = "002", sample_id = NULL)
maxMinScore(spe, gene_sets)
maxMinPermutationTest(x, gene_sets, mode = c("within", "pooled"), n_permutations = 5000L, legacy_alpha = 0.01, fdr_cutoff = 0.05, seed = 1L, chunk_size = 100L, output_dir = NULL, resume = TRUE, BPPARAM = BiocParallel::SerialParam())
runSEAC(input_dir, output_dir = NULL, bin_size = "002", permutation = FALSE, gene_sets = NULL, ...)
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
names or 10x gene symbols from \code{rowData(x)$Symbol}. For
\code{runSEAC}, these are optional unless \code{permutation = TRUE}; when
supplied they are also used to score the returned post-QC object.}
\item{x}{A post-QC object, directory, or named list of either.}
\item{mode}{Within-sample or pooled-sample label permutations.}
\item{n_permutations}{Number of Monte Carlo permutations.}
\item{legacy_alpha}{Original lower-tail hotspot threshold.}
\item{fdr_cutoff}{BH adjusted p-value cutoff.}
\item{seed}{Integer random seed.}
\item{chunk_size}{Permutations per checkpoint.}
\item{resume}{Whether a matching interrupted run may resume.}
\item{BPPARAM}{A BiocParallel backend.}
\item{permutation}{Whether to run the optional test.}
\item{barcodes}{Optional result barcode subset.}
\item{...}{Arguments passed to \code{maxMinPermutationTest}.}
}
\value{\code{bladeMaskHD} returns a mask data frame. \code{runBladeQC} and \code{maxMinScore} return a \code{SpatialExperiment}. \code{runSEAC} also persists supplied gene-set scores in \code{output_dir/post_qc/spe}. The permutation function returns a scored object or named list, while \code{loadSEACResults} returns a data frame.}
\details{The BLADE mask constructs a KD-tree from capture-border and non-tissue seeds, then marks every bin whose nearest seed is at Euclidean array-coordinate distance at most one. MaxMin is the maximum across input gene sets of the ceil-rounded nearest positive-bin distance. Lower values are more spatially co-localized. Gene sets resolve against both feature row names and 10x \code{Symbol} aliases; all rows sharing a requested symbol contribute.}
\examples{
positions <- data.frame(
    barcode = paste0("b", 1:4), in_tissue = c(0L, 1L, 1L, 1L),
    array_row = c(1L, 1L, 2L, 2L), array_col = c(1L, 2L, 1L, 2L)
)
bladeMaskHD(positions)
}
\references{Kummerfeld E, et al. (2025). Artifacts in spatial transcriptomics data: their detection, importance, prevalence, and prevention. \emph{Briefings in Bioinformatics}, 26(4), bbaf306. \doi{10.1093/bib/bbaf306}.}
