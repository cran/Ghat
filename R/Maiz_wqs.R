#' The Wisconsin Quality Synthetic (WQS) maize population datasets.
#'
#' @source Lorenz, A. J., Beissinger, T.M., Rodrigues, R., de Leon, N. 2015. 
#'	Selection for silage yield and composition did not affect genomic diversity 
#'	within the Wisconsin Quality Synthetic maize population. 
#'	Genes Genomes Genetics. DOI: 10.1534/g3.114.015263.
#'
#' @format A list of 4 data sets:
#' \describe{
#'  \item{Dataset-1}{Data frame including all SNP effects and changes in allele frequencies between cycle 2 and 5.}
#'  \item{Dataset-2}{BLUP breeding values for Acid Detergent Fiber (ADF) involving 5 generations of selection.}
#'  \item{Dataset-3}{Map file: Each line of the Map file describes a single marker and must contain at least three columns. 1: chromosome number; 2: SNP (snp id); 3: SNP position (in base-pairs (bp)).}
#'  \item{Dataset-4}{Maize Genotype (Illumina MaizeSNP50 BeadChip); an Infinium HD assay (Illumina, Inc. San Diego, CA). 
#'                   10,017 SNP markers (0,1 and 2) after filtration, distributed across the maize genome (Ganal et al.2011).}
#' }
#' @examples
#' maize <- Maize_wqs[[1]]
#' phe   <- Maize_wqs[[2]]
#' map   <- Maize_wqs[[3]]
#' gen   <- Maize_wqs[[4]]
"Maize_wqs"