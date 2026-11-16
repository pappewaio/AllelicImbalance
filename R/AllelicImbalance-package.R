#' A package meant to provide all basic functions for high-throughput allele
#' specific expression analysis
#' 
#' Package \code{AllelicImbalance} has functions for importing, filtering and
#' plotting high-throughput data to make an allele specific expression
#' analysis. A major aim of this package is to provide functions to collect as
#' much information as possible from regions of choice, and to be able to
#' explore the allelic expression of that region in detail.
#' 
#' \tabular{ll}{ Package: \tab AllelicImbalance\cr Type: \tab Package\cr
#' Version: \tab 1.2.0\cr Date: \tab 2014-08-24\cr License: \tab GPL-3\cr }
#' 
#' @name AllelicImbalance-package
#' @aliases AllelicImbalance-package AllelicImbalance
#' @docType package
#' @section Overview - standard procedure: Start out creating a \code{GRange}
#' object defining the region of interest. This can also be done using
#' \code{getAreaFromGeneNames} providing gene names as arguments. Then use
#' \code{BamImpGAList} to import reads from that reagion and find potential
#' SNPs using \code{scanForHeterozygotes}. Then retrieve the allele counts of
#' heterozygote sites by the function \code{getAlleleCount}. With this data
#' create an \code{ASEset}. At this point all pre-requisites for a 'basic'
#' allele specific expression analysis is available. Two ways to go on could be
#' to apply \code{\link{chisq.test}} or \code{\link{barplot}} on this ASEset
#' object.
#' @author Author: Jesper Robert Gadin Author: Lasse Folkersen
#' 
#' Maintainer: Jesper Robert Gadin <j.r.gadin@@gmail.com>
#' @seealso \itemize{ \item code?ASEset }
#' @references Reference to published application note (work in progress)
#' @keywords package
#'
#' @import methods
#' @import BiocGenerics
#' @import S4Vectors
#' @import grid
#' @import IRanges
#' @import GenomeInfoDb
#' @import GenomicRanges
#' @import SummarizedExperiment
#' 
#' @importFrom Biostrings subseq
#' @importFrom Biostrings DNAStringSet
#' 
#' @importFrom Rsamtools scanBamWhat
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools indexBam
#' @importFrom Rsamtools BamFileList
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBcfParam
#' @importFrom Rsamtools scanBcf
#' @importFrom Rsamtools BcfFile
#' @importFrom Rsamtools BcfFileList
#' @importFrom Rsamtools indexBcf
#' @importFrom Rsamtools indexFa
#' @importFrom Rsamtools FaFile
#' @importFrom Rsamtools scanFa
#' @importFrom Rsamtools scanFaIndex
#' @importFrom Rsamtools applyPileups
#' @importFrom Rsamtools ApplyPileupsParam
#' @importFrom Rsamtools PileupFiles
#' 
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importClassesFrom GenomicAlignments GAlignmentsList
#' 
#' @importFrom GenomicAlignments pileLettersAt
#' @importFrom GenomicAlignments cigar
#' @importFrom cigarillo cigars_as_RleList
#' @importFrom GenomicAlignments GAlignments
#' @importFrom GenomicAlignments GAlignmentPairs
#' @importFrom GenomicAlignments GAlignmentsList
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' 
#' @importFrom VariantAnnotation alt
#' @importFrom VariantAnnotation ref
#' @importFrom VariantAnnotation 'alt<-'
#' @importFrom VariantAnnotation 'ref<-'
#'
#' @importFrom lattice barchart
#' @importFrom lattice levelplot
#' @importFrom lattice xyplot
#' @importFrom lattice trellis.par.set
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.text
#' @importFrom lattice panel.barchart
#' @importFrom lattice panel.smoothScatter
#' @importFrom lattice panel.linejoin
#'
#' @importFrom seqinr write.fasta
#'
#' @importClassesFrom Gviz DataTrack
#' @importClassesFrom Gviz AnnotationTrack
#' @importClassesFrom Gviz DetailsAnnotationTrack
#' @importClassesFrom Gviz GeneRegionTrack
#' @importClassesFrom Gviz GenomeAxisTrack
#'
#' @importFrom Gviz DataTrack
#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz DetailsAnnotationTrack
#' @importFrom Gviz GeneRegionTrack
#' @importFrom Gviz GenomeAxisTrack
#' @importFrom Gviz plotTracks
#'
#' @importFrom gridExtra grid.arrange
#'
#' @importFrom latticeExtra doubleYScale
#'
#' @importClassesFrom GenomicFeatures TxDb
#' @importClassesFrom AnnotationDbi AnnotationDb
#' @importClassesFrom AnnotationDbi OrgDb
#' @importClassesFrom AnnotationDbi ChipDb
#' @importClassesFrom AnnotationDbi GODb
#' @importClassesFrom AnnotationDbi ReactomeDb
#'
#' @importMethodsFrom AnnotationDbi columns
#'
#' @importMethodsFrom BSgenome snpsByOverlaps
#' 
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keys
#' 
#' @importFrom GenomicFeatures exons
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicFeatures cds
#' @importFrom GenomicFeatures isActiveSeq
#' @importFrom GenomicFeatures 'isActiveSeq<-'
#'
#' @importFrom grDevices heat.colors
#'
#' @importFrom graphics abline 
#' @importFrom graphics lines
#' @importFrom graphics mtext
#' @importFrom graphics plot.default
#' @importFrom graphics symbols
#' @importFrom graphics text
#' @importFrom graphics title
#' @importFrom graphics points
#'
#' @importFrom stats anova
#' @importFrom stats dist
#' @importFrom stats lm
#' @importFrom stats setNames
#' @importFrom stats formula
#'
#' @importFrom nlme lme
#'
NULL

#' GRvariants object
#' 
#' this data is a \code{GRanges} object that contains the ranges for three
#' example SNPs.
#' 
#' 
#' @name GRvariants
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{reads}} which is another example
#' object }
#' @keywords data object example
#' @examples
#' 
#' #load example data
#' data(GRvariants)
#' 
#' 
NULL


#' reads object
#' 
#' This data set corresponds to the BAM-file data import illustrated in the
#' vignette. The data set consists of a chromosome 17 region from 20 RNA-seq
#' experiments of HapMap samples.
#' 
#' 
#' @name reads
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{GRvariants}} which is another
#' example object }
#' @references Montgomery SB et al. Transcriptome genetics using second
#' generation sequencing in a Caucasian population. Nature. 2010 Apr
#' 1;464(7289):773-7.
#' @keywords data object example
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(reads)
#' 
NULL 


#' ASEset.sim object
#' 
#' ASEset with simulated data with SNPs within the first 200bp of chromosome 17,
#' which is required to have example data for the refAllele function.
#' 
#' @name ASEset.sim
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords data object example
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(ASEset.sim)
#' 
NULL 

#' genomatrix object
#' 
#' genomatrix is an example of a matrix with genotypes
#' 
#' @name genomatrix
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords data genotype
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(genomatrix)
#' 
NULL 


#' ASEset.old object
#' 
#' old version of an ASEset which needs to be updated
#' 
#' @name ASEset.old
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords data object example
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(ASEset.old)
#' 
NULL 
