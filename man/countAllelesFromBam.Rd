% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxillary-methods.R
\docType{methods}
\name{countAllelesFromBam}
\alias{countAllelesFromBam}
\alias{countAllelesFromBam,GRanges-method}
\title{alleleCounts from bam file}
\usage{
countAllelesFromBam(gr, ...)

\S4method{countAllelesFromBam}{GRanges}(
  gr,
  pathToDir,
  flag = NULL,
  scanBamFlag = NULL,
  return.class = "array",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{gr}{GRanges that contains SNPs of interest}

\item{...}{arguments to pass on}

\item{pathToDir}{path to directory of bam files}

\item{flag}{specify one flag to use as filter, default is no filtering.
allowed flags are 99, 147, 83 and 163}

\item{scanBamFlag}{set a custom flag to use as filter}

\item{return.class}{type of class for the returned object}

\item{verbose}{makes funciton more talkative}
}
\description{
count alleles before creating ASEse.
}
\details{
counts the alleles in a bam file based on GRanges positions.

Important excerpt from the details section of the internal applyPileups
 function: Regardless of 'param' values, the algorithm follows samtools by
 excluding reads flagged as unmapped, secondary, duplicate, or
 failing quality control.
}
\examples{

data(GRvariants)
gr <- GRvariants

##not run at the moment
#pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#ar <- countAllelesFromBam(gr, pathToDir)

}
\author{
Jesper R. Gadin
}
\keyword{allelecount}
\keyword{counting}
