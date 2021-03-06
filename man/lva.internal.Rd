% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linkage-methods.R
\docType{methods}
\name{lva.internal}
\alias{lva.internal}
\alias{lva.internal,array-method}
\title{lva.internal}
\usage{
lva.internal(x, ...)

\S4method{lva.internal}{array}(
  x,
  grp,
  element = 3,
  type = "lm",
  subject = NULL,
  covariates = matrix(),
  ...
)
}
\arguments{
\item{x}{regionSummary array phased for maternal allele}

\item{...}{arguments to forward to internal functions}

\item{grp}{group 1-3 (1 for 0:0, 2 for 1:0 or 0:1, and 3 for 1:1)}

\item{element}{which column in x contains the values to use with lm.}

\item{type}{which column in x contains the values to use with lm.}

\item{subject}{which samples belongs to the same individual}

\item{covariates}{add data.frame with covariates (only integers and numeric)}
}
\description{
make an almlof regression for arrays (internal function)
}
\details{
internal method that takes one array with results from regionSummary
and one matrix with group information for each risk SNP (based on phase).
Input and output objects can change format slightly in future.
}
\examples{

data(ASEset) 
a <- ASEset
# Add phase
set.seed(1)
p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
	nrow=nrow(a), ncol(a))

phase(a) <- p

#add alternative allele information
mcols(a)[["alt"]] <- inferAltAllele(a)

# in this example two overlapping subsets of snps in the ASEset defines the region
region <- split(granges(a)[c(1,2,2,3)], c(1,1,2,2))
rs <- regionSummary(a, region, return.class="array", return.meta=FALSE)

# use  (change to generated riskSNP phase later)
phs <- array(c(phase(a,return.class="array")[1,,c(1, 2)], 
			 phase(a,return.class="array")[2,,c(1, 2)]), dim=c(20,2,2))
grp <- matrix(2, nrow=dim(phs)[1], ncol=dim(phs)[2])		 
grp[(phs[,,1] == 0) & (phs[,,2] == 0)] <- 1
grp[(phs[,,1] == 1) & (phs[,,2] == 1)] <- 3
#only use mean.fr at the moment, which is col 3
lva.internal(x=assays(rs)[["rs1"]],grp=grp, element=3)

}
\author{
Jesper R. Gadin, Lasse Folkersen
}
\keyword{phase}
