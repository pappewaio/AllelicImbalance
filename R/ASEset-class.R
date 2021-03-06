#'@include initialize-methods.R
NULL

#' ASEset objects
#' 
#' Object that holds allele counts, genomic positions and map-bias for a set of
#' SNPs
#' 
#' An ASEset object differs from a regular RangedSummarizedExperiment object
#' in that the assays contains an array instead of matrix. This array has
#' ranges on the rows, sampleNames on the columns and variants in the third
#' dimension.
#' 
#' It is possible to use the commands barplot and locationplot on an ASEset
#' object see more details in \code{\link{barplot}} and
#' \code{\link{locationplot}}.
#' 
#' Three different alleleCount options are available. The simples one is the
#' * option, and is for experiments where the strand information is not
#' known e.g. non-stranded data. The unknown strand could also be for strand 
#' specific data when the aligner could not find
#' any strand associated with the read, but this should normally not happen, 
#' and if it does probably having an extremely low mapping quality. 
#' Then there are an option too add plus
#' and minus stranded data. When using this, it is essential to make sure that
#' the RNA-seq experiment under analysis has in fact been created so that
#' correct strand information was obtained. The most functions will by default 
#' have their strand argument set to '*'.
#'
#' The phase information is stored by the convention of 
#' 'maternal chromosome|paternal chromosome', with 0 as reference allele and 1 
#' as alternative allele. '|' when the phase is known and '/' when the phase is
#' unknown. Internally the information will be stored as an three dimensional 
#' array, dim 1 for SNPs, dim 2 for Samples and dim 3 which is fixed and stores 
#' maternal chromosome, paternal chromosome and phased (1 equals TRUE).
#'
#' 
#' @name ASEset-class
#' @rdname ASEset-class
#' @aliases ASEset-class ASEset alleleCounts mapBias fraction arank 
#' frequency genotype genotype<- phase phase<- alleleCounts,ASEset-method mapBias,ASEset-method
#' fraction,ASEset-method arank,ASEset-method 
#' frequency,ASEset-method genotype,ASEset-method genotype<-,ASEset-method
#' alleleCounts<- alleleCounts<-,ASEset-method phase,ASEset-method phase<-,ASEset-method
#' 
#' @docType class
#' @param x ASEset object
#' @param strand which strand of '+', '-' or '*'
#' @param verbose makes function more talkative
#' @param return.type return 'names', rank or 'counts'
#' @param return.class return 'list' or 'array'
#' @param top.fraction.criteria 'maxcount', 'ref' or 'phase'
#' @param value replacement variable
#' @param ... additional arguments
#' @return An object of class ASEset containing location information and allele
#' counts for a number of SNPs measured in a number of samples on various
#' strand, as well as mapBias information. All data is stored in a manner
#' similar to the SummarizedExperiment class.
#' 
#' @section Table: table(...)
#' 
#' \describe{
#' Arguments: \item{...}{An \code{ASEset object} that contains the
#' variants of interest} 
#'
#' The generics for table does not easily allow more than one argument
#' so in respect to the different strand options, \code{table} will
#' return a SimpleList with length 3, one element for each strand.
#' }
#'
#' @section Frequency: frequency(x, 
#' return.class = "list", strand = "*",
#' threshold.count.sample = 15)
#' 
#' \describe{
#' Arguments: \item{x}{An \code{ASEset object} that contains the
#' variants of interest} 
#' 
#' \item{x}{threshold.count.samples} if sample has fewer counts the function
#' return NA.
#'
#' }
#' @section Constructor: ASEsetFromCountList(rowRanges, countListNonStranded =
#' NULL, countListPlus = NULL, countListMinus = NULL, countListUnknown = NULL,
#' colData = NULL, mapBiasExpMean = array(), verbose=FALSE, ...)
#' 
#' \describe{
#' 
#' Arguments: \item{rowRanges}{A \code{GenomicRanges object} that contains the
#' variants of interest} \item{countListNonStranded}{A \code{list} where each
#' entry is a matrix with allele counts as columns and sample counts as rows}
#' \item{countListPlus}{A \code{list} where each entry is a matrix with allele
#' counts as columns and sample counts as rows} \item{countListMinus}{A
#' \code{list} where each entry is a matrix with allele counts as columns and
#' sample counts as rows} \item{countListUnknown}{A \code{list} where each
#' entry is a matrix with allele counts as columns and sample counts as rows}
#' \item{colData}{A \code{DataFrame} object containing sample specific data}
#' \item{mapBiasExpMean}{A 3D \code{array} describing mapping bias. The SNPs
#' are in the 1st dimension, samples in the 2nd dimension and variants in the
#' 3rd dimension.} \item{verbose}{Makes function more talkative}
#' \item{...}{arguments passed on to SummarizedExperiment constructor} }
#'
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#' 
#' 
#' #make example countList
#' set.seed(42)
#' countListPlus <- list()
#' snps <- c('snp1','snp2','snp3','snp4','snp5')
#' for(snp in snps){
#'   count<-matrix(rep(0,16),ncol=4,dimnames=list(
#' c('sample1','sample2','sample3','sample4'),
#' c('A','T','G','C')))
#'   #insert random counts in two of the alleles 
#'   for(allele in sample(c('A','T','G','C'),2)){
#' count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
#'   }
#'   countListPlus[[snp]] <- count
#' }
#' 
#' #make example rowRanges
#' rowRanges <- GRanges(
#'   seqnames = Rle(c('chr1', 'chr2', 'chr1', 'chr3', 'chr1')),
#'   ranges = IRanges(1:5, width = 1, names = head(letters,5)),
#'   snp = paste('snp',1:5,sep='')
#' )
#'
#' #make example colData
#' colData <- DataFrame(Treatment=c('ChIP', 'Input','Input','ChIP'), 
#'  row.names=c('ind1','ind2','ind3','ind4'))
#' 
#' #make ASEset 
#' a <- ASEsetFromCountList(rowRanges, countListPlus=countListPlus, 
#' colData=colData)
#'
#'
#' #example phase matrix (simple form)
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#' 
#' phase(a) <- p
#'
#'
#' #generate ASEset from array
#' snps <- 999
#' samples <-5
#' ar <-array(rep(unlist(lapply(1:snps,
#' 			function(x){(sample(c(TRUE,FALSE,TRUE,FALSE), size = 4))})), samples), 
#' 			dim=c(4,snps,samples))
#' ar2 <- array(sample(50:300, 4*snps*samples,replace=TRUE), dim=c(4,snps,samples))
#' ar2[ar] <- 0
#' ar2 <- aperm(ar2, c(2, 3, 1))
#' dimnames(ar2) <- list(paste("snp",1:snps,sep=""),paste("sample",1:samples,sep=""),
#'							c("A","C","G","T"))
#' gr <- GRanges(seqnames=c("chr2"), ranges=IRanges(start=1:dim(ar2)[1], width=1), strand="*")
#' a <- ASEsetFromArrays(gr, countsUnknown=ar2)
#' 
#'
#' @exportClass ASEset
#' @exportMethod alleleCounts alleleCounts<- mapBias fraction arank
#' frequency genotype genotype<- phase phase<-
#' 
#' @export 

setClass("ASEset", contains = "RangedSummarizedExperiment", 
	representation(variants = "vector"))

#' @rdname ASEset-class
#' @export 
setGeneric("alleleCounts", function(x, strand = "*", return.class="list") {
    standardGeneric("alleleCounts")
})

#' @rdname ASEset-class
#' @export 
setMethod("alleleCounts", signature(x = "ASEset"), function(x, strand = "*",
	return.class="list") {

	#check that strand parameter is correct set
    if (!sum(strand %in% c("+", "-", "*", "both")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' or 'both' ")
    }

	#check if the assays are present
	if(strand=="+" ){
		if(!"countsPlus" %in% assayNames(x)){
			stop(paste("strand",strand,"is not present in ASEset object",sep=""))
		}
	}
	if(strand=="-" ){
		if(!"countsMinus" %in% assayNames(x)){
			stop(paste("strand",strand,"is not present in ASEset object",sep=""))
		}
	}
	if(strand=="*" ){
		if(!("countsMinus" %in% assayNames(x) | "countsPlus" %in% assayNames(x))){
			stop(paste("for strand '*' at least one of '+' or '-' is required to be",
					   " present in ASEset object. Old ASEsets set the '*' strand in",
					   " countsUnknown, but that has been deprecated. If so, move your counts",
					   " to the plus or minus slot instead ",
					   " assays(yourASEset, withDimnames=FALSE)[['countsPlus']] <- assays(yourASEset)[['countsUnknown']]",
					   sep=""))
		}
	}
	if(strand=="both" ){
		if(!"countsMinus" %in% assayNames(x) & "countsPlus" %in% assayNames(x)){
			stop(paste("for strand 'both' both '+' and '-' is required to be",
					   " present in the ASEset object",sep=""))
		}
	}

	#access the correct assays
    if (strand == "+") {
        ar <- assays(x)[["countsPlus"]]
    } else if (strand == "-") {
        ar <- assays(x)[["countsMinus"]]
    } else if (strand == "*") {

		if(!("countsMinus" %in% assayNames(x))){
			ar <- assays(x)[["countsPlus"]]
		}else if(!("countsPlus" %in% assayNames(x))){
			ar <- assays(x)[["countsMinus"]]
		}else{
			ar <- assays(x)[["countsMinus"]] + assays(x)[["countsPlus"]]
		}
    } else if (strand == "both") {
        ar <- array(c(assays(x)[["countsPlus"]], assays(x)[["countsMinus"]]), dim=c(dim(assays(x)[["countsPlus"]]), 2))
    } else {
        stop("not existing strand option")
    }
    
    
	if(return.class=="array"){
		dimnames(ar)[[3]]<- x@variants
		if (strand == "both") {
			dimnames(ar)[[4]]<- c("+","-")
		}
		ar

	}else if(return.class=="list"){


		if(strand=="both"){
			strands <- 2
			alleleCountList2 <- list()
		}else{
			strands <- 1
		}
		for( j in 1:strands){

			alleleCountList <- list()

			for (i in 1:nrow(ar)) {
				if(strand=="both"){
					mat <- ar[i, , ,j]
				}else{
					mat <- ar[i, , ]
				}
				if (class(mat)[1] == "integer") {
					mat <- t(as.matrix(mat))
				}
				if (class(mat)[1] == "numeric") {
					mat <- t(mat)
				}
				colnames(mat) <- x@variants
				rownames(mat) <- colnames(x)

				alleleCountList[[i]] <- mat
			}
			# add snp id
			names(alleleCountList) <- rownames(x)

			if(strand=="both"){
				alleleCountList2[[j]] <- alleleCountList
			}
		}

		if(strand=="both"){
			names(alleleCountList2) <- c("+","-")
			alleleCountList2
		}else{
			alleleCountList
		}
	}else{
		stop("return.class has to be 'list' or 'array'")
	}

    
})
### -------------------------------------------------------------------------
### helpers for alleleCounts
###
.plusStrandCountsExists <- function(x){
	"countsPlus" %in% assayNames(x)
}
.minusStrandCountsExists <- function(x){
	"countsMinus" %in% assayNames(x)
}
.unknownStrandCountsExists <- function(x){
	"countsUnknown" %in% assayNames(x)
}



#' @rdname ASEset-class
#' @export 
setGeneric("alleleCounts<-", function(x, ..., value) {
    standardGeneric("alleleCounts<-")
})

#' @rdname ASEset-class
#' @export 
setReplaceMethod("alleleCounts", "ASEset", 
	function(x, strand="*", return.class="array", ..., value) {

	#for dev
	#message("valueX: ", paste0(value, collapse=", "))

    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else {
        stop("not existing strand option")
    }
    


	#check that value has the right dimensions (make validObject method)
	#if(!all(dim(value)==dim(assays(x, withDimnames=FALSE)[["countsUnknown"]]))){
	#	stop("dimensions in replacement object is not correct")
	#}

	assays(x, withDimnames=FALSE)[[el]][] <- value

	#return object
	x
})

#' @rdname ASEset-class
#' @export 
setGeneric("mapBias", function(x, ...) {
    standardGeneric("mapBias")
})

#' @rdname ASEset-class
#' @export 
setMethod("mapBias", signature(x = "ASEset"), function(x,
	return.class="list") {
    # assume alleleCount information is stored as element 1

	if(return.class=="array"){
		return(assays(x)[["mapBias"]])

	}else if(return.class=="list"){
		mapBiasList <- list()
		for (i in 1:nrow(x)) {
			mat <- assays(x)[["mapBias"]][i, , ]
			if (class(mat)[1] == "numeric") {
				dim(mat) <- c(1, 4)
				rownames(mat) <- colnames(x)
			}
			colnames(mat) <- x@variants

			mapBiasList[[i]] <- mat
		}
		# add snp id
		names(mapBiasList) <- rownames(x)

		return(mapBiasList)
	}
    
})

#' @rdname ASEset-class
#' @export 
setGeneric("fraction", function(x, ...) {
    standardGeneric("fraction")
})

#' @rdname ASEset-class
#' @export 
setMethod("fraction", signature(x = "ASEset"), function(x, strand = "*", 
    top.fraction.criteria="maxcount", verbose = FALSE, ...) {
    
    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }

	#calculate frequency
	fr <- frequency(x, strand=strand, return.class="array", ...)

	#check and use top.fraction.criteria=="phase"
	if(top.fraction.criteria=="phase"){
		if(!"phase" %in% assayNames(x)){
			stop("the phase slot has not been initialized")
		}
		if(is.null(assays(x)[["phase"]])){
			stop("the phase slot cannot be empty if 'top.fraction.criteria=\"phase\"'")
		}

		#take out reference allele information
		ref <- ref(x)
		#select only ref rows (here is mcols() ref required)
		ar <- .arrayFromAlleleVector(x@variants, ref,nc=ncol(x))
		#use the array to subset the ref frequencies
		ret <- .subsetArrayToMatrix(fr, ar)
		#use maternal phase freq, by reverse values in mat that are not ref (0)
		ret <- .returnMaternalPhaseFrequency(phase(x,return.class="array")[,,1], ret)
	}else if(top.fraction.criteria=="ref"){
		#take out reference allele information
		ref <- ref(x)
		#select only ref rows (here is mcols() ref required)
		ar <- .arrayFromAlleleVector(x@variants, ref,nc=ncol(x))
		#use the array to subset the ref frequencies
		ret <- .subsetArrayToMatrix(fr, ar)
	}else if(top.fraction.criteria=="maxcount"){
		#use output from rank as maxcount, select only 1st rank 
		ar <- .arrayFromAlleleVector(x@variants, 
				arank(x, strand = strand, return.class="matrix")[,1], ncol(x))
		#use the array to subset the ref frequencies
		ret <- .subsetArrayToMatrix(fr, ar)
	}
    # return object matrix
	colnames(ret) <- colnames(x)
	rownames(ret) <- rownames(x)
    ret
})
### -------------------------------------------------------------------------
### helpers for fraction
###

#ph is an array, rf is the fraction matrix for ref frequencies
.returnMaternalPhaseFrequency <- function(ph,rf){
	rf[ph==1] <- 1 - rf[ph==1]
	rf
}

#' @rdname ASEset-class
#' @export 
setGeneric("arank", function(x, return.type = "names", 
	return.class = "list", strand = "*", ...) {
    standardGeneric("arank")
})

setMethod("arank", signature(x = "ASEset"), function(x, return.type = "names", 
	return.class = "list", strand = "*", ...) {

	if(return.class=="matrix"){
		if (return.type == "names") {
			ar <- t(apply(apply(alleleCounts(x,
						strand=strand,return.class="array"),c(1,3),sum),
					1, function(x){rank(x,ties.method="first")}))
			ar2 <- t(apply(ar,1,function(x){
						   x <- sort(x,index.return=TRUE,decreasing=TRUE)$ix
						   x
					}))

			mat <- matrix(x@variants[ar2],ncol=4, nrow=nrow(x), byrow=FALSE,
					  dimnames=list(dimnames(ar)[[1]],c(1,2,3,4)))

			return(mat)
		}else if (return.type == "rank") {
			return(t(apply(apply(alleleCounts(x,
						strand=strand,return.class="array"),c(1,3),sum),
					1, function(x){rank(x,ties.method="first")})))
		}else if (return.type == "counts") {
			return(apply(alleleCounts(x,strand=strand,return.class="array"),c(1,3),sum))
		}else{stop("return.type is not valid")}

	}else if(return.class=="list"){

		acounts <- alleleCounts(x, strand = strand)

		if (return.type == "names") {
			lapply(acounts, function(x) {
				names(sort(apply(x, 2, sum), decreasing = TRUE))
			})
		}else if (return.type == "counts") {
			lapply(acounts, function(x) {
				as.numeric(sort(apply(x, 2, sum), decreasing = TRUE))
			})
		}else if (return.type == "rank") {
			stop("rank for return.class 'list' is not implemented yet")
		}else{stop("return.type is not valid")}
	}
}) 

#setGeneric("table")
#setGeneric("table", function(x, strand = "*", sortBy="none", ...) {
#    standardGeneric("table")
#})

# @rdname ASEset-class
#setMethod("table", signature(... = "ASEset"), function(...) {
#
#	args <- list(...)
#	if (length(args) > 1)
#	  stop("Only one argument in '...' supported")
#	x <- args[[1L]]
#
#	#because the generis of table is rubbish we have to return a list for each strand
#	#retList <- list()
#
#	#for(strand in c("+","-","*")){
#	#	df <- data.frame(row.names=rownames(x))
#	#	df[,c("chromosome","position")] <- c(as.character(seqnames(x)),start(x))
#	#	df <- cbind(df,as.data.frame(arank(x, return.type="counts",
#	#				   return.class="matrix",strand=strand)))
#
#	#	#if only one sample add fraction to table()
#	#	if(ncol(x)==1){
#	#		df[,"fraction"] <- as.vector(fraction(x,strand=strand))
#	#	}
#
#	#	retList[[strand]] <- df
#
#	}
#	#return(SimpleList(retList))
#
#	df()
#
#})

#' @rdname ASEset-class
#' @export 
setGeneric("frequency")

setMethod("frequency", signature(x = "ASEset"), function(x, 
	return.class = "list", strand = "*",
	threshold.count.sample = 1) {
	#threshold.count.sample cannot be zero (further down f=1/0 woulf fail)
	if(threshold.count.sample<1){stop("threshold.count.sample needs to be >1")}
	#calculate frequency from Allele counts
	ar2 <- .calcFrequencyFromAlleleCounts(alleleCounts(x, strand=strand, return.class="array"), threshold.count.sample)
	if(return.class=="array"){
		return(ar2)
	}else if(return.class=="list"){
		return(.array2MatrixList(ar2))
	}else{
		stop("return.class has to be 'array' or 'list'")
	}
})
### -------------------------------------------------------------------------
### helpers for frequemcy
###
.calcFrequencyFromAlleleCounts <- function(ar, th){
	actot <- apply(ar, c(1,2), sum, na.rm=TRUE)
	actot[actot  < th] <- NaN
	ar * array(as.vector(1/actot),dim=dim(ar))
}



#' @rdname ASEset-class
#' @export 
setGeneric("genotype", function(x, ...){
    standardGeneric("genotype")
})

#' @rdname ASEset-class
#' @export 
setMethod("genotype", signature(x = "ASEset"), function(x,
			return.class="matrix"){

    if (!("phase" %in% assayNames(x))) {
		stop(paste("phase is not present in assays in",
				   " ASEset object, see '?inferGenotypes' "))
    }
	if(!("ref" %in% names(mcols(x)))){
		stop(paste("ref allele is not present in mcols in",
				   " ASEset object, see '?ASEset' "))
	}
	if(!("alt" %in% names(mcols(x)))){
		stop(paste("alt allele is not present in mcols in",
				   " ASEset object, see '?inferAltAllele' "))
	}

	phase2genotype(phase(x, return.class="array"), ref(x), alt(x), return.class=return.class)

})

#' @rdname ASEset-class
#' @export 
setGeneric("genotype<-", function(x, value){
    standardGeneric("genotype<-")
})

#' @rdname ASEset-class
#' @export 
setMethod("genotype<-", signature(x = "ASEset"), function(x, value){

	#check dimensions
	if(!nrow(x)==nrow(value)){
		stop("nrow(x) is not equal to nrow(value)")
	}
	if(!ncol(x)==ncol(value)){
		stop("ncol(x) is not equal to ncol(value)")
	}
	#check prexence of ref
	if(!"ref" %in% names(mcols(x))){
		stop("ref() is not set in ASEset")
	}
        #It would be better to not allow withDimnames=FALSE
        #I suspect that we need to change something in the init object to fix that.
	assays(x, withDimnames=FALSE)[["phase"]] <- genotype2phase(value, ref(x))	
	x
})


#' @rdname ASEset-class
#' @export 
setGeneric("countsPerSnp", function(x, ...){
    standardGeneric("countsPerSnp")
})

#' @rdname ASEset-class
#' @export 
setMethod("countsPerSnp", signature(x = "ASEset"), function(x, 
	return.class = "matrix", return.type="mean", strand = "*") {

	if(return.class=="matrix"){
		return(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum))
	}else if(return.class=="vector"){
		if(return.type=="all"){
			return(apply(alleleCounts(x, strand=strand, return.class="array"), 1, sum))
		}else if(return.type=="mean"){
			return(apply(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum), 1, mean))
		}
	}else{
		stop("return.class has to be 'vector' or 'matrix'")
	}
})

#' @rdname ASEset-class
#' @export 
setGeneric("countsPerSample", function(x, ...){
    standardGeneric("countsPerSample")
})

#' @rdname ASEset-class
#' @export 
setMethod("countsPerSample", signature(x = "ASEset"), function(x, 
	return.class = "matrix", return.type="mean", strand = "*") {

	if(return.class=="matrix"){
		return(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum))
	}else if(return.class=="vector"){
		if(return.type=="all"){
			return(apply(alleleCounts(x, strand=strand, return.class="array"), 2, sum))
		}else if(return.type=="mean"){
			return(apply(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum), 2, mean))
		}
	}else{
		stop("return.class has to be 'vector' or 'matrix'")
	}
})

#' @rdname ASEset-class
#' @export 
setGeneric("phase", function(x, ...){
    standardGeneric("phase")
})

#' @rdname ASEset-class
#' @export 
setMethod("phase", signature(x = "ASEset"), function(x, 
	return.class = "matrix" ) {

	if(return.class=="matrix"){
		mat <- phaseArray2phaseMatrix(assays(x)[["phase"]])
		colnames(mat) <- colnames(x)
		rownames(mat) <- rownames(x)
		mat
	}else if(return.class=="array"){
		assays(x)[["phase"]] 
	}
})

#' @rdname ASEset-class
#' @export 
setGeneric("phase<-", function(x, value){
    standardGeneric("phase<-")
})

#' @rdname ASEset-class
#' @export 
setMethod("phase<-", signature(x = "ASEset"), function(x,value) {

	if(class(value)[1]=="matrix") {

		if(!identical(dim(x),dim(value))){
			stop("dimension of value does not correspond to the values of object ASEset")
		}
	
		assays(x, withDimnames=FALSE)[["phase"]] <- phaseMatrix2Array(value,dimnames=dimnames(x))

	}else if(class(value)[1]=="array"){
		assays(x, withDimnames=FALSE)[["phase"]] <- value
	}

	x
})

#' @rdname ASEset-class
#' @export 
setGeneric("mapBias<-", function(x, value){
    standardGeneric("mapBias<-")
})

#' @rdname ASEset-class
#' @export 
setMethod("mapBias<-", signature(x = "ASEset"), function(x,value) {

	if(class(value)[1]=="array") {
		assays(x, withDimnames=FALSE)[["mapBias"]] <- value
	}else {
		stop("class has to be array")
	}
	x
})

#' @rdname ASEset-class
#' @export 
setGeneric("refExist", function(x){
    standardGeneric("refExist")
})

#' @rdname ASEset-class
#' @export 
setMethod("refExist", signature(x = "ASEset"), function(x) {
		!is.null(ref(x))
})

#' @rdname ASEset-class
#' @export 
setGeneric("ref")

#' @rdname ASEset-class
#' @export 
setMethod("ref", signature(x = "ASEset"), function(x) {

		mcols(x)[["ref"]]

})

#' @rdname ASEset-class
#' @export 
setGeneric("ref<-")

#' @rdname ASEset-class
#' @export 
setMethod("ref<-", signature(x = "ASEset"), function(x, value) {

	if(class(value)[1]=="character") {

		mcols(x)[["ref"]] <- value
	}else{

		stop("wrong class")
	}

	x
})

#' @rdname ASEset-class
#' @export 
setGeneric("altExist", function(x){
    standardGeneric("altExist")
})

#' @rdname ASEset-class
#' @export 
setMethod("altExist", signature(x = "ASEset"), function(x) {
		!is.null(ref(x))
})

#' @rdname ASEset-class
#' @export 
setGeneric("alt")

#' @rdname ASEset-class
#' @export 
setMethod("alt", signature(x = "ASEset"), function(x) {

		mcols(x)[["alt"]]

})

#' @rdname ASEset-class
#' @export 
setGeneric("alt<-")

#' @rdname ASEset-class
#' @export 
setMethod("alt<-", signature(x = "ASEset"), function(x, value) {

	if(class(value)[1]=="character") {

		mcols(x)[["alt"]] <- value
	}else{

		stop("wrong class")
	}

	x
})

#' @rdname ASEset-class
#' @export 
setGeneric("aquals", function(x, ...){
    standardGeneric("aquals")
})

#' @rdname ASEset-class
#' @export 
setMethod("aquals", signature(x = "ASEset"), function(x) {
		assays(x)[["aquals"]]
})

#' @rdname ASEset-class
#' @export 
setGeneric("aquals<-", function(x, value){
    standardGeneric("aquals<-")
})

#' @rdname ASEset-class
#' @export 
setMethod("aquals<-", signature(x = "ASEset"), function(x,value) {

	if(class(value)[1]=="array") {

		if(!identical(dim(x),dim(value)[1:2])){
			stop("dimension of value does not correspond to the values of object ASEset")
		}
		assays(x, withDimnames=FALSE)[["aquals"]] <- value

	}
	x
})


#' @rdname ASEset-class
#' @export 
setGeneric("maternalAllele", function(x, ...){
    standardGeneric("maternalAllele")
})

#' @rdname ASEset-class
#' @export 
setMethod("maternalAllele", signature(x = "ASEset"), 
		function(x) {

		mat <- phase(x,return.class="array")[,,1]
		ref <- mcols(x)[["ref"]]
		alt <- mcols(x)[["alt"]]

		apply(t(mat),1,function(y, ref, alt){

			vec <- rep(NA,length(y))
			if(any(y == 1)){
				vec[y == 1] <- alt[y == 1]
			}
			if(any(y == 0)){
				vec[y == 0] <- ref[y == 0]
			}
			vec
		}, ref=ref, alt=alt)
})

#' @rdname ASEset-class
#' @export 
setGeneric("paternalAllele", function(x, ...){
    standardGeneric("paternalAllele")
})

#' @rdname ASEset-class
#' @export 
setMethod("paternalAllele", signature(x = "ASEset"), 
		function(x) {

		mat <- phase(x,return.class="array")[,,2]
		ref <- mcols(x)[["ref"]]
		alt <- mcols(x)[["alt"]]

		apply(t(mat),1,function(y, ref, alt){

			vec <- rep(NA,length(y))
			if(any(y == 1)){
				vec[y == 1] <- alt[y == 1]
			}
			if(any(y == 0)){
				vec[y == 0] <- ref[y == 0]
			}
			vec
		}, ref=ref, alt=alt)
})


##' @rdname ASEset-class
##' @export 
#setGeneric("testet<-", signature="x", function(x, ..., value){
#    standardGeneric("testet<-")
#})
##we need to set a normal generic to be able to use the replacement generic
#setGeneric("testet", signature="x", function(x, ...){
#    standardGeneric("testet")
#})
#
##' @rdname ASEset-class
##' @export 
#setReplaceMethod("testet", c("ASEset"), 
#				 function(x, ... ,value){
#					message("value0: ", paste0(value, collapse=", "))
#					.set_ASEset_ref(x, value)
#				 }
#)
#
#setMethod("testet","ASEset",
#	function(x) ref(x)
#)
#
#.set_ASEset_ref <- function(x, value)
#{
#	message("value1: ", paste0(value, collapse=", "))
#    ## Use 'x@width[]' instead of 'x@width' so the right value is recycled.
#	ref(x)[] <- value
#	x
#}
#
#ref(x)
#testet(x) <- c("O","P","Q")
#testet(x)
#testet(x)[1:2] <- c("S","I")
#
#
##another test
#
#setGeneric("getLetters", function(x)standardGeneric("getLetters"))
#setGeneric("getLetters<-", function(x, value)standardGeneric("getLetters<-"))
#
#setMethod("getLetters", "ASEset", function(x) {
#		    x@variants
#})
#
#setReplaceMethod("getLetters", c("ASEset", "character"), function(x, value) {
#				   message("value: ", paste0(value, collapse=", "))
#				     x@variants <- value
#				     x
#})
#
## [1] "<0x3716b40>"
#getLetters(x)
##  [1] "A" "B" "C" "D" "E" "F" "G" "H" "I" "J"
#getLetters(x)[1:2] <- letters[1:2]
#
#
#
