#' @include initialize-methods.R
NULL

#' RegionSummary class
#' 
#' Object that holds results from the regionSummary method
#'
#' The RegionSummary-class objects contains summaries for specified regions
#'
#' @name RegionSummary-class
#' @rdname RegionSummary-class
#' @aliases RegionSummary-class RegionSummary RegionSummary-method
#' @docType class 
#' @param x RegionSummary object
#' @param ... pass arguments to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class RegionSummary
#' @examples
#'
#' #some code
#'
#' @exportClass RegionSummary
NULL

#' @rdname RegionSummary-class
#' @exportClass RegionSummary
setClass("RegionSummary", contains = "RangedSummarizedExperiment",
	representation(
		meta = "list",
		sumnames = "character"
	)
)

#' @rdname RegionSummary-class
#' @export 
setGeneric("sumnames", function(x, ...){
    standardGeneric("sumnames")
})

#' @rdname RegionSummary-class
#' @export 
setMethod("sumnames", signature(x = "RegionSummary"), function(x) {
	x@sumnames
})

#' @rdname RegionSummary-class
#' @export 
setGeneric("basic", function(x, ...){
    standardGeneric("basic")
})

#' @rdname RegionSummary-class
#' @export 
setMethod("basic", signature(x = "RegionSummary"), function(x) {
	ar <- assays(x)[["rs1"]]
	dimnames(ar)[[3]] <- x@sumnames
	lst <- .array2MatrixList(ar)
	names(lst) <- rownames(x)
	lst
})



