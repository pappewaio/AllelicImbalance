#'@include AllelicImbalance-package.R
NULL

### =========================================================================
### General helper functions not exported
### =========================================================================

#supposed to merge paths and file irrespective OS and presence of trailing slash
.mergeDirAndFilename <- function(dir, files){
	#check for presence of / in filename in that case remove
	files <- sub("/","",files)
	paste(normalizePath(dir),"/",files, sep="")
}

.matrixFromLmListCommonParam <- function(lst){
	mat <- t(sapply(lst, function(x){
		s <- summary(x)$coefficients	   
		c(s[1,1],s[2,1],s[1,2],s[2,2],s[1,3],s[2,3],s[1,4],s[2,4])
	}))
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	mat
}

#transform an index vector into a IRanges object
.IRangesFromIntegerList <- function(idx){
	IRanges(c(1,idx@partitioning@end[-length(idx@partitioning@end)]+1), idx@partitioning@end)
}

#first dimension of array will make up the length of the list
.array2MatrixList <- function(ar){
	lapply(setNames(1:dim(ar)[1], dimnames(ar)[[1]]), function(i, ar){
		ar[i,,]
	}, ar=ar)
}

#important helper to pick put the frequence of the reference allele
.arrayFromAlleleVector <- function(var, sel, nc){
	selmat <- array(sel, dim=c(length(sel), length(var)))
	varmat <- array(var, dim=c(length(var), length(sel)))
	tf <- t(selmat) == varmat

	#expand to array
	aperm(array(tf,dim=c(nrow(tf),ncol(tf),nc)),c(2,3,1) )
		

}

#this is eg. the follow up on arrayFromAlleleVector to extract frequency from that spcific variant
#It accepts all dimensions, as long as dim(fr)==dim(ar)
.subsetArrayToMatrix <- function(fr, ar){
	array(aperm(fr,c(3,1,2))[aperm(ar,c(3,1,2))],dim=c(nrow(fr),ncol(fr)))
}

#this is eg. the follow up on arrayFromAlleleVector to extract frequency from that spcific variant
#It accepts dim 1 or 2, and is the one
.expandMatrixToArray <- function(mat, len){
	array(mat, dim=c(nrow(mat),ncol(mat), len))
}

#when we need to make sure the input is character
#returns target as character and issues a warning if the input wasnt
.verboseCoerceToCharacter <- function(x){
	if(!class(x)[1]=="character"){
		warning("arg was not of type character, and will therefore be coerced to character")
		as.character(x)
	}else x
}

#convert NA to FALSE
.Na2False <- function(x){
	if(class(x[1])[1]=="numeric" | class(x[1])[1]=="integer") 
		stop("not possible to replace a numeric or integer with FALSE and get the expected result")
	x[is.na(x)] <- FALSE
	x
}

