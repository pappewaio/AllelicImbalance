#'@include locationplot-methods.R
NULL

#' ASEset-gviztrack ASEset objects
#' 
#' plotting ASE effects over a specific genomic region
#' 
#' For information of how to use these tracks in more ways, visit the Gviz
#' package manual.
#' 
#' @name ASEset-gviztrack
#' @rdname ASEset-gviztrack
#' @aliases ASEset-gviztrack CoverageDataTrack ASEDAnnotationTrack
#' CoverageDataTrack,ASEset-method ASEDAnnotationTrack,ASEset-method
#' @docType methods
#' @param x an ASEset object.
#' @param GR genomic range of plotting
#' @param type 'fraction' or 'count'
#' @param strand '+','-'. This argument determines which strand is plotted.
#' @param BamList GAlignmnentsList object of reads from the same genomic region
#' as the ASEset
#' @param start start position of reads to be plotted
#' @param trackName name of track (ASEDAnnotationTrack)
#' @param meanCoverage mean of coverage over samples (CoverageGataTrack)
#' @param trackNameVec names of tracks (CoverageDataTrack)
#' @param end end position of reads to be plotted
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run
#' @param ... arguments passed on to barplot function
#' @author Jesper R. Gadin
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the functions
#' can be called up on.}
#' @keywords ASEDAnnotationTrack CoverageDataTrack
#' @examples
#' 
#' data(ASEset)
#' x <- ASEset[,1:2]
#' r <- reads[1:2]
#' genome(x) <- 'hg19'
#' seqlevels(r) <- seqlevels(x)
#' 
#' GR <- GRanges(seqnames=seqlevels(x),
#'			ranges=IRanges(start=min(start(x)),end=max(end(x))),
#'			strand='+', genome=genome(x))
#' 
#' deTrack <- ASEDAnnotationTrack(x, GR=GR, type='fraction',strand='+')
#' covTracks <- CoverageDataTrack(x,BamList=r,strand='+') 
#' 
#' lst <- c(deTrack,covTracks)
#' 
#' sizes <- c(0.5,rep(0.5/length(covTracks),length(covTracks)))
#' #temporarily do not run this function 
#' #plotTracks(lst, from=min(start(x)), to=max(end(x)), 
#' #sizes=sizes, col.line = NULL, showId = FALSE, main='mainText', 
#' #cex.main=1, title.width=1, type='histogram')
#' 
#'
NULL


#' @rdname ASEset-gviztrack
#' @export
setGeneric("ASEDAnnotationTrack", function(x, GR = rowRanges(x), type = "fraction", 
    strand = "*", trackName = paste("deTrack", type), verbose = TRUE, 
    ...) {
    standardGeneric("ASEDAnnotationTrack")
})

#' @rdname ASEset-gviztrack
#' @export
setMethod("ASEDAnnotationTrack", signature(x = "ASEset"), function(x, GR = rowRanges(x), 
    type = "fraction", strand = "*",  trackName = paste("deTrack", 
        type), verbose = TRUE, ...) {
    
    # check genome
    if (is.null(genome(x)) | is.na(genome(x))) {
        stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))
    }
    
    # check seqnames has length=1
    if (!(length(seqlevels(x)) == 1)) {
        stop("This function can only use objects with one seqlevel")
    }
    
    if (!nrow(x) == 0) {
		if(strand %in% c("+","-","*")){
			GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
				end = max(end(x))), strand = strand, genome = genome(x))
		}else if (strand=="both"){
			GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
				end = max(end(x))), strand = "*", genome = genome(x))
		}else{
			stop("strand has to be +, -, * or 'both'")
		}
	}else{
		stop("there are no rows in ASEset")
	}

	#make an environment from ...
    if (length(list(...)) == 0) {
        e <- new.env(hash = TRUE)
    } else {
        e <- list2env(list(...))
    }

	#always TRUE
	e$deAnnoPlot <- TRUE
	e$x <- x
	
	#print(ls(envir=e))
	#print(e$mainvec)	

    if (!exists("mainvec", envir = e, inherits = FALSE)) {
		e$mainvec <- rep("",nrow(e$x))
	}

    if (!exists("cex.mainvec", envir = e, inherits = FALSE)) {
		e$cex.mainvec <- 1
	}

    if (!exists("ylab", envir = e, inherits = FALSE)) {
        e$ylab <- ""
    }
    if (!exists("xlab", envir = e, inherits = FALSE)) {
        e$xlab <- ""
    }
		
    if (!exists("middleLine", envir = e, inherits = FALSE)) {
        e$middleLine <- TRUE
    }

    if (!exists("top.fraction.criteria", envir = e, inherits = FALSE)) {
        e$top.fraction.criteria <- "maxcount"
    }

    ranges <- rowRanges(x)
    
    colnames(x) <- 1:ncol(x)
    
    details <- function(identifier, ...) {
        
		if (length(list(...)) == 0) {
			e2<- new.env(hash = TRUE)
		} else {
			e2<- list2env(list(...))
		}

		type<- e2$type 
        
        if (type == "fraction") {
            print(barplotLatticeFraction(identifier, 
                ...), newpage = FALSE, prefix = "plot")
            
        } else if (type == "count") {
            print(barplotLatticeCounts(identifier, 
                ...), newpage = FALSE, prefix = "plot")
        }
        
    }

    # plot the fraction
    deTrack <- AnnotationTrack(
		range = ranges,
		genome = genome(x),
		id = rownames(x),
		name = trackName,
		stacking = "squish",
		fun = details,
		detailsFunArgs = list(
			ylab = e$ylab,
			xlab = e$xlab,
			deAnnoPlot = e$deAnnoPlot,
			cex.mainvec = e$cex.mainvec,
			mainvec=list(list(e$mainvec)),
			type = type,
			x = x,
			astrand = strand,
			ids = list(list(rownames(x))),
			middleLine = e$middleLine,
			top.fraction.criteria = e$top.fraction.criteria))
    deTrack
})

#' @rdname ASEset-gviztrack
#' @export
setGeneric("CoverageDataTrack", function(x, GR = rowRanges(x), BamList = NULL, strand = NULL, 
    start = NULL, end = NULL, trackNameVec = NULL, meanCoverage=FALSE, verbose = TRUE, ...) {
    standardGeneric("CoverageDataTrack")
#' @export
})
setMethod("CoverageDataTrack", signature(x = "ASEset"), function(x, GR = rowRanges(x), 
	 BamList = NULL, strand = "*", start = NULL, end = NULL, trackNameVec = NULL,
	 meanCoverage=FALSE, verbose = TRUE, ...) {
    
     # GR is not in use atm. Missing is a subset of the return matrix based on the GR
     # values.
    
    if (!is.null(strand)) {
        if (strand == "+") {
            pstrand = TRUE
        } else if (strand == "-") {
            mstrand = TRUE
        } else if (strand == "*") {
            pstrand = TRUE
            mstrand = TRUE
        } else if (strand == "both") {
            pstrand = TRUE
            mstrand = TRUE
        } else {
            stop("strand has to be '+', '-', '*' or 'both' if not NULL\n")
        }
    } else {
         stop("strand has to be '+', '-', '*' or 'both' if not NULL\n")
    }
    
    # check genome
    if (is.null(genome(x)) | is.na(genome(x))) {
        stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))
    }
    
    # check seqnames has length=0
    if (!(length(seqlevels(x)) == 1)) {
        stop("This function can only use objects with one seqlevel")
    }
   
	if(strand=="both"){
		GR.p <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
			end = max(end(x))), strand = '+')
		GR.m <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
			end = max(end(x))), strand = '-')
	}else{
		GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
			end = max(end(x))), strand = strand)
	}
    
    trackList <- list()
	
    if (is.null(BamList)) {
        stop("must include GappedAlignmentsList as BamList ")
    } 
	
	# check that only one chromosome is present
	if (!length(seqlevels(BamList)) == 1) {
		stop("can only be one seq level\n")
	}

	if(strand=="both"){
		covMatList <- coverageMatrixListFromGAL(BamList, strand='both')

		mat.p <- covMatList[["mat"]][["plus"]]
		mat.m <- covMatList[["mat"]][["minus"]]
		start <- covMatList[["start"]]
		end <- covMatList[["end"]]

		if(meanCoverage){
			vec.p <- apply(mat.p,2,mean)
			mat.p <- matrix(vec.p,nrow=1)	
			dimnames(mat.p)[1] <- list("meancoverage")
			vec.m <- apply(mat.m,2,mean)
			mat.m <- matrix(vec.m,nrow=1)	
			dimnames(mat.m)[1] <- list("meancoverage")
		}



		if (is.null(trackNameVec)) {
			trackNameVec <- vector()

			if(!meanCoverage){
				for (i in 1:ncol(x)){
					trackNameVec <- c(trackNameVec,
									  paste(colnames(x)[i],"(+)",sep="-"),
									  paste(colnames(x)[i],"(-)",sep="-"))
				}
			 }else{
					trackNameVec <- c(trackNameVec,
									  paste(colnames(x)[i],"(+)",sep="-"),
									  paste(colnames(x)[i],"(-)",sep="-"))
			}
		} else {
			if (!((length(trackNameVec)*2) == nrow(mat.p))) {
				stop(paste("length of trackNameVec must be a two times the numer",
					 "of rows in out matrix"))
			}
		}
			
		# prepare Gviz dtracks
		for (j in 1:nrow(mat.p)) {

			#merge strands
			data <- matrix(0,nrow=2,ncol=ncol(mat.m))
			data[1,] <- mat.p[j,]
			data[2,] <- -mat.m[j,]
			rownames(data) <- c(paste(colnames(x)[i], "(+)", sep="-"),
								paste(colnames(x)[i], "(-)", sep="-"))

			trackList[[length(trackList) + 1]] <- DataTrack(data = data, start = start:end, 
				width = 1, chromosome = seqlevels(x), genome = genome(x),
				name = trackNameVec[j], groups=rownames(data),
				type = "s")
		}

	}else{

		covMatList <- coverageMatrixListFromGAL(BamList, strand)

		mat <- covMatList[["mat"]]
		start <- covMatList[["start"]]
		end <- covMatList[["end"]]
		
		if(meanCoverage){
			vec <- apply(mat,2,mean)
			mat <- matrix(vec,nrow=1)	
			dimnames(mat)[1] <- list("meancoverage")
		}


		if (is.null(trackNameVec)) {
			if(!meanCoverage){
				trackNameVec <- vector()
				trackNameVec[1:ncol(x)] <- colnames(x)
			}else{
				trackNameVec <- "mean coverage"
			}
		} else {
			if (!(length(trackNameVec) == nrow(mat))) {
				stop("length of trackNameVec must be equal to cols in (ASEset)")
			}
		}
		
		# prepare Gviz dtracks
		for (j in 1:nrow(mat)) {
			trackList[[length(trackList) + 1]] <- DataTrack(data = mat[j, ], start = start:end, 
				width = 1, chromosome = seqlevels(x), genome = genome(x), name = trackNameVec[j], 
				type = "s")
		}
	}
    
    trackList
}) 
