#'@include barplot-methods.R
NULL

#' locationplot ASEset objects
#' 
#' plotting ASE effects over a specific genomic region
#' 
#' The locationplot methods visualises how fractions are distributed over a
#' larger region of genes on one chromosome. It takes and ASEset object as well
#' as additional information on plot type (see \code{\link{barplot}}), strand
#' type (see \code{\link{getAlleleCounts}}), colouring, as well as annotation.
#' The annotation is taken either from the bioconductor OrgDb-sets, the TxDb
#' sets or both. It is obviously important to make sure that the genome build
#' used is the same as used in aligning the RNA-seq data.
#' 
#' @name ASEset-locationplot
#' @rdname locationplot
#' @aliases ASEset-locationplot locationplot locationplot,ASEset-method
#' @docType methods
#' @param x an ASEset object.
#' @param type 'fraction' or 'count'
#' @param strand '+','-','*' or 'both'. This argument determines
#' which strand is plotted. See \code{getAlleleCounts} for more information on
#' strand.
#' @param yaxis wheter the y-axis is to be displayed or not
#' @param xaxis wheter the x-axis is to be displayed or not
#' @param ylab showing labels for the tic marks
#' @param xlab showing labels for the tic marks
#' @param ylab.text ylab text
#' @param xlab.text xlab text
#' @param legend.colnames gives colnames to the legend matrix
#' @param size will give extra space in the margins of the inner plots
#' @param main text to use as main label
#' @param pValue Display p-value
#' @param cex.main set main label size
#' @param cex.ylab set ylab label size
#' @param cex.legend set legend label size
#' @param top.fraction.criteria 'maxcount', 'ref' or 'phase'
#' @param OrgDb an OrgDb object from which to plot a gene map. If given
#' together with argument TxDb this will only be used to extract genesymbols.
#' @param TxDb a TxDb object from which to plot an exon map.
#' @param allow.whole.chromosome logical, overrides 200kb region limit, defaults to FALSE
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run
#' @param ... arguments passed on to barplot function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the
#' locationplot function can be called up on.  }
#' @keywords locationplot
#' @examples
#' 
#' 
#' data(ASEset)
#' locationplot(ASEset)
#' 
#' #SNPs are plotted in the order in which they are found. 
#' #This can be sorted according to location as follows:
#' locationplot(ASEset[order(start(rowRanges(ASEset))),])
#' 
#' #for ASEsets with fewer SNPs the 'count' type plot is
#' # useful for detailed visualization.
#' locationplot(ASEset,type='count',strand='*')
#' 
#' @export
setGeneric("locationplot", function(x, ...) {
    standardGeneric("locationplot")
})

#' @rdname locationplot
#' @export
setMethod("locationplot", signature(x = "ASEset"), function(x, type = "fraction", 
    strand = "*", yaxis = TRUE, xaxis = FALSE, xlab = FALSE, ylab = TRUE, 
    xlab.text="", ylab.text="", legend.colnames = "", size = c(0.8, 1), main = NULL, pValue = FALSE, cex.main = 0.7, 
    cex.ylab = 0.6, cex.legend = 0.5, OrgDb = NULL, TxDb = NULL, verbose = TRUE, 
    top.fraction.criteria="maxcount",allow.whole.chromosome=FALSE, ...) {

    # check basic chromosome and region requirement
    if (!length(unique(seqnames(rowRanges(x)))) == 1) {
        stop("this plot only allows one chromosome")
    }
	if(!allow.whole.chromosome){
		if ((max(end(rowRanges(x))) - min(start(rowRanges(x)))) > 2e+05) {
			stop("this plot only allows a 200kb region")
		}
	}
    
    # check type
    if (class(type)[1] != "character") 
        stop(paste("type should be of class character, not", class(type)[1]))
    if (length(type) != 1) 
        stop(paste("type should be of length 1, not", length(type)))
    okPlotTypes <- c("fraction", "count")
    if (!type %in% okPlotTypes) 
        stop(paste("type can't be '", type, "' - it should be one of these: ", paste(okPlotTypes, 
            collapse = ", "), sep = ""))
    
    # check strand
    if (class(strand)[1] != "character") 
        stop(paste("strand should be of class character, not", class(strand)[1]))
    if (length(strand) != 1) 
        stop(paste("strand should be of length 1, not", length(strand)))
    okStrandTypes <- c("both", "+", "-", "*")
    if (!strand %in% okStrandTypes) 
        stop(paste("strand can't be '", strand, "' - it should be one of these: ", 
            paste(okStrandTypes, collapse = ", "), sep = ""))
    
    # check if strand type is present in object
    if (strand == "+") {
        el <- "countsPlus"
        if (!(el %in% assayNames(x))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else if (strand == "-") {
        el <- "countsMinus"
        if (!(el %in% assayNames(x))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else if (strand == "*") {
        el <- c("countsPlus","countsMinus")
        if (sum(el %in% assayNames(x))==0) {
            stop("neither plus or minus strand not present as assay in ASEset object")
        }
    } else if (strand == "both") {
        el <- "countsPlus"
        if (!(el %in% assayNames(x))) {
            stop("strand is not present as assay in ASEset object")
        }
        el <- "countsMinus"
        if (!(el %in% assayNames(x))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else {
        stop("unknown strand option")
    }
    
    
    # check if annotation is present, and if it is, then check if it is the right
    # class.
    if (!is.null(OrgDb)) {
        if (!class(OrgDb)[1] %in% c("OrgDb")) 
            stop(paste("if given, annotation should be of class OrgDb, not", class(OrgDb)[1]))
    }
    if (!is.null(TxDb)) {
        if (!class(TxDb)[1] %in% c("TxDb")) 
            stop(paste("if given, annotation should be of class TxDb, not", class(TxDb)[1]))
    }
    
    # check verbose argument is a logical with length 1
    if (class(verbose)[1] != "logical") 
        stop(paste("verbose should be of class logical, not", class(verbose)[1]))
    if (length(verbose) != 1) 
        stop(paste("verbose should be of length 1, not", length(verbose)))
    
    
    # get annotation
    if (!is.null(OrgDb)) {
        if (verbose) 
            message("extracting genes from OrgDb annotation")
        if (is.null(TxDb)) {
            genesInRegion <- getGenesFromAnnotation(OrgDb, rowRanges(x), verbose = verbose)
        } else {
            genesInRegion <- getGenesFromAnnotation(OrgDb, rowRanges(x), getUCSC = TRUE, 
                verbose = verbose)  #if TxDb is also given, then make sure to get the UCSC ids
            
        }
    }
    if (!is.null(TxDb)) {
        if (verbose) 
            message("extracting exons from TxDb annotation")
        exonsInRegion <- getExonsFromAnnotation(TxDb, rowRanges(x), verbose = verbose)
    }
    # if user gave both TxDb and OrgDb we can cross-reference and get gene names as
    # well as UCSC-IDs
    if (!is.null(TxDb) & !is.null(OrgDb)) {
        if ("UCSCKG" %in% colnames(mcols(genesInRegion))) {
            genesymbolKey <- mcols(genesInRegion)[, "Symbol"]
            names(genesymbolKey) <- mcols(genesInRegion)[, "UCSCKG"]
            # this line adds the genesymbol from OrgDb to each of the tx_name entries
            renamed <- as(lapply(mcols(exonsInRegion)[, "tx_name"], function(x, genesymbolKey) {
                paste(genesymbolKey[x], " (", x, ")", sep = "")
            }, genesymbolKey), "CharacterList")
            mcols(exonsInRegion)[, "tx_name"] <- renamed
        } else {
            if (verbose) 
                message("The OrgDb object did not contain cross-reference to TxDb object. No genesymbols were added")
        }
    }
    
    
    # check yaxis
    if (!is.logical(yaxis)) {
        stop("yaxis has to be logical, TRUE or FALSE")
    }
    if (!length(yaxis) == 1) {
        stop("yaxis has to be of length 1")
    }
    
    # check xaxis
    if (!is.logical(xaxis)) {
        stop("xaxis has to be logical, TRUE or FALSE")
    }
    if (!length(xaxis) == 1) {
        stop("xaxis has to be of length 1")
    }
    
    # check ylab
    if (!is.logical(ylab)) {
        stop("ylab has to be logical, TRUE or FALSE")
    }
    if (!length(ylab) == 1) {
        stop("ylab has to be of length 1")
    }
    
    # check xlab
    if (!is.logical(xlab)) {
        stop("xlab has to be logical, TRUE or FALSE")
    }
    if (!length(xlab) == 1) {
        stop("xlab has to be of length 1")
    }
    
    # check legend.colnames
    if (!class(legend.colnames)[1] == "character") {
        stop("legend.colnames has to be of class character")
    }
    
    # check pValue
    if (!is.logical(pValue)) {
        stop("pValue has to be logical, TRUE or FALSE")
    }
    if (!length(pValue) == 1) {
        stop("pValue has to be of length 1")
    }
    
    # check main
    if (!is.null(main)) {
        if (!class(main)[1] == "character") {
            stop("main has to be of class character")
        }
    }
    
    # check cex.main
    if (!class(cex.main)[1] == "numeric") {
        if (cex.main <= 0 | !length(cex.main) == 1) 
            stop("cex.main has to be of class numeric, and have length 1, and a value above 0")
    }
    # check cex.ylab
    if (!class(cex.ylab)[1] == "numeric") {
        if (cex.ylab <= 0 | !length(cex.ylab) == 1) 
            stop("cex.ylab has to be of class numeric, and have length 1, and a value above 0")
    }
    # check cex.legend
    if (!class(cex.legend)[1] == "numeric") {
        if (cex.legend <= 0 | !length(cex.legend) == 1) 
            stop("cex.legend has to be of class numeric, and have length 1, and a value above 0 ")
    }
    
    # check OrgDb
    if (!is.null(OrgDb)) {
        if (!class(OrgDb)[1] == "OrgDb") {
            stop("class of OrgDb has to be an OrgDb class")
        }
    }
    if (!is.null(TxDb)) {
        if (!class(TxDb)[1] == "TxDb") {
            stop("class of TxDb has to be a TxDb class")
        }
    }
    
    
    
    # begin calculations
    chromosome <- unique(as.character(seqnames(rowRanges(x))))
    
    xlim <- range(min(start(rowRanges(x))), max(end(rowRanges(x))))
    
    # increase xlim borders by 10%
    xlim[2] <- xlim[2] + (xlim[2] - xlim[1]) * 0.1
    xlim[1] <- max(c(0, xlim[1] - (xlim[2] - xlim[1]) * 0.1))
    
    
    # standardize ylim across all SNPs
    if (type == "count") {
        if (strand != "both") {
            ylim <- c(0, max(sapply(alleleCounts(x, strand), max)))
        } else {
            maxPlus <- max(sapply(alleleCounts(x, "+"), max))
            maxMinus <- max(sapply(alleleCounts(x, "-"), max))
            max <- max(c(maxPlus, maxMinus))
            ylim <- c(-max, max)
        }
    } else {
        ylim <- NULL
    }
    
    
    # start device
    plot.default(NULL, xlim = xlim, ylim = c(-0.5, 1.1), ylab = "", xlab = paste("genomic position on chromosome", 
        chromosome), yaxt = "n")
    
    
    for (i in 1:nrow(x)) {
        
        # retrieve genomic position
        genomicPos <- start(rowRanges(x[i]))
        # calculate on-plot position and barplot size (for evenly spaced barplots)
        sizeHere <- c((xlim[2] - xlim[1])/nrow(x), 1)
        lowerLeftCorner <- c(xlim[1] + (i - 1) * sizeHere[1], 0)
        sizeHere <- sizeHere * size
        
        # do bar plots
        barplot(x[i], type = type, strand = strand, lowerLeftCorner = lowerLeftCorner, 
            size = sizeHere, addHorizontalLine = 0.5, add = TRUE, pValue = pValue, 
            cex.ylab = cex.ylab, legend.colnames = legend.colnames, yaxis = yaxis, 
            xaxis = xaxis, ylab = ylab, xlab = xlab, main = main, cex.main = cex.main, 
            cex.legend = cex.legend, top.fraction.criteria=top.fraction.criteria,
			ylab.text=ylab.text, xlab.text=xlab.text)#, ...)
        
        # create lines indicating at which genomic position the Snp is found
        lines(x = rep(lowerLeftCorner[1] + size[1]/2, 2), y = c(0, -0.1), col = "dodgerblue")
        lines(x = c(lowerLeftCorner[1] + size[1]/2, genomicPos), y = c(-0.1, -0.2), 
            col = "dodgerblue")
        lines(x = c(genomicPos, genomicPos), y = c(-0.2, -0.25), col = "dodgerblue")
        points(x = genomicPos, y = -0.25, pch = 16)
        
    }
    
    # indicate genomic position on a horizontal line a bit below the plots (at
    # y=-0.25)
    abline(h = -0.25, col = "darkblue")
    
    # try to get the genes and put them on the plot
    if (!is.null(TxDb)) {
        if (!length(exonsInRegion) == 0) {
            decorateWithExons(x, exonsInRegion, xlim = xlim, ylim = c(-0.5, -0.3), 
                chromosome)
        }
	}
    if (!is.null(OrgDb)) {
        if (!length(genesInRegion) == 0) {
            decorateWithGenes(x, genesInRegion, xlim = xlim, ylim = c(-0.5, -0.3), 
                chromosome)
        }
    }
})


#' glocationplot ASEset objects
#' 
#' plotting ASE effects over a specific genomic region using Gviz functionality
#' 
#' The glocationplot methods visualises the distribution of ASE over a larger
#' region on one chromosome. It takes and ASEset object as well as additional
#' information on plot type (see \code{\link{gbarplot}}), strand type (see
#' \code{\link{getAlleleCounts}}), Annotation tracks are created from the Gviz
#' packageh. It is obviously important to make sure that the genome build used
#' is set correctly, e.g. 'hg19'.
#' 
#' sizes has to be of the same length as the number of tracks used.
#' 
#' @name ASEset-glocationplot
#' @aliases ASEset-glocationplot glocationplot glocationplot,ASEset-method
#' @docType methods
#' @param x an ASEset object.
#' @param type 'fraction' or 'count'
#' @param strand '+','-','*' or 'both'. This argument determines which strand is
#' plotted. See \code{getAlleleCounts} for more information of choice of strand.
#' @param BamGAL GAlignmentsList covering the same genomic region as the ASEset
#' @param GenomeAxisTrack include an genomic axis track
#' @param add add to existing plot
#' @param TxDb a TxDb object which provides annotation
#' @param sizes vector with the sum 1. Describes the size of the tracks 
#' @param trackNameDeAn  trackname for deAnnotation track
#' @param verbose if set to TRUE it makes function more talkative
#' @param ... arguments passed on to barplot function
#' @author Jesper R. Gadin
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the
#' glocationplot function can be called up on.  }
#' @keywords glocationplot
#' @examples
#' 
#' data(ASEset)
#' genome(ASEset) <- 'hg19'
#' 
#' glocationplot(ASEset,strand='+')
#' 
#' #for ASEsets with fewer SNPs the 'count' type plot is useful 
#' glocationplot(ASEset,type='count',strand='+')
#' 

#' @exportMethod glocationplot
setGeneric("glocationplot", function(x, type = "fraction", strand = "*", 
    BamGAL = NULL, GenomeAxisTrack = FALSE, trackNameDeAn = paste("deTrack", type), 
    TxDb=NULL, sizes=NULL, add = FALSE, verbose = FALSE, ...) {
    standardGeneric("glocationplot")
})

setMethod("glocationplot", signature(x = "ASEset"), function(x, type = "fraction", 
    strand = "*", BamGAL = NULL, GenomeAxisTrack = FALSE, trackNameDeAn = paste("deTrack", 
        type), TxDb=NULL, sizes=NULL, add = FALSE, verbose = FALSE, ...) {
   
	#tmp
	#if(!is.null(TxDb)){stop("the functionality with TxDb is not yet ready")}

    # check genome
    if (is.null(genome(x)) | is.na(genome(x))) {
        stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))
    }
    # check type argument
    if (!(type %in% c("fraction", "count"))) {
        stop(paste("type has to be", "fraction", "or", "count"))
    }
    
    # check seqnames has length=0
    if (!(length(seqlevels(x)) == 1)) {
        stop("This function can only use objects with one seqlevel")
    }
    
    #if (sum(strand == "+" | strand == "-") == 0) {
    #    stop("strand must be plus or minus at the moment")
    #}

    if (!nrow(x) == 1) {
       
		if(strand %in% c("+","-","*")){
			GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
				end = max(end(x))), strand = strand, genome = genome(x))
		}else if (strand=="both"){
			GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
				end = max(end(x))), strand = "*", genome = genome(x))
		}else{
			stop("strand has to be +, -, * or 'both'")
		}
    }

    
	#make an environment from ...
    if (length(list(...)) == 0) {
        e <- new.env(hash = TRUE)
    } else {
        e <- list2env(list(...))
    }

	e$x <- x


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

    # make deTrack 
    if (verbose) {
        (cat("preparing detailedAnnotationTrack\n"))
	}
    deTrack <- ASEDAnnotationTrack(x, GR = GR, type, strand, 
								   trackName = trackNameDeAn,
								   mainvec=e$mainvec,
								   cex.mainvec=e$cex.mainvec,
								   ylab=e$ylab,
								   xlab=e$xlab,
								   top.fraction.criteria=e$top.fraction.criteria
								   )
    lst <- list(deTrack)
   
	start <- min(start(GR))
	end <- max(end(GR))



    if (!is.null(BamGAL)) {
        if (verbose) {
            (cat("preparing coverageDataTrack\n"))
		}

        seqlevels(BamGAL) <- seqlevels(x)
        start <- min(start(x))
        end <- max(end(x))
        
        covTracks <- unlist(CoverageDataTrack(x, BamList = BamGAL, strand = strand, meanCoverage=TRUE))
        

        #lst[[length(lst) + 1]] <- covTracks
		lst <- c(lst,covTracks)
	}


	if(!is.null(TxDb)){
        if (verbose) {
            (cat("preparing transcriptDB track\n"))
		}
		txTrack <- GeneRegionTrack(TxDb, 
			start=start(GR), end=end(GR), 
			chromosome=seqlevels(GR)
		)	   

        lst[[length(lst) + 1]] <- txTrack
	}


    if (GenomeAxisTrack) {
        if (verbose) {
            (cat("preparing GenomeAxisTrack\n"))
		}
        axTrack <- GenomeAxisTrack()
        #lst[[length(lst) + 1]] <- axTrack
		lst <- c(lst,axTrack)
    }


	#set sizes
	parts <- 1/length(lst) #need mean coverage
	if(is.null(sizes)) {
		sizes <- c(rep(parts, length(lst)))
	}else{
		if(!length(sizes)==length(lst)){
			stop("sizes vector has to be same length as the number of tracks used")
		}
	}
	#sizes <- 1
	# plot
	plotTracks(lst, from = start, to = end, sizes = sizes, col.line = NULL, showId = FALSE, 
		title.width = 1, type = "histogram", 
		add = add)


}) 
