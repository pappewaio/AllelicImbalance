#######################
# Deprecated functions
#######################

#was deprecated 2015-03-26
implodeList <- function()
{
	    .Deprecated("implodeList.old",msg="no longer serving any purpose in the AllelicImbalance package and is deprecated")
}



#was deprecated 2015-03-25
realCigarPositions <- function()
{
	    .Deprecated("realCigarPositions.old",msg="no longer serving any purpose in the AllelicImbalance package and is deprecated")
}

#was deprecated 2015-03-25
realCigarPositionsList <- function()
{
	    .Deprecated("realCigarPositionsList.old",msg="no longer serving any purpose in the AllelicImbalance package and is deprecated")
}

#was deprecated 2015-03-25
realCigarPosition <- function()
{
	    .Deprecated("realCigarPosition",msg="no longer serving any purpose in the AllelicImbalance package and is deprecated")
}

#was deprecated 2015-03-24
impBamGRL <- function()
{
	    .Deprecated("impBamGRL.old",msg="no longer serving any purpose in the AllelicImbalance package and is deprecated")
}

#in 2014
getAlleleCount <- function() {
    .Deprecated("getAlleleCounts")
}

#' Import Bam-2
#'
#' Imports bla bal bal a specified genomic region from a bam file using a GenomicRanges
#' object as search area.
#'
#' These functions are right  on tahea wrappers to import bam files into R and store them into
#' either GRanges, GAlignments or GappedAlignmentpairs objects.
#'
#' It is recommended to use the impBamGAL() which takes information of gaps
#' into account. It is also possible to use the other variants as well, but
#' then pre-filtering becomes important keps to understand because gapped, intron-spanning reads
#' will cause problems. This is because the GRanges objects can not handle if
#' gaps are present and will then give a wrong result when calculating the
#' allele (SNP) count table.
#'
#' @name import-bam-2
#' @rdname import-bam-2
#' @aliases import-bam-2 impBamGRL
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges object} that contains the regions of
#' interest
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run.
#' @return \code{impBamGRL} returns a GRangesList object containing the RNA-seq
#' reads in the region defined by the \code{searchArea} argument.
#' \code{impBamGAL} returns a list with GAlignments objects containing the
#' RNA-seq reads in the region defined by the \code{searchArea} argument.
#' \code{funImpBamGAPL} returns a list with GappedAlignmentPairs object
#' containing the RNA-seq reads in the region defined by the \code{searchArea}
#' argument.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bam import
#' @examples
#'
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c('17'), ranges=IRanges(79478301,79478361))
#'
#' #Relative or full path
#' pathToFiles <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
#'
#'
#' @export impBamGRL
NULL


#' @rdname import-bam-2
impBamGRL.old <- function(UserDir, searchArea, verbose = TRUE) {
    # Set parameters
    which <- searchArea  #A GRanges, IntegerRangesList, or missing object, from which a IRangesList instance will be constructed.
    what <- scanBamWhat()  #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    param <- ScanBamParam(flag = flag, which = which, what = what)  #store ScanBamParam in param.

    # Point to correct directory and create a BamFileList object
    bamDir <- normalizePath(UserDir)
    allFiles <- list.files(bamDir, full.names = TRUE)
    bamFiles <- allFiles[grep(".bam$", allFiles)]
    if (length(bamFiles) == 0) {
        stop(paste("No bam files found in", bamDir))
    }
    if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
        if (verbose) {
            cat(paste("The bam files in UserDir are required to also have", ".bam.bai index files.",
                " Trying to run indexBam function on each", "\n"), )
        }
        indexBam(bamFiles)
        if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
            stop("The bam files in UserDir are required to also have", ".bam.bai index files.")
        } else {
            if (verbose) {
                cat(paste("Succesfully indexed all bamFiles in UserDir", UserDir,
                  "\n"))
            }
        }
    }
    # store all the .bam paths in a BamFile.
    bamFilesList <- BamFileList(bamFiles)

    # check that sequences in searchArea are actually found in the bam files
    header <- scanBamHeader(bamFiles)
    checkSeqNameExists <- function(bamHeader, requestedSeqNames) {
        as.character(requestedSeqNames) %in% names(bamHeader[["targets"]])
    }
    if (!all(unlist(lapply(header, checkSeqNameExists, seqnames(searchArea))))) {
        # not all searchArea requested seq-names found in bam files. Create nice error
        # report and stop
        seqNotFoundErrors <- lapply(header, checkSeqNameExists, seqnames(searchArea))
        seqNotFounds <- vector()
        for (sampleName in names(seqNotFoundErrors)) {
            seqNotFounds <- c(seqNotFounds, as.character(seqnames(searchArea)[!seqNotFoundErrors[[sampleName]]]))
        }
        stop(paste("The following seq name(s) not found in the bam files:", paste(sort(unique(seqNotFounds)),
            collapse = ", ")))
    }

    # Loop through, open scanBam, store in GRList and then close each object in the
    # BamFileList object.
    i <- 1
    BamGRL <- GRangesList()
    for (bamName in names(bamFilesList)) {
        # Description
        bf <- bamFilesList[[bamName]]
        open(bf)
        if (verbose) {
            cat(paste("Reading bam file", i, "with filename", basename(bamName)),
                "\n")
        }
        bam <- scanBam(bf, param = param)
        # Description
        for (rangeName in names(bam)) {

            # if NA values your in trouble. That means the read didnt map
            ranges <- IRanges(start = bam[[rangeName]][["pos"]], width = cigarWidthAlongReferenceSpace(bam[[rangeName]][["cigar"]]))
            GRangeBam <- GRanges(seqnames = as.character(bam[[rangeName]][["rname"]]),
                ranges = ranges, strand = bam[[rangeName]][["strand"]], names = bam[[rangeName]][["qname"]],
                flag = bam[[rangeName]][["flag"]], cigar = bam[[rangeName]][["cigar"]],
                mapq = bam[[rangeName]][["mapq"]], mpos = bam[[rangeName]][["mpos"]],
                isize = bam[[rangeName]][["isize"]], seq = bam[[rangeName]][["seq"]],
                qual = bam[[rangeName]][["qual"]])
            # This way of merging the different chromosomes to the same GRangeObject is maybe
            # not the best way. Later try to store them in a separate list, and then unlist
            # before importing to GrangeBam Store GRangeBam in BamGRL (which is the GRange
            # List object)
            if (basename(bamName) %in% names(BamGRL)) {
                BamGRL[[basename(bamName)]] <- c(BamGRL[[basename(bamName)]], GRangeBam)
            } else {
                BamGRL[[basename(bamName)]] <- GRangeBam
            }
        }
        if (verbose) {
            cat(paste("stored", basename(bamName), "in BamGRL"), "\n")
        }
        i <- 1 + i
        gc()
        close(bf)
    }
    return(BamGRL)
}



#'@include ASEset-class.R
NULL

#' realCigarPosition
#'
#' From a GAlignments calculate the real corresponding position for each read
#' based on its cigar.
#'
#' The main intention for these functions are to be the internal functions for
#' \code{scanForHeterozygotes} and \code{getAlleleCount}.
#'
#' @name cigar-utilities
#' @rdname cigar-utilities
#' @aliases cigar-utilities realCigarPosition realCigarPositions
#' realCigarPositionsList
#' @param RleCigar An \code{Rle} containing cigar information
#' @param RleCigarList An \code{RleList} containing cigar information
#' @param BpPos the absolute position on the chromosome of interest
#' @return realCigarPosition returns the new position
#' realCigarPositions returns a vector with the corrected positions to
#' be subsetted from a read.  \code{realCigarPositionsList} returns a list
#' where each element i a vector with the corrected positions to be subsetted
#' from a read.
#' @author Jesper R. Gadin
#' @keywords internal
#' @examples
#'
#'   RleCigarList <-  cigarillo::cigars_as_RleList('3M4I93M')
#'   BpPos <- 5
#'
#'   newPos <- realCigarPosition.old(RleCigar=RleCigarList[[1]], BpPos)
#'   newPositions <- realCigarPositions.old(RleCigar=RleCigarList[[1]])
#'   newPositionsList <- realCigarPositionsList.old(RleCigarList=RleCigarList)
NULL

#' @rdname cigar-utilities
#' @export
realCigarPosition.old <- function(RleCigar, BpPos) {

    # because of speed issues, checks are best performed outside this function.
    if (!is(RleCigar, "Rle")) {
        stop("RleCigar must an Rle object")
    }

    e <- as.character(RleCigar)

    # changeVector
    v <- rep(0, length = length(e))
    names(v) <- e

    v[e == "M"] <- 1
    v[e == "I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar) == "I"],
        function(x) {
            c(x + 1, rep(1, x - 1))
        }))
    # v[e=='D'] <- 0 #already zero v[e=='N'] <- 0 #already zero

    # sum all until interesting position
    cs <- cumsum(v)

    if (names(cs[BpPos]) == "D") {
        retPos <- 0
    } else if (names(cs[BpPos]) == "N") {
        retPos <- -1
    } else {
        retPos <- cs[BpPos]
        names(retPos) <- NULL
        if (retPos > sum(e == "M" | e == "I")) {
            retPos <- -1  # the position went outside the read
        }
    }
    retPos
}



#' @rdname cigar-utilities
#' @export
realCigarPositions.old <- function(RleCigar) {
    # returns a vector that have order all positions to match with the cigar

    # because of speed issues, checks are best performed outside this function.
    if (!is(RleCigar, "Rle")) {
        stop("RleCigar must be an Rle object")
    }


    e <- as.character(RleCigar)
    # make a new representation vector
    v <- rep(0, length = length(e))
    names(v) <- e

    v[e == "M"] <- 1
    v[e == "I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar) == "I"],
        function(x) {
            c(x + 1, rep(1, x - 1))
        }))
    # v[e=='D'] <- 0 v[e=='N'] <- 0

    # sum all until interesting position
    cs <- cumsum(v)

    cs <- cs[!names(cs) == "D"]
    cs <- cs[!names(cs) == "N"]

	#remove matches overexeeding the length of the read
	cs <- cs[!cs>length(v)]

    cs
}


#' @rdname cigar-utilities
#' @export
realCigarPositionsList.old <- function(RleCigarList) {

    # because of speed issues, checks are best performed outside this function.
    if (!is(RleCigarList, "RleList")) {
        stop("RleCigarList must be an RleList object")
    }

    lapply(RleCigarList, function(RleCigar) {
        e <- as.character(RleCigar)
        # make a new representation vector
        v <- rep(0, length = length(e))
        names(v) <- e

        v[e == "M"] <- 1
        v[e == "I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar) == "I"],
            function(x) {
                c(x + 1, rep(1, x - 1))
            }))
        # v[e=='D'] <- 0 v[e=='N'] <- 0

        # sum all until interesting position
        cs <- cumsum(v)

        cs <- cs[!names(cs) == "D"]
        cs <- cs[!names(cs) == "N"]


		#remove matches overexeeding the length of the read
		cs <- cs[!cs>length(v)]

        cs
    })
}

#####
# Kept for performance comparions
# will be removed in short
#####
#' scanForHeterozygotes-old
#'
#' Identifies the positions of SNPs found in BamGR reads.
#'
#' This function scans all reads stored in a \code{GAlignmentsList} for
#' possible heterozygote positions. The user can balance the sensitivity of the
#' search by modifying the minimumReadsAtPos, maximumMajorAlleleFrequency and
#' minimumBiAllelicFrequency arguments.
#'
#' @rdname scanForHeterozygotes-old
#' @param BamList A \code{GAlignmentsList object}
#' @param minimumReadsAtPos minimum number of reads required to call a SNP at a
#' given position
#' @param maximumMajorAlleleFrequency maximum frequency allowed for the most
#' common allele. Setting this parameter lower will minimise the SNP calls
#' resulting from technical read errors, at the cost of missing loci with
#' potential strong ASE
#' @param minimumBiAllelicFrequency minimum frequency allowed for the first and
#' second most common allele. Setting a Lower value for this parameter will
#' minimise the identification of loci with three or more alleles in one
#' sample. This is useful if sequencing errors are suspected to be common.
#' @param maxReads max number of reads of one list-element allowed
#' @param verbose logical indicating if process information should be displayed
#' @return \code{scanForHeterozygotes.old} returns a GRanges object with the SNPs
#' for the BamList object that was used as input.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{getAlleleCounts}} which is a
#' function that count the number of reads overlapping a site.  }
#' @keywords scan SNP heterozygote
#' @examples
#'
#' data(reads)
#' s <- scanForHeterozygotes.old(reads,verbose=FALSE)
#'
#' @export scanForHeterozygotes.old
scanForHeterozygotes.old <- function(BamList, minimumReadsAtPos = 20, maximumMajorAlleleFrequency = 0.9,
    minimumBiAllelicFrequency = 0.9, maxReads = 15000, verbose = TRUE) {

    # if just one element of, make list (which is a convenient way of handling this
    # input type)
    if (is(BamList, "GAlignments")) {
        BamList <- GAlignmentsList(BamList)
    }
    if (is(BamList, "GRanges")) {
        BamList <- GRangesList(BamList)
    }

    if (!is.numeric(minimumReadsAtPos))
        stop(paste("minimumReadsAtPos must be a numeric vector, not a", class(minimumReadsAtPos)))
    if (length(minimumReadsAtPos) != 1)
        stop(paste("minimumReadsAtPos must be of length 1, not", length(minimumReadsAtPos)))
    if (!is.numeric(maximumMajorAlleleFrequency))
        stop(paste("maximumMajorAlleleFrequency must be a numeric vector, not a", class(maximumMajorAlleleFrequency)))
    if (length(maximumMajorAlleleFrequency) != 1)
        stop(paste("maximumMajorAlleleFrequency must be of length 1, not", length(maximumMajorAlleleFrequency)))
    if (maximumMajorAlleleFrequency < 0 | maximumMajorAlleleFrequency > 1)
        stop("maximumMajorAlleleFrequency must be between 0 and 1")
    if (!is.numeric(minimumBiAllelicFrequency))
        stop(paste("minimumBiAllelicFrequency must be a numeric vector, not a", class(minimumBiAllelicFrequency)))
    if (length(minimumBiAllelicFrequency) != 1)
        stop(paste("minimumBiAllelicFrequency must be of length 1, not", length(minimumBiAllelicFrequency)))
    if (minimumBiAllelicFrequency < 0 | maximumMajorAlleleFrequency > 1)
        stop("minimumBiAllelicFrequency must be between 0 and 1")
    if (!is.logical(verbose))
        stop(paste("verbose must be a logical, not a", class(verbose)))
    if (length(verbose) != 1)
        stop(paste("verbose must be of length 1, not", length(verbose)))

    # checking that BamList format is ok (has to be done quite carefully for the
    # list-class gappedAlignments
    if (is(BamList, "GRangesList"))
        stop("The use of GRangesList is not recommended. Use BamImpGAList or BamImpGAPList")
    if (!is(BamList, "GAlignmentsList"))
        stop("BamList has to be a GAlignmentsList object")

    # checks specific for GAlignments
    if (all(vapply(BamList, is, logical(1), "GAlignments"))) {
        if (!all(unlist(lapply(BamList, function(x) {
            all(c("cigar", "qwidth") %in% colnames(mcols(x)))
        })))) {
            stop("BamList given as GAlignmentsLists of GAlignments objects must contain mcols named qwidth and cigar")
        }
    }

    # check that we dont create a too large and memory-consuming matrix
    if (sum(unlist(lapply(BamList, length)) > maxReads) > 0) {
        stop("you may consume too much memory. If there is plenty of memory, then increase maxReads to allow more reads")
    }

    ans <- GRanges()
    chromosomeLevels <- unique(unlist(lapply(BamList, function(x) {
        levels(droplevels(runValue(seqnames(x))))
    })))

    for (chr in chromosomeLevels) {
        if (verbose)
            cat(paste("Investigating chromosome", chr), "\n")

        BamListchr <- GAlignmentsList(mapply(function(x, y) {
            x[y]
        }, BamList, seqnames(BamList) == chr))


        for (sample in 1:length(BamListchr)) {
            if (verbose)
                cat(paste("Investigating sample", sample), "out of", length(BamListchr),
                  "\n")

            # extract samples
            BamListHere <- BamListchr[[sample]]

            if (!(length(BamListHere) == 0)) {

                # then iterate over all reads
                cigarRle <- cigars_as_RleList(mcols(BamListHere)[, "cigar"])
                toKeep <- realCigarPositionsList.old(cigarRle)

                seq <- mcols(BamListHere)[, "seq"]
                charList <- strsplit(as.character(seq), "")

                start <- start(BamListHere) - min(start(BamListHere)) + 1

                nw <- max(end(BamListHere)) - min(start(BamListHere)) + 2

                # populate matrix
                new <- matrix(NA, nrow = max(end(BamListHere)) - min(start(BamListHere)) +
                  1, ncol = length(BamListHere))
                for (i in 1:ncol(new)) {
                  new[start[i]:(start[i] + length(toKeep[[i]]) - 1), i] <- charList[[i]][toKeep[[i]]]

                }

                # set rownames
                rownames(new) <- as.character(1:(nw - 1))

                new <- new[apply(!is.na(new), 1, sum) > minimumReadsAtPos, ]

                if (!nrow(new) == 0) {

                  # tabulate countsPerPosition (cpp)
                  cpp <- apply(new, 1, table)

                  TFl <- unlist(lapply(cpp, function(x) {
                    if (length(x) > 1) {

                      MajorAlleleFrequency <- x[order(x, decreasing = TRUE)[1]]/sum(x)
                      if (MajorAlleleFrequency < maximumMajorAlleleFrequency) {
                        MinorAlleleFrequency <- x[order(x, decreasing = TRUE)[2]]/sum(x)
                        if ((MinorAlleleFrequency + MajorAlleleFrequency) > minimumBiAllelicFrequency) {
                          TRUE
                        } else {
                          FALSE
                        }
                      } else {
                        FALSE
                      }
                    } else {
                      FALSE
                    }
                  }))
                  if (!all(!TFl)) {
                    GR <- GRanges(ranges = IRanges(start = (as.numeric(names(cpp[TFl])) +
                      min(start(BamListHere)) - 1), width = 1), seqnames = chr)
                    ans <- c(ans, GR)
                  }
                }
            }
        }
    }
    # merge from all individuals
    ans <- unique(ans)

    # Add a Snp name based on position
    if (!(length(ans) == 0)) {
        names(ans) <- paste("chr", seqnames(ans), "_", start(ans),
            sep = "")

    }
    return(ans)
}


#' implode list of arguments into environment
#'
#' apply on list of variables to be put in the local environment
#'
#' help the propagation of e.g. graphical paramters
#'
#' @rdname implodeList-old
#' @param x list of variables
#' @author Jesper R. Gadin
#' @keywords implode
#' @examples
#'
#' lst <- list(hungry='yes', thirsty='no')
#' implodeList.old(lst)
#' #the check ls()
#'  ls()
#' @export implodeList.old
implodeList.old <- function(x) {
    oname <- deparse(substitute(x))
    eval(parse(text = paste0("for(i in 1:length(", oname, ")){assign(names(", oname,
        ")[i],", oname, "[[i]])}")), parent.frame())
}
