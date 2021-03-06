#'@include ASEset-class.R
NULL


#' lva
#' 
#' make an almlof regression for arrays
#' 
#' internal method that takes one array with results from regionSummary
#' and one matrix with group information for each risk SNP (based on phase)
#'
#' @name lva
#' @rdname lva
#' @aliases lva,array-method
#' @docType methods
#' @param x ASEset object with phase and 'ref'/'alt' allele information
#' @param rv RiskVariant object with phase and 'ref'/'alt' allele information
#' @param region RiskVariant object with phase and alternative allele information
#' @param settings RiskVariant object with phase and alternative allele information
#' @param return.class 'LinkVariantAlmlof' (more options in future)
#' @param verbose logical, if set TRUE, then function will be more talkative
#' @param type "lm" or "nlme", "nlme" needs subject information
#' @param covariates add data.frame with covariates (only integers and numeric)
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#' 
#' data(ASEset) 
#' a <- ASEset
#' # Add phase
#' set.seed(1)
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#' 
#' phase(a) <- p
#' 
#' #add alternative allele information
#' mcols(a)[["alt"]] <- inferAltAllele(a)
#' 
#' #init risk variants
#' p.ar <- phaseMatrix2Array(p)
#' rv <- RiskVariantFromGRangesAndPhaseArray(x=GRvariants, phase=p.ar)
#'
#' #colnames has to be samea and same order in ASEset and RiskVariant
#' colnames(a) <- colnames(rv)
#'
#' # in this example each and every snp in the ASEset defines a region
#' r1 <- granges(a)
#' 
#' #use GRangesList to merge and use regions defined by each element of the
#' #GRangesList
#' r1b <- GRangesList(r1)
#' r1c <- GRangesList(r1, r1)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' r2 <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#'
#' # link variant almlof (lva)
#' lva(a, rv, r1)
#' lva(a, rv, r1b)
#' lva(a, rv, r1c)
#' lva(a, rv, r2)
#'
#' # Use covariates (integers or nuemric)
#' cov <- data.frame(age=sample(20:70, ncol(a)), sex=rep(c(1,2), each=ncol(a)/2),  
#' row.names=colnames(a))
#' lva(a, rv, r1, covariates=cov)
#' lva(a, rv, r1b, covariates=cov)
#' lva(a, rv, r1c, covariates=cov)
#' lva(a, rv, r2, covariates=cov)
#'
#' # link variant almlof (lva), using nlme
#' a2 <- a
#' ac <- assays(a2)[["countsPlus"]]
#' jit <- sample(c(seq(-0.10,0,length=5), seq(0,0.10,length=5)), size=length(ac) , replace=TRUE)
#' assays(a2, withDimnames=FALSE)[["countsPlus"]] <- round(ac * (1+jit),0)
#' ab <- cbind(a, a2)
#' colData(ab)[["subject.group"]] <- c(1:ncol(a),1:ncol(a))
#' rv2 <- rv[,c(1:ncol(a),1:ncol(a))]
#' colnames(ab) <- colnames(rv2)
#'
#' lva(ab, rv2, r1, type="nlme")
#' lva(ab, rv2, r1b, type="nlme")
#' lva(ab, rv2, r1c, type="nlme")
#' lva(ab, rv2, r2, type="nlme")
#'
#' 
NULL

#' @rdname lva
#' @export
setGeneric("lva", function(x, ... 
	){
    standardGeneric("lva")
})

#' @rdname lva
#' @export
setMethod("lva", signature(x = "ASEset"),
		function(x, rv, region, settings=list(),
				 return.class="LinkVariantAlmlof", type="lm",
				 verbose=FALSE, covariates=matrix(),  ...
	){

		#safety check
		if(any(!colnames(x) %in% colnames(rv)) | any(!colnames(rv) %in% colnames(x)))
				stop("missmatch of colnames for x and rv")

		if("threshold.distance" %in% names(settings)){
			distance <- settings[["threshold.distance"]]
		}else{
			distance <- 200000
		}
		#if(nrow(covariates)==ncol(x)){stop("nrow(covariates) must match ncol(x)")}

		#region summary
		rs <- regionSummary(x, region)


		#match riskVariant to rs granges
		hits <- findOverlaps(rv, granges(rs) + distance)
		#stop if no overlap
		if(length(hits)==0){stop(paste("no rs and rv are close with distance: ",distance, sep=""))}
		#if overlap, then subset hits
		rs2 <- rs[subjectHits(hits)]
		rv2 <- rv[queryHits(hits),, drop=FALSE]
		#make groups for regression based on (het hom het)
		grp <- .groupBasedOnPhaseAndAlleleCombination(phase(rv2, return.class="array")[,,c(1, 2), drop=FALSE])
		plotGroups <- .lvaGroups(mcols(rv2)[["ref"]], mcols(rv2)[["alt"]])
		#call internal regression function	
		if(type=="lm"){
			mat <- lva.internal(x = assays(rs2)[["rs1"]], grp = grp, element = 3, 
								type=type, covariates=covariates)
		}else if(type=="nlme"){
			#covariates have not been implemented completely
			mat <- lva.internal(x = assays(rs2)[["rs1"]], grp = grp, element = 3, type=type, 
								subject=colData(rs2)[["subject.group"]], covariates=covariates)
		}

		#make txSNP specific lva test
		rs2 <- .addLva2ASEset(rs2, grp, type=type, covariates=covariates)

		#create return object
		if(return.class=="LinkVariantAlmlof"){
			sset <- SummarizedExperiment(
						assays = SimpleList(rs1=assays(rs2)[["rs1"]], lvagroup=grp), 
						colData = colData(rs2),
						rowRanges = granges(rs2))

			rownames(sset) <- rownames(rs2)
			mcols(sset)[["RiskVariantMeta"]] <- DataFrame(GR=granges(rv2), rsid=rownames(rv2))
			mcols(sset)[["RiskVariantMetaFull"]] <- rv2
			mcols(sset)[["LMCommonParam"]] <- DataFrame(mat, row.names=NULL)
			mcols(sset)[["LvaPlotGroups"]] <- DataFrame(plotGroups, row.names=NULL)

			#create an object with results
			new("LinkVariantAlmlof", sset,
				meta = list()
			)
		}
})


### -------------------------------------------------------------------------
### helpers for lva
###

#send in an array rows=SNPs, cols=sampels, 3dim=phase with maternal as el. 1 and paternal as el. 2)
#returns a matrix with three groups, het A/B =1, hom (A/A, B/B) =2 and het B/A =3
.groupBasedOnPhaseAndAlleleCombination <- function(ar){
		grp <- matrix(2, nrow=dim(ar)[1], ncol=dim(ar)[2])
		grp[(ar[,,1] == 1) & (ar[,,2] == 0)] <- 1                             	
		grp[(ar[,,1] == 0) & (ar[,,2] == 1)] <- 3
		grp
}


.lvaGroups <- function(ref, alt){
		fir <- paste(alt,"|", ref, sep="")
		sec <- paste(ref,"|", ref, " AND ",alt,"|", alt, sep="")
		thi <- paste(ref,"|", alt, sep="")

		matrix(c(fir, sec, thi), ncol=3)
}

.addLva2ASEset <- function(rs2, grp, type, covariates=matrix()){
		lst <- mcols(rs2)[["ASEsetMeta"]][[1]]
		for(i in 1:length(lst)){
		  fr <- assays(lst[[i]])[["matfreq"]]
		  grp2 <- grp[i,]
		  if(type=="lm"){ 
				  lmcomparam <- .lvaRegressionReturnCommonParamMatrixTxSNPspecific(fr,grp2, covariates)
		  }else if(type=="nlme"){
				  subj <- colData(rs2)[["subject.group"]]
				  lmcomparam <- .lvaRegressionReturnCommonParamMatrixTxSNPspecific.nlme(fr,grp2, subj)
		  }
		  mcols(mcols(rs2)[["ASEsetMeta"]][[1]][[i]])[["lmcomparam"]] <- DataFrame(lmcomparam)
		}
		rs2
}


#' lva.internal
#' 
#' make an almlof regression for arrays (internal function)
#' 
#' internal method that takes one array with results from regionSummary
#' and one matrix with group information for each risk SNP (based on phase).
#' Input and output objects can change format slightly in future.
#'
#' @name lva.internal
#' @rdname lva.internal
#' @aliases lva.internal,array-method
#' @docType methods
#' @param x regionSummary array phased for maternal allele
#' @param grp group 1-3 (1 for 0:0, 2 for 1:0 or 0:1, and 3 for 1:1)
#' @param element which column in x contains the values to use with lm.
#' @param type which column in x contains the values to use with lm.
#' @param subject which samples belongs to the same individual
#' @param covariates add data.frame with covariates (only integers and numeric)
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#' 
#' data(ASEset) 
#' a <- ASEset
#' # Add phase
#' set.seed(1)
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#' 
#' phase(a) <- p
#' 
#' #add alternative allele information
#' mcols(a)[["alt"]] <- inferAltAllele(a)
#' 
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' region <- split(granges(a)[c(1,2,2,3)], c(1,1,2,2))
#' rs <- regionSummary(a, region, return.class="array", return.meta=FALSE)
#'
#' # use  (change to generated riskSNP phase later)
#' phs <- array(c(phase(a,return.class="array")[1,,c(1, 2)], 
#'				 phase(a,return.class="array")[2,,c(1, 2)]), dim=c(20,2,2))
#' grp <- matrix(2, nrow=dim(phs)[1], ncol=dim(phs)[2])		 
#' grp[(phs[,,1] == 0) & (phs[,,2] == 0)] <- 1
#' grp[(phs[,,1] == 1) & (phs[,,2] == 1)] <- 3

#' #only use mean.fr at the moment, which is col 3
#' lva.internal(x=assays(rs)[["rs1"]],grp=grp, element=3)
#' 
NULL

#' @rdname lva.internal
#' @export
setGeneric("lva.internal", function(x, ... 
	){
    standardGeneric("lva.internal")
})

#' @rdname lva.internal
#' @export
setMethod("lva.internal", signature(x = "array"),
		function(x, grp, element=3, type="lm", subject=NULL, covariates=matrix(), ...
	){
		
		#unlist(.lvaRegressionPvalue(x, grp, element))
		#normal regression
		if(type=="lm"){
			.lvaRegressionReturnCommonParamMatrix(ar=x, grp, element, covariates=covariates)
		#mixed models lme4 regression
		} else if(type=="nlme"){
			if(!is.null(subject)){
			.lvaRegressionReturnCommonParamMatrix.nlme(ar=x, grp, subject, element)
			#.lvaRegressionReturnCommonParamMatrix.nlme(ar, grp, subject, element)
			}else{
				stop("subject cannot be null when using nlme method")
			}
		}else{
				stop("type version not specified")
		}

})

### -------------------------------------------------------------------------
### helpers for lva.internal
###

# input ar array 1d=SNP, 2d=samples, 3d=variable
# input grp matrix with group 1 2 3. 2d=samples, 1d=SNP
# lapply over each snp and make regression over variable element based on grp
# output is a list with one result for each SNP
.lvaRegressionPvalue <- function(ar, grp, element){
	lapply(1:dim(ar)[1], function(i, y, x){
				summary(lm(y[i, ,element]~x[, i]))$coefficients[2, 4]
		}, y=ar, x=grp)
}



.lvaRegressionReturnCommonParamMatrix <- function(ar, grp, element, covariates){

	mat <- matrix(NA, ncol=dim(ar)[3], nrow=nrow(ar))
	nocalc <- apply(ar[,,3, drop=FALSE], 1, function(x){sum(!(is.na(x)))==0})

	#use covariates if they exist
	cov2 <- covariates
	if(!length(covariates)==1){
	  if(!ncol(grp)==nrow(covariates)) stop("grp and cov has to be same length")
	  #cov2 <-covariates[!nocalc,drop=FALSE,] #nocalc is for SNPs not samples
	}

	#only make regression if there is at least one row possible to compute
	if(any(!nocalc)){
		mat[!nocalc,] <- t(sapply(which(!nocalc), function(i, y, x, c){
						mat2 <- matrix(NA, ncol=2, nrow=4)
						covform <- paste(colnames(c), collapse="+")
						if(!length(c)==1) form <- formula(paste("y2~x2",covform, sep="+"))
						if(length(c)==1) form <- formula("y2~x2")
						df <- cbind(c, data.frame(y2=y[i, ,element], x2=x[i, ]))
						s <-summary(lm(form, data=df))$coefficients
						mat2[,1:2] <- s[rownames(s) %in% c("(Intercept)","x2"),]
						c(mat2)
					}, y=ar[!nocalc,,,drop=FALSE], x=grp[!nocalc,,drop=FALSE], c=cov2))
	}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	mat
}

#mixed model variant lme4 (gives no straight forward p-value)
#.lvaRegressionReturnCommonParamMatrix.lme4 <- function(ar, grp, subject, element){
#
#	mat <- matrix(NA, ncol=dim(ar)[3], nrow=nrow(ar))
#	nocalc <- apply(ar[,,3, drop=FALSE], 1, function(x){sum(!(is.na(x)))==0})
#
#	#only make regression if there is at least one row possible to compute
#	if(any(!nocalc)){
#		mat[!nocalc,] <- t(sapply(which(!nocalc), function(i, y, x, s){
#						mat2 <- matrix(NA, ncol=2, nrow=4)
#						nas <- !(is.na(y[i, ,element]) | is.na(x[, i]) | is.na(s))
#						l <- lmer(y[i,nas ,element]~x[nas, i]+(1|s[nas]))
#
#						s <-l$coefficients
#						mat2[,1:nrow(s)] <- s
#						c(mat2)
#					}, y=ar[!nocalc,,,drop=FALSE], x=grp[,!nocalc,drop=FALSE], s=subject))
#	}
#	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
#	mat
#}
#
#inspired by ben bolkers answer here
#http://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
.lvaRegressionReturnCommonParamMatrix.nlme <- function(ar, grp, subject, element){

	mat <- matrix(NA, ncol=dim(ar)[3], nrow=nrow(ar))
	nocalc <- apply(ar[,,3, drop=FALSE], 1, function(x){sum(!(is.na(x)))==0})

	#only make regression if there is at least one row possible to compute
	if(any(!nocalc)){
		mat[!nocalc,] <- t(sapply(which(!nocalc), function(i, y, x, s){
						#for(i in 1:5){
							mat2 <- matrix(NA, ncol=2, nrow=4)
							nas <- (is.na(y[i, ,3]) | is.na(x[, i]) | is.na(s))
							few <- length(unique(x[!nas,i])) == 1
							few2 <- length(x[!nas,i]) == 2
							
							#zero counts (which logically should not exist, but does in one instance)
							zc <- all(y[i,!nas ,3]==0)

							#less than four we do not calculate
							few3 <- length(s[!nas]) <=5

							if(!few & !few2 & !few3 & !zc){
								df <- data.frame(res=y[i,!nas ,3], exp=x[!nas, i], ran=s[!nas])
								m1 <- lme(res~exp, random=~1|ran, data=df)
								mat2[7:8] <- anova(m1)$'p-value'
								c(mat2)
							}else{
								c(mat2)
							}
						#}
					}, y=ar[!nocalc,,,drop=FALSE], x=t(grp)[,!nocalc,drop=FALSE], s=subject))
	}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	mat
}


.lvaRegressionReturnCommonParamMatrixTxSNPspecific <- function(fr, grp, covariates=matrix()){
	fr2 <- t(fr)
	grp2 <- grp
	mat <- matrix(NA, ncol=8, nrow=ncol(fr2))
	nocalc <- apply(fr2[,, drop=FALSE], 2, function(x){sum(!(is.na(x)))==0})

	#use covariates if they exist
	cov2 <- covariates
	if(!length(covariates)==1){
	  if(!length(grp2)==nrow(covariates)) stop("grp and cov has to be same length")
	  #cov2 <-covariates[!nocalc,drop=FALSE,] #nocalc is for SNPs not samples
	}

		#y <- fr2[!nocalc,,drop=FALSE]
		#x <- grp[,!nocalc,drop=FALSE]
	    x <- grp2
		for(i in 1:ncol(fr2)){
			y <- fr2[,i]
			mat2 <- matrix(NA, ncol=2, nrow=4)
			if(!(nocalc[i])){
				covform <- paste(colnames(cov2), collapse="+")
				if(!length(cov2)==1) form <- formula(paste("y2~x2",covform, sep="+"))
				if(length(cov2)==1) form <- formula("y2~x2")
				df <- cbind(cov2, data.frame(y2=y, x2=x))
				s <-summary(lm(form, data=df))$coefficients
				mat2[,1:2] <- s[rownames(s) %in% c("(Intercept)","x2"),]
			  mat[i,] <- c(mat2)
			}
		}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")

	mat
}

.lvaRegressionReturnCommonParamMatrixTxSNPspecific.nlme <- function(fr, grp, subj){
	fr2 <- t(fr)
	grp2 <- grp
	s <- subj
	mat <- matrix(NA, ncol=8, nrow=ncol(fr2))
	nocalc <- apply(fr2[,, drop=FALSE], 2, function(x){sum(!(is.na(x)))==0})
	#y <- fr2[!nocalc,,drop=FALSE]
	#x <- grp[,!nocalc,drop=FALSE]
	x <- grp2
	if(!(any(nocalc))){
      for(i in 1:ncol(fr2)){
		y <- fr2[,i,drop=TRUE]
		mat2 <- matrix(NA, ncol=2, nrow=4)
		nas <- (is.na(y) | is.na(x) | is.na(s))
		few <- length(unique(x[!nas])) == 1
		if(!few){
		 df <- data.frame(res=y[!nas], exp=x[!nas], ran=s[!nas])
		 m1 <- lme(res~exp, random=~1|ran, data=df)
		 mat2[7:8] <- anova(m1)$'p-value'
		 mat[i,] <- c(mat2)
		}else{
		 mat[i,] <- c(mat2)
	    }
	  }
	}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")

	mat
}


