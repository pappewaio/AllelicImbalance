#' detectAI
#' 
#' detection of AllelicImbalance 
#' 
#' threshold.frequency is the least fraction needed to classify as bi tri or
#' quad allelic SNPs. If 'all' then all of bi tri and quad allelic SNPs will use the same
#' threshold. Everything under the treshold will be regarded as noise. 'all' will return 
#' a matrix with snps as rows and uni bi tri and quad will be columns. For this function
#' Anything that will return TRUE for tri-allelicwill also return TRUE for uni and bi-allelic
#' for the same SNP an Sample. 
#'
#' return.type 'ref' return only AI when reference allele is more expressed. 'alt' return only
#' AI when alternative allele is more expressed or 'all' for both 'ref' and 'alt' alleles.
#' Reference allele is the one present in the reference genome on the forward strand.
#'
#' threshold.delta.frequency and function.test will use the value in mapBias(x) as expected value. 
#' 
#' function.test will use the two most expressed alleles for testing. Make therefore sure there
#' are no tri-allelic SNPs or somatic mutations among the SNPs in the ASEset. 
#' 
#' inferGenotype(), set TRUE it should be used with as much samples as possible. If you split up the
#' samples and run detectAI() on each sample separately, please make sure you have inferred 
#' the genotypes in before hand, alternatively used the genotypes detected by another variantCaller
#' or chip-genotypes. Use ONLY biallelic genotypes.
#'
#' @name detectAI
#' @rdname detectAI
#' @aliases detectAI
#' detectAI,ASEset-method 
#' @docType methods
#' @param x ASEset
#' @param strand strand to infer from
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @param return.class class to return (atm only class 'logical')
#' @param threshold.delta.frequency minimum of frequency difference 
#' from 0.5 (or mapbias adjusted value)
#' @param threshold.pvalue pvalue over this number will be filtered out
#' @param function.test At the moment the only available option is 'binomial.test'
#' @param inferGenotype infer genotypes based on count data in ASEset object
#' @param random.ref set the reference as random if you dont know. Affects interpretation of results.
#' @param verbose makes function more talkative
#' @param gc use garbage collection when possible to save space
#' @param biasMatrix use biasMatrix in ASEset, or use default expected frequency of 0.5 for all sites
#' @param ... internal arguments 
#' @author Jesper R. Gadin
#' @keywords detection
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' dai <- detectAI(a)
#' 
#'
#' @export
setGeneric("detectAI", function(x, ...){
    standardGeneric("detectAI")
})

#' @rdname detectAI
#' @export
setMethod("detectAI", signature(x = "ASEset"), function(x, 
	return.class = "DetectedAI", strand = "*",
	threshold.frequency=0, threshold.count.sample=1,
	threshold.delta.frequency=0, threshold.pvalue=0.05,
	inferGenotype=FALSE,
	random.ref=FALSE,
	function.test="binom.test",
	verbose=TRUE,
	gc=FALSE, biasMatrix=FALSE) {

	if(inferGenotype){
		genotype(x) <- inferGenotypes(x,threshold.frequency = 0.05, return.allele.allowed ="bi")
	}
	#check for presence of genotype data
    if (!("genotype" %in% assayNames(x))) {
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' ",
				   "or set 'inferGenotypes=TRUE'"))
    }

	if(!random.ref){
		if(!("ref" %in% colnames(mcols(x)))){
			stop("column name 'ref' in mcols(x) is required")
		}
	}else{
		mcols(x)[,"ref"] <- randomRef(x)
	}

	#1) make refFreq array (third dim has length 1 and softest condition)
	if(verbose){cat("calculating reference fractions\n")}
	fr <- fraction(x, strand=strand,
				top.fraction.criteria="ref",
				threshold.count.sample = 1)

	#2) make t.c.s array
	if(verbose){cat("checking count thresholds\n")}
	acounts <- alleleCounts(x, strand=strand, return.class="array")
	#acounts[array(!hetFilt(x), dim=c(nrow(x), ncol(x), 4))] <- 0  # unnecessary het is checked from fraction later
	allele.count.tot <- apply(acounts, c(1,2), sum)
	t.c.s  <- !array(allele.count.tot,dim=c(nrow(x),ncol(x),length(threshold.count.sample))) < 
					aperm(array(threshold.count.sample,
						dim=c(length(threshold.count.sample),nrow(x),ncol(x) )),
					c(2,3,1))
	dimnames(t.c.s) <- list(rownames(fr),colnames(fr),paste(threshold.count.sample))

	if(gc){
		if(verbose){cat("removing and garbage collect temporary variables\n")}
		rm(acounts,allele.count.tot)
		gc()
	}

	#3) make thr.req array (if TRUE it fullfills the condition)
	if(verbose){cat("checking frequency thresholds\n")}
	newDims <- length(threshold.frequency)
	thr.fr.freq <- array(fr,dim=c(nrow(x),ncol(x),newDims))
	thr.freq <- aperm(array(threshold.frequency, dim=c(newDims,nrow(x),ncol(x))),c(2,3,1))

	thr.freq.ret <- !(thr.fr.freq < thr.freq | thr.fr.freq > (1-thr.freq))
	dimnames(thr.freq.ret) <- list(rownames(fr),colnames(fr),paste(threshold.frequency))

	if(gc){
		if(verbose){cat("removing and garbage collect temporary variables\n")}
		rm(newDims, thr.fr.freq, thr.freq)
		gc()
	}

	#4) delta freq array 
	if(verbose){cat("checking delta frequency thresholds\n")}

	if(biasMatrix){
		biasmatRef <- mapBiasRef(x)
	}else{
		biasmatRef <- matrix(0.5, ncol=ncol(x), nrow=nrow(x))
	}

	fr2 <- fr
	fr2[is.na(fr2)] <- biasmatRef[is.na(fr2)] 

	#make fr2 and biasmatRef to array dim 3
	newDims <- length(threshold.delta.frequency)

	t.d.f <- aperm(array(threshold.delta.frequency, dim=c(newDims,ncol(x),nrow(x))),dim=c(3,2,1))

	#to keep1 (fullfills cond min.delta.freq)
	#if(any(fr2<biasmatRef)){
	#	#fr3 <- fr2
	#	#fr3[fr2<biasmatRef] <- biasmatRef[fr2<biasmatRef] - fr2[fr2<biasmatRef]
	#	#fr3[fr2<=biasmatRef] <- 1	
	#	#reset the NA again
	#	#fr3[is.na(fr)] <- NA
	#	#make 3d array return object 
	#	#fr3 <- array(fr3,dim=c(nrow(x),ncol(x),newDims))
	#	#tf.keep1 <- !(fr3 < t.d.f)
	#	tf.not.keep1 <- !array(abs(fr2-biasmatRef),dim=c(nrow(x),ncol(x),newDims)) > t.d.f

	#}else{
	#	#set all FALSE (NAs will not be kept )
	#	tf.keep1 <- matrix(FALSE,ncol=ncol(fr),nrow=nrow(fr))
	#	#tf.keep1[is.na(fr)] <- NA
	#	#make 3d array return object 
	#	tf.not.keep1 <- !array(tf.keep1,dim=c(nrow(x),ncol(x),newDims))
	#}

	##to keep2 (fullfills cond min.delta.freq)
	#if(any(fr2>biasmatRef)){
	#	#fr3 <- fr2
	#	##fr3[fr2>biasmatRef] <- fr[fr2>biasmatRef] - biasmatRef[fr2>biasmatRef] 
	#	#fr3[fr2>biasmatRef] <- fr[fr2>biasmatRef] - biasmatRef[fr2>biasmatRef] 
	#	#fr3[fr2<=biasmatRef] <- 1	
	#	##fr3[is.na(fr)] <- NA
	#	##make 3d array return object 
	#	#fr3 <- array(fr3,dim=c(nrow(x),ncol(x),newDims))
	#	#tf.keep2 <- !(fr3 < t.d.f)
	#	tf.not.keep1 <- !array((fr2-biasmatRef),dim=c(nrow(x),ncol(x),newDims)) > t.d.f
	#}else{
	#	#set all FALSE (NAs will not be kept )
	#	tf.keep2 <- matrix(TRUE,ncol=ncol(fr),nrow=nrow(fr))
	#	#tf.keep2[is.na(fr)] <- NA
	#	#make 3d array return object 
	#	tf.not.keep2 <- !array(tf.keep2,dim=c(nrow(x),ncol(x),newDims))
	#}
	
	tf.keep <- array(abs(fr2-biasmatRef),dim=c(nrow(x),ncol(x),newDims)) >= t.d.f

	#delta.freq <- tf.not.keep1 | tf.not.keep2
	delta.freq <- tf.keep
	dimnames(delta.freq) <- list(rownames(fr),colnames(fr),1:length(threshold.delta.frequency))
	#delta.freq[is.na(fr)] <- NA

	if(gc){
		if(verbose){cat("removing and garbage collect temporary variables\n")}
		rm(fr2, biasmatRef, tf.keep, newDims)
		gc()
	}

	#5)  p-value array
	if(verbose){cat("checking p-value thresholds\n")}
	if(function.test=="binom.test"){
		idx <- which(!is.na(fr))
		arn <- arank(x,return.type="names",return.class="matrix") 
		ac <- alleleCounts(x, strand=strand, return.class="array")

		mat1 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		mat2 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))


		allele1 <- mat1[idx] 
		allele2 <- mat2[idx]

		#test the two most expressed alleles
		biasmat1 <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
			x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))


		biasAllele1 <- biasmat1[idx] 
		
		#maybe put a check here later that bias1 and bias2 sum to 1

		ml <- mapply(allele1,allele2,biasAllele1,FUN=function(x,y,z){
			binom.test(c(x,y), p = z)[[3]]
		})

		#cerate matrix in same size as fr	
		pv <- fr	
		pv[idx] <- ml
		#replace na with 1
		pv[is.na(fr)] <- 1
		

		#make pv array
		pv <- array(pv,dim=c(nrow(fr), ncol(fr),length(threshold.pvalue)),
					dimnames=list(rownames(fr),colnames(fr),paste(threshold.pvalue)))

		#filter
		newDims <- length(threshold.pvalue)
		pv.thr <- aperm(array(threshold.pvalue, dim=c(newDims,nrow(x),ncol(x))), c(2,3,1))
		pv.thr.ret <- pv <= pv.thr
	
		if(gc){
			if(verbose){cat("removing and garbage collect temporary variables\n")}
			rm(idx,arn,ac,mat2,mat1,allele1,allele2,biasmat1,biasAllele1,ml,pv,newDims,pv.thr)
			gc()
		}

	}else{stop("function.test must be binom.test")}

	if(return.class=="DetectedAI"){

		#make DetectedAI object
		if(verbose){cat("creating DetectedAI object\n")}
		DetectedAIFromArray(
			x, 
			strand=strand,
			reference.frequency=fr,
			threshold.frequency=thr.freq.ret,
			threshold.count.sample=t.c.s,
			threshold.delta.frequency=delta.freq,
			threshold.pvalue=pv.thr.ret,

			#reference.frequency.names=reference.frequency,
			threshold.frequency.names=paste(threshold.frequency),
			threshold.count.sample.names=paste(threshold.count.sample),
			threshold.delta.frequency.names=paste(threshold.delta.frequency),
			threshold.pvalue.names=paste(threshold.pvalue)

		)

	}else if(return.class=="list"){
		stop("return.class list as option is not available atm")
	}else{
		stop(paste("return.class",return.class,"as option is not available"))
	}

})

