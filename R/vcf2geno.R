vcf2geno <- function(vcf, ped, none="0/0", one=c("0/1"), both="1/1", use.rownames=FALSE, allowDifference=FALSE,
		removeMonomorphic=TRUE, changeMinor=FALSE){
	if(!is(vcf, "CollapsedVCF") & !is.matrix(vcf))
		stop("vcf must be either a matrix or an object of class collapsedVCF")
	if(is(vcf, "CollapsedVCF")){
		require(VariantAnnotation) || stop("The package VariantAnnotation is required.")
		vcf <- geno(vcf)$GT
		if(is.null(vcf))
			stop("vcf does not seem to contain the genotype data.")
	}
	tmp <- c(none, one, both)
	if(length(tmp) != length(unique(tmp)))
		stop("The genotype codings specified by none, one, and both must be unique.")
	if(!is.data.frame(ped))
		stop("ped must be a data frame.")
	if(!all(c("famid", "pid", "fatid", "motid") %in% colnames(ped)))
		stop("ped must contain columns called famid, pid, fatid, and motid comprising\n",
			"the family ID, the personal ID as well as the IDs of the father and the mother.")
	cn <- colnames(vcf)
	if(is.null(cn))
		stop("The genotype matrix in vcf must have column names specifying the IDs for the samples.")
	ids.kid1 <- ped$fatid != 0
	ids.kid2 <- ped$motid != 0
	if(any(ids.kid1 != ids.kid2))
		stop("fatid and motid must both be either zero or non-zero.")
	if(use.rownames){
		cnped <- rownames(ped)
		if(is.null(cnped))
			stop("ped must have rownames if use.rownames = TRUE.")
		if(length(cnped) != unique(length(cnped)))
			stop("The rownames of ped must be unique.")
	}
	else{
		if(any(duplicated(ped$pid))){
			ped$pid <- paste(ped$famid, ped$pid, sep="_")
			ped$fatid[ids.kid1] <- paste(ped$famid[ids.kid1], ped$fatid[ids.kid1], sep="_")
			ped$motid[ids.kid2] <- paste(ped$famid[ids.kid2], ped$motid[ids.kid2], sep="_")
			warning("Since the individual IDs in pid are not unique,\n", 
				"they are made unique by combining famid with pid.")
		}
		cnped <- ped$pid
		if(any(duplicated(cnped)))
			stop("Even after combining famid and pid, the individual IDs are not unique.\n",
				"So either there exists at least one subject being more than once in ped,\n",
				"or there exist at least two subjects with the same famid and pid.") 
	}
	if(mean(cn %in% cnped) < 0.1 | mean(cnped %in% cn) < 0.1)
		stop("Less than 10% of the samples in vcf are also in ped, or vice versa.")
	if(!all(cn %in% cnped)){
		if(allowDifference){
			vcf <- vcf[ , cn %in% cnped]
			cn <- cn[cn %in% cnped]
			warning("For some subjects in vcf, no information is available in ped.\n",
				"These subjects are removed from vcf.")
		}
		else
			stop("All subjects in vcf must also appear (with the same ID) in ped.")
	}
	if(!all(cnped %in% cn)){
		if(allowDifference){
			ped <- ped[cnped %in% cn, ]
			cnped <- cnped[cnped %in% cn]
			warning("For some subjects in ped, no genotypes are available in vcf.\n",
				"These subjects are removed from ped.")
		}
		else
			stop("All subjects in ped must also appear (with the same ID) in vcf.")
	}
	nr <- min(nrow(vcf), 100)
	if(!any(vcf =="0/0"))
		stop("None of the genotypes in the first ", nr, " rows of the genotype matrix\n",
			"contains a genotype coded by none.")
	geno <- matrix(NA, nrow=nrow(vcf), ncol=ncol(vcf), dimnames=dimnames(vcf))
	if(length(none) == 1)
		geno[vcf == none] <- 0
	else
		geno[vcf %in% none] <- 0
	if(length(one) == 1)
		geno[vcf == one] <- 1
	else
		geno[vcf %in% one] <- 1
	if(length(both) == 1)
		geno[vcf == both] <- 2
	else
		geno[vcf %in% both] <- 2
	mat.kid <- as.matrix(ped[ids.kid1, c("fatid", "motid", "pid")])
	vec.ids <- as.vector(t(mat.kid))
	if(use.rownames){
		if(all(vec.ids == ped$pid))
			vec.ids <- rownames(ped)
		else{
			m <- match(vec.ids, ped$pid)
			vec.ids <- rownames(ped)[m]
		}
	
	}
	#return(list(geno=geno, vec=vec.ids, vcf=vcf, ped=ped))		
	mat.trio <- t(geno[,vec.ids])
	if(is.null(colnames(mat.trio)))
		colnames(mat.trio) <- paste("SNV", 1:ncol(mat.trio), sep="")
	if(changeMinor){
		idsMAF <- colMAFtrio(mat.trio) > 0.5
		mat.trio[,idsMAF] <- 2 - mat.trio[,idsMAF]
	}
	if(removeMonomorphic){
		idsMono <- colMeans(mat.trio, na.rm=TRUE) == 0
		mat.trio <- mat.trio[,!idsMono]
		warning("Since ", sum(idsMono), " of the SNVs were monomorphic, these SNVs were removed.")
	}
	mat.trio
}
	
