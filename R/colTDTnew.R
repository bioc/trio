# New version of colTDT in which now also the output of colNtrios can be used.

# New argument: 
# matNumber: Matrix containing the numbers of trios for the different genotypes, i.e. the output
#	of colNtrios (or maybe in future also of similar functions)
 

colTDT <- function(mat.snp, matNumber = NULL, model = c("additive", "dominant", "recessive"), 
		size = 50){
	if(missing(mat.snp) & is.null(matNumber))
		stop("Either mat.snp or matNumber must be specified.")
	if(!missing(mat.snp) & !is.null(matNumber))
		stop("Only one of mat.snp and matNumber can be specified.")
	type <- match.arg(model)
	if(!is.null(matNumber)){
		out.ntrios <- ntriosTDT(matNumber, type)
		return(out.ntrios)
	}
	checkMatSNP(mat.snp, size = size)
	fastTDT(mat.snp, type, size = size)
}

ntriosTDT <- function(matNumber, type){
	checkMatN(matNumber)
	fun <- match.fun(switch(type, "additive" = fastTDTsplit, "dominant" = fastTDTdomSplit, 
		"recessive" = fastTDTrecSplit))
	fun(matNumber = matNumber)
}


fastTDTsplit <- function(geno, matNumber = NULL, size = 50){
	if(!is.null(matNumber)){
		num <- rowSums(matNumber[ , c(2, 4, 6, 7, 7)])
		used <- rowSums(matNumber[ , 1:7])
		denom <- used + rowSums(matNumber[ , c(5:7)])
		cn <- rownames(matNumber)
		n.snp <- nrow(matNumber)
	}
	else{
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
		num <- denom <- used <- numeric(n.snp)
		for(i in 1:(length(int)-1)){
			tmp <- fastTDTchunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
			num[int[i]:(int[i+1]-1)] <- tmp$num
			denom[int[i]:(int[i+1]-1)] <- tmp$denom
			used[int[i]:(int[i+1]-1)] <- tmp$used
		}
		cn <- colnames(geno)
	}
	beta <- logit(num/denom)
	se <- sqrt(denom / ((denom-num)*num))
	stat <- beta/se
	stat <- stat*stat
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(is.null(cn))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:n.snp, sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- cn
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="additive", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}


fastTDTdomSplit <- function(geno, matNumber = NULL, size = 50){
	if(!is.null(matNumber)){
		dmat <- ntrios2Dom(matNumber, check = FALSE, quiet = TRUE)[ , 1:4]
		n.snp <- nrow(matNumber)
		cn <- rownames(matNumber)
	}
	else{
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
		dmat <- matrix(0, n.snp, 4)
		for(i in 1:(length(int)-1))
			dmat[int[i]:(int[i+1]-1),] <- fastTDTdomChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		cn <- colnames(geno)
	}
	rownames(dmat) <- cn
	h <- (dmat[,1]/3 - dmat[,2] + dmat[,3] - dmat[,4]/3) / (2*(dmat[,1]+dmat[,3]))
	tmp <- (dmat[,2]+dmat[,4])/(3*(dmat[,1]+dmat[,3])) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	tmp <- (dmat[,2]+dmat[,1])*or/(or+1)^2 + (dmat[,3]+dmat[,4])*or/(3*(or+1/3)^2)
	se <- sqrt(1/tmp)
	stat <- beta/se
	stat <- stat * stat
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	used <- rowSums(dmat)
	if(is.null(cn))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:n.snp, sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- cn
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="dominant", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}


fastTDTrecSplit <- function(geno, matNumber = NULL, size = 50){
	if(!is.null(matNumber)){
		rmat <- ntrios2Rec(matNumber, check = FALSE, quiet = TRUE)[ , 1:4]
		n.snp <- nrow(matNumber)
		cn <- rownames(matNumber)
	}
	else{
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
		rmat <- matrix(0, n.snp, 4)
		for(i in 1:(length(int)-1))
			rmat[int[i]:(int[i+1]-1),] <- fastTDTrecChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		cn <- colnames(geno)
	}
	rownames(rmat) <- cn
	h <- (3*rmat[,1] - rmat[,2] + rmat[,3] - 3*rmat[,4]) / (2 * (rmat[,1]+rmat[,3]))
	tmp <- 3 * (rmat[,2] + rmat[,4]) / (rmat[,1] + rmat[,3]) + h*h
	or <- sqrt(tmp) - h		 
	beta <- log(or)
	tmp <- (rmat[,1] + rmat[,2]) * or/(or+1)^2 + 3 *(rmat[,3] + rmat[,4]) * or / (or+3)^2
	se <- sqrt(1/tmp)
	stat <- (beta/se)^2
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	used <- rowSums(rmat)
	if(is.null(cn))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:n.snp, sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- cn
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="recessive", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}	



