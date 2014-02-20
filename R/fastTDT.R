logit <- function(x) log(x/(1-x))

fastTDT <- function(mat.geno, type, size=50){
	if(size < 1)
		stop("size should be at least 1.")
	fun <- match.fun(switch(type, "additive"=fastTDTsplit, "dominant"=fastTDTdomSplit, 
		"recessive"=fastTDTrecSplit))
	fun(mat.geno, size=size)
}



fastTDTsplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
	num <- denom <- used <- numeric(n.snp)
	for(i in 1:(length(int)-1)){
		tmp <- fastTDTchunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		num[int[i]:(int[i+1]-1)] <- tmp$num
		denom[int[i]:(int[i+1]-1)] <- tmp$denom
		used[int[i]:(int[i+1]-1)] <- tmp$used
	}
	beta <- logit(num/denom)
	se <- sqrt(denom / ((denom-num)*num))
	stat <- beta/se
	stat <- stat*stat
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(is.null(colnames(geno)))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:ncol(geno), sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- colnames(geno)
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="additive", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}
	

fastTDTchunk <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),,drop=FALSE]
 	mom <- geno[seq.int(2, n.row, 3),,drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),,drop=FALSE]
	het1 <- (dad==1) & (mom==1)
	a567 <- colSums(het1 & !is.na(kid), na.rm=TRUE)
	a67 <- colSums(het1 * kid, na.rm=TRUE)
	mom <- mom + dad
	het1 <- mom == 1
	a1 <- colSums(het1 & kid == 0, na.rm=TRUE)
	a2 <- colSums(het1 & kid == 1, na.rm=TRUE)
	het1 <- mom == 3
	a3 <- colSums(het1 & kid == 1, na.rm=TRUE)
	a4 <- colSums(het1 & kid == 2, na.rm=TRUE)
	num <- a2 + a4 + a67
	used <- a1 + a2 + a3 + a4 + a567
	# denom <- a1 + a2 + a3 + a4 + 2 * a567
	return(list(num=num, denom=used+a567, used=used))
}	
	

fastTDTdomSplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	dmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		dmat[int[i]:(int[i+1]-1),] <- fastTDTdomChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
	rownames(dmat) <- colnames(geno)
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
	if(is.null(colnames(geno)))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:ncol(geno), sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- colnames(geno)
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="dominant", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}


fastTDTdomChunk <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	dmat <- matrix(0, ncol(dad), 4)
	het1 <- mom+dad == 1
	dmat[,1] <- colSums(het1 & kid==0, na.rm=TRUE)
	dmat[,2] <- colSums(het1 & kid==1, na.rm=TRUE)
	het1 <- (dad==1) & (mom==1)
	dmat[,3] <- colSums(het1 & kid==0, na.rm=TRUE)
	dmat[,4] <- colSums(het1 & kid!=0, na.rm=TRUE)
	dmat
}


fastTDTrecSplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	rmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		rmat[int[i]:(int[i+1]-1),] <- fastTDTrecChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
	rownames(rmat) <- colnames(geno)
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
	if(is.null(colnames(geno)))
		names(beta) <- names(stat) <- names(pval) <- names(used) <- paste("SNP", 1:ncol(geno), sep="")
	else
		names(beta) <- names(stat) <- names(pval) <- names(used) <- colnames(geno)
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		ia=FALSE, type="recessive", usedTrios=used, add=FALSE)
	class(out) <- "colTDT"
	out
}

fastTDTrecChunk <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	rmat <- matrix(0, ncol(dad), 4)
	het3 <- mom + dad == 3
	rmat[,1] <- colSums(het3 & kid==1, na.rm=TRUE)
	rmat[,2] <- colSums(het3 & kid==2, na.rm=TRUE)
	het3 <- (dad==1) & (mom==1)
	rmat[,3] <- colSums(het3 & kid!=2, na.rm=TRUE)
	rmat[,4] <- colSums(het3 & kid==2, na.rm=TRUE)
	rmat
}

	





