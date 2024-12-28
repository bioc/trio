# Removed fastTDTsplit, fastTDTdomSplit, fastTDTrecSplit.
# Updated version is in colTDTnew.R


logit <- function(x) log(x/(1-x))

fastTDT <- function(mat.geno, type, size=50){
	if(size < 1)
		stop("size should be at least 1.")
	fun <- match.fun(switch(type, "additive"=fastTDTsplit, "dominant"=fastTDTdomSplit, 
		"recessive"=fastTDTrecSplit))
	fun(mat.geno, size=size)
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

