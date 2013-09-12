colTAT <- function(mat.snp, stratified=FALSE, size=50, bothHet=0){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0,1,2,NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	if(bothHet>1 | bothHet<0)
		stop("bothHet must be between 0 and 1.")
	n.snp <- ncol(mat.snp)
	int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
	stat <- ntrio <- numeric(n.snp)
	if(stratified){
		sp <- sm <- numeric(n.snp)
		matStrat <- matrix(nrow=n.snp, ncol=4)
	}
	for(i in 1:(length(int) - 1)){
		tmp <- tatChunk(mat.snp[, int[i]:(int[i+1]-1), drop=FALSE], bothHet=bothHet, 
			stratified=stratified)
		stat[int[i]:(int[i+1]-1)] <- tmp$stat
		ntrio[int[i]:(int[i+1]-1)] <- tmp$n
		if(stratified){
			sp[int[i]:(int[i+1]-1)] <- tmp$sp
			sm[int[i]:(int[i+1]-1)] <- tmp$sm
			print(i)
			print(dim(matStrat[int[i]:(int[i+1]-1),,drop=FALSE]))
			print(dim(tmp$matObs))
			print(tmp$matObs)
			matStrat[int[i]:(int[i+1]-1),] <- tmp$matObs 
		}
	}
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(is.null(colnames(mat.snp)))
		names(stat) <- names(pval) <- paste("SNP", 1:n.snp, sep="")
	else
		names(stat) <- names(pval) <- colnames(mat.snp)
	if(!stratified)
		out <- list(stat=stat, pval=pval, usedTrios=ntrio, strat=FALSE)
	else{
		pp <- pchisq(sp, 1, lower.tail=FALSE)
		pm <- pchisq(sm, 1, lower.tail=FALSE)
		names(sp) <- names(sm) <- names(pp) <- names(pm) <- names(stat)
		rownames(matStrat) <- names(stat)
		colnames(matStrat) <- c("PaternalTrans", "PaternalNonT", "MaternalTrans", "MaternalNonT")
		out <- list(stat=stat, pval=pval, usedTrios=ntrio, matStrat=matStrat, statPaternal=sp,
			pvalPaternal=pp, statMaternal=sm, pvalMaternal=pm, strat=TRUE)
	}
	class(out) <- "tat"
	out
}

tatChunk <- function(geno, bothHet=0, stratified=FALSE){
	matObs <- getTATnumbers(geno, bothHet=bothHet)
	ttp <- rowSums(matObs[,c(1,2)])
	ttm <- rowSums(matObs[,c(3,4)])
	tpm <- rowSums(matObs[,c(1,3)])
	ntpm <- rowSums(matObs[,c(2,4)])
	n <- ttp + ttm
	matExp <- matrix(nrow = ncol(geno), ncol=4)
	matExp[,1] <- ttp * tpm 
	matExp[,2] <- ttp * ntpm
	matExp[,3] <- ttm * tpm 
	matExp[,4] <- ttm * ntpm
	matExp <- matExp / n 
	stat <- rowSums(matObs * matObs / matExp) - n
	if(!stratified)
		return(list(stat=stat, n=n))
	sp <- 2 * rowSums(matObs[,1:2] * matObs[,1:2]) / ttp - ttp
	sm <- 2 * rowSums(matObs[,3:4] * matObs[,3:4]) / ttm - ttm
	list(stat=stat, n=n, sp=sp, sm=sm, matObs=matObs)	
}
	
getTATnumbers <- function(geno, bothHet=0){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	het <- (dad == 1)
	hethom <- het & (mom == 2)
	n212 <- colSums(hethom & (kid == 2), na.rm=TRUE)
	n211 <- colSums(hethom & (kid == 1), na.rm=TRUE)
	hethom <- het & (mom == 0)
	n011 <- colSums(hethom & (kid == 1), na.rm=TRUE)
	n010 <- colSums(hethom & (kid == 0), na.rm=TRUE)
	het <- (mom == 1)
	hethom <- het & (dad == 2)
	n122 <- colSums(hethom & (kid == 2), na.rm=TRUE)
	n121 <- colSums(hethom & (kid == 1), na.rm=TRUE)
	hethom <- het & (dad == 0)
	n101 <- colSums(hethom & (kid == 1), na.rm=TRUE)
	n100 <- colSums(hethom & (kid == 0), na.rm=TRUE)
	matObs <- matrix(nrow = ncol(geno), ncol = 4)
	matObs[,1] <- n212 + n011   # tp
	matObs[,2] <- n211 + n010   # ntp
	matObs[,3] <- n122 + n101   # tm
	matObs[,4] <- n121 + n100   # ntm
	if(bothHet > 0){
		het <- (mom == 1) & (dad == 1)
		n112 <- colSums(het & (kid == 2), na.rm=TRUE)
		matObs[,c(1,3)] <- matObs[,c(1,3)] + bothHet * n112
		n110 <- colSums(het & (kid == 0), na.rm=TRUE)
		matObs[,c(2,4)] <- matObs[,c(2,4)] + bothHet * n110
	}
	matObs
}

print.tat <- function(x, top = 5, digits = 4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Statistic=x$stat, "p-value"=pval, Trios=x$usedTrios,
		check.names=FALSE, stringsAsFactors=FALSE)
	if(x$strat){
		pp <- format.pval(x$pvalParental, digits=digits)
		pm <- format.pval(x$pvalMaternal, digits=digits)
		out <- data.frame(out, x$matStrat[,1:2], 
			"Paternal p-value"=x$pvalPaternal,
			x$matStrat[,3:4], 
			"Maternal p-value"=x$pvalMaternal,
			check.names=FALSE, stringsAsFactors=FALSE)
		space <- "               "
	}
	else
		space <- "   "
	cat(space, "Transmission Asymmetry Test\n\n", sep="")
	if(!is.na(top) && top > 0 && top <= length(pval)){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top ", top, " SNPs:\n", sep="")
	}
	print(format(out, digits=digits))
}

   
	 	 