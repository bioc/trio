colTAT <- function(mat.snp, size=50, bothHet=0){
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
	for(i in 1:(length(int) - 1)){
		tmp <- tatChunk(mat.snp[, int[i]:(int[i+1]-1), drop=FALSE], bothHet=bothHet)
		stat[int[i]:(int[i+1]-1)] <- tmp$stat
		ntrio[int[i]:(int[i+1]-1)] <- tmp$n
	}
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(is.null(colnames(mat.snp)))
		names(stat) <- names(pval) <- paste("SNP", 1:n.snp, sep="")
	else
		names(stat) <- names(pval) <- colnames(mat.snp)
	out <- list(stat=stat, pval=pval, usedTrios=ntrio)
	class(out) <- "tat"
	out
}

tatChunk <- function(geno, bothHet=0){
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
	list(stat=stat, n=n)
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
	cat("   Transmission Asymmetry Test\n\n", sep="")
	if(!is.na(top) && top > 0 && top <= length(pval)){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top ", top, " SNPs:\n", sep="")
	}
	print(format(out, digits=digits))
}

   
	 	 