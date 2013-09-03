colPOlrt <- function(mat.snp, size=20){
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
	n.snp <- ncol(mat.snp)
	int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
	full <- red <- numeric(n.snp)
	for(i in 1:(length(int) - 1)){
		tmp <- polrtChunk(mat.snp[, int[i]:(int[i+1]-1), drop=FALSE])
		full[int[i]:(int[i+1]-1)] <- tmp$full
		red[int[i]:(int[i+1]-1)] <- tmp$red
	}
	stat <- 2 * (red - full)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	names(stat) <- names(pval) <- colnames(mat.snp)
	out <- list(stat=stat, pval=pval, full=-full, red=-red)
	class(out) <- "polrt"
	out
} 

polrtChunk <- function(geno){
	matCounts <- getPOlrtCounts(geno)
	n.snp <- ncol(geno)
	matRed <-matrix(nrow=n.snp, ncol=6)
	matRed[,1:3] <- matCounts[,2:4]
	matRed[,4] <- matCounts[,5] + matCounts[,8]
	matRed[,5] <- matCounts[,6] + matCounts[,9]
	matRed[,6] <- matCounts[,7] + matCounts[,10]
	full <- red <- numeric(n.snp)
	for(i in seq(along=full)){
		param <- optim(par=c(0,0,0), fn=negPOfullLike, vec=matCounts[i,])$par
		full[i] <- negPOfullLike(param, matCounts[i,])
		param2 <- optim(par=c(0,0), fn=negPOredLike, vec=matRed[i,])$par
		red[i] <- negPOredLike(param2, matRed[i,])
	}
	list(red=red, full=full)
}

getPOlrtCounts <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	kid1 <- kid == 1
	kid0 <- kid != 1
	matCounts <- matrix(nrow=ncol(geno), ncol=10)
	matCounts[,1] <- colSums(mom > dad & kid1, na.rm=TRUE)
	matCounts[,2] <- colSums(mom == 1 & dad == 0, na.rm=TRUE)
	matCounts[,3] <- colSums(mom == 2 & dad == 0, na.rm=TRUE)
	matCounts[,4] <- colSums(mom == 2 & dad == 1, na.rm=TRUE)
	mp <- mom + dad
	matCounts[,5] <- colSums(mp == 1 & kid1, na.rm=TRUE)
	matCounts[,6] <- colSums(mp == 2 & mom != dad & kid1, na.rm=TRUE)
	matCounts[,7] <- colSums(mp == 3 & kid1, na.rm=TRUE)
	matCounts[,8] <- colSums(mp == 1 & kid0, na.rm=TRUE)
	matCounts[,9] <- colSums(mp == 2 & mom != dad & kid0, na.rm=TRUE)
	matCounts[,10] <- colSums(mp == 3 & kid0, na.rm=TRUE)
	matCounts
}

negPOfullLike <- function(beta, vec){
	-vec[1] * beta[1] - vec[2] * beta[3] - vec[3] * beta[2] - 
		vec[4] * (beta[2] - beta[3]) +
		vec[5] * log(1 + exp(beta[1] + beta[3])) + 
		vec[6] * log(1 + exp(beta[1] + beta[2])) + 
		vec[7] * log(1 + exp(beta[1] + beta[2] - beta[3])) + 
		vec[8] * log(1 + exp(beta[3])) +
		vec[9] * log(1 + exp(beta[2])) +
		vec[10] * log(1 + exp(beta[2] - beta[3]))
}

negPOredLike <- function(beta, vec){
	- vec[1] * beta[2] - vec[2] * beta[1] - vec[3] * (beta[1] - beta[2]) +
		vec[4] * log(1 + exp(beta[2])) + 
		vec[5] * log(1 + exp(beta[1])) + 
		vec[6] * log(1 + exp(beta[1] - beta[2]))
}

print.polrt <- function(x, top = 5, digits = 4, ...){
	pval <- format.pval(x$pval, digits = digits)
	out <- data.frame("LL (with)" = x$full, "LL (without)"=x$red, Statistic=x$stat, "P-Value"=pval,
		check.names=FALSE)
	cat("         Parent-Of-Origin Likelihood Ratio Test\n\n", sep="")
	if(!is.na(top) && top > 0 && top <= length(pval)){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top ", top, " SNPs:\n", sep="")
	}
	print(format(out, digits=digits))
}

 
	





	