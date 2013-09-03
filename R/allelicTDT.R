allelicTDT <- function(mat.snp, size=50, correct=FALSE){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	n.snp <- ncol(mat.snp)
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- numeric(n.snp)
	for(i in 1:(length(int)-1))
		stat[int[i]:(int[i+1]-1)] <- aTDTchunk(mat.snp[,int[i]:(int[i+1]-1), drop=FALSE], correct=correct)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	out <- list(stat=stat, pval=pval)
	class(out) <- "aTDT"
	out
}

aTDTchunk <- function(geno, correct=FALSE){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	dad <- 100 * dad + 10 * mom + kid
	tmp1 <- colSums(dad==11 | dad==101, na.rm=TRUE)
	tmp2 <- colSums(dad==122 | dad==212, na.rm=TRUE)
	tmp3 <- colSums(dad==111, na.rm=TRUE)
	tmp4 <- 2 * colSums(dad==112, na.rm=TRUE)
	transMinor <- tmp1 + tmp2 + tmp3 + tmp4
	tmp1 <- colSums(dad==10 | dad==100, na.rm=TRUE)
	tmp2 <- colSums(dad==121 | dad==211, na.rm=TRUE)
	tmp4 <- 2 * colSums(dad==110, na.rm=TRUE)
	transMajor <- tmp1 + tmp2 + tmp3 + tmp4
	tmp1 <- transMinor - transMajor
	if(correct)
		tmp1 <- abs(tmp1) - 1
	tmp1 * tmp1 / (transMinor + transMajor)
}

print.aTDT <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Statistic=x$stat, "p-value"=pval, check.names=FALSE, stringsAsFactors=FALSE)
	cat("      Allelic TDT", "\n\n", sep="")
	if(!is.na(top) && top > 0 && top <= length(pval)){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top ", top, " SNPs:\n", sep="")
	}
	print(format(out, digits=digits))
}
	