colTDTmaxStat <- function(geno, size=50){
	if(!is.matrix(geno))
		stop("geno must be a matrix.")
	if(nrow(geno) %%3 != 0)
		stop("geno does not seem to contain trio data, as is number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(geno)))
		stop("geno does not seem to be a matrix in genotype format,\n",
			"as the row names are missing.")
	if(any(!geno %in% c(0,1,2,NA)))
		stop("The values in geno must be 0, 1, or 2.")
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	if(is.null(colnames(geno)))
		rn <- paste("SNP", 1:n.snp, sep="")
	else
		rn <- colnames(geno)
	matStats <- matrix(0, n.snp, 8, dimnames=list(rn, 
		c("a12","a34","a567","max","denom","Additive","Dominant","Recessive")))
	for(i in 1:(length(int)-1))
		matStats[int[i]:(int[i+1]-1),] <- getStatsChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
	out <- list(matA=matStats[,1:3], stat=matStats[,4], denom=matStats[,5], mat.stat=matStats[,6:8])
	class(out) <- "maxStatTrio"
	out
}

print.maxStatTrio <- function(x, top=5, digits=4, ...){
	out <- data.frame("Max-Statistic"=x$stat, x$mat.stat, check.names=FALSE)
	cat("          Maximum Genotypic TDT Statistic\n\n")
	if(length(x$stat) > top){
		ord <- order(x$stat, decreasing=TRUE)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNPs:\n")
	}
	print(format(out, digits=digits))
}
	

getStatsChunk <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	tmp.mat <- matrix(0, ncol(dad), 8)
	het1 <- (dad==1) & (mom==1)
	tmp.mat[,3:4] <- colSums(het1 & !is.na(kid), na.rm=TRUE)
	a5 <- colSums(het1 & kid==0, na.rm=TRUE)
	a6 <- colSums(het1 & kid==1, na.rm=TRUE)
	a7 <- colSums(het1 & kid==2, na.rm=TRUE)
	mom <- mom + dad
	het1 <- mom==1
	a1 <- colSums(het1 & kid == 0, na.rm=TRUE)
	a2 <- colSums(het1 & kid == 1, na.rm=TRUE)
	tmp.mat[,1] <- a1 + a2
	het1 <- mom==3
	a3 <- colSums(het1 & kid == 1, na.rm=TRUE)
	a4 <- colSums(het1 & kid == 2, na.rm=TRUE)
	tmp.mat[,2] <- a3+a4
	h <- a2 + a4 + a6 + 2*a7
	tmp <- rowSums(tmp.mat[,1:4, drop=FALSE])
	beta <- logit(h/tmp)
	v <- tmp / ((tmp - h) * h)
	add <- beta * beta / v
	tmp.mat[,5] <- tmp
	d <- a6 + a7
	h <- (1/3*a1 - a2 + a5 - 1/3*d) / (2 * (a1+a5))
	tmp <- (a2 + d) / (3 * (a1+a5)) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	tmp <- tmp.mat[,1] * or/(or+1)^2 + tmp.mat[,3] * or/(3*(or+1/3)^2)
	v <- 1/tmp
	dom <- beta * beta / v
	d <- a5 +a6
	h <- (3*a3 - a4 + d - 3*a7) / (2*(a3+d))
	tmp <- 3 * (a4+a7) / (a3+d) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	tmp <- tmp.mat[,2] * or/(or+1)^2 + 3 * tmp.mat[,3] * or/(or+3)^2
	v <- 1/tmp
	rec <- beta * beta / v
	tmp.mat[,4] <- pmax(add,dom,rec, na.rm=TRUE)
	tmp.mat[,6] <- add
	tmp.mat[,7] <- dom
	tmp.mat[,8] <- rec	
	tmp.mat
}


colTDTmaxTest <- function(geno, perm=10000, size=50, chunk=10000, minimum=0.001, verbose=FALSE){
	if(verbose)
		cat("Starting with computing MAX-Statistics at ", date(), ".\n")
	out <- colTDTmaxStat(geno, size=size)
	if(verbose)
		cat("Finished at ", date(), ".\n", 
			"Starting with computing permutation-based p-values.\n", sep="")
	matA <- out$matA
	stat <- out$stat
	denom <- out$denom
	mat.stat <- out$mat.stat
	rm(out)
	pval <- rep.int(NA, length(stat))
	vec.perm <- diff(unique(c(seq(1, perm, chunk), perm+1)))
	le.perm <- length(vec.perm)
	ids.in <- stat > minimum
	pval[!ids.in] <- perm
	ids.in <- !is.na(stat) & ids.in	
	ids.het <- matA[,3]==0 & ids.in
	if(verbose)
		cat("Considering SNPs without trios with two heterozygous parents.\n")
	tmpPval <- numeric(sum(ids.het))
	for(i in 1:le.perm){
		if(verbose)
			cat(i, "of", le.perm, "chunks. Considering", vec.perm[i], "permutations.\n")
		tmpPval <- tmpPval + noBothHet(matA[ids.het,1], matA[ids.het,2], vec.perm[i], stat[ids.het])
	}
	pval[ids.het] <- tmpPval
	if(verbose)
		cat("Finished at ", date(), ".\n")
	ids.het <- !ids.het & ids.in
	if(verbose)
		cat("Considering ", sum(ids.het), " SNPs with trios composed of two heterozygous parents.\n",
			sep="")
	matA <- matA[ids.het,]
	uni.val <- unique(matA[,3])
	if(verbose)
		cat("These SNPs show ", length(uni.val), " different numbers of trios with two ",
			"heterozygous parents.\n", sep="")
	tmpPval <- numeric(sum(ids.het))
	for(i in 1:le.perm){
		if(verbose)
			cat(i, "of", le.perm, "chunks. Considering", vec.perm[i], "permutations.\n")
		tmpPval <- tmpPval + bothHet(matA[,1], matA[,2], matA[,3], vec.perm[i], uni.val,
			stat[ids.het], denom[ids.het], verbose=verbose)
	}
	pval[ids.het] <- tmpPval
	names(pval) <- names(stat)
	out <- list(pval=pval/perm, stat=stat, mat.stat=mat.stat)
	class(out) <- "maxTestTrio"
	out
} 

print.maxTestTrio <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("Max-Statistic"=x$stat, x$mat.stat, "p-Value"=pval,
		check.names=FALSE)
	cat("                  Maximum Genotypic TDT\n\n")
	if(length(x$stat) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNPs:\n")
	}
	print(format(out, digits=digits))
}
	

bothHet <- function(a12, a34, a567, perm, uni.val, stat, denom, verbose=FALSE){
	vec <- numeric(length(stat))
	count <- 0
	if(verbose)
		n.val <- length(uni.val)
	for(i in uni.val){
		if(verbose){
			count <- count+1
			cat("Computing p-values for SNPs with ", i, " trios showing two heterozygous parents (",
				count, " of ", n.val, ").\n", sep="")
		}	 
		tmpids <- a567==i
		vec[tmpids] <- compBothHet(a12[tmpids], a34[tmpids], i, perm, stat[tmpids], denom[tmpids])
		if(verbose)
			cat("Finished at ", date(), ".\n", sep="") 
	}
	vec
} 


compBothHet <- function(a12, a34, val, perm, stat, denom){
	matA567 <- rmultinom(perm, val, c(0.25, 0.5, 0.25))
	a5 <- matA567[1,]
	a6 <- matA567[2,]
	a7 <- matA567[3,]
	rm(matA567)
	tmpvec <- length(stat)
	for(i in 1:length(stat)){
		a2 <- rbinom(perm, a12[i], 0.5)
		a1 <- a12[i] - a2
		a4 <- rbinom(perm, a34[i], 0.5)
		a3 <- a34[i] - a4
		h <- a2 + a4 + a6 + 2*a7
		beta <- logit(h/denom[i])
		v <- denom[i] / ((denom[i] - h) * h)
		add <- beta * beta / v
		d <- a6 + a7
		h <- (1/3*a1 - a2 + a5 - 1/3*d) / (2 * (a1+a5))
		tmp <- (a2 + d) / (3 * (a1+a5)) + h*h
		or <- sqrt(tmp) - h
		beta <- log(or)
		tmp <- a12[i] * or/(or+1)^2 + val * or/(3*(or+1/3)^2)
		v <- 1/tmp
		dom <- beta * beta / v
		d <- a5 +a6
		h <- (3*a3 - a4 + d - 3*a7) / (2*(a3+d))
		tmp <- 3 * (a4+a7) / (a3+d) + h*h
		or <- sqrt(tmp) - h
		beta <- log(or)
		tmp <- a34[i] * or/(or+1)^2 + 3 * val * or/(or+3)^2
		v <- 1/tmp
		rec <- beta * beta / v
		tmpvec[i] <- sum(stat[i] <= pmax(add,dom,rec, na.rm=TRUE), na.rm=TRUE)
	}
	tmpvec
}		 
		

noBothHet <- function(a12, a34, perm, stat){
	vec <- numeric(length(a12))
	tmpids <- a34==0
	tmpvec <- numeric(sum(tmpids))
	tmpa12 <- a12[tmpids]
	for(i in unique(tmpa12)){
		a2 <- rbinom(perm, i, 0.5)
		beta <- logit(a2/i)
		v <- i / ((i-a2) * a2)
		add <- beta*beta/v
		a1 <- i - a2
		h <- (1/3*a1 - a2) / (2*a1)
		or <- sqrt(a2/(3*a1) + h*h) - h
		beta <- log(or)
		v <- (or+1)^2 / (i * or)
		dom <- beta*beta/v
		tmpstat <- stat[tmpids & a12==i]
		for(k in 1:length(tmpstat))
			tmpstat[k] <- sum(tmpstat[k] <= pmax(add, dom, na.rm=TRUE), na.rm=TRUE)
		tmpvec[tmpa12==i] <- tmpstat
	}
	vec[tmpids] <- tmpvec
	tmpvec <- numeric(sum(!tmpids))
	tmpa12 <- a12[!tmpids]
	tmpa34 <- a34[!tmpids]
	for(i in unique(tmpa12)){
		a2 <- rbinom(perm, i, 0.5)
		a1 <- i - a2
		for(j in tmpa34[tmpa12==i]){
			a4 <- rbinom(perm, j, 0.5)
			a3 <- i - a4
			num <- a2+a4
			denom <- i+j
			beta <- logit(num/denom)
			v <- denom/((denom-num)*num)
			add <- beta*beta/v
			num <- (1/3*a1-a2) / (2*a1)
			or <- sqrt(a2/(3*a1) + num*num) - num
			beta <- log(or)
			v <- (or+1)^2 / (i*or)
			dom <- beta*beta/v
			num <- (3*a3 - a4) / (2*(a3))
			or <- sqrt(3*a4/a3 + num*num) - num 
			beta <- log(or)
			v <- (or+1)^2 / (j*or)
			rec <- beta*beta/v
			ids2 <- which(a12==i & a34==j)
			for(k in ids2)
				vec[k] <- sum(stat[k] <= pmax(add, dom, rec, na.rm=TRUE), na.rm=TRUE)
		}
	}
	vec
}

			