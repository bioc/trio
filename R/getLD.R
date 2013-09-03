getLD <- function(x, which=c("both", "rSquare", "Dprime"), parentsOnly=FALSE, 
		iter=50, snp.in.col=TRUE, asMatrix=FALSE, addVarN=FALSE){
	if(!is.matrix(x))
		stop("x must be a matrix.")
	if(parentsOnly){
		if(!snp.in.col)
			stop("The columns of x must represent the SNPs if parentsOnly=TRUE.")
		x <- removeKids(x)
	}
	if(snp.in.col)
		x <- t(x)
	typeLD <- match.arg(which)
	n.snp <- nrow(x)
	if(is.null(rownames(x)))
		rownames(x) <- paste("S", 1:n.snp, sep="")
	D.out <- compD(x, iter=iter, addVar=addVarN)
	D <- D.out$D
	pA1A2 <- D.out$p1 * (1 - D.out$p2)
	pA2B1 <- (1 - D.out$p1) * D.out$p2
	n <- if(addVarN) D.out$n else NULL
	if(typeLD %in% c("rSquare", "both")){
		r2 <- D * D / (pA1A2 * pA2B1)
	}
	else
		r2 <- NULL
	if(typeLD %in% c("Dprime", "both")){
		Dmax <- rep.int(NA, length(D))
		ids <- which(D > 0)
		Dmax[ids] <- pmin(pA1A2[ids], pA2B1[ids])
		ids2 <- which(D <= 0)
		if(length(ids2) > 0){
			pA1B1 <- D.out$p1[ids2] * D.out$p2[ids2]
			pA2B2 <- (1-D.out$p1[ids2]) * (1-D.out$p2[ids2])
			Dmax[ids2] <- pmax(-pA1B1, -pA2B2)
		}
		if(addVarN){
			outDprime <- compVarDprime(D, Dmax, D.out$p1, D.out$p2, 
				D.out$n, D.out$varD)
			D <- outDprime$Dprime
			varDprime <- outDprime$varDprime
		}
		else{
			varDprime <- NULL
			D <- D/Dmax
		}
	
		if(any(D>1))
			D[D>1] <- 1
	}
	else
		D <- varDprime <- NULL
	if(!asMatrix){
		names(D) <- names(r2) <- NULL
		LDout <- list(Dprime=D, rSquare=r2, n=n, varDprime=varDprime, 
			rn=rownames(x))
		class(LDout) <- "getLD"
		return(LDout)
	}
	rn <- rownames(x)
	if(!is.null(r2)){
		mat.r2 <- matrix(nrow=n.snp, ncol=n.snp, dimnames=list(rn, rn))
		mat.r2[lower.tri(mat.r2)] <- r2
	}
	else
		mat.r2 <- NULL
	if(!is.null(D)){
		mat.D <- matrix(nrow=n.snp, ncol=n.snp, dimnames=list(rn, rn))
		mat.D[lower.tri(mat.D)] <- D
	}
	else
		mat.D <- NULL
	LDout <- list(Dprime=mat.D, rSquare=mat.r2, n=n, varDprime=varDprime)
	class(LDout) <- "getLD"
	LDout
}

compD <- function(x, iter=50, addVar=FALSE){
	xrange <- range(x, na.rm=TRUE)
	if(xrange[1]!=0 | xrange[2]!=2)
		stop("x must contain values between 0 and 2 specifying the numbers of",
			"minor alleles.", call.=FALSE)
	out <- est.pAB(x, iter=iter)
	if(any(out$n<10))
		stop("For at least one pair of SNPs, the genotypes of less than 5 individuals\n",
			"are available at both SNPs.", call.=FALSE)
	D <- out$f11 - out$p1 * out$p2
	if(addVar){
		tmp <- D * (1-2*out$p1) * (1-2*out$p2) - D^2
		varD <- (out$p1 * (1-out$p1) * out$p2 * (1-out$p2) + tmp) / out$n
	}
	else
		varD <- NULL
	return(list(D=D, p1=out$p1, p2=out$p2, n=out$n, varD=varD))
}	

est.pAB <- function(x, iter=50){
	out <- compX11(x)
	n <- 2 * out$n
	x11 <- out$x11 / 2
	n22 <- out$n22
	if(is.null(out$p1)){
		combs <- allCombs(nrow(x))
		pA <- 1- 0.5 * rowMeans(x, na.rm=TRUE)
		p1 <- pA[combs[,1]]
		p2 <- pA[combs[,2]]
	}
	else{
		p1 <- out$p1 / n
		p2 <- out$p2 / n
	}
	f11 <- (2*x11 + n22) / n - p1 * p2
	for(i in 1:iter){
		fac <- f11 * (1 - p1 - p2 + f11)
		f11 <- (x11 + n22 * fac / (fac + (p1-f11) * (p2-f11))) / n
	}
	return(list(f11=f11, p1=p1, p2=p2, n=n))
}

compX11 <- function(x){
	x <- 2-x
	naIdentifier <- !is.na(x)
	if(all(naIdentifier)){
		n <- ncol(x)
		p1 <- p2 <- NULL
	}
	else{
		tmp.mat <- naIdentifier %*% t(naIdentifier)
		n <- tmp.mat[lower.tri(tmp.mat)]
		x[!naIdentifier] <- 0
		tmp.mat <- naIdentifier %*% t(x)
		p1 <- tmp.mat[lower.tri(tmp.mat)]
		p2 <- tmp.mat[allCombs(nrow(x), upperID=TRUE)]
		
	}
	tmp.mat <- x %*% t(x)
	x11 <- tmp.mat[lower.tri(tmp.mat)]
	heteIdentifier <- x==1
	tmp.mat <- heteIdentifier %*% t(heteIdentifier)
	n22 <- tmp.mat[lower.tri(tmp.mat)]
	return(list(x11=x11-n22, n22=n22, n=n, p1=p1, p2=p2))
}
	
allCombs <- function(m, upperID=FALSE){
	vec1 <- rep.int(1:(m-1), (m-1):1)
	txt <- paste("c(", paste(2:m, "m", sep=":", collapse=","), ")", sep="")
	vec2 <- eval(parse(text=txt))
	if(upperID)
		return(vec1 + m * (vec2 - 1))
	cbind(vec1, vec2)
}

plot.getLD <- function(x, y="rSquare", start=1, end=NA, squared=TRUE, col=NULL, 
		xlab="", ylab="", cexAxis=0.8, alpha=0.1, ciLD=c(0.7,0.98), 
		cuRecomb=0.9, ...){
	if(!y %in% c("rSquare", "Dprime", "gabriel"))
		stop("y must be either rSquare, Dprime, or gabriel.")
	if(y=="gabriel"){
		if(is.null(x$varDprime))
			stop("varDprime not available in x.")
		y <- "Dprime"
		gab <- TRUE
	}
	else
		gab <- FALSE
	mat <- x[[y]]
	if(is.null(mat))
		stop(y, " is not available in x.")
	if(gab)
		mat <- getCalls4LD(mat, x$varDprime, alpha=alpha, ciLD=ciLD, 
			cuRecomb=cuRecomb)
	if(!is.matrix(mat)){
		tmp <- mat
		mat <- matrix(nrow=length(x$rn), ncol=length(x$rn),
			dimnames=list(x$rn, x$rn))
		mat[lower.tri(mat)] <- tmp
	} 
	if(is.character(start)){
		start <- which(rownames(mat) == start)
		if(length(start)==0)
			stop("start does not specify the name of a SNP.")
	}
	if(is.na(end))
		end <- nrow(mat)
	else{
		if(is.character(end)){
			end <- which(rownames(mat) == end)
			if(length(end)==0)
				stop("end does not specify the name of a SNP.")
		}
	}
	if(start<1)
		stop("start must be a positiv integer.")
	if(end > nrow(mat))
		stop("end is not allowed to be larger than the number of SNPs.")
	if(start>=end)
		stop("start must be larger than end.")
	if(y=="rSquare" && !squared)
		mat <- sqrt(mat)
	mat <- as.matrix(rev(as.data.frame(t(mat[start:end, start:end]))))
	mat <- mat[-nrow(mat), -ncol(mat)]
	if(gab){
		if(!is.null(col) && length(col)!=3)
			stop("For y='gabriel', col must be a vector of length 3.")
		colSel <- sort(unique(mat[!is.na(mat)])) + 2
	}
	par(mar=c(5,5,1,1))
	if(is.null(col)){
		col <- if(gab) c("white", "yellow", "blue")[colSel]
			else col <- gray(c(200:0)/200)
	}
	image(1:nrow(mat), 1:ncol(mat), mat, col=col, xlab=xlab, ylab=ylab, xaxt="n",
		yaxt="n", ...)		
	axis(2, at=1:ncol(mat), labels=colnames(mat), cex.axis=cexAxis, las=1)
	axis(1, at=1:nrow(mat), labels=rownames(mat), cex.axis=cexAxis, las=3)		
	if(gab)
		legend("topright", legend=c("Recombination", "Other", "Strong LD")[colSel], 
			col=col, bg="grey85", pch=15, cex=0.9) 
}


compVarDprime <- function(D, Dmax, p1, p2, n, varD){
	a <- ifelse(Dmax>0, p2, 1-p2)
	b <- 1-a
	pAB <- D + p1 * p2
	x <- rep(NA, length(D))
	ids <- which(Dmax == -p1*p2)
	if(length(ids) > 0)
		x[ids] <- pAB[ids]
	ids <- which(Dmax == p1*(1-p2))
	if(length(ids) > 0)
		x[ids] <- (p1-pAB)[ids]
	ids <- which(Dmax == (1-p1)*p2)
	if(length(ids) > 0)
		x[ids] <- (p2-pAB)[ids]
	ids <- which(Dmax == -(1-p1)*(1-p2))
	if(length(ids) > 0)
		x[ids] <- (1-p1-p2+pAB)[ids]
	Dprime <- D/Dmax
	Dmax <- abs(Dmax)
	D <- abs(D)
	varDprime <- (1-Dprime) * (n * varD - Dprime * Dmax * (a*p1 + b*(1-p1) - 2*D))
	num <- n * Dmax^2
	varDprime <- (varDprime + Dprime * x* (1-x)) / num
	if(any(varDprime<0, na.rm=TRUE)){
		ids2 <- which(varDprime < 0)
		if(any(round(Dprime[ids2],8)<1, na.rm=TRUE))
			warning("At least one of the variance estimates for a D'<1 is negative.",
				call.=FALSE)
		varDprime[ids2] <- 0
	}
	return(list(Dprime=Dprime, varDprime=varDprime))
}


getCalls4LD <- function(Dprime, varDprime, alpha=0.1, ciLD=c(0.7, 0.98), cuRecomb=0.9){
	if(alpha>=0.5 | alpha<=0)
		stop("alpha must be larger than zero and smaller than 0.5")
	if(length(ciLD)!=2)
		stop("ciLD must contain exactly two values.")
	if(any(ciLD>=1 | ciLD<0))
		stop("The values in ciLD must between 0 and 1.")
	if(ciLD[1]>=ciLD[2])
		stop("The first element in ciLD must be larger than the second.")
	if(length(cuRecomb)!=1)
		stop("cuRecomb must be a numeric value.")
	if(cuRecomb>=ciLD[2] | cuRecomb<0)
		stop("cuRecomb must be larger than zero and smaller than ciLD[2].")
	if(is.matrix(Dprime))
		Dprime <- Dprime[lower.tri(Dprime)]
	zsd <- qnorm(1-alpha/2) * sqrt(varDprime)
	upper <- Dprime + zsd
	lower <- Dprime - zsd
	calls <- numeric(length(Dprime))
	calls[upper>=ciLD[2] & lower>=ciLD[1]] <- 1
	calls[upper<cuRecomb] <- -1
	calls
}

	
findLDblocks <- function(x, alpha=0.1, ciLD=c(0.7, 0.98), cuRecomb=0.9, ratio=9,
		alsoOthers=FALSE, parentsOnly=FALSE, iter=50, snp.in.col=TRUE){
	if(!is(x, "getLD") && !is(x, "getLDlarge") && !is.matrix(x))
		stop("x must either be a matrix or the output of getLD or getLDlarge.")
	if(is.matrix(x)){
		if((snp.in.col & ncol(x) > 500) | (!snp.in.col & nrow(x) > 500))
			stop("If x is a matrix, x must contain at most 500 SNPs. If x contains more SNPs,\n",
				"please consider to use getLDlarge.")	
		x <- getLD(x, which="Dprime", parentsOnly=parentsOnly, iter=iter, 
			snp.in.col=snp.in.col, addVarN=TRUE)
		large <- FALSE
		}
	else
		large <- is(x, "getLDlarge")
	calls <- getCalls4LD(x$Dprime, x$varDprime, alpha=alpha, ciLD=ciLD,
		cuRecomb=cuRecomb)
	n.snps <- length(x$rn)
	if(large)
		combs <- allNeighborCombs(n.snps, x$neighbors)[,c(2,1)]
	else
		combs <- allCombs(n.snps)
	i <- b <- 1
	endLD <- ifelse(alsoOthers, 0, 1)
	idsEnd <- calls >= endLD
	blocks <- numeric(n.snps)
	while(i<=n.snps){
		posj <- combs[combs[,1]==i & idsEnd, 2] 
		if(length(posj)==0){
			blocks[i] <- b
			i <- i+1
		}
		else{
			vec.ratio <- numeric(length(posj))
			for(j in 1:length(vec.ratio)){
				tmpCalls <- calls[combs[,1]>=i & combs[,2]<=posj[j]]
				vec.ratio[j] <- sum(tmpCalls==1) / max(1/ratio, sum(tmpCalls==-1))
			}
			idsRatio <- vec.ratio>=ratio
			if(any(idsRatio)){
				selj <- max(posj[idsRatio])
				blocks[i:selj] <- b
				i <- selj+1
			}
			else{
				blocks[i] <- b
				i <- i+1
			}
		}
		b <- b+1		
	}
	names(blocks) <- x$rn
	vec.blocks <- split(names(blocks), blocks)
	idsLength <- sapply(vec.blocks, length)
	vec.blocks <- vec.blocks[idsLength>1]
	param <- list(alpha=alpha, ciLD=ciLD, cuRecomb=cuRecomb, ratio=ratio)
	out <- list(ld=x, blocks=blocks, vec.blocks=vec.blocks, param=param, large=large)
	class(out) <- "LDblocks"
	out
}

print.LDblocks <- function(x, ...){
	le <- sapply(x$vec.blocks, length)
	cat("Found", length(x$vec.blocks), "LD blocks containing between", min(le), "and", 
		max(le), "SNPs.\n")
	cat(length(x$blocks)-sum(le),"of the", length(x$blocks), "SNPs do not belong to a LD block.",
		"\n\n")
	cat("Used Parameter:\n")
	cat("Strong LD:    ", "C_L >=", x$param$ciLD[1], "and C_U >=", x$param$ciLD[2],"\n")
	cat("Recombination:", "C_U <", x$param$cuRecomb, "\n")
	cat("(C_L and C_U are the lower and upper bound of\n")
	cat("the ", (1-x$param$alpha)*100, "%-confidence intervals for D')\n", sep="")
	cat("LD blocks: Ratio >=", x$param$ratio, "\n\n")
}

plot.LDblocks <- function(x, y="gabriel", col=NULL, start=1, end=NA, xlab="", ylab="", 
		cexAxis=0.8, block.col=2, block.lwd=3, ...){
	if(!y %in% c("gabriel", "Dprime"))
		stop("y must be either gabriel or Dprime.")
	if(is.character(start)){
		start <- which(names(x$blocks) == start)
		if(length(start)==0)
			stop("start does not specify the name of a SNP.")
	}
	if(is.na(end))
		end <- length(x$blocks)
	else{
		if(is.character(end)){
			end <- which(names(x$blocks) == end)
			if(length(end)==0)
				stop("end does not specify the name of a SNP.")
		}
	}
	namesBlock <- names(x$blocks)[start:end]
	n.snps <- length(namesBlock)
	getRangeBlock <- function(z, nm=namesBlock){
		ids <- which(nm %in% z)
		if(length(ids)==0)
			return(NULL)
		range(ids)
	}
	rangeBlock <- lapply(x$vec.blocks, getRangeBlock)
	if(start!=1 | end!=length(x$block)){
		nullBlock <- sapply(rangeBlock, is.null)
		rangeBlock <- rangeBlock[!nullBlock]
		vecBlock <- x$vec.blocks[!nullBlock]
		if(length(rangeBlock)>0){
			if(start==1 || rangeBlock[[1]][1]!=1)
				showSeg1 <- TRUE
			else{
				tmpids <- which(names(x$blocks)==namesBlock[1])
				showSeg1 <- ifelse(x$blocks[tmpids] == x$blocks[tmpids-1],
					FALSE, TRUE)
			}
			if(end==length(x$block) || rangeBlock[[length(rangeBlock)]][2]!=n.snps)
				showSegLast <- TRUE
			else{
				tmpids <- which(names(x$blocks) == namesBlock[n.snps])
				showSegLast <- ifelse(x$blocks[tmpids] == x$blocks[tmpids+1],
					FALSE, TRUE)
			}
		}
	}
	else
		showSeg1 <- showSegLast <- TRUE
	plot(x$ld, y=y, col=col, start=start, end=end, xlab=xlab, ylab=ylab, cexAxis=cexAxis,
		alpha=x$param$alpha, ciLD=x$param$ciLD, cuRecomb=x$param$cuRecomb,  ...)
	if(length(rangeBlock)>0){
		for(i in 1:length(rangeBlock)){
			tmpMin <- rangeBlock[[i]][1] - 0.5
			tmpMax <- rangeBlock[[i]][2] - 0.5
			if(showSegLast || i!=length(rangeBlock))
				segments(tmpMin, n.snps-tmpMax, tmpMax, n.snps-tmpMax, 
					lwd=block.lwd, col=block.col, xpd=NA)
			if(showSeg1 || i!=1)
				segments(tmpMin, n.snps-tmpMax, tmpMin, n.snps-tmpMin, 
					lwd=block.lwd, col=block.col, xpd=NA)
		}
	}
}


removeKids <- function(x){
	if(nrow(x) %%3 != 0)
		stop("x does not seem to be in genotype format, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(x)))
		stop("x does not seem to be in genotype format, as the row names are missing.")
	x <- x[-seq(3, nrow(x), 3), ]
	x <- x[!duplicated(rownames(x)),]
	x
}

			





