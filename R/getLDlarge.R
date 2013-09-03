getLDlarge <- function(x, neighbors=25, which=c("both", "rSquare", "Dprime"), parentsOnly=FALSE, iter=50,
		snp.in.col=TRUE, addVarN=FALSE){
	if(!is.matrix(x))
		stop("x must be a matrix.")
	if(neighbors < 1)
		stop("neighbors should be at least 1.")
	if(parentsOnly){
		if(!snp.in.col)
			stop("The columns of x must represent the SNPs if parentsOnly=TRUE.")
		x <- removeKids(x)
	}
	if(any(!x %in% c(0,1,2,NA)))
		stop("The values in x must be 0, 1, and 2.")
	if(snp.in.col)
		x <- t(x)
	typeLD <- match.arg(which)
	if(ncol(x) < 5)
		stop("There must be data available for at least 5 individuals.")
	n.snp <- nrow(x)
	if(n.snp/neighbors < 2)
		stop("The number of SNPs should be at least 2 * neighbors.")
	if(is.null(rownames(x)))
		rownames(x) <- paste("S", 1:n.snp, sep="")
	if(any(is.na(x)))
		n <- p1 <- p2 <- NULL
	else{
		n <- 2*ncol(x)
		pA <- 1 - 0.5*rowMeans(x)
		combs <- neighborCombs(neighbors)
	}
	x <- 2-x
	intS <- c(seq.int(1, n.snp, neighbors), n.snp+1)
	nD <- neighbors * (n.snp-neighbors) + neighbors*(neighbors-1)/2
	intD <- c(seq.int(1, nD, neighbors^2), nD+1)
	if(typeLD %in% c("rSquare", "both")){
		r2 <- rep.int(NA, nD)
		isr2 <- TRUE
	}
	else{
		isr2 <- FALSE
		r2 <- NULL
	}
	if(typeLD %in% c("Dprime", "both")){
		Dprime <- rep.int(NA, nD)
		isDp <- TRUE
	}
	else{
		isDp <- addVarN <- FALSE
		Dprime <- NULL
	}
	if(addVarN)
		varDprime <- vecN <- rep.int(NA, nD)
	else
		varDprime <- vecN <- NULL
	leS <- length(intS)
	idsStart <- 1	
	for(i in 1:(leS-2)){
		if(!is.null(n)){
			tmppa <- pA[intS[i]:(intS[i+2]-1)]
			if(i==leS-2)
				combs <- combs[combs[,2]<=length(tmppa),]
			p1 <- tmppa[combs[,1]]
			p2 <- tmppa[combs[,2]]
		}
		out <- compLDsplit(x[intS[i]:(intS[i+2]-1),], p1=p1, p2=p2, n=n, iter=iter, 
			type=typeLD, addVarN=addVarN, neighbors=neighbors)
		idsEnd <- idsStart + out$leD - 1
		if(isr2)
			r2[idsStart:idsEnd] <- out$r2
		if(isDp)
			Dprime[idsStart:idsEnd] <- out$Dprime
		if(addVarN){
			varDprime[idsStart:idsEnd] <- out$varDprime
			vecN[idsStart:idsEnd] <- out$n
		}
		idsStart <- idsEnd + 1
	}
	if(intS[leS]-intS[leS-1]>1){
		out <- getLD(x[intS[leS-1]:(intS[leS]-1),], which=typeLD, iter=iter, snp.in.col=FALSE,
			addVarN=addVarN)
		nLast <- max(length(out$rSquare), length(out$Dprime))
		startLast <- nD - nLast + 1
		if(isr2)
			r2[startLast:nD] <- out$rSquare
		if(isDp)
			Dprime[startLast:nD] <- out$Dprime
		if(addVarN){
			varDprime[startLast:nD] <- out$varDprime
			vecN[startLast:nD] <- out$n
		}
	}
	out <- list(Dprime=Dprime, rSquare=r2, n=vecN, varDprime=varDprime, rn=rownames(x), neighbors=neighbors)
	class(out) <- "getLDlarge"
	out 
}



compLDsplit <- function(x, p1=NULL, p2=NULL, n=NULL, iter=50, type="both", addVarN=FALSE, neighbors=25){
	out <- est.pAB2Mat(x, x[1:neighbors,], p1=p1, p2=p2, n=n, iter=iter)
	D <- out$f11 - out$p1 * out$p2
	pA1A2 <- out$p1 * (1-out$p2)
	pA2B1 <- (1-out$p1) * out$p2
	D2 <- D*D
	leD <- length(D)
	if(type %in% c("rSquare", "both"))
		r2 <- D2 / (pA1A2 * pA2B1)
	else
		r2 <- NULL
	if(type %in% c("Dprime", "both")){
		Dmax <- rep.int(NA, length(D))
		idsL <- which(D > 0)
		Dmax[idsL] <- pmin(pA1A2[idsL], pA2B1[idsL])
		idsSE <- which(D <= 0)
		if(length(idsSE) > 0){
			pA1B1 <- out$p1[idsSE] * out$p2[idsSE]
			pA2B2 <- (1-out$p1[idsSE]) * (1-out$p2[idsSE])
			Dmax[idsSE] <- pmax(-pA1B1, -pA2B2)
		}
		if(addVarN){
			tmp <- D * (1 - 2 * out$p1) * (1 - 2 * out$p2) - D2
			varD <- (pA1A2 * pA2B1 + tmp) / out$n
			outD <- compVarDprime(D, Dmax, out$p1, out$p2, out$n, varD)
			D <- outD$Dprime
			varDprime <- outD$varDprime
		}
		else{
			varDprime <- NULL
			D <- D / Dmax
		}
		idsLarger1 <- which(D>1)
		if(length(idsLarger1) > 0)
			D[idsLarger1] <- 1
	}
	else
		D <- varDprime <- NULL
	list(r2=r2, Dprime=D, varDprime=varDprime, n=out$n, leD=leD)
}


est.pAB2Mat <- function(x1, x2, p1=NULL, p2=NULL, n=NULL, iter=50){
	withna <- is.null(n)
	out <- comp2MatX11(x1, x2, withna=withna)
	x11 <- out$x11/2
	n22 <- out$n22
	if(withna){
		n <- out$n
		p1 <- out$p1/n
		p2 <- out$p2/n
	}
	f11 <- (2 * x11 + n22)/n - p1 * p2
	for(i in 1:iter){
		fac <- f11 * (1 - p1 - p2 + f11)
		f11 <- (x11 + n22 * fac / (fac + (p1 - f11) * (p2 - f11))) / n
	}
	list(f11=f11, p1=p1, p2=p2, n=n)
}


neighborCombs <- function(neighbors){
	vec1 <- rep(1:neighbors, each=neighbors)
	txt <- paste("c(", paste(2:(neighbors+1), (neighbors+1):(2*neighbors), sep=":", collapse=","),
		")", sep="")
	vec2 <- eval(parse(text=txt))
	cbind(vec1,vec2)
}



comp2MatX11 <- function(x1, x2, withna=FALSE){
	if(withna){
		naId1 <- !is.na(x1)
		naId2 <- !is.na(x2)
		n <- naId1 %*% t(naId2)
		if(any(n<5))
			stop("For at least one pair of SNPs, the genotypes of less than 5 individuals\n",
				"are available for both SNPs.", call.=FALSE)
		x1[!naId1] <- 0
		x2[!naId2] <- 0
		p1 <- naId1 %*% t(x2)
		p2 <- x1 %*% t(naId2)
	}
	tmp.mat <- x1 %*% t(x2)
	idsMat <- neighborIds(tmp.mat, nrow(x2))
	x11 <- tmp.mat[idsMat]
	tmp.mat <- (x1==1) %*% t(x2==1)
	n22 <- tmp.mat[idsMat]
	if(withna)
		return(list(x11=x11-n22, n22=n22, p1=p1[idsMat], p2=p2[idsMat], n=2*n[idsMat]))
	list(x11=x11-n22, n22=n22)
}

neighborIds <- function(x, neighbors){
	tmp <- row(x) - col(x)
	tmp>0 & tmp <= neighbors
}

	
plot.getLDlarge <- function(x, y="rSquare", start=NA, end=NA, squared=TRUE, col=NULL, xlab="", ylab="", cexAxis=0.8,
		alpha=0.1, ciLD=c(0.7,0.98), cuRecomb=0.9, ...){
	if(!y %in% c("rSquare", "Dprime", "gabriel"))
		stop("y must be either 'rSquare', 'Dprime', or 'gabriel'.")
	if(y == "gabriel"){
		if(is.null(x$varDprime))
			stop("varDprime not available in x, but required when using y='gabriel'.")
		if(!is.null(col) && length(col) != 3)
			stop("For y='gabriel', col must be a vector of length 3.")
		y <- "Dprime"
		gab <- TRUE
	}
	else
		gab <- FALSE
	vecLD <- x[[y]]
	if(is.null(vecLD))
		stop(y, " is not available in x.")
	if(is.na(start) | is.na(end))
		stop("Both start and end must be specified.")
	if(is.character(start)){
		start <- which(x$rn==start)
		if(length(start)!=1)
			stop("start does not specify the name of a SNP.")
	}
	if(is.character(end)){
		end <- which(x$rn==end)
		if(length(end)!=1)
			stop("end does not specify the name of a SNP.")
	}
	if(start<1)
		stop("start must be a positive integer.")
	if(start > length(x$rn)-2)
		stop("Since x contains data for ", length(x$rn), " SNPs, start must be at most ",
			length(x$rn)-2, ".")
	if(end > length(x$rn))
		stop("end must be smaller than or equal to the number of considered SNPs.") 
	if(start >= end)
		stop("start must be larger than end.")
	idsLD <- allNeighborCombs(length(x$rn), x$neighbors, start=start, end=end)
	vecLD <- vecLD[idsLD]
	if(gab){
		vecLD <- getCalls4LD(vecLD, x$varDprime[idsLD], alpha=alpha, ciLD=ciLD, 
			cuRecomb=cuRecomb)
		colSel <- sort(unique(vecLD)) + 2
	}		
	if(y == "rSquare" && !squared)
		vecLD <- sqrt(vecLD)
	n.snps <- end - start + 1
	mat <- matrix(NA, n.snps, n.snps)
	idsPos <- (row(mat) - col(mat) <= x$neighbors) & (row(mat) - col(mat) > 0)
	mat[idsPos] <- vecLD
	colnames(mat) <- rownames(mat) <- x$rn[start:end]
	mat <- as.matrix(rev(as.data.frame(t(mat))))
	mat <- mat[-nrow(mat), -ncol(mat)]
	if(is.null(col))
		col <- if(gab) c("white", "yellow", "blue")[colSel]  else gray(c(200:0)/200)	 
	par(mar=c(5,5,1,1))
	image(1:nrow(mat), 1:ncol(mat), mat, col=col, xlab=xlab, ylab=ylab, xaxt="n", yaxt="n", ...)
	axis(2, at=1:ncol(mat), labels=colnames(mat), cex.axis=cexAxis, las=1)
	axis(1, at=1:nrow(mat), labels=rownames(mat), cex.axis=cexAxis, las=3)
	if(gab)
		legend("topright", legend=c("Recombination", "Other", "Strong LD")[colSel], col=col,
			bg="grey85", pch=15, cex=0.9)
}
	

allNeighborCombs <- function(m, n, start=NA, end=NA){
	vec2 <- c(rep(1:(m-n), each=n), rep.int((m-n+1):(m-1), (n-1):1))
	txt <- paste("c(", paste(2:m, c((n+1):m, rep(m, n-1)), sep=":", collapse=","), ")", sep="")
	vec1 <- eval(parse(text=txt))
	if(is.na(start))
		return(cbind(vec1=vec1, vec2=vec2))
	vec2 >= start & vec1 <= end
} 
	

