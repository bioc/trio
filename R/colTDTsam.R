colTDTsam <- function(mat.snp, model=c("additive", "dominant", "recessive", "max"),
		approx=NULL, B=1000, size=10, chunk=100, rand=NA){
	requireNamespace("siggenes")
	if(!is.matrix(mat.snp))
		stop("mat.snp must be a matrix.")
	if(nrow(mat.snp) %%3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0:2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	type <- match.arg(model)
	if(is.null(approx))
		approx <- type != "max"
	if(!is.logical(approx))
		stop("approx must be a logic value.")
	sam(mat.snp, type, method=gtdt.stat, gene.names=colnames(mat.snp), approx=approx, B=B, 
		size=size, chunk=chunk, rand=rand)
}


gtdt.stat <- function(data, cl, approx=TRUE, B=1000, size=10, chunk=100, rand=NA){
	if(cl == "max"){
		if(approx)
			stop("If a MAX-statistic should be used, then approx must be FALSE.")
		return(maxStatSAM(data, size=size, B=B, chunk=chunk, rand=rand))
	}
	out <- fastTDT(data, cl, size=size)
	if(approx)
		null.out <- gtdt.null.approx(out$pval)
	else
		null.out <- gtdt.null.samp(data, cl, out$stat, size=size, B=B, chunk=chunk, rand=rand)
	cap <- gsub("\\b(\\w)", "\\U\\1", cl, perl=TRUE)
	msg <- paste("SAM for Trio Data (Based on", cap, "gTDT)\n\n")
	list(d=out$stat, d.bar=null.out$d.bar, p.value=null.out$p.value,
		vec.false=null.out$vec.false, s=out$se, s0=numeric(0),
		mat.samp=null.out$mat.samp, msg=msg, fold=numeric(0))
}


gtdt.null.approx <- function(pval){
	m <- length(na.exclude(pval))
	d.bar <- qchisq(((1:m)-0.5)/m, 1)
	list(d.bar=d.bar, p.value=pval, vec.false=m*pval, mat.samp=matrix(numeric(0)))
}

		
maxStatSAM <- function(geno, size=10, B=1000, chunk=100, rand=NA){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- numeric(n.snp)
	matA <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1)){
		tmpMat <- getStatsChunk(geno[, int[i]:(int[i+1]-1), drop=FALSE])
		stat[int[i]:(int[i+1]-1)] <- tmpMat[,4]
		matA[int[i]:(int[i+1]-1),] <- tmpMat[,c(1:3,5)]
	}
	vec.perm <- diff(unique(c(seq.int(1, B, chunk), B + 1)))
	le.perm <- length(vec.perm)
	ids.in <- !is.na(stat)
	matA <- matA[ids.in,]
	statnotna <- stat[ids.in]
	ids.het <- matA[,3] == 0 
	d.rank <- rank(-statnotna, ties.method="first")
	n.notna <- length(statnotna)
	vecf <- -d.rank * le.perm
	mat.dperm <- matrix(0, n.notna, le.perm)
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:le.perm){
		tmpMat <- matrix(NA, n.notna, vec.perm[i])
		tmpMat[ids.het,] <- noBothHetMat(matA[ids.het, 1], matA[ids.het, 2], vec.perm[i])
		tmpMat[!ids.het,] <- bothHetMat(matA[!ids.het, 1], matA[!ids.het, 2], matA[!ids.het, 3],
			vec.perm[i], matA[!ids.het, 4])
		tmpMat <- apply(tmpMat, 2, sort, na.last=FALSE)
		mat.dperm[,i] <- rowSums(tmpMat, na.rm=TRUE)
		tmple <- sum(!is.na(tmpMat))
		vecf <- vecf + rank(-c(as.vector(tmpMat), statnotna), ties.method="first", 
			na.last=NA)[tmple + (1:n.notna)]
	}
	d.bar <- rowSums(mat.dperm, na.rm=TRUE) / B
	vec.false <- rep(NA, n.snp)
	vec.false[ids.in] <- vecf/B
	msg <- "SAM for Trio Data (Based on MAX gTDT)\n\n"
	structure(list(d=stat, d.bar=d.bar, vec.false=vec.false, p.value=vec.false/n.notna, 
		s=numeric(0), s0=numeric(0), mat.samp=matrix(numeric(0)), msg=msg, fold=numeric(0)))
}


noBothHetMat <- function(a12, a34, perm){
	mat <- matrix(NA, length(a12), perm)
	tmpids <- a34==0
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
		mat[tmpids & a12==i,] <- pmax(add, dom, na.rm=TRUE)
	}
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
			mat[a12==i & a34==j,] <- pmax(add, dom, rec, na.rm=TRUE)
		}
	}
	mat
}

bothHetMat <- function(a12, a34, a567, perm, denom){
	uni.val <- unique(a567)
	mat <- matrix(NA, length(a567), perm)
	for(i in uni.val){
		tmpids <- which(a567 == i)
		matA567 <- rmultinom(perm, i, c(0.25, 0.5, 0.25))
		a5 <- matA567[1,]
		a6 <- matA567[2,]
		a7 <- matA567[3,]
		for(j in tmpids){
			a2 <- rbinom(perm, a12[j], 0.5)
			a1 <- a12[j] - a2
			a4 <- rbinom(perm, a34[j], 0.5)
			a3 <- a34[j] - a4
			h <- a2 + a4 + a6 + 2*a7
			beta <- logit(h/denom[j])
			v <- denom[j] / ((denom[j] - h) * h)
			add <- beta * beta / v
			d <- a6 + a7
			h <- (1/3*a1 - a2 + a5 - 1/3*d) / (2 * (a1+a5))
			tmp <- (a2 + d) / (3 * (a1+a5)) + h*h
			or <- sqrt(tmp) - h
			beta <- log(or)
			v <- a12[j] * or/(or+1)^2 + i * or/(3*(or+1/3)^2)
			#v <- 1/tmp
			dom <- beta * beta * v
			d <- a5 + a6
			h <- (3*a3 - a4 + d - 3*a7) / (2*(a3+d))
			tmp <- 3 * (a4+a7) / (a3+d) + h*h
			or <- sqrt(tmp) - h
			beta <- log(or)
			v <- a34[j] * or/(or+1)^2 + 3 * i * or/(or+3)^2
			rec <- beta * beta * v
			mat[j,] <- pmax(add, dom, rec, na.rm=TRUE)
		}
	}
	mat
}


gtdt.null.samp <- function(geno, type, stat, size=10, B=1000, chunk=100, rand=NA){
	ids.in <- !is.na(stat)
	statnotna <- stat[ids.in]
	geno <- geno[,ids.in]
	n.notna <- ncol(geno)
	int <- unique(c(seq.int(1, n.notna, size), n.notna + 1))
	matA <- matrix(NA, n.notna, 3 + as.numeric(type=="additive"))
	for(i in 1:(length(int) - 1))
		matA[int[i]:(int[i+1]-1),] <- getNullStats(geno[,int[i]:(int[i+1]-1), drop=FALSE], type)
	vec.perm <- diff(unique(c(seq(1, B, chunk), B+1)))
	n.perm <- length(vec.perm)
	ids.het <- matA[,3] == 0
	nullfun <- match.fun(paste("gtdt.null", substring(type, 1, 3), sep="."))
	mat.dperm <- matrix(0, n.notna, n.perm)
	d.rank <- rank(-statnotna, ties.method="first")
	vecf <- -d.rank * n.perm
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:n.perm){
		tmpmat <- nullfun(matA, vec.perm[i], ids.het)
		tmpmat <- apply(tmpmat, 2, sort, na.last=FALSE)
		mat.dperm[,i] <- rowSums(tmpmat, na.rm=TRUE)
		tmple <- sum(!is.na(tmpmat))
		vecf <- vecf + rank(-c(as.vector(tmpmat), statnotna), ties.method="first",
			na.last=NA)[tmple + (1:n.notna)]
	}		
	d.bar <- rowSums(mat.dperm, na.rm=TRUE) / B
	vec.false <- rep(NA, length(stat))
	vec.false[ids.in] <- vecf / B
	structure(list(d.bar=d.bar, vec.false=vec.false, p.value=vec.false/n.notna,
		mat.samp=matrix(numeric(0))))
}	


getNullStats <- function(mat, type){
	n.row <- nrow(mat)
	tmp <- mat[seq.int(1, n.row, 3),, drop=FALSE] + 10 * mat[seq.int(2, n.row, 3),, drop=FALSE]
	outmat <- matrix(NA, ncol(mat), 4)
	outmat[,3:4] <- colSums(tmp==11, na.rm=TRUE)
	if(type != "recessive")
		outmat[,1] <- colSums(tmp==1, na.rm=TRUE) + colSums(tmp==10, na.rm=TRUE)
	if(type != "dominant")
		outmat[,2] <- colSums(tmp==12, na.rm=TRUE) + colSums(tmp==21, na.rm=TRUE)
	if(type != "additive")
		return(outmat[,1:3])
	outmat[,4] <- rowSums(outmat, na.rm=TRUE)
	outmat
}
  
gtdt.null.add <- function(matA, perm, idsHet){
	mat <- matrix(NA, nrow(matA), perm)
	tmpids <- idsHet & matA[,2]==0
	tmpa12 <- matA[tmpids, 1]
	for(i in unique(tmpa12)){
		a2 <- rbinom(perm, i, 0.5)
		beta <- logit(a2/i)
		v <- i / ((i-a2) * a2)
		mat[tmpids & matA[,1]==i,] <- beta*beta/v
	}
	tmpa12 <- matA[idsHet & matA[,2]!=0, 1]
	tmpa34 <- matA[idsHet & matA[,2]!=0, 2]
	for(i in unique(tmpa12)){
		a2 <- rbinom(perm, i, 0.5)
		for(j in tmpa34[tmpa12==i]){
			a4 <- rbinom(perm, j, 0.5)
			num <- a2 + a4
			denom <- i + j
			beta <- logit(num/denom)
			v <- denom/((denom-num)*num)
			mat[idsHet & matA[,1]==i & matA[,2]==j,] <- beta*beta/v
		}
	}
	matA <- matA[!idsHet,]
	mat2 <- matrix(NA, nrow(matA), perm)
	uni.val <- unique(matA[, 3])
	for(i in uni.val){
		tmpids <- which(matA[,3]==i)
		matMulti <- rmultinom(perm, i, c(0.25, 0.5, 0.25))
		for(j in tmpids){
			a2 <- rbinom(perm, matA[j,1], 0.5)
			a4 <- rbinom(perm, matA[j,2], 0.5)
			num <- a2 + a4 + matMulti[2,] + 2*matMulti[3,]
			beta <- logit(num/matA[j,4])
			v <- matA[j,4] / ((matA[j,4] - num) * num)
			mat2[j,] <- beta * beta / v
		}
	}
	mat[!idsHet,] <- mat2
	mat
}

gtdt.null.dom <- function(matA, perm, idsHet){
	mat <- matrix(NA, nrow(matA), perm)
	tmpids <- idsHet & matA[,1]>0
	tmpa12 <- matA[tmpids, 1]
	for(i in unique(tmpa12)){
		d2 <- rbinom(perm, i, 0.5)
		d1 <- i - d2
		h <- (1/3*d1 - d2) / (2*d1)
		or <- sqrt(d2 / (3*d1) + h*h) - h
		beta <- log(or)
		v <- (or+1)^2 / (i * or)
		mat[tmpids & matA[,1]==i,] <- beta*beta/v
	}
	matA <- matA[!idsHet,]
	mat2 <- matrix(NA, nrow(matA), perm)
	uni.val <- unique(matA[,3])
	for(i in uni.val){
		tmpids <- which(matA[,3]==i)
		d4 <- rbinom(perm, i, 0.75)
		d3 <- i - d4
		for(j in tmpids){
			d2 <- rbinom(perm, matA[j,1], 0.5)
			d1 <- matA[j,1] - d2
			h <- (1/3 * d1 - d2 + d3 - 1/3 * d4) / (2 * (d1 + d3))
			or <- sqrt((d2 + d4) / (3 * (d1 + d3)) + h * h) - h
			beta <- log(or)
			v <- matA[j,1] * or/(or+1)^2 + i * or / (3 * (or + 1/3)^2)
			mat2[j,] <- beta*beta*v
		}
	}
	mat[!idsHet,] <- mat2
	mat
}

gtdt.null.rec <- function(matA, perm, idsHet){
	mat <- matrix(NA, nrow(matA), perm)
	tmpids <- idsHet & matA[,2] > 0
	tmpa34 <- matA[tmpids, 2]
	for(i in unique(tmpa34)){
		r2 <- rbinom(perm, i, 0.5)
		r1 <- i - r2
		h <- (3*r1 - r2) / (2*r1)
		or <- sqrt(3*r2 / r1 + h*h) - h
		beta <- log(or)
		v <- (or + 1)^2/(i * or)
		mat[tmpids & matA[,2]==i,] <- beta*beta/v
	}
	matA <- matA[!idsHet,]
	mat2 <- matrix(NA, nrow(matA), perm)
	uni.val <- unique(matA[,3])
	for(i in uni.val){
		tmpids <- which(matA[,3]==i)
		r4 <- rbinom(perm, i, 0.25)
		r3 <- i - r4
		for(j in tmpids){
			r2 <- rbinom(perm, matA[j,2], 0.5)
			r1 <- matA[j,2] - r2
			h <- (3*r1 - r2 + r3 - 3*r4) / (2 * (r1 + r3))
			or <- sqrt(3 * (r2+r4) / (r1+r3) + h*h) - h
			beta <- log(or)
			v <- matA[j,2] * or/(or + 1)^2 + 3*i * or/(or + 3)^2
			mat2[j,] <- beta*beta*v
		}
	}
	mat[!idsHet,] <- mat2
	mat
}
	



colTDTebam <- function(mat.snp, model=c("additive", "dominant", "recessive", "max"), approx=NULL,
		B=1000, size=10, chunk=100, n.interval=NULL, df.ratio=3, df.dens=3, 
		knots.mode=TRUE, type.nclass=c("wand", "FD", "scott"), fast=FALSE, rand=NA){
	requireNamespace("siggenes")
	if(!is.matrix(mat.snp))
		stop("mat.snp must be a matrix.")
	if(nrow(mat.snp) %%3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in  genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0:2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	type <- match.arg(model)
	if(is.null(approx))
		approx <- type != "max"
	if(!is.logical(approx))
		stop("approx must be a logic variable.")		
	ebam(mat.snp, numeric(ncol(mat.snp)), method=gtdt.ebam, model=type, delta=0.9, gene.names=NULL, approx=approx,
		B=B, size=size, chunk=chunk, n.interval=n.interval, df.ratio=df.ratio, df.dens=df.dens, 
		knots.mode=knots.mode, type.nclass=type.nclass, fast=fast, rand=rand)
}


gtdt.ebam <- function(data, cl, model, approx=TRUE, B=1000, size=10, chunk=100, n.interval=NULL, df.ratio=3, 
		df.dens=3, knots.mode=TRUE, type.nclass=c("wand", "FD", "scott"), fast=FALSE, rand=NA){
	cl <- model
	type.nclass <- match.arg(type.nclass)
	if(cl=="max"){	
		if(approx)
			stop("If a MAX-statistic should be used, then approx must be set to FALSE.")
		return(maxStatEBAM(data, size=size, B=B, chunk=chunk, n.interval=n.interval,
			df.ratio=df.ratio, type.nclass=type.nclass, fast=fast,rand=rand))
	}
	out <- fastTDT(data, cl, size=size)
	cap <- gsub("\\b(\\w)", "\\U\\1", cl, perl=TRUE)
        msg <- paste("EBAM for Trio Data (Based on", cap, "gTDT)\n\n")
	ids.notna <- !is.na(out$stat)
	if(any(!ids.notna))
		warning("For ", sum(!ids.notna), " SNPs, the gTDT statistic cannot be computed.\n",
			"These SNPs are removed from the analysis and the results.")
	statnotna <- out$stat[ids.notna]
	names(statnotna) <- colnames(data[,ids.notna])
	if(approx)
		return(gtdt.null.approx2(statnotna, out$pval[ids.notna], n.interval=n.interval, df.dens=df.dens,
			knots.mode=knots.mode, type.nclass=type.nclass, msg=msg))
	if(is.null(n.interval)){
		tmpfun <- match.fun(paste("nclass", type.nclass, sep="."))
		n.interval <- max(tmpfun(statnotna), 139)
	}
	if(n.interval < 20)
		stop("n.interval should be at least 20.")
	interval <- seq(floor(100 * min(statnotna))/100, ceiling(100*max(statnotna))/100, length=n.interval+1)
	center <- (interval[2] - interval[1])/2 + interval[-length(interval)]
	success <- tabulate(cut(statnotna, interval, include.lowest=TRUE), n.interval)
	fail.out <- compFailTDT(data[,ids.notna], statnotna, interval, cl, B=B, size=size, chunk=chunk, 
		fast=fast, rand=rand)
	ratio <- estimateRatioTDT(center, success, fail.out$vec.fail, statnotna, df.ratio=df.ratio)
	if(fast)
		return(list(z=statnotna, ratio=ratio, success=success, failure=fail.out$vec.fail,
			center=center, mat.samp=matrix(rep(NA, B)), msg=msg))
	structure(list(z=statnotna, ratio=ratio, vec.pos=fail.out$vec.pos/B, vec.neg=fail.out$vec.neg/B, 
		mat.samp=matrix(numeric(0)), msg=msg))
}
	

gtdt.null.approx2 <- function(stat, pval, n.interval=NULL, df.dens=3, knots.mode=TRUE, type.nclass="wand",
		msg=""){
	requireNamespace("siggenes", quietly = TRUE)
	ids.notna <- !is.na(stat)
	statnotna <- stat[ids.notna]
	z.dens <- vec.neg <- rep(NA, length(stat))
	z.dens[ids.notna] <- denspr(statnotna, n.interval=n.interval, df=df.dens, knots.mode=knots.mode,
		type.nclass=type.nclass)$y
	vec.pos <- length(statnotna) * pval
	vec.neg[ids.notna] <- 0
	z.null <- dchisq(stat, 1)
	structure(list(z=stat, ratio=z.null/z.dens, vec.pos=vec.pos, vec.neg=vec.neg, msg=msg))
}

compFailTDT <- function(geno, z, interval, type, B=1000, size=10, chunk=100, fast=FALSE, rand=NA){
	n.snp <- length(z)
	int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
	matA <- matrix(NA, n.snp, 3 + as.numeric(type=="additive"))
	for(i in 1:(length(int)-1))
		matA[int[i]:(int[i+1]-1),] <- getNullStats(geno[,int[i]:(int[i+1]-1), drop=FALSE], type)
	vec.perm <- diff(unique(c(seq(1, B, chunk), B+1)))
	n.perm <- length(vec.perm)
	ids.het <- matA[,3] == 0
	nullfun <- match.fun(paste("gtdt.null", substring(type, 1, 3), sep="."))
	z.min <- min(z)
	z.max <- max(z)
	n.int <- length(interval) - 1
	vec.fail <- numeric(n.int)
	if(!fast){
		vec.neg <- numeric(n.snp)
		vec.pos <- -n.perm * rank(-z, ties.method="first")
	}
	else
		vec.pos <- vec.neg <- NULL
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:n.perm){
		tmpmat <- nullfun(matA, vec.perm[i], ids.het)
		z.perm <- tmpmat[!is.na(tmpmat)]
		z.perm[z.perm < z.min] <- z.min
		z.perm[z.perm > z.max] <- z.max
		tmpbin <- cut(z.perm, interval, include.lowest=TRUE)
		vec.fail <- vec.fail + tabulate(tmpbin, n.int)
		if(!fast)
			vec.pos <- vec.pos + rank(-c(z.perm, z), ties.method="first")[length(z.perm) + (1:n.snp)]
	}
	structure(list(vec.fail=vec.fail, vec.pos=vec.pos, vec.neg=vec.neg))
}

estimateRatioTDT <- function(center, succ, fail, stat, df.ratio=3){
	# requireNamespace("splines")
	n.obs <- succ + fail
	ids <- n.obs > 0
	tmpmat <- ns.out <- splines::ns(center[ids], df.ratio)
	class(tmpmat) <- "matrix"
	dat <- data.frame(s=succ[ids], f=fail[ids], tmpmat)
	glm.out <- glm(cbind(s, f) ~ ., data=dat, family=binomial)
	pred <- predict(ns.out, stat)
	class(pred) <- "matrix"
	probs <- predict(glm.out, data.frame(pred), type="response")
	B <- sum(fail) / sum(succ)
	(1-probs) / (B * probs)
}


maxStatEBAM <- function(geno, size=10, B=1000, chunk=100, n.interval=NULL, df.ratio=3, type.nclass="wand",
		fast=FALSE, rand=NA){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- numeric(n.snp)
	matA <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1)){
		tmpMat <- getStatsChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		stat[int[i]:(int[i+1]-1)] <- tmpMat[,4]
		matA[int[i]:(int[i+1]-1), ] <- tmpMat[,c(1:3,5)]
	}
	vec.perm <- diff(unique(c(seq.int(1, B, chunk), B+1)))
	le.perm <- length(vec.perm)
	ids.in <- !is.na(stat)
	matA <- matA[ids.in,]
	stat <- stat[ids.in]
	names(stat) <- colnames(geno)[ids.in]
	if(any(!ids.in))
		warning("For ", sum(!ids.in), " SNPs, the gTDT statistic cannot be computed.\n",
			"These SNPs are removed from the analysis and the results.")
	if(is.null(n.interval)){
		tmpfun <- match.fun(paste("nclass", type.nclass, sep="."))
		n.interval <- max(tmpfun(stat), 139)
	}
	if(n.interval < 20)
		stop("n.interval should be at least 20.")
	interval <- seq(floor(100 * min(stat))/100, ceiling(100*max(stat))/100, length=n.interval+1)
	center <- (interval[2] - interval[1])/2 + interval[-length(interval)]
	success <- tabulate(cut(stat, interval, include.lowest=TRUE), n.interval)
	fail.out <- compFailTDTmax(geno[,ids.in], stat, matA, interval, B=B, size=size, chunk=chunk,
		fast=fast, rand=rand)
	ratio <- estimateRatioTDT(center, success, fail.out$vec.fail, stat, df.ratio=df.ratio)
	msg <- "EBAM for Trio Data (Based on MAX gTDT)\n\n"
	if(fast)
		out <- list(z=stat, ratio=ratio, success=success, failure=fail.out$vec.fail, center=center,
			mat.samp=matrix(rep(NA,B)), msg=msg)
	else
		out <- list(z=stat, ratio=ratio, vec.pos=fail.out$vec.pos/B, vec.neg=fail.out$vec.neg/B,
			mat.samp=matrix(numeric(0)), msg=msg)
	out
}



compFailTDTmax <- function(geno, stat, matA, interval, B=1000, size=10, chunk=100, fast=FALSE, rand=NA){
	n.snp <- length(stat)
	vec.perm <- diff(unique(c(seq(1, B, chunk), B+1)))
	n.perm <- length(vec.perm)
	ids.het <- matA[,3]==0
	z.min <- min(stat)
	z.max <- max(stat)
	n.int <- length(interval) - 1
	vec.fail <- numeric(n.int)
	if(!fast){
		vec.neg <- numeric(n.snp)
		vec.pos <- -n.perm * rank(-stat, ties.method="first")
	}
	else
		vec.pos <- vec.neg <- NULL
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:n.perm){
		tmpmat <- matrix(NA, n.snp, vec.perm[i])
		tmpmat[ids.het,] <- noBothHetMat(matA[ids.het, 1], matA[ids.het, 2], vec.perm[i])
		tmpmat[!ids.het,] <- bothHetMat(matA[!ids.het, 1], matA[!ids.het, 2], matA[!ids.het, 3],
			vec.perm[i], matA[!ids.het, 4])
		z.perm <- tmpmat[!is.na(tmpmat)]
		z.perm[z.perm < z.min] <- z.min
		z.perm[z.perm > z.max] <- z.max
		tmpbin <- cut(z.perm, interval, include.lowest=TRUE)
		vec.fail <- vec.fail + tabulate(tmpbin, n.int)
		if(!fast)
			vec.pos <- vec.pos + rank(-c(z.perm, stat), ties.method="first")[length(z.perm)+(1:n.snp)]
	}
	structure(list(vec.fail=vec.fail, vec.pos=vec.pos, vec.neg=vec.neg))
}

