fastGxG <- function(mat.snp, model=c("additive", "dominant", "recessive"), genes=NULL, 
		interval=c(-10,10), tol=10^-8, maxiter=1000, size=20){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	modeltype <- match.arg(model) 
	type <- switch(modeltype, additive="Add", dominant="Dom", recessive="Rec")
	out <- getGenoComb(mat.snp, size=size, type=type)
	kid <- out$kid
	mat.ids <- out$mat.comb
	cn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:ncol(mat.snp), sep="") else colnames(mat.snp)
	if(is.null(genes))
		combs <- allCombs(ncol(kid))
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != ncol(mat.snp))
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
	}
	rm(mat.snp, out)
	n.combs <- nrow(combs)
	int <- unique(c(seq.int(1, n.combs, size), n.combs + 1))
	betaFun <- match.fun(paste("getBeta", type, sep=""))
	beta <- neg2ndInv <- rep.int(NA, n.combs)
	for(i in 1:(length(int)-1)){
		b.out <- betaFun(mat.ids, kid, combs[int[i]:(int[i+1]-1),], interval=interval, 
			tol=tol, maxiter=maxiter)
		beta[int[i]:(int[i+1]-1)] <- b.out$beta
		neg2ndInv[int[i]:(int[i+1]-1)] <- b.out$negInv2
	}
	stat <- beta * beta * neg2ndInv
	se <- neg2ndInv^-0.5
	se[is.infinite(se)] <- 0
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	names(beta) <- names(stat) <- names(pval) <- paste(cn[combs[,1]], cn[combs[,2]], sep=" : ")
	if(!is.null(genes))
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	out <- list(coef=beta, se=se, stat=stat, pval=pval, RR=exp(coef), lowerRR=lower, upperRR=upper, ia=TRUE,
		type=modeltype, add=FALSE, genes=genes, maf=NULL, matMAF=NULL)
	class(out) <- "colTDT"
	out
}

fastGxGrec <- function(mat.snp, genes=NULL, size=20){
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the names of the rows are missing.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	out <- getGenoComb(mat.snp, size=size, type="Rec")
	kid <- out$kid
	mat.ids <- out$mat.comb
	cn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:ncol(mat.snp), sep="") else colnames(mat.snp)
	rm(mat.snp, out)
	if(is.null(genes))
		combs <- allCombs(ncol(kid))
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != ncol(mat.snp))
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
	}
	n.combs <- nrow(combs)
	int <- unique(c(seq.int(1, n.combs, size), n.combs + 1))
	beta <- neg2ndInv <- rep.int(NA, n.combs)
	for(i in 1:(length(int)-1)){
		b.out <- getBetaRec2(mat.ids, kid, combs[int[i]:(int[i+1]-1),])
		beta[int[i]:(int[i+1]-1)] <- b.out$beta
		neg2ndInv[int[i]:(int[i+1]-1)] <- b.out$negInv2
	}
	stat <- beta*beta *neg2ndInv
	se <- neg2ndInv^-0.5
	se[is.infinite(se)] <- 0
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	names(beta) <- names(stat) <- names(pval) <- paste(cn[combs[,1]], cn[combs[,2]], sep=" : ")
	if(!is.null(genes))
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	out <- list(coef=beta, se=se, stat=stat, pval=pval, OR=exp(coef), lowerOR=lower, upperOR=upper, ia=TRUE,
		type="recessive", add=FALSE, genes=genes, maf=NULL, matMAF=NULL)
	class(out) <- "colTDT"
	out
}

getBetaRec <- function(mat.ids, kid, combs, interval=c(-10,10), tol=10^-8, maxiter=1000){
	n.combs <- nrow(combs)
	facBeta <- colSums(kid[,combs[,1]] * kid[,combs[,2]], na.rm=TRUE)
	mat.comb <- 10 * mat.ids[,combs[,1]] + mat.ids[,combs[,2]]
	mat.num <- matrix(0, n.combs, 6)
	mat.num[,1] <- colSums(mat.comb==11, na.rm=TRUE)
	mat.num[,2] <- colSums(mat.comb==21 | mat.comb==12, na.rm=TRUE)
	mat.num[,3] <- colSums(mat.comb==31 | mat.comb==13, na.rm=TRUE)
	mat.num[,4] <- colSums(mat.comb==22, na.rm=TRUE)
	mat.num[,5] <- colSums(mat.comb==32 | mat.comb==23, na.rm=TRUE)
	mat.num[,6] <- colSums(mat.comb==33, na.rm=TRUE)
	beta <- rep.int(NA, n.combs)
	for(i in 1:n.combs){
		uni.out <- try(uniroot(derivLLrec, interval=interval, n=mat.num[i,], k=facBeta[i], 
			tol=tol, maxiter=maxiter), silent=TRUE)
		if(!is(uni.out, "try-error"))
			beta[i] <- uni.out$root
	}
	beta[beta==min(interval)] <- 0
	or <- exp(beta)
	tmp <- derivExpit(or, 3) * (mat.num[,1] + mat.num[,5]) + derivExpit(or, 7) * mat.num[,2] +
		derivExpit(or, 1) * mat.num[,3] + derivExpit(or, 15) * mat.num[,4] 
	list(beta=beta, negInv2=tmp)
}

getBetaRec2 <- function(mat.ids, kid, combs){
	n.combs <- nrow(combs)
	matKid <- kid[,combs[,1]] * kid[,combs[,2]]
	mat.comb <- 10 * mat.ids[,combs[,1]] + mat.ids[,combs[,2]]
	mat.num <- matrix(0, n.combs, 5)
	mat.num[,1] <- colSums(mat.comb==11 | mat.comb==23 | mat.comb==32, na.rm=TRUE)
	mat.num[,2] <- colSums(mat.comb==12 | mat.comb==21, na.rm=TRUE)
	mat.num[,3] <- colSums(mat.comb==13 | mat.comb==31, na.rm=TRUE)
	mat.num[,4] <- colSums(mat.comb==22, na.rm=TRUE)
	mat.num[,5] <- colSums(matKid, na.rm=TRUE) - colSums(mat.comb==33, na.rm=TRUE)
	matFac <- matrix(c(23,19,25,11,-26, 127,63,171,31,-196, 105,45,315,21,-486), 5)
	abcde <- matrix(0, n.combs, 5)
	abcde[,1] <- colSums(matKid==0, na.rm=TRUE)
	abcde[,2:4] <- mat.num %*% matFac
	abcde[,5] <- -315*mat.num[,5]
	beta <- rep.int(NA, n.combs)
	if(any(abcde[,1]==0)){
		ids <- which(abcde[,1]!=0)
		rs <- rowSums(abcde==0)==5
		if(any(rs))
			beta[rs] <- 0
	}
	else
		ids <- 1:n.combs
	beta[ids] <- getLogRoot(abcde[ids,,drop=FALSE])
	or <- exp(beta)	
	tmp <- derivExpit(or, 3) * (mat.num[,1]) + derivExpit(or, 7) * mat.num[,2] +
		derivExpit(or, 1) * mat.num[,3] + derivExpit(or, 15) * mat.num[,4] 
	list(beta=beta, negInv2=tmp)
}	

getLogRoot <- function(mat){
	out <- poly4rootMat(mat)
	out[out<=10^-8] <- NA
	if(any(rowSums(!is.na(out))>1, na.rm=TRUE))
		stop("There is more than one root. Please contact the package maintainer.")
	log(rowSums(out, na.rm=TRUE))
}


getBetaDom <- function(mat.ids, kid, combs, interval=c(-10,10), tol=10^-8, maxiter=1000){
	n.combs <- nrow(combs)
	facBeta <- colSums(kid[,combs[,1]] * kid[,combs[,2]], na.rm=TRUE)
	mat.comb <- 10 * mat.ids[,combs[,1]] + mat.ids[,combs[,2]]
	mat.num <- matrix(0, n.combs, 6)
	mat.num[,1] <- colSums(mat.comb==11, na.rm=TRUE)
	mat.num[,2] <- colSums(mat.comb==21 | mat.comb==12, na.rm=TRUE)
	mat.num[,3] <- colSums(mat.comb==31 | mat.comb==13, na.rm=TRUE)
	mat.num[,4] <- colSums(mat.comb==22, na.rm=TRUE)
	mat.num[,5] <- colSums(mat.comb==32 | mat.comb==23, na.rm=TRUE)
	mat.num[,6] <- colSums(mat.comb==33, na.rm=TRUE)
	beta <- rep.int(NA, n.combs)
	for(i in 1:n.combs){
		uni.out <- try(uniroot(derivLLdom, interval=interval, n=mat.num[i,], k=facBeta[i],
			tol=tol, maxiter=maxiter), silent=TRUE)
		if(!is(uni.out, "try-error"))
			beta[i] <- uni.out$root
	}
	beta[beta==min(interval)] <- 0
	or <- exp(beta)
	tmp <- derivExpit(or, 3) * mat.num[,1] + derivExpit(or, 5/3) * mat.num[,2] + 
		derivExpit(or, 1) * mat.num[,3] + derivExpit(or, 7/9) * mat.num[,4] +
		derivExpit(or, 1/3) * mat.num[,5]
	list(beta=beta, negInv2=tmp)
}

getBetaAdd <- function(mat.ids, kid, combs, interval=c(-10,10), tol=10^-8, maxiter=1000){
	n.combs <- nrow(combs)
	facBeta <- colSums(kid[,combs[,1]] * kid[,combs[,2]], na.rm=TRUE)
	mat.comb <- 10 * mat.ids[,combs[,1]] + mat.ids[,combs[,2]]
	mat.num <- matrix(0, n.combs, 12)
	mat.num[,1] <- colSums(mat.comb==11, na.rm=TRUE)
	mat.num[,2] <- colSums(mat.comb==12 | mat.comb==21, na.rm=TRUE)
	mat.num[,3] <- colSums(mat.comb==13 | mat.comb==31, na.rm=TRUE)
	mat.num[,4] <- colSums(mat.comb==14 | mat.comb==41, na.rm=TRUE)
	mat.num[,5] <- colSums(mat.comb==15 | mat.comb==51, na.rm=TRUE)
	mat.num[,6] <- colSums(mat.comb==22, na.rm=TRUE)
	mat.num[,7] <- colSums(mat.comb==23 | mat.comb==32, na.rm=TRUE)
	mat.num[,8] <- colSums(mat.comb==24 | mat.comb==42, na.rm=TRUE)
	mat.num[,9] <- colSums(mat.comb==25 | mat.comb==52, na.rm=TRUE)
	mat.num[,10] <- colSums(mat.comb==33, na.rm=TRUE)
	mat.num[,11] <- colSums(mat.comb==34 | mat.comb==43, na.rm=TRUE)
	mat.num[,12] <- colSums(mat.comb==35 | mat.comb==53, na.rm=TRUE)
	facBeta <- facBeta - colSums(mat.comb==44, na.rm=TRUE) - 
		2 * colSums(mat.comb==45 | mat.comb==54, na.rm=TRUE) -
		4 * colSums(mat.comb==55, na.rm=TRUE)
	beta <- rep.int(NA, n.combs)
	for(i in 1:n.combs){
		uni.out <- try(uniroot(derivLLadd, interval=interval, n=mat.num[i,], k=facBeta[i],
			tol=tol, maxiter=maxiter), silent=TRUE)
		if(!is(uni.out, "try-error"))
			beta[i] <- uni.out$root
	}
	beta[beta==min(interval)] <- 0
	tmp <- valDeriv2LLadd(exp(beta), mat.num)
	list(beta=beta, negInv2=tmp)
}

derivLLrec <- function(beta, n, k){
	x <- exp(beta)
	k - n[6] - x/(3+x) * (n[1] + n[5]) - x/(7+x) * n[2] -
		x/(1+x) * n[3] - x/(15+x) * n[4]
}

derivLLdom <- function(beta, n, k){
	x <- exp(beta)
	k - n[6] - x/(3+x) * n[1] - x/(5/3 + x) * n[2] - x/(1+x) * n[3] - 
		x/(7/9 + x) * n[4] - x/(1/3 + x) * n[5]
}

derivLLadd <- function(beta, n, k){
	x <- exp(beta)
	x2 <- x*x
	x3 <- x2*x
	x4 <- x3*x
	k - x/(3+x) * n[1] - (x+2*x2)/(2+x+x2) * n[2] - (x+x2)/(2.5+x+0.5*x2) * n[3] -
		x/(1+x) * n[4] - 2*x2/(1+x2) * n[5] - (1+4*x+4*x3)/(1+2*x+x3) *n[6] -
		(x+3*x2+2*x4)/(1+x+1.5*x2+.5*x4) * n[7] - (1+2*x)/(1+x) * n[8] -
		(2+4*x2)/(1+x2) * n[9] - (x+2*x2+x4)/(1.75+x+x2+0.25*x4) * n[10] -
		(x+x2)/(0.5+x+0.5*x2) * n[11] - (x2+x4)/(0.25+0.5*x2+0.25*x4) * n[12]
}

valDeriv2LLadd <- function(x, n){
	mat <- matrix(NA, length(x), 12)
	mat[,1] <- derivExpit(x, 3)
	x2 <- x*x
	x3 <- x2*x
	h <- 2 + x + x2
	mat[,2] <- (2*x+8*x2+x3)/(h*h)
	h <- 2.5 + x + 0.5*x2
	mat[,3] <- (2.5*x + 5*x2 + 0.5*x3)/(h*h)
	h <- derivExpit(x, 1)
	mat[,4] <- h 
	mat[,8] <- h 
	h <- 4*x2/(1+x2)^2
	mat[,5] <- h 
	mat[,9] <- h 
	x4 <- x2*x2
	h <- 1+2*x+x3
	g <- 2*x + 9*x3 + 8*x4
	mat[,6] <- g/(h*h) 
	h <- 1+x+1.5*x2+0.5*x4
	g <- x + 3*x2+2*x4
	gp <- x + 6*x2 + 8*x4
	mat[,7] <- (gp*h-g*g)/(h*h) 
	g <- x + 2*x2 + x4
	h <- 1.75 + x + x2 + 0.25*x4
	gp <- x + 4*x2 + 4*x4
	mat[,10] <- (gp*h-g*g)/(h*h) 
	g <- 0.5*x+x2+0.5*x3
	h <- 0.5 + x + 0.5*x2
	mat[,11] <- g/(h*h) 
	g <- 0.5*x2 + x4 + 0.5 * x4*x2
	h <- 0.25 + 0.5*x2 + 0.25*x4
	mat[,12] <- g/(h*h) 
	rowSums(mat*n, na.rm=TRUE)
}
	

derivExpit <- function(expbeta, a) a*expbeta / (a + expbeta)^2


getGenoComb <- function(geno, size=20, type=""){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
	mat.comb <- mat.kid <- matrix(NA, nrow(geno)/3, n.snp)
	chunkFun <- match.fun(paste("chunkGenoComb", type, sep=""))
	for(i in 1:(length(int)-1)){
		tmp <- chunkFun(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		mat.comb[,int[i]:(int[i+1]-1)] <- tmp$mat.comb
		mat.kid[,int[i]:(int[i+1]-1)] <- tmp$kid
	}
	list(mat.comb=mat.comb, kid=mat.kid)
}

chunkGenoCombRec <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	mat <- matrix(NA, n.row/3, ncol(geno))
	het3 <- mom + dad == 3
	mat[het3 & kid>0] <- 1 
	mat[mom==1 & dad==1 & !is.na(kid)] <- 2 
	mat[mom==2 & dad==2 & kid==2] <- 3
	kid[is.na(mat)] <- NA
	kid <- (kid==2) * 1
	list(mat.comb=mat, kid=kid)
}

chunkGenoCombDom <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	mat <- matrix(NA, n.row/3, ncol(geno))
	het1 <- mom + dad
	mat[het1==1 & kid<2] <- 1
	mat[dad==1 & mom==1 & !is.na(kid)] <- 2
	mat[het1==3 & kid!=0] <- 3
	mat[het1==4 & kid==2] <- 3
	mat[het1==2 & mom!=1 & kid==1] <- 3
	kid[is.na(mat)] <- NA
	kid <- (kid>0) * 1
	list(mat.comb=mat, kid=kid)
}


chunkGenoCombAdd <- function(geno){
	n.row <- nrow(geno)
	dad <- geno[seq.int(1, n.row, 3),, drop=FALSE]
	mom <- geno[seq.int(2, n.row, 3),, drop=FALSE]
	kid <- geno[seq.int(3, n.row, 3),, drop=FALSE]
	mat <- matrix(NA, n.row/3, ncol(geno))
	mat[dad==1 & mom==1 & !is.na(kid)] <- 3
	mom <- mom + dad
	mat[mom==1 & kid<2] <- 1
	mat[mom==3 & kid>0] <- 2
	mat[mom==4 & kid==2] <- 5
	mat[mom==2 & dad!=1 & kid==1] <- 4
	kid[is.na(mat)] <- NA
	list(mat.comb=mat, kid=kid)
}
