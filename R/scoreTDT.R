scoreTDT <- function(mat.snp, model=c("additive", "dominant", "recessive"), size=20){
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
	type <- match.arg(model)	
	if(size < 1)
		stop("size should be at least 1.")
	fun <- match.fun(switch(type, "additive"=scoreTDTsplit, "dominant"=scoreTDTdomSplit, 
		"recessive"=scoreTDTrecSplit))
	fun(mat.snp, size=size)
}

print.scoreTDT <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Score=x$score, Statistic=x$stat, "p-value"=pval, check.names=FALSE, stringsAsFactors=FALSE)
	if(x$ia)
		cat("   Score Test for Two-Way Interactions\n\n")
	else
		cat("   Score Test for Individual SNPs\n\n")
	cat("Model Type: ", switch(x$type, "additive"="Additive", "dominant"="Dominant",
		"recessive"="Recessive"), "\n\n", sep="")
	if(length(x$score) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, ifelse(x$ia, "SNP interactions:\n", "SNPs:\n"))
	}
	print(format(out, digits=digits))
}


scoreTDTsplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
	num <- denom <- rep.int(NA, n.snp)
	for(i in 1:(length(int)-1)){
		tmp <- scoreTDTchunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
		num[int[i]:(int[i+1]-1)] <- tmp$num
		denom[int[i]:(int[i+1]-1)] <- tmp$denom
	}
	beta <- num - 0.5 * denom
	v <- denom/4
	stat <- beta*beta/v
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	names(beta) <- names(stat) <- names(pval) <- colnames(geno)
	out <- list(score=beta, info=v, stat=stat, pval=pval, ia=FALSE, type="additive")
	class(out) <- "scoreTDT"
	out
}

scoreTDTchunk <- function(geno){
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
	denom <- a1 + a2 + a3 + a4 + 2*a567
	return(list(num=num, denom=denom))
}	
	

scoreTDTdomSplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	dmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		dmat[int[i]:(int[i+1]-1),] <- scoreTDTdomChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
	beta <- -0.5*dmat[,1] + 0.5*dmat[,2] - 0.75*dmat[,3] + 0.25*dmat[,4]
	v <- (dmat[,1]+dmat[,2])/4 + 3/16*(dmat[,3] + dmat[,4]) 
	stat <- beta*beta/v
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	names(beta) <- names(stat) <- names(pval) <- colnames(geno)
	out <- list(score=beta, info=v, stat=stat, pval=pval, ia=FALSE, type="dominant")
	class(out) <- "scoreTDT"
	out
}

scoreTDTdomChunk <- function(geno){
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


scoreTDTrecSplit <- function(geno, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	rmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		rmat[int[i]:(int[i+1]-1),] <- scoreTDTrecChunk(geno[,int[i]:(int[i+1]-1), drop=FALSE])
	beta <- -0.5*rmat[,1] + 0.5*rmat[,2] - 0.25*rmat[,3] + 0.75*rmat[,4]
	v <- 0.25*(rmat[,1]+rmat[,2]) + 3/16 * (rmat[,3]+rmat[,4])
	stat <- beta*beta/v
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	names(beta) <- names(stat) <- names(pval) <- colnames(geno)
	out <- list(score=beta, info=v, stat=stat, pval=pval, ia=FALSE, type="recessive")
	class(out) <- "scoreTDT"
	out
}

scoreTDTrecChunk <- function(geno){
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




scoreGxE <- function(mat.snp, env, model=c("additive", "dominant", "recessive"), size=20, famid=NULL){
	if (!is.matrix(mat.snp)) 
        	stop("mat.snp has to be a matrix.")
    	if (nrow(mat.snp)%%3 != 0) 
        	stop("mat.snp does not seem to contain trio data, as its number of rows is\n", 
            		"not dividable by 3.")
    	if (is.null(rownames(mat.snp))) 
        	stop("mat.snp does not seem to be a matrix in genotype format,\n", 
            		"as the names of the rows are missing.")
    	if (any(!mat.snp %in% c(0, 1, 2, NA))) 
        	stop("The values in mat.snp must be 0, 1, and 2.")
	if(size<1)
		stop("size should be at least 1.")
	if(length(env)!=nrow(mat.snp)/3)
		stop("The length of env must be equal to the number of trios in mat.snp.")
	if(!is.null(famid))
		env <- reorderEnv(env, famid, rownames(mat.snp))
	if(any(is.na(env))){
		tmpEnv <- rep(env, e=3)
		mat.snp <- mat.snp[!is.na(tmpEnv),]
		env <- env[!is.na(env)]
	}
	if(any(!env %in% 0:1))
		stop("The values in env must be 0 and 1.")
	if(sum(env)<5 | sum(1-env)<5)
		stop("Each of the two groups specified by env must contain at least 5 trios.") 
	type <- match.arg(model)
	env2 <- rep(env, e=3)
	fun <- match.fun(switch(type, additive=scoreGxEsplit, dominant=scoreGxEdomSplit,
		recessive=scoreGxErecSplit))
	tmp1 <- fun(mat.snp, env2==0, size=size)
	tmp2 <- fun(mat.snp, env2==1, size=size)
	beta <- cbind(SNP=tmp1$coef, GxE=tmp2$coef-tmp1$coef)
	v <- cbind(SNP=tmp1$v, GXE=tmp2$v + tmp1$v)
	stat <- beta*beta/v
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	rownames(beta) <- rownames(stat) <- rownames(pval) <- colnames(mat.snp)
	out <- list(score=beta, info=v, stat=stat, pval=pval, type=type)
	class(out) <- "scoreGxE"
	out
}

print.scoreGxE <- function(x, top=5, digits=4, onlyGxE=FALSE, ...){
	if(!onlyGxE){
		pvalG <- format.pval(x$pval[,1], digits=digits)
		outG <- data.frame(Score=x$score[,1], Statistic=x$stat[,1], "p-value"=pvalG, check.names=FALSE,
			stringsAsFactors=FALSE)
	}
	pvalGE <- format.pval(x$pval[,2], digits=digits)
	outGE <- data.frame(Score=x$score[,2], Statistic=x$stat[,2], "p-value"=pvalGE, check.names=FALSE,
		stringsAsFactors=FALSE)
	cat("   Score Test for GxE Interactions with Binary E\n\n", "Model Type: ",
		switch(x$type, "additive"="Additive", "dominant"="Dominant", "recessive"="Recessive"),
		"\n\n", sep="")
	if(nrow(x$score) > top){
		ord <- order(x$pval[,2])[1:top]
		if(!onlyGxE)
			outG <- outG[ord,]
		outGE <- outGE[ord,]
		cat("Top", top, "GxE Interactions:\n")
	}
	else
		cat("Statistics for the GxE Interactions:\n")
	print(format(outGE, digits=digits))
	if(!onlyGxE){
		cat("\n\n", "Statistics for the SNPs in the Corresponding GxE Models:\n", sep="")
		print(format(outG, digits=digits))
	}
}
	

scoreGxEsplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
	num <- denom <- rep.int(NA, n.snp)
	for(i in 1:(length(int)-1)){
		tmp <- scoreTDTchunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
		num[int[i]:(int[i+1]-1)] <- tmp$num
		denom[int[i]:(int[i+1]-1)] <- tmp$denom
	}
	beta <- num - 0.5 * denom
	v <- denom/4
	list(coef=beta, v=v)
}

scoreGxEdomSplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	dmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		dmat[int[i]:(int[i+1]-1),] <- scoreTDTdomChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	beta <- -0.5*dmat[,1] + 0.5*dmat[,2] - 0.75*dmat[,3] + 0.25*dmat[,4]
	v <- (dmat[,1]+dmat[,2])/4 + 3/16*(dmat[,3] + dmat[,4]) 
	list(coef=beta, v=v)
}

scoreGxErecSplit <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	rmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		rmat[int[i]:(int[i+1]-1),] <- scoreTDTrecChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	beta <- -0.5*rmat[,1] + 0.5*rmat[,2] - 0.25*rmat[,3] + 0.75*rmat[,4]
	v <- 0.25*(rmat[,1]+rmat[,2]) + 3/16 * (rmat[,3]+rmat[,4])
	list(coef=beta, v=v)
}



scoreMaxStat <- function(mat.snp, size=20){
	if(!is.matrix(mat.snp))
		stop("mat.snp must be a matrix.")
	if(nrow(mat.snp) %%3 != 0)
		stop("mat.snp does not seem to contain trio data, as is number of rows is\n",
			"not dividable by 3.")
	if(is.null(rownames(mat.snp)))
		stop("mat.snp does not seem to be a matrix in genotype format,\n",
			"as the row names are missing.")
	if(any(!mat.snp %in% c(0,1,2,NA)))
		stop("The values in mat.snp must be 0, 1, or 2.")
	n.snp <- ncol(mat.snp)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- rep.int(NA, n.snp)
	mat.stat <- matrix(NA, n.snp, 3)
	for(i in 1:(length(int)-1)){
		tmp <- getScoresChunk(mat.snp[,int[i]:(int[i+1]-1), drop=FALSE])
		stat[int[i]:(int[i+1]-1)] <- tmp$stat
		mat.stat[int[i]:(int[i+1]-1),] <- tmp$mat.stat
	}
	names(stat) <- colnames(mat.snp)
	colnames(mat.stat) <- c("Additive", "Dominant", "Recessive")
	out <- list(stat=stat, mat.stat=mat.stat)
	class(out) <- "maxScoreTrio"
	out
}

print.maxScoreTrio <- function(x, top=5, digits=4, ...){
	out <- data.frame("Max-Statistic"=x$stat, x$mat.stat, check.names=FALSE)
	cat("      Maximum Score Statistic\n\n")
	if(length(x$stat) > top){
		ord <- order(x$stat, decreasing=TRUE)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNPs:\n")
	}
	print(format(out, digits=digits))
}

getScoresChunk <- function(geno){
	n.row <- nrow(geno)
    	dad <- geno[seq.int(1, n.row, 3), , drop = FALSE]
    	mom <- geno[seq.int(2, n.row, 3), , drop = FALSE]
    	kid <- geno[seq.int(3, n.row, 3), , drop = FALSE]
    	het1 <- (dad == 1) & (mom == 1)
	a567 <- colSums(het1 & !is.na(kid), na.rm=TRUE)
    	a5 <- colSums(het1 & kid == 0, na.rm = TRUE)
    	a6 <- colSums(het1 & kid == 1, na.rm = TRUE)
    	a7 <- colSums(het1 & kid == 2, na.rm = TRUE)
    	mom <- mom + dad
    	het1 <- mom == 1
    	a1 <- colSums(het1 & kid == 0, na.rm = TRUE)
    	a2 <- colSums(het1 & kid == 1, na.rm = TRUE)
    	het1 <- mom == 3
    	a3 <- colSums(het1 & kid == 1, na.rm = TRUE)
    	a4 <- colSums(het1 & kid == 2, na.rm = TRUE)
	mat.stat <- matrix(0, ncol(geno), 3)
	num <- a2 + a4 + a6 + 2*a7
	denom <- a1 + a2 + a3 + a4 + 2 * a567
	beta <- num - 0.5 * denom
	v <- denom/4
	mat.stat[,1] <- beta*beta/v
	num <- a6+a7
	beta <- -0.5*a1 + 0.5*a2 - 0.75*a5 + 0.25*num
	v <- (a1+a2)/4 + 3/16*a567 
	mat.stat[,2] <- beta*beta/v
	beta <- -0.5*a3 + 0.5*a4 - 0.25*(a5+a6) + 0.75*a7
	v <- 0.25*(a3+a4) + 3/16 * a567
	mat.stat[,3] <- beta*beta/v
	stat <- pmax(mat.stat[,1], mat.stat[,2], mat.stat[,3], na.rm=TRUE)
	list(stat=stat, mat.stat=mat.stat)
}



scoreGxG <- function(mat.snp, model=c("additive", "dominant", "recessive"), genes=NULL, size=20){
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
	n.combs <- nrow(combs)
	int <- unique(c(seq.int(1, n.combs, size), n.combs + 1))
	if(type == "Add")
		betaFun <- scoreBetaAdd
	if(type == "Dom")
		betaFun <- scoreBetaDom
	if(type == "Rec")
		betaFun <- scoreBetaRec
	beta <- v <- rep.int(NA, n.combs)
	for(i in 1:(length(int)-1)){
		tmp <- betaFun(mat.ids, kid, combs[int[i]:(int[i+1]-1),])
		beta[int[i]:(int[i+1]-1)] <- tmp$beta
		v[int[i]:(int[i+1]-1)] <- tmp$v
	}
	stat <- beta * beta / v
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(!is.null(genes))
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	names(beta) <- names(stat) <-  names(pval) <- paste(cn[combs[,1]], cn[combs[,2]], sep=" : ")
	out <- list(score=beta, info=v, stat=stat, pval=pval, genes=genes, type=modeltype, ia=TRUE)
	class(out) <- "scoreTDT"
	out
}



scoreBetaRec <- function(mat.ids, kid, combs){
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
	beta <- facBeta - mat.num[,6] - 0.25 * (mat.num[,1] + mat.num[,5]) - 0.125 * mat.num[,2] - 
		0.5 * mat.num[,3] - 1/16 * mat.num[,4]
	v <- derivExpit(1, 3) * (mat.num[,1] + mat.num[,5]) + derivExpit(1, 7) * mat.num[,2] +
		derivExpit(1, 1) * mat.num[,3] + derivExpit(1, 15) * mat.num[,4] 
	list(beta=beta, v=v)
}


scoreBetaDom <- function(mat.ids, kid, combs){
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
	beta <- facBeta - mat.num[,6] - 0.25*mat.num[,1] - 0.375*mat.num[,2] - 0.5*mat.num[,3] - 
		0.5625*mat.num[,4] - 0.75*mat.num[,5]
	v <- derivExpit(1, 3) * mat.num[,1] + derivExpit(1, 5/3) * mat.num[,2] + 
		derivExpit(1, 1) * mat.num[,3] + derivExpit(1, 7/9) * mat.num[,4] +
		derivExpit(1, 1/3) * mat.num[,5]
	list(beta=beta, v=v)
}

scoreBetaAdd <- function(mat.ids, kid, combs){
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
	beta <- facBeta - 0.25*mat.num[,1] - 0.75*mat.num[,2] - 0.5*mat.num[,3] - 0.5*mat.num[,4] -
		mat.num[,5] - 2.25*mat.num[,6] - 1.5*mat.num[,7] - 1.5*mat.num[,8] - 3*mat.num[,9] - 
		mat.num[,10] - mat.num[,11] - 2*mat.num[,12]
	v <- 0.1875*mat.num[,1] + 0.6875*mat.num[,2] + 0.5*mat.num[,3] + 0.25*mat.num[,4] + 
		mat.num[,5] + 1.1875*mat.num[,6] + 1.5*mat.num[,7] + 0.25*mat.num[,8] + mat.num[,9] +
		1.25*mat.num[,10] + 0.5*mat.num[,11] + 2*mat.num[,12]
	list(beta=beta, v=v)
}

