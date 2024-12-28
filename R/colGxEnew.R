colGxE <- function(mat.snp, env = NULL, listNumber = NULL, model = c("additive", "dominant", "recessive"), 
		alpha1=1, size=50, addGandE=TRUE, whichLRT=c("both", "2df", "1df", "none"), add2df=TRUE, 
		addCov=FALSE, famid=NULL, unstructured=FALSE, step1Stats = c("alpha", "add", "only")){
	if(missing(mat.snp) & is.null(listNumber))
		stop("Either mat.snp or listNumber must be specified.")
	if(!missing(mat.snp) & !is.null(listNumber))
		stop("Only one of mat.snp and listNumber can be specified.")
	if(!is.null(listNumber)){
		isMatSNP <- FALSE
		checkMatN(listNumber, usedName = "listNumber")
		n.snp <- nrow(listNumber[[1]])
		cn <- rownames(listNumber[[1]])
	}
	else{
		isMatSNP <- TRUE
		checkMatSNP(mat.snp, size = size)
		n.trio <- nrow(mat.snp) / 3
		env2 <- checkEnv(env, n.trio, rownames(mat.snp), famid = famid)
		if(any(is.na(env))){
			mat.snp <- mat.snp[!is.na(env2),]
			env2 <- env2[!is.na(env2)]
		}
		n.snp <- ncol(mat.snp)
		cn <- colnames(mat.snp)		
	}
	if(unstructured){
		if(!isMatSNP)
			stop("A fully parameterized model (unstructured = TRUE) can currently only be fitted,\n",
				" when mat.snp is specified.")
		return(colGxEunstructured(mat.snp, env))
	}
	if(alpha1 > 1 | alpha1 <= 0)
		stop("alpha1 must be larger than 0 and smaller than or equal to 1.") 
	typeStep1 <- match.arg(step1Stats)	
	compStep1 <- alpha1 < 1 | typeStep1 != "alpha" 
	typeLRT <- match.arg(whichLRT)
	addLRT1 <- typeLRT %in% c("both", "1df") 
	addLRT2 <- typeLRT %in% c("both", "2df")
	type <- match.arg(model)
	if(compStep1){
		if(!isMatSNP)
			stop("The first step statistics of the Gauderman procedure can currently\n",
				" only be computed and used, when mat.snp is specified.")
		if(typeStep1 != "alpha" & alpha1 < 1){
			warning("Since step1Stats is set to \"", typeStep1, 
				"\", the specification of alpha1 is ignored.")
			alpha1 <- 1
		}	
		if(type=="additive")
			datStats1 <- EvsG4(mat.snp, env2, size = size)
		else
			datStats1 <- EvsG2(mat.snp, env2, size = size, type=type)
		rownames(datStats1) <- cn
		if(typeStep1 == "only")
			return(datStats1)
		step1pval <- datStats1$pval
		names(step1pval) <- cn
		if(all(datStats1$pval > alpha1)){
			out <- list(statsStep1 = datStats1, step1pval = step1pval, n.select=0)
			class(out) <- "colGxE2"
			return(out)
		}
		if(alpha1 < 1){
			mat.snp <- mat.snp[ , datStats1$pval <= alpha1, drop = FALSE]
			n.snp <- ncol(mat.snp)
			cn <- colnames(mat.snp)	
		}
	}   
	fun <- match.fun(switch(type, additive = fastGxEsplit, dominant = fastGxEdomSplit,
		recessive = fastGxErecSplit))
	if(isMatSNP){
		tmp1 <- fun(mat.snp, env2==0, size=size, addLRT1=addLRT1, addLRT2=addLRT2)
		tmp2 <- fun(mat.snp, env2==1, size=size, addLRT1=addLRT1, addLRT2=addLRT2)
	}
	else{
		tmp1 <- fun(matNumber = listNumber[[1]], addLRT1 = addLRT1, addLRT2 = addLRT2)
		tmp2 <- fun(matNumber = listNumber[[2]], addLRT1 = addLRT1, addLRT2 = addLRT2)
	}
	beta <- se <- matrix(NA, n.snp, 2 + addGandE)
	beta[,1] <- tmp1$beta
	beta[,2] <- tmp2$beta - tmp1$beta
	se[,1] <- sqrt(tmp1$var)
	tmpVar2 <- tmp2$var + tmp1$var
	se[,2] <- sqrt(tmpVar2)
	used <- cbind(NotExposed=tmp1$used, Exposed=tmp2$used)
	if(addGandE){
		beta[,3] <- tmp2$beta
		se[,3] <- sqrt(tmp2$var)
		colnames(beta) <- colnames(se) <- c("SNP", "GxE", "GandE")
	}
	else
		colnames(beta) <- colnames(se) <- c("SNP", "GxE")		
	rownames(beta) <- rownames(se) <- rownames(used) <- cn
	stat <- beta[,1:2, drop = FALSE] / se[,1:2, drop = FALSE]
	stat <- stat * stat
	if(typeLRT != "none")
		lrtFull <-  tmp1$lBeta + tmp2$lBeta
	if(addLRT2){
		lrtNull <- log(0.25) * (tmp1$used + tmp2$used)
		tmpStat <- -2 * (lrtNull - lrtFull)
		tmpP <- pchisq(tmpStat, 2, lower.tail=FALSE)
		lrt2df <- cbind(Statistic=tmpStat, pValue= tmpP)
	}
	if(addLRT1){
		funG <- match.fun(switch(type, additive = lrtGxE, dominant = lrtGxEdom, 
			recessive = lrtGxErec))
		lrtG <- funG(tmp1$mat + tmp2$mat)
		tmpStat <- -2 * (lrtG - lrtFull)
		tmpP <- pchisq(tmpStat, 1, lower.tail=FALSE)
		lrt1df <- cbind(Statistic=tmpStat, pValue=tmpP)
	}	
	if(add2df){
		tmpFac <- beta[,2]/tmp2$var
		tmpW <- beta[,1] * beta[,1] * tmpVar2 / (tmp1$var * tmp2$var)
		tmpStat <- tmpW + (2 * beta[,1] + beta[,2]) * tmpFac
		tmpP <- pchisq(tmpStat, 2, lower.tail=FALSE)
		wald2df <- cbind(Statistic=tmpStat, pValue= tmpP)
	}
	lower <- exp(beta - qnorm(0.975) * se)
	upper <- exp(beta + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	out <- list(coef=beta[,1:2, drop = FALSE], se=se[,1:2, drop = FALSE], stat=stat, pval=pval, RR=exp(beta), 
		lowerRR=lower, upperRR=upper, usedTrios=used, env=env, type=type, addGandE=addGandE, 
		addOther=c(addLRT2, add2df, addLRT1), n.select=n.snp, typeStep1 = typeStep1,
		alpha1 = alpha1)
	if(compStep1){
		out$statsStep1 <- datStats1
		out$step1pval <- step1pval
	}
	if(addCov)
		out$cov <- -tmp1$var
	if(addLRT2)
		out$lrt2df <- lrt2df
	if(add2df)
		out$wald2df <- wald2df
	if(addLRT1)
		out$lrt1df <- lrt1df
	class(out) <- "colGxE"
	out
}

checkMatSNP <- function(mat.snp, size = 50){
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
	if(size < 1)
		stop("size should be at least 1.")
}

checkEnv <- function(env, n.trio, rn, famid = NULL){
	if(is.null(env))
		stop("env has to be specified, when mat.snp is specified.")
	if(length(env) != n.trio)
		stop("The length of env must be equal to the number of trios in mat.snp.")
	if(!is.null(famid))
		env <- reorderEnv(env, famid, rn)
	if(any(!env %in% 0:1, na.rm = TRUE))
		stop("The values in env must be 0 and 1.")
	if(sum(env, na.rm = TRUE) < 5 | sum(1 - env, na.rm = TRUE) < 5)
		stop("Each of the two groups specified by env must contain at least 5 trios.")
	rep(env, each = 3)
}

checkMatN <- function(matNumber, usedName = "matNumber"){
	if(is.list(matNumber)){
		if(length(matNumber) != 2)
			stop(usedName, " does not seem to be the output of colNtrios.")
		lnames <- names(matNumber)
		if(is.null(lnames) | any(lnames != c("ntriosEnv0", "ntriosEnv1")))
			stop(usedName, " does not seem to be the output of colNtrios.")
		usedName <- paste("An element of", usedName)
		lapply(matNumber, checkMatN, usedName = usedName)
	}
	else{
		if(!is.matrix(matNumber))
			stop(usedName, " must be a matrix.") 
		n.col <- ncol(matNumber)
		if(!n.col %in% c(7, 10))
			stop(usedName, " must be the output of colNtrios.")
		cn <- c("G010", "G011", "G121", "G122", "G110", "G111", "G112")
		if(n.col == 10)
			cn <- c(cn, "G021", "G222", "G000")
		if(any(cn != colnames(matNumber)))
			stop(usedName, " does not seem to be the output of colNtrios.")	
	}
}


fastGxEsplit <- function(geno, env2, matNumber = NULL, size=50, addLRT1=TRUE, addLRT2=TRUE){
	if(!is.null(matNumber)){
		num <- rowSums(matNumber[ , c(2, 4, 6, 7, 7)])
		used <- rowSums(matNumber[ , 1:7])
		denom <- used + rowSums(matNumber[ , c(5:7)])
	}
	else{
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
		num <- denom <- used <- numeric(n.snp)
		for(i in 1:(length(int)-1)){
			tmp <- fastTDTchunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
			num[int[i]:(int[i+1]-1)] <- tmp$num
			denom[int[i]:(int[i+1]-1)] <- tmp$denom
			used[int[i]:(int[i+1]-1)] <- tmp$used
		}
	}
	beta <- logit(num/denom)
	if(any(is.infinite(beta)))
		beta[is.infinite(beta)] <- NA
	tmpFac <- denom/(denom-num)
	se <- tmpFac / num
	if(!addLRT1 & !addLRT2)
		return(list(beta=beta, var=se, used=used))
	oneHet <- 2 * used - denom
	lnTerm <- log(tmpFac) * denom
	lBeta <- beta * num - log(2) * oneHet - lnTerm
	if(addLRT1)
		return(list(beta=beta, var=se, used=used, lBeta=lBeta, mat=cbind(num, denom, used)))
	list(beta=beta, var=se, used=used, lBeta=lBeta)	
}
	


fastGxEdomSplit <- function(geno, env2, matNumber = NULL, size=50, addLRT1=TRUE, addLRT2=TRUE){
	if(!is.null(matNumber))
		dmat <- ntrios2Dom(matNumber, check = FALSE, quiet = TRUE)[ , 1:4]
	else{	
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
		dmat <- matrix(0, n.snp, 4)
		for(i in 1:(length(int)-1))
			dmat[int[i]:(int[i+1]-1),] <- fastTDTdomChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	}
	h <- (dmat[,1]/3 - dmat[,2] + dmat[,3] - dmat[,4]/3) / (2*(dmat[,1]+dmat[,3]))
	d24 <- dmat[,2] + dmat[,4]
	tmp <- d24/(3*(dmat[,1]+dmat[,3])) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	if(any(is.infinite(beta)))
		beta[is.infinite(beta)] <- NA
	d12 <- dmat[,2] + dmat[,1]
	d34 <- dmat[,3] + dmat[,4]
	tmp <- d12*or/(or+1)^2 + d34*or/(3*(or+1/3)^2)
	used <- rowSums(dmat)
	if(!addLRT1 & !addLRT2)
		return(list(beta=beta, var=1/tmp, used=used))
	lBeta <- d24 * beta - d12 * log(2 + 2*or) - d34 * log(1+3*or)
	if(addLRT1)
		return(list(beta=beta, var=1/tmp, used=used, lBeta=lBeta, mat=dmat))
	list(beta=beta, var=1/tmp, used=used, lBeta=lBeta)
}

fastGxErecSplit <- function(geno, env2, matNumber = NULL, size=50, addLRT1=TRUE, addLRT2=TRUE){
	if(!is.null(matNumber))
		rmat <- ntrios2Rec(matNumber, check = FALSE, quiet = TRUE)[ , 1:4]
	else{
		n.snp <- ncol(geno)
		int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
		rmat <- matrix(0, n.snp, 4)
		for(i in 1:(length(int)-1))
			rmat[int[i]:(int[i+1]-1),] <- fastTDTrecChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	}
	h <- (3*rmat[,1] - rmat[,2] + rmat[,3] - 3*rmat[,4]) / (2 * (rmat[,1]+rmat[,3]))
	r24 <- rmat[,2] + rmat[,4]
	tmp <- 3 * r24 / (rmat[,1] + rmat[,3]) + h*h
	or <- sqrt(tmp) - h	 
	beta <- log(or)
	if(any(is.infinite(beta)))
		beta[is.infinite(beta)] <- NA
	r12 <- rmat[,1] + rmat[,2]
	r34 <- rmat[,3] + rmat[,4]
	tmp <- r12 * or/(or+1)^2 + 3 * r34 * or / (or+3)^2
	used <- rowSums(rmat)
	if(!addLRT1 & !addLRT2)
		return(list(beta=beta,var=1/tmp, used=used))
	lBeta <- r24 * beta - r12 * log(2+2*or) - r34 * log(3+or)
	if(addLRT1)
		return(list(beta=beta, var=1/tmp, used=used,lBeta=lBeta, mat=rmat))
	list(beta=beta, var=1/tmp, used=used, lBeta=lBeta)	
}


EvsG4 <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	matStat <- matrix(nrow = n.snp, ncol = 3)
	for(i in 1:(length(int)-1))
		matStat[int[i]:(int[i+1]-1), ] <- EvsG4split(geno[,int[i]:(int[i+1]-1), drop=FALSE], 
			env2==1, size=size)
	colnames(matStat) <- c("coef", "se", "stat")
	pval <- pchisq(matStat[ , 3], 1, lower.tail=FALSE)
	data.frame(matStat, pval = pval)
}	

EvsG4split <- function(geno, env3, size=50){ 
	matA <- EvsG4chunk(geno)
	matS <- EvsG4chunk(geno[env3,,drop=FALSE], onlySum=TRUE)
	coef <- matrix(NA, ncol(geno), 2)
	for(i in 1:ncol(geno))
		coef[i,] <- optim(c(0,0), llEvsG4, a=matA[i,], s=matS[i,], method="BFGS")$par
	varbeta1 <- compVarEvsG4(coef, matA)	
	stat <- coef[,2] * coef[,2] / varbeta1
	cbind(coef[,2], sqrt(varbeta1), stat)	
}	

EvsG2 <- function(geno, env2, size=50, type="dominant"){
	if(type=="dominant")
		geno <- (geno > 0) * 1
	else
		geno <- (geno == 2) * 1
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	matStat <- matrix(nrow = n.snp, ncol = 3)
	for(i in 1:(length(int)-1))
		matStat[int[i]:(int[i+1]-1), ] <- EvsG2split(geno[,int[i]:(int[i+1]-1), drop=FALSE], 
			env2==1, size=size)
	colnames(matStat) <- c("coef", "se", "stat")
	pval <- pchisq(matStat[ , 3], 1, lower.tail=FALSE)
	data.frame(matStat, pval = pval)
}

EvsG2split <- function(geno, env3, size=50){
	matA <- EvsG2chunk(geno)
	matS <- EvsG2chunk(geno[env3,,drop=FALSE], onlySum=TRUE)
	coef <- matrix(NA, ncol(geno), 2)
	for(i in 1:ncol(geno))
		coef[i,] <- optim(c(0,0), llEvsG2, a=matA[i,], s=matS[i,], method="BFGS")$par
	varbeta1 <- compVarEvsG2(coef, matA)
	stat <- coef[,2] * coef[,2] / varbeta1
	cbind(coef[,2], sqrt(varbeta1), stat)
}


print.colGxE <- function(x, top=5, digits=4, onlyGxE=FALSE, ...){
	if(x$n.select>0){
		if(!onlyGxE){
			pvalG <- format.pval(x$pval[,1], digits=digits)
			outG <- data.frame(Coef=x$coef[,1], RR=x$RR[,1], Lower=x$lowerRR[,1], Upper=x$upperRR[,1],
				SE=x$se[,1], Statistic=x$stat[,1], "p-value"=pvalG, check.names=FALSE, stringsAsFactors=FALSE)
			if(x$addGandE)
				outOR <- data.frame(RR=x$RR[,3], Lower=x$lowerRR[,3], Upper=x$upperRR[,3], check.names=FALSE,
					stringsAsFactors=FALSE)
			if(any(x$addOther)){
				outMore <- data.frame(row.names=rownames(outG))
				if(x$addOther[1])
					outMore <- data.frame(outMore, "2df Stat"=x$lrt2df[,1], 
						"2df p-Value"=format.pval(x$lrt2df[,2], digits=digits), check.names=FALSE,
						stringsAsFactors=FALSE)
				if(x$addOther[2])
					outMore <- data.frame(outMore, "Wald Stat"=x$wald2df[,1],
						"Wald p-value"=format.pval(x$wald2df[,2], digits=digits), check.names=FALSE,
						stringsAsFactors=FALSE)
				if(x$addOther[3])
					outMore <- data.frame(outMore, "1df Stat"=x$lrt1df[,1], 
						"1df p-Value"=format.pval(x$lrt1df[,2], digits=digits), check.names=FALSE,
						stringsAsFactors=FALSE)
			}
			
		}
		pvalGE <- format.pval(x$pval[,2], digits=digits)
		outGE <- data.frame(Coef=x$coef[,2], RR=x$RR[,2], Lower=x$lowerRR[,2], Upper=x$upperRR[,2],
			SE=x$se[,2], Statistic=x$stat[,2], "p-value"=pvalGE, Trios0=x$usedTrios[,1],
			Trios1=x$usedTrios[,2], check.names=FALSE, stringsAsFactors=FALSE)
		if(x$n.select == 1){
			rn <- rownames(x$coef)
			rownames(outGE) <- rn
			if(!onlyGxE){
				rownames(outG) <- rn
				if(x$addGandE)
					rownames(outOR) <- rn
				if(any(x$addOther))
					rownames(outMore) <- rn
			}
		}		
	}
	if(x$alpha1 < 1)
		cat("         Gauderman's Two-Step Procedure for Testing GxE Interactions\n\n", 
			"First Step: Logistic Regression of E vs. Sum of Genotypes of Parents\n\n",
			ifelse(x$n.select==0, "None", x$n.select), " of ", length(x$step1pval),
			" SNPs selected for the second step  (using alpha1 = ", x$alpha1, ").\n\n\n", sep="")
	if(x$n.select>0){
		if(x$alpha1 < 1)
			cat("Second Step: Genotypic TDT for GxE Interactions with Binary E\n\n", sep="")
		else
			cat("          Genotypic TDT for GxE Interactions with Binary E\n\n", sep="")
		cat("Model Type: ", switch(x$type, "additive"="Additive", "dominant"="Dominant","recessive"="Recessive"), 
			"\n\n",sep="")
		if(!is.na(top) && top>0 && top <= nrow(x$coef)){
			ord <- order(x$pval[,2])[1:top]
			if(!onlyGxE){
				outG <- outG[ord,]
				if(x$addGandE)
					outOR <- outOR[ord,]
				if(any(x$addOther))
					outMore <- outMore[ord,]
			}
			outGE <- outGE[ord,]
			cat("Top", top, "GxE Interactions:\n")
		}
		else
			cat("Effects of the GxE Interactions:\n")
		print(format(outGE, digits=digits))
		if(!onlyGxE){
			cat("\n\n", "Effects of the SNPs in the Corresponding GxE Models:\n", sep="")
			print(format(outG, digits=digits))
			if(x$addGandE){
				cat("\n\n", "RRs for Exposed Cases:\n", sep="")
				print(format(outOR, digits=digits))
			}
			if(any(x$addOther)){
				txt <- paste(c("2 df Likelihood Ratio Test", "2 df Wald Test", 
					"1 df Likelihood Ratio Test")[x$addOther], collapse=", ")
				cat("\n\n", txt, ":\n", sep="")
				print(format(outMore, digits=digits))
			}
		}
	}
}


