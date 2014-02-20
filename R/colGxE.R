colGxE <- function(mat.snp, env, model=c("additive", "dominant", "recessive"), alpha1=1, size=50, addGandE=TRUE,
		whichLRT=c("both", "2df", "1df", "none"), add2df=TRUE, addCov=FALSE, famid=NULL){
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
	if(alpha1 > 1 | alpha1 <= 0)
		stop("alpha1 must be larger than 0 and smaller than or equal to 1.") 
	typeLRT <- match.arg(whichLRT)
	addLRT1 <- typeLRT %in% c("both", "1df") 
	addLRT2 <- typeLRT %in% c("both", "2df")
	type <- match.arg(model)
	env2 <- rep(env, e=3)
	if(alpha1 < 1){
		if(type=="additive")
			step1pval <- EvsG4(mat.snp, env2, size=size)
		else
			step1pval <- EvsG2(mat.snp, env2, size=size, type=type)
		names(step1pval) <- colnames(mat.snp)
		if(all(step1pval > alpha1)){
			out <- list(step1pval=step1pval, n.select=0)
			class(out) <- "colGxE"
			return(out)
		}
		mat.snp <- mat.snp[,step1pval <= alpha1, drop=FALSE]
	}   
	fun <- match.fun(switch(type, additive=fastGxEsplit, dominant=fastGxEdomSplit,
		recessive=fastGxErecSplit))
	tmp1 <- fun(mat.snp, env2==0, size=size, addLRT1=addLRT1, addLRT2=addLRT2)
	tmp2 <- fun(mat.snp, env2==1, size=size, addLRT1=addLRT1, addLRT2=addLRT2)
	beta <- se <- matrix(NA, ncol(mat.snp), 2+addGandE)
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
	rownames(beta) <- rownames(se) <- rownames(used) <- colnames(mat.snp)
	stat <- beta[,1:2] / se[,1:2]
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
		funG <- match.fun(switch(type, additive=lrtGxE, dominant=lrtGxEdom, recessive=lrtGxErec))
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
	out <- list(coef=beta[,1:2], se=se[,1:2], stat=stat, pval=pval, RR=exp(beta), lowerRR=lower, upperRR=upper,
		usedTrios=used, env=env, type=type, addGandE=addGandE, addOther=c(addLRT2, add2df, addLRT1), 
		n.select=ncol(mat.snp))
	if(alpha1<1)
		out$step1pval <- step1pval
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

reorderEnv <- function(env, famid, rn){
	if(length(famid) != length(env))
		stop("famid has another length as env.")
	if(!is.null(names(env))){
		if(any(famid != names(env)))
			stop("If the names of env are specified, they must be equal to famid.")
	}
	fid <- sapply(strsplit(rn[seq.int(3, length(rn), 3)], "_"), function(x) x[1])
	names(env) <- famid
	if(any(!famid %in% fid))
		stop("At least one famid is not a family ID used in mat.snp.")
	if(all(famid==fid))
		return(env)
	env <- env[fid]
	env
} 


fastGxEsplit <- function(geno, env2, size=50, addLRT1=TRUE, addLRT2=TRUE){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))	
	num <- denom <- used <- numeric(n.snp)
	for(i in 1:(length(int)-1)){
		tmp <- fastTDTchunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
		num[int[i]:(int[i+1]-1)] <- tmp$num
		denom[int[i]:(int[i+1]-1)] <- tmp$denom
		used[int[i]:(int[i+1]-1)] <- tmp$used
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
	


fastGxEdomSplit <- function(geno, env2, size=50, addLRT1=TRUE, addLRT2=TRUE){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	dmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		dmat[int[i]:(int[i+1]-1),] <- fastTDTdomChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
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

fastGxErecSplit <- function(geno, env2, size=50, addLRT1=TRUE, addLRT2=TRUE){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	rmat <- matrix(0, n.snp, 4)
	for(i in 1:(length(int)-1))
		rmat[int[i]:(int[i+1]-1),] <- fastTDTrecChunk(geno[env2,int[i]:(int[i+1]-1), drop=FALSE])
	# rownames(rmat) <- colnames(geno)
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


print.colGxE <- function(x, top=5, digits=4, onlyGxE=FALSE, ...){
	if(x$n.select>0){
		if(!onlyGxE){
			pvalG <- format.pval(x$pval[,1], digits=digits)
			outG <- data.frame(Coef=x$coef[,1], RR=x$RR[,1], Lower=x$lowerRR[,1], Upper=x$upperRR[,1],
				SE=x$se[,1], Statistic=x$stat[,1], "p-value"=pvalG, check.names=FALSE, stringsAsFactors=FALSE)
			if(any(x$addGandE))
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
	}
	if(!is.null(x$step1pval))
		cat("         Gauderman's Two-Step Procedure for Testing GxE Interactions\n\n", 
			"First Step: Logistic Regression of E vs. Sum of Genotypes of Parents\n\n",
			ifelse(x$n.select==0, "None", x$n.select), " of ", length(x$step1pval),
			" SNPs selected for the second step.\n\n\n", sep="")
	if(x$n.select>0){
		if(!is.null(x$step1pval))
			cat("Second Step: Genotypic TDT for GxE Interactions with Binary E\n\n", sep="")
		else
			cat("          Genotypic TDT for GxE Interactions with Binary E\n\n", sep="")
		cat("Model Type: ", switch(x$type, "additive"="Additive", "dominant"="Dominant","recessive"="Recessive"), 
			"\n\n",sep="")
		if(!is.na(top) && top>0 && top <= nrow(x$coef)){
			ord <- order(x$pval[,2])[1:top]
			if(!onlyGxE){
				outG <- outG[ord,]
				outOR <- outOR[ord,]
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

getGxEstats <- function(x, top=NA, sortBy=c("none", "gxe", "lrt2df", "wald2df", "lrt1df", "g")){
	if(!is(x, "colGxE"))
		stop("x must be the output of colGxE.")
	if(x$n.select==0)
		stop("No SNP has been chosen for the second step of Gauderman's procedure for testing GxE interactions.")
	type <- match.arg(sortBy)
	if(!is.na(top)){
		if(type=="none")
			stop("If sortBy is set to none, top must be set to NA.")
		if(top <= 0 | top>nrow(x$pval))
			top <- NA
	}		
	dat <- data.frame("Coef GxE"=x$coef[,2], "RR GxE"=x$RR[,2], "Lower GxE"=x$lowerRR[,2], "Upper GxE"=x$upperRR[,2],
		"SE GxE"=x$se[,2], "Stat GxE"=x$stat[,2], "pval GxE"=x$pval[,2], "Trios0"=x$usedTrios[,1],
		Trios1=x$usedTrios[,2], "Coef G"=x$coef[,1], "RR G"=x$RR[,1], "Lower G"=x$lowerRR[,1], "Upper G"=x$upperRR[,1],
		"SE G"=x$se[,1], "Stat G"=x$stat[,1], "pval G"=x$pval[,1], check.names=FALSE, stringsAsFactors=FALSE)
	if(x$addGandE)
		dat <- data.frame(dat, "RR G&E"=x$RR[,3], "Lower G&E"=x$lowerRR[,3], "Upper G&E"=x$upperRR[,3], 
			check.names=FALSE, stringsAsFactors=FALSE)
	if(x$addOther[1])
		dat <- data.frame(dat, "Stat LRT 2df"=x$lrt2df[,1], "pval LRT 2df"=x$lrt2df[,2], check.names=FALSE, 
			stringsAsFactors=FALSE)
	if(x$addOther[2])
		dat <- data.frame(dat, "Stat Wald 2df"=x$wald2df[,1], "pval Wald 2df"=x$wald2df[,2], check.names=FALSE, 
			stringsAsFactors=FALSE)
	if(x$addOther[3])
		dat <- data.frame(dat, "Stat LRT 1df"=x$lrt1df[,1], "pval LRT 1df"=x$lrt1df[,2], check.names=FALSE, 
			stringsAsFactors=FALSE)
	if(type=="none")
		return(dat)
	if(type=="gxe")
		ord <- order(x$pval[,2])
	if(type=="lrt2df"){
		if(!x$addOther[1])
			stop("x does not contain results from the 2 df Likelihood Ratio Test.")
		ord <- order(x$lrt2df[,2])
	}
	if(type=="wald2df"){
		if(!x$addOther[2])
			stop("x does not contain results from the 2 df Wald Test.")
		ord <- order(x$wald2df[,2])
	}
	if(type=="lrt1df"){
		if(!x$addOther[3])
			stop("x does not contain results from the 1 df Likelihood Ratio Test.")
		ord <- order(x$lrt1df[,2])
	}
	if(type=="g")
		ord <- order(x$pval[,1])
	if(!is.na(top))
		ord <- ord[1:top]
	dat[ord,]
}

lrtGxE <- function(stats){
	beta <- logit(stats[,1]/stats[,2])
	if(any(is.infinite(beta)))
		beta[is.infinite(beta)] <- NA
	tmpFac <- stats[,2]/(stats[,2]-stats[,1])
	lnTerm <- log(tmpFac) * stats[,2]
	oneHet <- 2*stats[,3] - stats[,2]
	beta * stats[,1] - log(2) * oneHet - lnTerm
}

lrtGxEdom <- function(stats){
	h <- (stats[,1]/3 - stats[,2] + stats[,3] - stats[,4]/3) / (2*(stats[,1]+stats[,3]))
	d24 <- stats[,2] + stats[,4]
	tmp <- d24/(3*(stats[,1]+stats[,3])) + h*h
	or <- sqrt(tmp) - h
	beta <- log(or)
	if(any(is.infinite(beta)))
		beta[is.infinite(beta)] <- NA
	d12 <- stats[,2] + stats[,1]
	d34 <- stats[,3] + stats[,4]
	d24 * beta - d12 * log(2 + 2*or) - d34 * log(1+3*or)
}

lrtGxErec <- function(stats){
	h <- (3*stats[,1] - stats[,2] + stats[,3] - 3*stats[,4]) / (2 * (stats[,1]+stats[,3]))
	r24 <- stats[,2] + stats[,4]
	tmp <- 3 * r24 / (stats[,1] + stats[,3]) + h*h
	or <- sqrt(tmp) - h		 
	beta <- log(or)
	r12 <- stats[,1] + stats[,2]
	r34 <- stats[,3] + stats[,4]
	r24 * beta - r12 * log(2+2*or) - r34 * log(3+or)
}

EvsG4 <- function(geno, env2, size=50){
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- rep.int(NA, n.snp)
	for(i in 1:(length(int)-1))
		stat[int[i]:(int[i+1]-1)] <- EvsG4split(geno[,int[i]:(int[i+1]-1), drop=FALSE], env2==1, size=size)
	pchisq(stat, 1, lower.tail=FALSE)
}	

EvsG4split <- function(geno, env3, size=50){ 
	matA <- EvsG4chunk(geno)
	matS <- EvsG4chunk(geno[env3,,drop=FALSE], onlySum=TRUE)
	coef <- matrix(NA, ncol(geno), 2)
	for(i in 1:ncol(geno))
		coef[i,] <- optim(c(0,0), llEvsG4, a=matA[i,], s=matS[i,], method="BFGS")$par
	varbeta1 <- compVarEvsG4(coef, matA)	
	coef[,2] * coef[,2] / varbeta1
	
}	

EvsG4chunk <- function(geno, onlySum=FALSE){
	n.row <- nrow(geno)
	x <- geno[seq.int(1, n.row, 3),, drop=FALSE] + geno[seq.int(2, n.row, 3),, drop=FALSE]
	mat <- matrix(NA, ncol(geno), 5)
	for(i in 0:4)
		mat[,i+1] <- colSums(x==i, na.rm=TRUE)
	if(!onlySum)
		return(mat)	
	vecPos <- mat[,2:5] %*% (1:4)
	vecAll <- rowSums(mat, na.rm=TRUE)
	cbind(vecAll, vecPos)
}	
	
llEvsG4 <- function(beta, a, s){
	a[1] * log(1+exp(beta[1])) + a[2] * log(1+exp(beta[1]+beta[2])) + a[3] * log(1+exp(beta[1]+2*beta[2])) +
		a[4] * log(1+exp(beta[1] + 3*beta[2])) + a[5] * log(1+exp(beta[1]+4*beta[2])) - 
		beta[1] * s[1] - beta[2] * s[2] 	
}

compVarEvsG4 <- function(coef, matA){
	expb <- exp(coef[,1])
	fac0 <- matA[,1] * expb / (1+expb)^2
	expb1 <- exp(coef[,2])
	mat.fac <- matrix(NA, nrow(coef), 4)
	for(i in 1:4){
		expb <- expb * expb1
		mat.fac[,i] <- matA[,i+1] * expb / (1+expb)^2
	}
	b0b0 <- fac0 + rowSums(mat.fac)
	b0b1 <- mat.fac %*% (1:4)
	b1b1 <- mat.fac %*% c(1,4,9,16)
	detInfo <- b0b0 * b1b1 - b0b1 * b0b1
	b0b0/detInfo
}

EvsG2 <- function(geno, env2, size=50, type="dominant"){
	if(type=="dominant")
		geno <- (geno > 0) * 1
	else
		geno <- (geno == 2) * 1
	n.snp <- ncol(geno)
	int <- unique(c(seq.int(1, n.snp, size), n.snp+1))
	stat <- rep.int(NA, n.snp)
	for(i in 1:(length(int)-1))
		stat[int[i]:(int[i+1]-1)] <- EvsG2split(geno[,int[i]:(int[i+1]-1), drop=FALSE], env2==1, size=size)
	pchisq(stat, 1, lower.tail=FALSE)
}

EvsG2split <- function(geno, env3, size=50){
	matA <- EvsG2chunk(geno)
	matS <- EvsG2chunk(geno[env3,,drop=FALSE], onlySum=TRUE)
	coef <- matrix(NA, ncol(geno), 2)
	for(i in 1:ncol(geno))
		coef[i,] <- optim(c(0,0), llEvsG2, a=matA[i,], s=matS[i,], method="BFGS")$par
	varbeta1 <- compVarEvsG2(coef, matA)
	coef[,2] * coef[,2] / varbeta1
}

EvsG2chunk <- function(geno, onlySum=FALSE){
	n.row <- nrow(geno)
	x <- geno[seq.int(1, n.row, 3),, drop=FALSE] + geno[seq.int(2, n.row, 3),, drop=FALSE]
	mat <- matrix(NA, ncol(geno), 3)
	for(i in 0:2)
		mat[,i+1] <- colSums(x==i, na.rm=TRUE)
	if(!onlySum)
		return(mat)
	vecPos <- mat[,2] + 2*mat[,3]
	vecAll <- rowSums(mat, na.rm=TRUE)
	cbind(vecAll, vecPos)
}

llEvsG2 <- function(beta, a, s){
	a[1] * log(1+exp(beta[1])) + a[2] * log(1+exp(beta[1]+beta[2])) + a[3] * log(1+exp(beta[1]+2*beta[2])) - 
		beta[1] * s[1] - beta[2] * s[2]
}

compVarEvsG2 <- function(coef, matA){
	expb <- exp(coef[,1])
	fac0 <- matA[,1] * expb / (1+expb)^2
	expb1 <- exp(coef[,2])
	expb <- expb * expb1
	fac1 <- matA[,2] * expb / (1+expb)^2
	expb <- expb * expb1
	fac2 <- matA[,3] * expb / (1+expb)^2
	b0b0 <- fac0 + fac1 + fac2
	b0b1 <- fac1 + 2*fac2
	b1b1 <- fac1 + 4*fac2
	detInfo <- b0b0 * b1b1 - b0b1 * b0b1
	b0b0/detInfo
}



