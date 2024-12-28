# Removed colGxE, fastGxEsplit, fastGxEdomSplit, fastGxErecSplit, EvsG4, EvsG4split,
# EvsG2, EvsG2split, print.colGxE.
# Updated in versions in colGxEnew.R.

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


colGxEunstructured <- function(mat.snp, env){
	requireNamespace("survival", quietly=TRUE)
	mat.code <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 
		1,1,1,1, 2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 
		1,0,1,2, 2,0,1,1, NA,NA,NA,NA), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112", "NANANA")
	colnames(mat.code) <- cn
	n.snp <- ncol(mat.snp)
	mat.pseudo <- matrix(NA, 4/3 * nrow(mat.snp), n.snp)
	for(i in 1:n.snp){
		mat.trio <- matrix(mat.snp[,i], ncol=3, byrow=TRUE)
		mat.trio[rowSums(is.na(mat.trio)) > 0, ] <- NA
		code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
		if(any(!code %in% cn)){
			tmp.ids <- !code %in% cn
			warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
				call.=FALSE)
			code[tmp.ids] <- "NANANA"
		}
		mat.pseudo[,i] <- as.vector(mat.code[,code])
	}
	n.trio <- nrow(mat.snp) / 3
	env2 <- rep(env, e=4)
	strat <- rep(1:n.trio, e=4)
	y <- rep.int(c(1, rep.int(0, 3)), n.trio) 
	ll.main <- ll.full <- rep.int(NA, n.snp)
	wa <- options()$warn
	options(warn=2)
	for(i in 1:n.snp){
		x1 <- mat.pseudo[,i] == 1
		x2 <- mat.pseudo[,i] == 2
		x1e <- x1 * env2
		x2e <- x2 * env2
		woIA <- try(clogit(y ~ x1 + x2 + strata(strat)), silent=TRUE)
		full <- try(clogit(y ~ x1 + x2 + x1e + x2e + strata(strat)), silent=TRUE)
		ll.main[i] <- if(is(woIA, "try-error")) NA else woIA$loglik[2]
		ll.full[i] <- if(is(full, "try-error")) NA else full$loglik[2]
	}
	options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("For some interactions, the fitting of at least one of the models has failed.\n",
			"Therefore, the corresponding test statistic and the p-value are thus set to NA.",
			call.=FALSE)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 2, lower.tail=FALSE)
	cn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:n.snp, sep="") else colnames(mat.snp)	
	names(ll.full) <- names(ll.main) <- names(stat) <- names(pval) <- cn
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval)
	class(out) <- "colGxEunstruct"
	out
}
	
print.colGxEunstruct <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("LL (with IA)"=x$ll.full, "LL (w/o IAs)"=x$ll.main, Statistic=x$stat,
		"P-Value"=x$pval, check.names=FALSE)
	cat("   Genotypic TDT for GxE Interactions Using a Fully Parameterized Model\n\n")
	if(top > 0 && length(x$ll.main) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "GxE Interactions (Likelihood Ratio Test):\n")
	}
	else
		cat("Likelihood Ratio Test:\n")
	print(format(out, digits=digits))
}


