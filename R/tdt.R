# Removed colTDT. Updated version in colTDTnew.R.


tdt <- function(snp, model=c("additive", "dominant", "recessive")){
	requireNamespace("survival")
	n <- length(snp)
	if(n%%3 != 0)
		stop("snp has a length not dividable by 3.")
	if(any(!snp %in% c(0:2, NA)))
		stop("The values in snp must be 0, 1, and 2.")
	type <- match.arg(model)
	mat.trio <- matrix(snp, ncol=3, byrow=TRUE)
	mat.trio <- mat.trio[rowSums(is.na(mat.trio))==0,]
	mat <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 1,1,1,1,
		2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 1,0,1,2, 
		2,0,1,1), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112")
	colnames(mat) <- cn
	code <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
	if(any(!code %in% cn)){
		tmp.ids <- !code %in% cn
		warning(sum(tmp.ids), " trios show Mendelian errors. These are removed.",
			call.=FALSE)
		code <- code[!tmp.ids]
	}
	x <- as.vector(mat[,code])
	if(type=="dominant")
		x <- (x >= 1) * 1
	if(type=="recessive")
		x <- (x > 1) * 1
	y <- rep.int(c(1,0,0,0), length(code))
	strat <- rep(1:length(code), e=4)
	c.out <- clogit(y ~ x + strata(strat))
	coef <- c.out$coefficients
	names(coef) <- NULL
	se <- sqrt(diag(c.out$var))
	stat <- (coef/se)^2
	conf <- exp(coef + c(-1, 1) * qnorm(0.975) * se)
	out <- list(coef=coef, se=se, stat=stat, pval= 1 - pchisq(stat, 1),
		RR=exp(coef), lowerRR=conf[1], upperRR=conf[2], ia=FALSE, type=type)  
	class(out) <- "tdt"
	out
}

print.tdt <- function(x, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame(Coef=x$coef, RR=x$RR, Lower=x$lowerRR, Upper=x$upperRR, 
		SE=x$se, Statistic=x$stat, "p-Value"=pval, check.names=FALSE)
	if(length(x$coef)==1)
		rownames(out) <- ""
	if(x$ia)
		cat("Genotypic TDT for Two-Way Interaction (Using 15 Pseudo Controls)",
			"\n\n")
	else
		cat("        Genotypic TDT Based on 3 Pseudo Controls","\n\n")
	cat("Model Type:", switch(x$type, "additive"="Additive", "dominant"="Dominant",
		"recessive"="Recessive"), "\n\n")
	if(x$ia && is.na(x$stat))
		cat("Fitting failed. Possible reason: SNPs might be in (strong) LD.\n\n")
	
	else
		print(format(out, digits=digits))
	if(!is.na(x$pval) && x$pval <= .Machine$double.eps)
		warning("The results might be misleading, as the very small p-value might be caused\n",
			"by non-biological artefacts (e.g., sparseness of the data).", call.=FALSE)
}

tdt2way <- function(...)
	cat("tdt2way has been renamed to tdtGxG. So please use tdtGxG instead.")

tdtGxG <- function(snp1, snp2, test=c("epistatic", "lrt", "full", "screen"), model=c("additive", "dominant", "recessive")){
	requireNamespace("survival")
	n1 <- length(snp1)
	n2 <- length(snp2)
	if(n1 != n2)
		stop("snp1 must have the same length as snp2.")
	if(n1%%3 != 0)
		stop("The SNPs do not seem to contain trio data, as their length\n",
			"is not dividable by 3.")
	if(any(!snp1 %in% c(0, 1, 2, NA)))
		stop("The values in snp1 must be 0, 1, and 2.")
	if(any(!snp2 %in% c(0, 1, 2, NA)))
		stop("The values in snp2 must be 0, 1, and 2.")
	type <- match.arg(model)
	testtype <- match.arg(test)
	mat.trio <- cbind(matrix(snp1, ncol=3, byrow=TRUE), matrix(snp2, ncol=3, byrow=TRUE))
	mat.trio <- mat.trio[rowSums(is.na(mat.trio))==0,]
	mat <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 1,1,1,1,
		2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 1,0,1,2, 
		2,0,1,1), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112")
	colnames(mat) <- cn
	code1 <- paste(mat.trio[,1], mat.trio[,2], mat.trio[,3], sep="")
	code2 <- paste(mat.trio[,4], mat.trio[,5], mat.trio[,6], sep="")
	if(any(!code1 %in% cn)){
		tmp.ids <- !code1 %in% cn
		warning(sum(tmp.ids), " trios in snp1 show Mendelian errors. These are removed.",
			call.=FALSE)
		code1 <- code1[!tmp.ids]
		code2 <- code2[!tmp.ids]
	}
	if(any(!code2 %in% cn)){
		tmp.ids <- !code2 %in% cn
		warning(sum(tmp.ids), " trios show Mendelian errors in snp2 (but not in snp1).\n",
			"These trios are removed.", call.=FALSE)
		code2 <- code2[!tmp.ids]
		code1 <- code1[!tmp.ids]
	}  
	n.trio <- length(code1)		
	x1 <- as.vector(mat[,code1])
	x2 <- as.vector(mat[,code2])
	if(testtype=="epistatic")
		return(tdtEpistatic(x1, x2, n.trio))
	if(type=="dominant"){
		x1 <- (x1 >= 1) * 1
		x2 <- (x2 >= 1) * 1
	}
	if(type=="recessive"){
		x1 <- (x1 > 1) * 1
		x2 <- (x2 > 1) * 1
	} 
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	IA <- x1[mat.ids[,1]] * x2[mat.ids[,2]]
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	wa <- options()$warn
	options(warn=2)
	if(testtype=="screen")
		c.out <- try(clogit(y ~ IA + strata(strat)), silent=TRUE)
	else{
		SNP1 <- x1[mat.ids[,1]]
		SNP2 <- x2[mat.ids[,2]]
		c.out <- try(clogit(y ~ SNP1 + SNP2 + IA + strata(strat)), silent=TRUE)
		if(testtype=="lrt")
			c.out2 <- try(clogit(y ~ SNP1 + SNP2 + strata(strat)), silent=TRUE)
	}
	options(warn=wa)
	if(testtype=="lrt"){
		ll.main <- if(is(c.out2, "try-error")) NA else c.out2$loglik[2]
		ll.full <- if(is(c.out, "try-error")) NA else c.out$loglik[2]
		if(!is.na(ll.full) && !is.na(ll.main)){
			stat <- -2 * (ll.main - ll.full)
			pval <- pchisq(stat, 4, lower.tail=FALSE)
		}
		else
			stat <- pval <- NA
		out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval, 
			full.model=c.out, testtype="lrt")
		return(out)
	}
	if(is(c.out, "try-error"))
		coef <- se <- stat <- pval <- lower <- upper <- NA
	else{
		coef <- c.out$coefficients
		se <- sqrt(diag(c.out$var))
		if(type=="screen")
			names(coef) <- NULL
		stat <- (coef/se)^2
		pval <- pchisq(stat, 1, lower.tail=FALSE)
		lower <- exp(coef - qnorm(0.975) * se)
		upper <- exp(coef + qnorm(0.975) * se)
	}
	out <- list(coef=coef, se=se, stat=stat, pval= pval,
		RR=exp(coef), lowerRR=lower, upperRR=upper, ia=TRUE, type=type)  
	class(out) <- "tdt"
	out

}


colTDT2way <- function(...){
	cat("colTDT2way has been renamed to colGxG. So please use colGxG instead.\n")
}

colGxG <- function(mat.snp, test=c("epistatic", "lrt", "full", "screen"), genes=NULL, maf=FALSE,
		model=c("additive", "dominant", "recessive")){
	requireNamespace("survival")
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
	testtype <- match.arg(test)
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
	if(maf){
		mat.snp <- mat.snp[-seq(3, nrow(mat.snp), 3), ]
		mat.snp <- mat.snp[!duplicated(rownames(mat.snp)), ]
		valMAF <- colSums(mat.snp, na.rm=TRUE) / (2 * colSums(!is.na(mat.snp)))
	}
	else
		valMAF <- NULL
	n.trio <- length(code)
	if(testtype=="epistatic")
		return(colTDTepistatic(mat.pseudo, n.trio, colnames(mat.snp), genes=genes, valMAF=valMAF))
	if(type=="dominant")
		mat.pseudo <- (mat.pseudo >= 1) * 1
	if(type=="recessive")
		mat.pseudo <- (mat.pseudo > 1) * 1
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	if(is.null(genes)){
		coef <- se <- numeric(n.snp * (n.snp - 1) / 2)
		combs <- allCombs(n.snp)
	}
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
		coef <- se <- numeric(nrow(combs))
	}
	if(testtype=="lrt")
		return(colGxGlrt(mat.pseudo, mat.ids, combs, y, strat, rn=colnames(mat.snp), genes=genes,
			valMAF=valMAF))
	add <- testtype=="full"
	wa <- options()$warn
	options(warn=2)
	for(i in 1:nrow(combs)){
		x1 <- mat.pseudo[mat.ids[,1], combs[i,1]] 
		x2 <- mat.pseudo[mat.ids[,2], combs[i,2]]
		IA <- x1 * x2
		if(add)
			c.out <- try(clogit(y ~ IA + x1 + x2 + strata(strat)), silent=TRUE)
		else
			c.out <- try(clogit(y ~ IA + strata(strat)), silent=TRUE)
		coef[i] <- if(is(c.out, "try-error")) NA else c.out$coefficients[1]
		se[i] <- if(is(c.out, "try-error")) NA else sqrt(diag(c.out$var))[1]
	}
	options(warn=wa)
	if(any(is.na(coef)))
		warning("For at least one interaction,the fitting of the model has failed.\n",
			"Therefore, all statistics for these interactions are set to NA.", call.=FALSE)
	stat <- (coef / se)^2
	lower <- exp(coef - qnorm(0.975) * se)
	upper <- exp(coef + qnorm(0.975) * se)
	pval <- pchisq(stat, 1, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	rn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:n.snp, sep="") else colnames(mat.snp)
	names(coef) <- names(stat) <- names(pval) <- paste(rn[combs[,1]], rn[combs[,2]], sep=" : ")
	if(!is.null(genes))
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	if(maf){
		mat.maf <- cbind(round(valMAF, 4)[combs[,1]], round(valMAF, 4)[combs[,2]])
		rownames(mat.maf) <- names(coef) 
		colnames(mat.maf) <- c("First SNP", "Second SNP")
	}
	else
		mat.maf <- NULL
	out <- list(coef=coef, se=se, stat=stat, pval=pval, RR=exp(coef), 
		lowerRR=lower, upperRR=upper, ia=TRUE, type=type, add=add, genes=genes, maf=valMAF,
		matMAF=mat.maf) 
	class(out) <- "colTDT"
	out	

}


print.colTDT <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	if(x$ia){
		out <- data.frame(Coef=x$coef, RR=x$RR, Lower=x$lowerRR, Upper=x$upperRR, 
			SE=x$se, Statistic=x$stat, "p-Value"=pval,
			check.names=FALSE, stringsAsFactors=FALSE)
		cat("      Genotypic TDT for Two-Way Interaction (Using 15 Pseudo Controls)",
			"\n\n")
	}
	else{
		out <- data.frame(Coef=x$coef, RR=x$RR, Lower=x$lowerRR, Upper=x$upperRR, 
			SE=x$se, Statistic=x$stat, "p-Value"=pval, Trios=x$usedTrios,
			check.names=FALSE, stringsAsFactors=FALSE)
		if(!is.null(x$pMendelErr))
			out <- data.frame(out, "P(Mendel Error)"=x$pMendelErr, 
				check.names=FALSE, stringsAsFactors=FALSE)
		cat("        Genotypic TDT Based on 3 Pseudo Controls","\n\n")
	}
	cat("Model Type: ", switch(x$type, "additive"="Additive", "dominant"="Dominant",
		"recessive"="Recessive"), "\n", 
		if(x$add) "Model also contains the two respective individual SNPs.\n",
		"\n", sep="")
	if(!is.na(top) && top>0 && top <= length(x$coef)){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, ifelse(x$ia, "SNP Interactions:\n", "SNPs:\n"))
	}
	print(format(out, digits=digits))
}

tdtEpistatic <- function(x1, x2, n.trio){
	x1 <- x1 - 1
	z1 <- ifelse(x1==0, 0.5, -0.5)
	x2 <- x2 - 1
	z2 <- ifelse(x2==0, 0.5, -0.5)

	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	x1 <- x1[mat.ids[,1]]
	x2 <- x2[mat.ids[,2]]
	z1 <- z1[mat.ids[,1]]
	z2 <- z2[mat.ids[,2]]
	x1x2 <- x1 * x2
	x1z2 <- x1 * z2
	z1x2 <- z1 * x2
	z1z2 <- z1 * z2
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	wa <- options()$warn
	options(warn=2)
	woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
	full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1x2 + x1z2 + z1x2 + z1z2 + 
		strata(strat)), silent=TRUE)
	options(warn=wa)
	ll.main <- if(is(woIA, "try-error")) NA else woIA$loglik[2]	
	ll.full <- if(is(full, "try-error")) NA else full$loglik[2]
	if(!is.na(ll.full) && !is.na(ll.main)){
		stat <- -2 * (ll.main - ll.full)
		pval <- pchisq(stat, 4, lower.tail=FALSE)
	}
	else
		stat <- pval <- NA
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval,
		full.model=full, testtype="epistatic") 
	class(out) <- "tdtEpi"
	out
}

colTDTepistatic <- function(mat.pseudo, n.trio, rn, genes=NULL, valMAF=NULL){
	x <- mat.pseudo - 1 
	z <- (x == 0) - 0.5
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	n.snp <- ncol(mat.pseudo)
	if(is.null(genes)){
		ll.main <- ll.full <- numeric(n.snp * (n.snp - 1) / 2)
		combs <- allCombs(n.snp)
	}
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
		ll.main <- ll.full <- numeric(nrow(combs))
	}
	wa <- options()$warn
	options(warn=2)
	for(i in 1:nrow(combs)){
		x1 <- x[mat.ids[,1], combs[i,1]]
		z1 <- z[mat.ids[,1], combs[i,1]] 
		x2 <- x[mat.ids[,2], combs[i,2]]
		z2 <- z[mat.ids[,2], combs[i,2]]
		woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
		full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + 
			z1*z2 + strata(strat)), silent=TRUE)
		ll.main[i] <- if(is(woIA, "try-error")) NA else woIA$loglik[2]
		ll.full[i] <- if(is(full, "try-error")) NA else full$loglik[2]
	}
	options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("For some interactions, the fitting of at least one of the models has failed.\n",
			"Therefore, the corresponding test statistic and the p-value are thus set to NA.",
			call.=FALSE)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	if(is.null(rn))
		rn <- paste("SNP", 1:n.snp, sep="")
	names(ll.full) <- names(stat) <- names(pval) <- paste(rn[combs[,1]], rn[combs[,2]], sep=" : ")
	if(!is.null(genes)){
		ind.genes <- genes
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	}
	else
		ind.genes <- NULL
	if(!is.null(valMAF)){
		mat.maf <- cbind(round(valMAF, 4)[combs[,1]], round(valMAF, 4)[combs[,2]])
		rownames(mat.maf) <- names(ll.full)
		colnames(mat.maf) <- c("First MAF", "Second MAF") 
	}
	else
		mat.maf <- NULL
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval, 
		maf=valMAF, matMAF=mat.maf, genes=genes, ind.genes=ind.genes, testtype="epistatic")
	class(out) <- "colTDTepi"
	out
}


colGxGlrt <- function(mat.pseudo, mat.ids, combs, y, strat, rn=NULL, genes=NULL, valMAF=NULL){
	ll.main <- ll.full <- rep.int(NA, nrow(combs))
	wa <- options()$warn
	options(warn=2)
	for(i in 1:nrow(combs)){
		x1 <- mat.pseudo[mat.ids[,1], combs[i,1]]
		x2 <- mat.pseudo[mat.ids[,2], combs[i,2]]
		IA <- x1*x2
		full <- try(clogit(y ~ x1 + x2 + IA + strata(strat)), silent=TRUE)
		main <- try(clogit(y ~ x1 + x2 + strata(strat)), silent=TRUE)
		ll.full[i] <- if(is(full, "try-error")) NA else full$loglik[2]
		ll.main[i] <- if(is(main, "try-error")) NA else main$loglik[2]
	}
	options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("For some interactions, the fitting of at least one of the models has failed.\n",
			"Therefore, the corresponding test statistic and the p-value are thus set to NA.",
			call.=FALSE)
	stat <- -2 * (ll.main - ll.full)
	pval <- pchisq(stat, 4, lower.tail=FALSE)
	if(any(pval <= .Machine$double.eps, na.rm=TRUE))
		warning("Some of the p-values are very small (< 2.2e-16). These results might be\n",
			"misleading, as the small p-values might be caused by non-biological artefacts\n",
			"(e.g., sparseness of the data).", call.=FALSE)
	if(is.null(rn))
		rn <- paste("SNP", 1:ncol(mat.pseudo), sep="")
	names(ll.full) <- names(stat) <- names(pval) <- paste(rn[combs[,1]], rn[combs[,2]], sep=" : ")
	if(!is.null(genes)){
		ind.genes <- genes
		genes <- paste(genes[combs[,1]], genes[combs[,2]], sep=" : ")
	}
	else
		ind.genes <- NULL
	if(!is.null(valMAF)){
		mat.maf <- cbind(round(valMAF, 4)[combs[,1]], round(valMAF, 4)[combs[,2]])
		rownames(mat.maf) <- names(ll.full)
		colnames(mat.maf) <- c("First MAF", "Second MAF") 
	}
	else
		mat.maf <- NULL
	out <- list(ll.main=ll.main, ll.full=ll.full, stat=stat, pval=pval, 
		maf=valMAF, matMAF=mat.maf, genes=genes, ind.genes=ind.genes, testtype="lrt")
	class(out) <- "colTDTepi"
	out	
}	

print.tdtEpi <- function(x, digits=3, ...){
	pval <- format.pval(x$pval, digits=digits-1)
	if(x$testtype=="epistatic")
		cat("Likelihood Ratio Test for Epistatic Interactions Based on Genotypic TDTs", 
			"\n\n", sep="")
	else
		cat("Likelihood Ratio Test Based on Genotypic TDTs\n", sep="")
	if(is.na(x$stat))
		cat("Failed. At least one of the models did not fit properly",
			"\n\n")
	else{
		cat("Loglikelihood (with Interactions): ", round(x$ll.full, digits), "\n", sep="")
		cat("Loglikelihood (without IAs):       ", round(x$ll.main, digits), "\n", sep="")
		cat("Test Statistic: ", round(x$stat, digits), "\n", sep="")
		cat("P-Value:        ", pval, "\n\n", sep="")
	}
}

print.colTDTepi <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
	out <- data.frame("LL (with IAs)"=x$ll.full, "LL (w/o IAs)"=x$ll.main,
		Statistic=x$stat, "P-Value"=pval, check.names=FALSE)
	if(!is.null(x$genes))
		out <- data.frame(out, Genes=x$genes, check.names=FALSE, stringsAsFactors=FALSE)
	if(!is.null(x$matMAF))
		out <- data.frame(out, x$matMAF, check.names=FALSE, stringsAsFactors=FALSE) 
	cat("      Genotypic TDT for", ifelse(x$testtype=="epistatic", "Epistatic", "SNP-SNP"), 
		"Interactions (Using 15 Pseudo Controls)", "\n\n")
	if(length(x$ll.main) > top){
		ord <- order(x$pval)[1:top]
		out <- out[ord,]
		cat("Top", top, "SNP Interactions (Likelihood Ratio Test):\n")
	}	
	else
		cat("Likelihood Ratio Test:\n")
	print(format(out, digits=digits))
}
	


colTDTinter2way <- function(...){
	cat("colTDTinter2way has been removed from trio. Please use colGxG instead.\n",
		"In colGxG, the argument genes can be used to specify the different sets of SNPs.\n", sep="")
}
	
allBetweenCombs <- function(gene){
	ids <- 1:length(gene)
	mat.comb <- NULL
	for(i in 1:(max(gene) - 1)){
		tmp1 <- ids[gene == i]
		tmp2 <- ids[gene > i]
		tmpmat <- cbind(rep(tmp1, e=length(tmp2)), rep.int(tmp2, length(tmp1)))
		mat.comb <- rbind(mat.comb, tmpmat)
	}
	mat.comb
}

getMatPseudo <- function(mat.snp){
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
	colnames(mat.pseudo) <- colnames(mat.snp)
	mat.pseudo
}
