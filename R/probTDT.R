# mat.geno should be nSNPs rows by 9*nFams columns (three probs for each family member, complete trios only, arranged parent 1, parent 2, proband)
probTDT<-function (mat.geno, model = c("additive", "dominant", "recessive"), size = 50) 
{
    if (!is.matrix(mat.geno)) 
        stop("mat.snp has to be a matrix.")
    if (ncol(mat.geno)%%9 != 0) 
        stop("mat.snp does not seem to contain trio genotype probability data, as its number of columns is\n", 
            "not dividable by 9.")
	# this probably takes too long, but it's a decent check
	if (!all(round(rowSums(mat.geno)) == ncol(mat.geno)/3))
		stop("mat.snp does not seem to contain genotype probability data as its rows do not sum to the\n",
		"number of subjects")
    type <- match.arg(model)    
	if (size < 1) 
        stop("size should be at least 1.")
    fun <- match.fun(switch(type, additive = probTDTsplit, dominant = probTDTdomSplit, 
        recessive = probTDTrecSplit))
    fun(mat.geno, size = size)
}

################################################
# Modified versions of functions from trioPaperCode.R
# Changes include: incorporation of usedTrios info for df adjustment
# 				   modification to prevent error in case of Mendelian inconsistency
#					inclusion of sum of mendelian error probs across fams for each marker
################################################

# basically a copy of fastTDTsplit with dimensions of geno changed, plus a couple more
# fields in return value from probTDTchunk
probTDTsplit<-function (geno, size = 50) 
{
    n.snp <- nrow(geno)
    int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
    num <- denom <- usedTrios <- pMendelErr <- numeric(n.snp)
    for (i in 1:(length(int) - 1)) {
        tmp <- probTDTchunk(geno[int[i]:(int[i + 1] - 1), , drop = FALSE])
        num[int[i]:(int[i + 1] - 1)] <- tmp$num
        denom[int[i]:(int[i + 1] - 1)] <- tmp$denom
        usedTrios[int[i]:(int[i + 1] - 1)] <- tmp$usedTrios
        pMendelErr[int[i]:(int[i + 1] - 1)] <- tmp$pMendelErr
		# remove printing to improve timing accuracy
		# if (i%%10 == 0) print(i)
    }
    logit <- function(x) log(x/(1 - x))
    beta <- logit(num/denom)
    se <- sqrt(denom/((denom - num) * num))
    stat <- beta/se
    stat <- stat * stat
    lower <- exp(beta - qnorm(0.975) * se)
    upper <- exp(beta + qnorm(0.975) * se)
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    if (is.null(rownames(geno))) 
        names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- names(pMendelErr) <- paste("SNP", 
            1:nrow(geno), sep = "")
    else names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- names(pMendelErr) <- rownames(geno)
    out <- list(coef = beta, se = se, stat = stat, pval = pval, 
        OR = exp(beta), lowerOR = lower, upperOR = upper, ia = FALSE, 
        type = "additive", usedTrios=usedTrios, add = FALSE, pMendelErr=pMendelErr)
	class(out)<-"colTDT"
    out
}

# also returns num informative families
# and prob of mendel error at each SNP
probTDTchunk<-function (geno) 
{
    n.col <- ncol(geno)
    dadAA <- geno[,seq.int(1, n.col, 9), drop = FALSE]
    dadAB <- geno[,seq.int(2, n.col, 9), drop = FALSE]
    dadBB <- geno[,seq.int(3, n.col, 9), drop = FALSE]
    momAA <- geno[,seq.int(4, n.col, 9), drop = FALSE]
    momAB <- geno[,seq.int(5, n.col, 9), drop = FALSE]
    momBB <- geno[,seq.int(6, n.col, 9), drop = FALSE]
    kidAA <- geno[,seq.int(7, n.col, 9), drop = FALSE]
    kidAB <- geno[,seq.int(8, n.col, 9), drop = FALSE]
    kidBB <- geno[,seq.int(9, n.col, 9), drop = FALSE]

	het1<-dadAA*momAB + dadAB*momAA
	a1<-het1*kidAA
	a2<-het1*kidAB
	
	hetHom<-dadAB*momBB + dadBB*momAB
	a3<-hetHom*kidAB
	a4<-hetHom*kidBB
	
	a567<-dadAB*momAB
	
	a67<-a567*(kidAB + 2*kidBB)

	a8910<-dadAA*momAA*kidAA + (dadAA*momBB + dadBB*momAA)*kidAB + dadBB*momBB*kidBB
	
	normFact<-a1+a2+a3+a4+a567+a8910
	
	usedTrios=rowSums((a1+a2+a3+a4+a567)/normFact, na.rm=TRUE)
	pMendelErr<-rowSums(1-normFact)
	
    num <- rowSums((a2 + a4 + a67)/normFact, na.rm=TRUE)
    denom <- rowSums((a1 + a2 + a3 + a4 + 2 * a567)/normFact, na.rm=TRUE)
    return(list(num = num, denom = denom, usedTrios=usedTrios, pMendelErr=pMendelErr))
}

# basically just a copy of fastTDTdomSplit with dimensions of geno changed
probTDTdomSplit<-function (geno, size = 50) 
{
    n.snp <- nrow(geno)
    int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
    dmat <- matrix(0, n.snp, 4)
	usedTrios <- pMendelErr <- numeric(n.snp)
    for (i in 1:(length(int) - 1)){
		tmp<-probTDTdomChunk(geno[int[i]:(int[i + 1] - 1), , drop = FALSE])
		dmat[int[i]:(int[i + 1] - 1), ] <- tmp$dmat
		usedTrios[int[i]:(int[i + 1] - 1)] <- tmp$usedTrios
		pMendelErr[int[i]:(int[i + 1] - 1)] <- tmp$pMendelErr
	}
    rownames(dmat) <- rownames(geno)
    h <- (dmat[, 1]/3 - dmat[, 2] + dmat[, 3] - dmat[, 4]/3)/(2 * 
        (dmat[, 1] + dmat[, 3]))
    tmp <- (dmat[, 2] + dmat[, 4])/(3 * (dmat[, 1] + dmat[, 3])) + 
        h * h
    or <- sqrt(tmp) - h
    beta <- log(or)
    tmp <- (dmat[, 2] + dmat[, 1]) * or/(or + 1)^2 + (dmat[, 
        3] + dmat[, 4]) * or/(3 * (or + 1/3)^2)
    se <- sqrt(1/tmp)
    stat <- beta/se
    stat <- stat * stat
    lower <- exp(beta - qnorm(0.975) * se)
    upper <- exp(beta + qnorm(0.975) * se)
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    if (is.null(rownames(geno))) 
        names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- names(pMendelErr) <- paste("SNP", 
            1:nrow(geno), sep = "")
    else names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- names(pMendelErr) <- rownames(geno)
    out <- list(coef = beta, se = se, stat = stat, pval = pval, 
        OR = exp(beta), lowerOR = lower, upperOR = upper, ia = FALSE, 
        type = "dominant", usedTrios = usedTrios, add = FALSE, pMendelErr = pMendelErr)
    class(out) <- "colTDT"
    out
}

# very similar to fastTDTdomChunk, but add a column to record additional Mendel-correct possibilities
probTDTdomChunk<-function (geno) 
{
    n.col <- ncol(geno)
    dadAA <- geno[,seq.int(1, n.col, 9), drop = FALSE]
    dadAB <- geno[,seq.int(2, n.col, 9), drop = FALSE]
    dadBB <- geno[,seq.int(3, n.col, 9), drop = FALSE]
    momAA <- geno[,seq.int(4, n.col, 9), drop = FALSE]
    momAB <- geno[,seq.int(5, n.col, 9), drop = FALSE]
    momBB <- geno[,seq.int(6, n.col, 9), drop = FALSE]
    kidAA <- geno[,seq.int(7, n.col, 9), drop = FALSE]
    kidAB <- geno[,seq.int(8, n.col, 9), drop = FALSE]
    kidBB <- geno[,seq.int(9, n.col, 9), drop = FALSE]
    dmat <- matrix(0, nrow(dadAA), 4)

	het1<-dadAA*momAB + dadAB*momAA
	a1<-het1*kidAA
	a2<-het1*kidAB
	
	hetHom<-dadAB*momBB + dadBB*momAB
	a3<-hetHom*kidAB
	a4<-hetHom*kidBB
	
	a567<-dadAB*momAB

	a5<-a567*kidAA

	# here, really want a6 + a7, not a6+2a7 like in additive case
	a67<-a567*(kidAB + kidBB)

	a8910<-dadAA*momAA*kidAA + (dadAA*momBB + dadBB*momAA)*kidAB + dadBB*momBB*kidBB

	normFact<-a1+a2+a3+a4+a567+a8910

    dmat[, 1] <- rowSums(a1/normFact, na.rm = TRUE)
    dmat[, 2] <- rowSums(a2/normFact, na.rm = TRUE)
    dmat[, 3] <- rowSums(a5/normFact, na.rm = TRUE)
    dmat[, 4] <- rowSums(a67/normFact, na.rm = TRUE)
	
	usedTrios<-rowSums((a1+a2+a5+a67)/normFact, na.rm=TRUE)
	pMendelErr<-rowSums(1-normFact)

    return(list(dmat=dmat, usedTrios=usedTrios, pMendelErr=pMendelErr))
}

probTDTrecSplit<-function (geno, size = 50) 
{
    n.snp <- nrow(geno)
    int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))
    rmat <- matrix(0, n.snp, 4)
	usedTrios <- pMendelErr <- numeric(n.snp)
    for (i in 1:(length(int) - 1)){
		tmp<-probTDTrecChunk(geno[int[i]:(int[i + 1] - 1), ,drop = FALSE])
		rmat[int[i]:(int[i + 1] - 1), ] <- tmp$rmat
		usedTrios[int[i]:(int[i + 1] - 1)] <- tmp$usedTrios
		pMendelErr[int[i]:(int[i + 1] - 1)] <- tmp$pMendelErr
    }
	rownames(rmat) <- colnames(geno)
    h <- (3 * rmat[, 1] - rmat[, 2] + rmat[, 3] - 3 * rmat[, 
        4])/(2 * (rmat[, 1] + rmat[, 3]))
    tmp <- 3 * (rmat[, 2] + rmat[, 4])/(rmat[, 1] + rmat[, 3]) + 
        h * h
    or <- sqrt(tmp) - h
    beta <- log(or)
    tmp <- (rmat[, 1] + rmat[, 2]) * or/(or + 1)^2 + 3 * (rmat[, 
        3] + rmat[, 4]) * or/(or + 3)^2
    se <- sqrt(1/tmp)
    stat <- (beta/se)^2
    lower <- exp(beta - qnorm(0.975) * se)
    upper <- exp(beta + qnorm(0.975) * se)
    pval <- pchisq(stat, 1, lower.tail = FALSE)
    if (is.null(rownames(geno))) 
        names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- paste("SNP", 
            1:nrow(geno), sep = "")
    else names(beta) <- names(stat) <- names(pval) <- names(usedTrios) <- rownames(geno)
    out <- list(coef = beta, se = se, stat = stat, pval = pval, 
        OR = exp(beta), lowerOR = lower, upperOR = upper, ia = FALSE, 
        type = "recessive", usedTrios = usedTrios, add = FALSE, pMendelErr=pMendelErr)
    class(out) <- "colTDT"
    out
}

probTDTrecChunk<-function (geno) 
{
    n.col <- ncol(geno)
    dadAA <- geno[,seq.int(1, n.col, 9), drop = FALSE]
    dadAB <- geno[,seq.int(2, n.col, 9), drop = FALSE]
    dadBB <- geno[,seq.int(3, n.col, 9), drop = FALSE]
    momAA <- geno[,seq.int(4, n.col, 9), drop = FALSE]
    momAB <- geno[,seq.int(5, n.col, 9), drop = FALSE]
    momBB <- geno[,seq.int(6, n.col, 9), drop = FALSE]
    kidAA <- geno[,seq.int(7, n.col, 9), drop = FALSE]
    kidAB <- geno[,seq.int(8, n.col, 9), drop = FALSE]
    kidBB <- geno[,seq.int(9, n.col, 9), drop = FALSE]
    rmat <- matrix(0, nrow(dadAA), 4)


	het1<-dadAA*momAB + dadAB*momAA
	a1<-het1*kidAA
	a2<-het1*kidAB
	
	hetHom<-dadAB*momBB + dadBB*momAB
	a3<-hetHom*kidAB
	a4<-hetHom*kidBB
	
	
	a56<-dadAB*momAB*(kidAB +kidAA)
	
	a7<-dadAB*momAB*(kidBB)

	a8910<-dadAA*momAA*kidAA + (dadAA*momBB + dadBB*momAA)*kidAB + dadBB*momBB*kidBB

	normFact<-a1+a2+a3+a4+a56+a7+a8910
	
	rmat[,1]<-rowSums(a3/normFact, na.rm=TRUE)
	rmat[,2]<-rowSums(a4/normFact, na.rm=TRUE)
	rmat[,3]<-rowSums(a56/normFact, na.rm=TRUE)
	rmat[,4]<-rowSums(a7/normFact, na.rm=TRUE)

	usedTrios<-rowSums((a3+a4+a56+a7)/normFact, na.rm=TRUE)
	pMendelErr<-rowSums(1-normFact)

    return(list(rmat=rmat, usedTrios=usedTrios, pMendelErr=pMendelErr))
}
