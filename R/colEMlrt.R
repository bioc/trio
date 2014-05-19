colEMlrt <- function(mat.snp, model = c("general", "dominant", "recessive"), maternal=FALSE,
		parentMissing = c("father", "mother", "either"),  iter=40, eps=10^-16){
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
  	if(iter < 5)
    		stop("iter should be at least 5.")
  	if(eps <= 0)
    		stop("eps must be larger than zero.")
  	if(eps >= 10^-4)
    		stop("eps must be smaller than 0.0001.")
  
  	type <- match.arg(model)
  	miss <- match.arg(parentMissing)
  
  	em.step <- match.fun(switch(miss,
		"either"   = em.step.both,
		"father" = em.step.father,
		"mother" = em.step.mother))
  
  	# number of SNPs, initialized loglikelihood vector
  	n.snp <- ncol(mat.snp)
  	ll.red <- ll.full <- rep.int(NA, n.snp)
  
  	# Indicator matrix for stratified regression model 
  	if(maternal) 
    		X.null <- matrix(c(rep(c(0,1,0), c(0,1,14)), c(0,1,1,0,1,0,1,rep(0,8)),
			c(0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), rep(c(0,1,0,1,0), c(5,1,3,1,5)), 
			c(0,0,0,0,0,0,0,0,1,0,1,0,1,1,0), c(0,0,1,1,0,0,1,1,1,0,0,1,1,0,0),
			c(0,0,0,0,0,0,0,0,0,1,1,0,0,1,1)), 15)
  	else
    		X.null <- matrix(c(rep(c(0,1,0), c(0,1,14)), c(0,1,1,0,1,0,1,rep(0,8)),
			c(0,0,0,1,0,0,0,1,0,0,0,1,0,0,0), rep(c(0,1,0,1,0), c(5,1,3,1,5)), 
			c(0,0,0,0,0,0,0,0,1,0,1,0,1,1,0)), 15)
  
  
  	# type specific design matrix 
  	if(type == "general"){
    		Z.full <- matrix(c(rep(0:1, c(11,4)), rep(c(0,1,0), c(4,7,4))), 15)
    		df <- 2
  	}
  	else if(type == "recessive"){
    		Z.full <- rep(0:1, c(11,4))
    		df <- 1
  	}
  	else if(type == "dominant"){
    		Z.full <- rep(0:1, c(4,11))
    		df <- 1
  	}
  
  	X.full <- cbind(X.null, Z.full)
  
  	# offset
  	modelOffset <- log(2)*rep(c(0,1,0), c(7,1,7))
  
  
  	##  snp-wise
  	for(i in 1:n.snp){
    		# in matrix
    		dat <- matrix(mat.snp[,i], 3)
    		dat[is.na(dat)] <- 27
    
    		## every number belongs to one combination
    		## impossible or unwanted combinations are discarded automatically (without warning)
    		tab <- table(c(c(0:1,3:4,10:16,22:23,25:27,30,36,39,42,48,51,81,82,90,91,92,100,101), colSums(dat*c(1,3,9)))) - 1  
    		tab <- tab[as.character(c(0:1,3:4,10:16,22:23,25:26,  # complete trios
			27,30,36,39,42,48,51,       # missing fathers
			81,82,90,91,92,100,101))]   # missing mothers
    
    		# complete trios
    		# 000, 100, 010, 110, 101, 201, 011, 111, 211, 021, 121, 112, 212, 122, 222
    		N <- as.vector(tab[1:15])
    
    		# hier muessten wir noch was einbauen (oder besser ausserhalb der for-Schleife), das den
    		# EM-Algorithmus erst gar nicht startet, wenn es keine fehlenden Werte gibt
    
    		# father missing
    		# M00, M00, M10, M10, M01, M01, M11, M11, M11, M21, M21, M12, M12, M22, M22
    		M <- c(rep(tab[16:22], c(2,2,2,3,2,2,2)))
    
    		# mother missing
    		# 0M0, 1M0, 0M0, 1M0, 1M1, 2M1, 0M1, 1M1, 2M1, 0M1, 1M1, 1M2, 2M2, 1M2, 2M2
    		L <- as.vector(tab[c(23,24,23,24,26,27,25,26,27,25,26,28,29,28,29)])   # possible without as.vector
    
    		# working copy, correcting for null cells
    		P <- switch(miss, "either" = L+M, "father" = M, "mother" = L)
    
    		Y <- N + P*rep(c(0.5, 1/3, 0.5), c(6,3,6))  # expected cell count, sums up to number of used trios
    		Y[Y == 0] <- 0.01  # weglassen? Appendix 2 aus Weinberg Paper
    
    		# EM algorithm
    		for(j in 1:iter){
      			# Update
      			Ynew <- em.step(Y=Y, N=N, M=M, L=L)
      
      			# MLE step omitted, cancels out in update step
      			# Y estimated cell count
      
     	 		tmpDiff <- sum((Y-Ynew)^2, na.rm=TRUE)
      			Y <- Ynew
      			if(tmpDiff < eps)
        			break
      
    		}
    		Y[is.na(Y)] <- 0
    
    		# Poisson regression
    		ll.full[i] <- glm(Y ~ X.full + offset(modelOffset), family=quasipoisson(link = log))$deviance
    		ll.red[i] <- glm(Y ~ X.null + offset(modelOffset), family=quasipoisson(link = log))$deviance
  	}
  	stat <- ll.red - ll.full
  	pval <- pchisq(stat, df, lower.tail=FALSE)
  
  	names(ll.red) <- names(ll.full) <- names(pval) <- names(stat) <- colnames(mat.snp)
  
  	out <- list(ll.red=ll.red/-2, ll.full=ll.full/-2, stat=stat, pval=pval, model=type, missing=miss)
  	class(out) <- "colEMlrt"
 	out
}

print.colEMlrt <- function(x, top=5, digits=4, ...){
	pval <- format.pval(x$pval, digits=digits)
  	out <- data.frame("LL (Full Model)" = x$ll.full, "LL (Reduced Model)" = x$ll.red, Statistic=x$stat,
		"P-Value"=pval, check.names=FALSE)
  	type <- switch(x$model, general="General", dominant="Dominant", recessive="Recessive")
  	miss <- switch(x$missing, either="Either Mother or Father", mother="Mother", father="Father")
	cat("       EM Likelihood Ratio Test", "\n\n", "Model: ", type, "\n", "Allowed Missing: ",
		miss, "\n\n", sep="")
  	if(top > 0 && length(x$ll.full) > top){
    		ord <- order(x$pval)[1:top]
    		out <- out[ord,]
    		cat("Top", top, "SNPs:\n")
  	}
 	 print(format(out, digits=digits))
}

em.step.both <- function(Y, N, M, L){
	N + M*(Y/( rep(c(Y[1]+Y[2], Y[3]+Y[4], Y[5]+Y[6], sum(Y[7:9]), Y[10]+Y[11], Y[12]+Y[13], Y[14]+Y[15]),
		c(2,2,2,3,2,2,2)))) + L*(Y/( c(Y[1]+Y[3], Y[2]+Y[4], Y[1]+Y[3], Y[2]+Y[4],
		Y[5]+Y[8]+Y[11], Y[6]+Y[9], Y[7]+Y[10], Y[5]+Y[8]+Y[11], Y[6]+Y[9], Y[7]+Y[10],
		Y[5]+Y[8]+Y[11], Y[12]+Y[14], Y[13]+Y[15], Y[12]+Y[14], Y[13]+Y[15])))
}


em.step.father <- function(Y, N, M, L){
	N + M*(Y/( rep(c(Y[1]+Y[2], Y[3]+Y[4], Y[5]+Y[6], sum(Y[7:9]), 
		Y[10]+Y[11], Y[12]+Y[13], Y[14]+Y[15]) ,c(2,2,2,3,2,2,2))))
}

em.step.mother <- function(Y, N, M, L){
  	N + L*(Y/( c(Y[1]+Y[3], Y[2]+Y[4], Y[1]+Y[3], Y[2]+Y[4], Y[5]+Y[8]+Y[11], Y[6]+Y[9], Y[7]+Y[10], 
		Y[5]+Y[8]+Y[11], Y[6]+Y[9], Y[7]+Y[10], Y[5]+Y[8]+Y[11], Y[12]+Y[14], Y[13]+Y[15],
		Y[12]+Y[14], Y[13]+Y[15])))
}
