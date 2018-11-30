getYmat <- function(n.trio, n.perm){
	mat <- matrix(0, 16*n.trio, n.perm)
	begin <- (0:(n.trio-1)) * 16
	for(i in 1:n.perm){
		ids <- sample(16, n.trio, TRUE)
		mat[begin + ids, i] <- 1
	}
	mat
}

getDummyX <- function(mat.snp, n.snp){
	mat.code <- matrix(c(0,0,0,0, 0,0,1,1, 0,0,1,1, 1,0,0,1, 1,0,0,1, 1,1,1,1, 
		1,1,1,1, 2,2,2,2, 1,1,2,2, 1,1,2,2, 2,1,1,2, 2,1,1,2, 0,1,1,2, 
		1,0,1,2, 2,0,1,1, NA,NA,NA,NA), nrow=4)
	cn <- c("000", "010", "100", "011", "101", "021", "201", "222", "121", "211",
		"122", "212", "110", "111", "112", "NANANA")
	colnames(mat.code) <- cn
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
	mat.pseudo - 1
}	

compPermTDT2way <- function(...)
	cat("compPermTDT2way has been renamed to colGxGPerms. So please use colGxGPerms instead.\n")

colGxGPerms <- function(mat.snp, n.perm=1000, genes=NULL, col.out=NULL,
		warnError=TRUE, verbose=TRUE, rand=NA){
	requireNamespace("survival", quietly=TRUE)
	if(!is.null(col.out)){
		if(!is(col.out, "colTDTepi"))
			stop("col.out must be the output of colTDT2way with epistatic=TRUE.")
		if(!is.null(genes))
			stop("genes should not be specified if col.out is specified.")
		genes <- col.out$ind.genes
	}
	if(!is.matrix(mat.snp))
		stop("mat.snp has to be a matrix.")
	if(nrow(mat.snp) %% 3 != 0)
		stop("mat.snp does not seem to contain trio data, as its number of rows is\n",
			"not dividable by 3.")
	if(any(!mat.snp %in% c(0, 1, 2, NA)))
		stop("The values in mat.snp must be 0, 1, and 2.")
	n.snp <- ncol(mat.snp)
	n.trio <- nrow(mat.snp) / 3
	x <- getDummyX(mat.snp, n.snp)
	z <- (x==0) - 0.5
	mat.ids <- cbind(rep.int(rep(1:4, e=4), n.trio), rep.int(1:4, 4*n.trio))
	mat.ids <- mat.ids + 4 * rep(0:(n.trio-1), e=16)
	if(is.null(genes))
		combs <- allCombs(n.snp)
	else{
		if(!is.character(genes))
			stop("genes must be a vector of character strings.")
		if(length(genes) != n.snp)
			stop("The length of genes must be equal to the number of columns of mat.snp.")
		ids.genes <- as.numeric(as.factor(genes))
		combs <- allBetweenCombs(ids.genes)
	}
	sn <- if(is.null(colnames(mat.snp))) paste("SNP", 1:n.snp, sep="") else colnames(mat.snp)
	ianames <- paste(sn[combs[,1]], sn[combs[,2]], sep=" : ")
	if(is.null(col.out)){
		if(verbose)
			cat("Computing the test statistics for the", nrow(combs), "interactions.\n",
				"Computation started at", date(),"\n")
		stat <- getOriginalStat(x, z, mat.ids, combs, n.trio, warnError=warnError)
		names(stat) <- ianames
		if(verbose)
			cat("and ended at ", date(), ".\n", sep="") 
	}
	else{
		stat <- col.out$stat
		if(any(names(stat) != ianames))
			stop("The names of the interactions in col.out differ from the interaction names\n",
				"generated from mat.snp.")
		if(any(is.na(stat)))
			warning("Some of the original values of the test statistic are missing values.\n",
				"The permuted values of the test statistic will not be computed for the corresponding interactions.",
				call.=FALSE)
	}
	if(!is.na(rand))
		set.seed(rand)
	mat.y <- getYmat(n.trio, n.perm)
	mat.perm <- getPermStat(x, z, mat.ids, combs, n.trio, mat.y, stat, warnError=warnError, verbose=verbose)
	rownames(mat.perm) <- ianames
	return(list(stat=stat, permStat=mat.perm, y.perm=mat.y))
}		
		
	
getOriginalStat <- function(x, z, mat.ids, combs, n.trio, warnError=TRUE){
	ll.main <- ll.full <- numeric(nrow(combs))
	y <- rep.int(c(1, rep.int(0, 15)), n.trio)
	strat <- rep(1:n.trio, e=16)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
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
	if(warnError)
		options(warn=wa)
	if(any(is.na(ll.full)) | any(is.na(ll.main)))
		warning("The fitting of some of the models has failed. A possible reason\n",
			"is that the two respective SNPs might be in (strong) LD.\n",
			"The corresponding test statistics are thus set to NA\n",
			"and the permuted test statistics will not be computed for this SNP interaction.",
			call.=FALSE)
	-2 * (ll.main - ll.full)
}
	
getPermStat <- function(x, z, mat.ids, combs, n.trio, mat.y, stat, warnError=TRUE, verbose=TRUE){
	n.combs <- nrow(combs)
	matLLmain <- matLLfull <- matrix(NA, n.combs, ncol(mat.y))
	strat <- rep(1:n.trio, e=16)
	if(warnError){
		wa <- options()$warn
		on.exit(options(warn=wa))
		options(warn=2)
	}
	for(i in 1:n.combs){
		if(!is.na(stat[i])){
			if(verbose)
				cat("Computing permuted test scores for ", names(stat)[i], "  (", i, " of ",
					n.combs, " interactions).\n", "Started at ", date(), ".\n", sep="")
			x1 <- x[mat.ids[,1], combs[i,1]]
			z1 <- z[mat.ids[,1], combs[i,1]] 
			x2 <- x[mat.ids[,2], combs[i,2]]
			z2 <- z[mat.ids[,2], combs[i,2]]
			for(j in 1:ncol(mat.y)){
				y <- mat.y[,j]
				woIA <- try(clogit(y ~ x1 + z1 + x2 + z2 + strata(strat)), silent=TRUE)
				full <- try(clogit(y ~ x1 + z1 + x2 + z2 + x1*x2 + x1*z2 + z1*x2 + 
					z1*z2 + strata(strat)), silent=TRUE)
				matLLmain[i, j] <- if(is(woIA, "try-error")) NA else woIA$loglik[2]
				matLLfull[i, j] <- if(is(full, "try-error")) NA else full$loglik[2]
			}
		}
	}
	if(warnError)
		options(warn=wa)
	-2 * (matLLmain - matLLfull)
}
