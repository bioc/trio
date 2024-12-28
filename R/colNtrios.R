colNtrios <- function(mat.snp, env = NULL, onlyContributing = FALSE, famid = NULL, size = 50){
	checkMatSNP(mat.snp, size = size)
	if(!is.null(env))
		out <- ntrioGxE(mat.snp, env, only = onlyContributing, famid = famid, size = size)
	else
		out <- ntrioTDT(mat.snp, only = onlyContributing, size = size)
	out
}


ntrios2Dom <- function(matNumber, check = TRUE, quiet = FALSE){
	if(is.list(matNumber)){
		if(length(matNumber) != 2)
			stop("matNumber does not seem to be the output of colNtrios.")
		lnames <- names(matNumber)
		if(is.null(lnames) | any(lnames != c("ntriosEnv0", "ntriosEnv1")))
			stop("matNumber does not seem to be the output of colNtrios.")
		out <- lapply(matNumber, ntrios2Dom, check = check, quiet = quiet)
		return(out)
	}
	if(check)
		checkMatN(matNumber)
	matN <- cbind(D0011 = matNumber[ , 1], D1001 = matNumber[ , 2], D0111 = matNumber[ , 5],
		D1011 = matNumber[ , 6] + matNumber[ , 7])
  	if(ncol(matNumber) == 10){
		n1111 <- rowSums(matNumber[ , c(3, 4, 8, 9)])
		matN <- cbind(matN, D0000 = matNumber[ , 10], D1111 = n1111)
	}
	else{
		if(!quiet)
			message("Since only the numbers of trios contributing to the likelihood estimation ",
				"were contributed,\n", "also for the dominant model only the numbers of ",
				"trios are computed contributing\n", "to the estimation in this model.\n")
	}
	matN
}   	

ntrios2Rec <- function(matNumber, check = TRUE, quiet = FALSE){
	if(is.list(matNumber)){
		if(length(matNumber) != 2)
			stop("matNumber does not seem to be the output of colNtrios.")
		lnames <- names(matNumber)
		if(is.null(lnames) | any(lnames != c("ntriosEnv0", "ntriosEnv1")))
			stop("matNumber does not seem to be the output of colNtrios.")
		out <- lapply(matNumber, ntrios2Rec, check = check, quiet = quiet)
		return(out)
	}
	if(check)
		checkMatN(matNumber)
	matN <- cbind(R0011 = matNumber[ , 3], R1001 = matNumber[ , 4], 
		R0001 = matNumber[ , 5] + matNumber[ , 6], R1000 = matNumber[ , 7])
  	if(ncol(matNumber) == 10){
		n0000 <- rowSums(matNumber[ , c(1, 2, 8, 10)])
		matN <- cbind(matN, R1111 = matNumber[ , 9], R0000 = n0000)
	}
	else{
		if(!quiet)
			message("Since only the numbers of trios contributing to the likelihood estimation ",
				"were contributed,\n", "also for the recessive model only the numbers of ",
				"trios are computed contributing\n", "to the estimation in this model.\n")
	}
	matN
} 


ntrioTDT <- function(mat.snp, only = FALSE, size = 50){
	n.snp <- ncol(mat.snp)
	int <- unique(c(seq.int(1, n.snp, size), n.snp + 1))	
	matNumber <- matrix(nrow = n.snp, ncol = ifelse(only, 7, 10))
	for(i in 1:(length(int) - 1)){
		matNumber[int[i]:(int[i + 1] - 1), ] <- ntrioChunk(mat.snp[ , int[i]:(int[i+1]-1), drop = FALSE],
			only = only)
	}
	rownames(matNumber) <- colnames(mat.snp)
	cn <- c("G010", "G011", "G121", "G122", "G110", "G111", "G112")
	if(!only)
		cn <- c(cn, "G021", "G222", "G000")
	colnames(matNumber) <- cn
	matNumber
}


ntrioChunk <- function(geno, only = FALSE){
	n.row <- nrow(geno)
	matN <- matrix(NA, nrow = ncol(geno), ncol = ifelse(only, 7, 10))
	dad <- geno[seq.int(1, n.row, 3), , drop = FALSE]
 	mom <- geno[seq.int(2, n.row, 3), , drop = FALSE]
	kid <- geno[seq.int(3, n.row, 3), , drop = FALSE]
	het1 <- (dad == 1) & (mom == 1)
	matN[ , 5] <- colSums(het1 & kid == 0, na.rm = TRUE)
	matN[ , 6] <- colSums(het1 & kid == 1, na.rm = TRUE)
	matN[ , 7] <- colSums(het1 & kid == 2, na.rm = TRUE)
	mom <- mom + dad
	het1 <- mom == 1
	matN[ , 1] <- colSums(het1 & kid == 0, na.rm = TRUE)
	matN[ , 2] <- colSums(het1 & kid == 1, na.rm = TRUE)
	het1 <- mom == 3
	matN[ , 3] <- colSums(het1 & kid == 1, na.rm = TRUE)
	matN[ , 4] <- colSums(het1 & kid == 2, na.rm = TRUE)
	if(!only){
		matN[ , 9] <- colSums(mom == 4 & kid == 2, na.rm = TRUE)
		matN[ , 10] <- colSums(mom == 0 & kid == 0, na.rm = TRUE)
		het1 <- mom == 2 & dad != 1
		matN[ , 8] <- colSums(het1 & kid == 1, na.rm = TRUE)
	}		 
	matN
}

ntrioGxE <- function(mat.snp, env, only = FALSE, famid = NULL, size = 50){
	n.trio <- nrow(mat.snp) / 3
	env2 <- checkEnv(env, n.trio, rownames(mat.snp), famid = famid)
	if(any(is.na(env))){
		mat.snp <- mat.snp[!is.na(env2),]
		env2 <- env2[!is.na(env2)]
	}
	ntrios0 <- ntrioTDT(mat.snp[env2 == 0, ], only = only, size = size)
	ntrios1 <- ntrioTDT(mat.snp[env2 == 1, ], only = only, size = size)
	list(ntriosEnv0 = ntrios0, ntriosEnv1 = ntrios1)
}









