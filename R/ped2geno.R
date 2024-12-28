ped2geno <- function(ped, snpnames = NULL, coded = c("12", "AB", "ATCG", "1234"), 
		naVal = 0, cols4ID = FALSE){
	cn <- colnames(ped)
	if(any(tolower(cn[1:6]) != c("famid", "pid", "fatid", "motid", "sex","affected")))
		stop("The names of the first six columns must be\n",
			"famid, pid, fatid, motid, sex, affected.")
	if(any(duplicated(ped$pid)))
		stop("The IDs of the individuals in the second column of ped must be unique.")
	if(any(ped$fatid == 0 & ped$motid != 0) | any(ped$fatid != 0 & ped$motid == 0))
		stop("In some of the rows, fatid is equal to zero and matid differs from zero,\n",
			"or vice versa.")
	mat.snp <- as.matrix(ped[,-(1:6)])
	rownames(ped) <- ped$pid
	if(ncol(mat.snp) %% 2 == 1)
		stop("ped must contain six columns for the information on the trios\n",
			"and two columns for each SNP.")
	n.snp <- ncol(mat.snp) / 2
	if(!is.null(snpnames) && n.snp != length(snpnames))
		stop("The length of snpnames is not equal to the number of SNPs in ped.")
	if(length(naVal)!=1)
		stop("naVal must consist of exactly one value.")
	coded <- match.arg(coded)
	allele <- unlist(strsplit(coded, ""))
	if(any(!mat.snp %in% c(allele, naVal)))
		stop("The SNP columns of ped contain other values than ", 
			paste(allele, collapse=", "), " and ", naVal, ".") 
	mat.snp[mat.snp==naVal] <- NA
	mat.allele <- matrix(0, ncol(mat.snp), length(allele))
	if(coded %in% c("12", "1234"))
		allele <- as.numeric(allele)
	if(!is.na(naVal) && any(naVal==allele))
		stop("naVal cannot be one of the letters/numbers coding for the alleles.")
	for(i in 1:length(allele))
		mat.allele[,i] <- colSums(mat.snp==allele[i], na.rm=TRUE)
	mat.count <- mat.allele[seq(1, 2*n.snp, 2),, drop=FALSE] + 
		mat.allele[seq(2, 2*n.snp, 2),, drop=FALSE]
	rsmc0 <- rowSums(mat.count != 0)
	if(any(rsmc0 > 2))
		stop("Each SNP must show at most two alleles.")
	if(any(rsmc0 < 2)){
		ids.rs2 <- rsmc0 == 2
		mat.count <- mat.count[ids.rs2,]
		seq.in <- rep(ids.rs2, e=2)
		mat.allele <- mat.allele[seq.in,]
		mat.snp <- mat.snp[,seq.in]
		n.snp <- ncol(mat.snp) / 2
		snpnames <- snpnames[ids.rs2]
		warning(sum(!ids.rs2), " of the SNPs are monomorph. These SNPs are removed.")
	}
	mat.recoded <- matrix(NA, nrow(mat.snp), ncol(mat.snp))
	if(length(allele) == 2){
		ids.major <- rep(max.col(mat.count), e=2)
		if(any(ids.major==1)){
			mat.recoded[,ids.major==1][mat.snp[,ids.major==1] == allele[1]] <- 0
			mat.recoded[,ids.major==1][mat.snp[,ids.major==1] == allele[2]] <- 1
		}
		if(any(ids.major==2)){
			mat.recoded[,ids.major==2][mat.snp[,ids.major==2] == allele[2]] <- 0
			mat.recoded[,ids.major==2][mat.snp[,ids.major==2] == allele[1]] <- 1
		}
	}
	else{
		ids.max <- max.col(mat.count)
		ids.major <- rep(ids.max, e=2)
		for(i in 1:4){
			if(any(ids.major==i)){
				mat.recoded[,ids.major==i][mat.snp[,ids.major==i] == allele[i]] <- 0
				mat.count[ids.max==i, i] <- 0
			}
		}
		ids.major <- rep(max.col(mat.count), e=2)
		for(i in 1:4){
			if(any(ids.major==i))
				mat.recoded[,ids.major==i][mat.snp[,ids.major==i] == allele[i]] <- 1
		}
	}	
	rownames(mat.recoded) <- ped$pid			  
	
	mat.kid <- ped[ped$fatid != 0 & ped$motid != 0, c(3, 4, 2)]
	mat.ids <- as.matrix(mat.kid)
	vec.ids <- as.vector(t(mat.ids))
	if(any(!vec.ids %in% ped$pid))
		stop("Data for some of the moms or dads seem to be missing.")
	mat.trio <- mat.recoded[vec.ids, seq(1, 2*n.snp, 2)] + 
		mat.recoded[vec.ids, seq(2, 2*n.snp, 2)]
	if(cols4ID){
		if(is.data.table(ped)){
			ped <- as.data.frame(ped)
			rownames(ped) <- ped$pid
		}
		mat.trio <- data.frame(ped[vec.ids, 1:2], mat.trio, stringsAsFactors = FALSE)
		colnames(mat.trio)[-(1:2)] <- if(!is.null(snpnames)) snpnames
			else paste("SNP", 1:n.snp, sep="")
	}
	else{
		if(is.matrix(mat.trio))
			colnames(mat.trio) <- if(!is.null(snpnames)) snpnames
				else paste("SNP", 1:n.snp, sep="")
	} 
	mat.trio
}
			 
		