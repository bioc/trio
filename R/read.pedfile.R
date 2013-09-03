read.pedfile <- function(file, first.row=NA, coded=NULL, naVal=0, sep=" ", 
		p2g=FALSE, non.rs.IDs=FALSE, cols4ID=FALSE){
	if(!is.null(coded) && !coded %in% c("12", "AB", "1234", "ATCG"))
		stop("coded must be either '12', or 'AB', or '1234', or 'ATCG'.")
	if(is.na(first.row)){
		if(non.rs.IDs){
			rs <- readLines(file, n=2)
			tmp1 <- unlist(strsplit(rs[1], sep))
			tmp2 <- unlist(strsplit(rs[2], sep))
			first.row <- length(tmp1) == length(tmp2)
			rs <- rs[1]
		}
		else{
			rs <- readLines(file, n=1)
			first.row <- tolower(substring(rs, 1, 2)) != "rs"
		}
		read <- TRUE
		cat("NOTE: first.row has not been specified. Since the first row ",
			ifelse(first.row, "does not seem", "seems"), "\n", 
			"to contain the rs-IDs, first.row is set to ", first.row,
			".\n\n", sep="")
	}
	else
		read <- FALSE
	if(!first.row){
		if(!read)
			rs <- readLines(file, n=1)
		snpnames <- unlist(strsplit(tolower(rs), sep))
		if(!non.rs.IDs && any(substring(snpnames, 1, 2) != "rs"))
			stop("All SNP names must be rs-IDs, if non.rs.IDs=FALSE.")
		ped <- read.table(file, stringsAsFactors=FALSE, skip=1)
	}
	else{
		ped <- read.table(file, stringsAsFactors=FALSE)
		snpnames <- paste("SNP", 1:((ncol(ped) - 6) / 2), sep="")
	}
	if(any(sapply(ped, is.logical))){
		idsLogical <- which(sapply(ped, is.logical))
		for(i in idsLogical){
			if(any(!ped[,i]))
				stop("ped contains a logic variable with values TRUE and FALSE.")
			ped[,i] <- "T"
		}
	}
	n.snp <- length(snpnames)
	colnames(ped) <- c("famid", "pid", "fatid", "motid", "sex", "affected",
		paste(rep(snpnames, e=2), rep(1:2, n.snp), sep="."))
	ids.kid1 <- ped[,3] != 0
	ids.kid2 <- ped[,4] != 0
	if(any(ids.kid1 != ids.kid2))
		stop("The third and fourth column of file (containing fatid and motid)",
			"\n", "must both be either zero or non-zero.")
	if(any(duplicated(ped[,2]))){
		ped[,2] <- paste(ped[,1], ped[,2], sep="_")
		ped[ids.kid1, 3] <- paste(ped[ids.kid1, 1], ped[ids.kid1, 3], sep="_")
		ped[ids.kid2, 4] <- paste(ped[ids.kid2, 1], ped[ids.kid2, 4], sep="_")
		warning("Since the individual IDs in the second column are not unique,\n",
			"they are made unique by combining the first and second column.")
	}
	if(any(duplicated(ped[,2])))
		stop("Even after combining the first and second column, the individual IDs\n",
			"in the second column are not unique. Please make them unique in file.")
	if(!p2g)
		return(ped)
	if(is.null(coded)){
		ids.select <- sample(n.snp, min(n.snp, 20)) * 2 + 6
		tmpmat <- as.matrix(ped[,ids.select])
		tabnames <- names(table(tmpmat))
		if(all(tabnames %in% c(naVal, 1:4)))
			coded <- "1234"
		else if(all(tabnames %in% c(naVal, 1:2)))
			coded <- "12"
		else if(all(tabnames %in% c(naVal, "A", "B")))
			coded <- "AB"
		else if(all(tabnames %in% c(naVal, "A", "T", "C", "G")))
			coded <- "ATCG"
		else stop("It is not clear how the SNPs and how missing values are coded.\n",
			"Please specify coded and naVal.")
		cat("NOTE: Since coded has not been specified, it is set to \"", coded, "\".\n\n", sep="")
	}
	ped2geno(ped, snpnames=snpnames, coded=coded, naVal=naVal, cols4ID=cols4ID)
}

		   
