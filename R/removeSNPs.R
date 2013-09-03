removeSNPs <- function(geno, maf=NA, perc.na=NA){
	if(!is.na(maf) && (maf<0 | maf>0.2))
		stop("maf must be between 0 and 0.2.")
	if(!is.na(perc.na) && (perc.na<=0 | perc.na>1))
		stop("perc.na must be between 0 and 1 (excluding 0).")
	if(!is.matrix(geno) & !is.data.frame(geno))
		stop("geno must be a matrix or a data frame.")
	if(is.null(rownames(geno)))
		stop("geno must have row names.")
	if(nrow(geno)%%3 != 0)
		stop("The number of rows of geno is not dividable by 3.")
	if(all(colnames(geno)[1:2] == c("famid", "pid"))){
		fam.info <- geno[,1:2]
		geno <- as.matrix(geno[,-(1:2)])
		info <- TRUE
	}
	else{
		if(!is.matrix(geno))
			stop("geno must be a matrix.")
		info <- FALSE
	}
	if(any(!geno %in% c(0:2, NA)))
		stop("The genotypes must be coded by 0, 1, 2; missing by NA.")
	if(!is.na(maf)){
		tmp.geno <- geno[-seq(3,nrow(geno),3),]
		tmp.geno <- tmp.geno[!duplicated(rownames(tmp.geno)),]
		obs.maf <- colSums(tmp.geno, na.rm=TRUE) / (2*colSums(!is.na(tmp.geno)))
		if(maf==0)
			geno <- geno[,obs.maf>0]
		else
			geno <- geno[,obs.maf>=maf]
	}
	if(!is.na(perc.na)){
		num.na <- colMeans(is.na(geno))
		geno <- geno[,num.na < perc.na]
	}
	if(info)
		geno <- data.frame(fam.info, geno, stringsAsFactors=FALSE)
	geno
}

orderSNPs <- function(geno, map, snp="SNP", orderBy=c("Chr", "Position")){
	if(!is.data.frame(map))
		stop("map must be a data.frame.")
	if(any(!c(snp,orderBy) %in% colnames(map)))
		stop("Some of the columns specified by snp and orderBy are not available in map.")
	if(is.factor(map[,snp]))
		stop("map contains factors. Please set stringsAsFactors=FALSE when reading map into R.")
	if(!is.numeric(map[,orderBy[1]]) | !is.numeric(map[,orderBy[2]]))
		stop("The columns specified by orderby must be numeric.")
	cn <- tolower(gsub("IL", "", colnames(geno)))
	marker <- tolower(map[,snp])
	if(any(!cn %in% c(marker, "famid", "pid")))
		stop("For some of the SNPs in geno, no annotations is available in map.")
	map <- map[marker %in% cn, ]
	map <- map[order(map[,orderBy[1]], map[,orderBy[2]]), ]
	cn.ordered <- tolower(map[,snp])
	if(any(cn %in% c("famid", "pid")))
		cn.ordered <- c("famid", "pid", cn.ordered)
	geno[,cn.ordered]
}

removeTrios <- function(geno, perc.na=1){
	if(!is.na(perc.na) & (perc.na<=0 | perc.na>1))
		stop("perc.na must be between 0 and 1 (excluding 0).")
	if(!is.matrix(geno) & !is.data.frame(geno))
		stop("geno must be a matrix or a data frame.")
	if(is.null(rownames(geno)))
		stop("geno must have row names.")
	if(nrow(geno)%%3 != 0)
		stop("The number of rows of geno is not dividable by 3.")
	ids.snps <- !colnames(geno) %in% c("famid", "pid")
	rs.na <- rowMeans(is.na(geno[,ids.snps]))
	trio <- rep(1:(nrow(geno)/3), e=3)
	out.trio <- trio[rs.na >= perc.na]
	geno[!trio %in% out.trio,]
}

	
	
	
	