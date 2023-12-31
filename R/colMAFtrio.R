colMAFtrio <- function(geno, changeMinor=FALSE){
	if(!is.matrix(geno) & !is.data.frame(geno))
		stop("geno must be a matrix or a data frame.")
	if(is.null(rownames(geno)))
		stop("geno does not seem to be a genotype matrix, as it\n",
			"does not have row names.")
	if(nrow(geno)%%3 != 0)
		stop("The number of rows of geno is not dividable by 3.")
	mat <- as.matrix(geno[-seq(3,nrow(geno),3) ,!colnames(geno) %in% c("famid", "pid")])
	if(any(!mat %in% c(0:2, NA)))
		stop("The genotypes in geno must be coded by 0, 1, 2; and missing values by NA.") 
	mat <- mat[!duplicated(rownames(mat)),]
	maf <- colSums(mat, na.rm=TRUE) / (2 * colSums(!is.na(mat)))
	if(changeMinor && any(maf > 0.5)){
		ids <- which(maf > 0.5)
		maf[ids] <- 1 - maf[ids]
		warning("Some of the MAFs were larger than 0.5. These MAFs are replaced by 1 - MAF.")
	}
	maf
}


	