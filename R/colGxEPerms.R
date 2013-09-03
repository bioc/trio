colGxEPerms <- function(mat.snp, env, model=c("additive", "dominant", "recessive"), B=10000, size=20,
		addPerms=TRUE, famid=NULL, rand=NA){
	if (length(env) != nrow(mat.snp)/3) 
        	stop("The length of env must be equal to the number of trios in mat.snp.")
    	if (!is.null(famid)) 
        	env <- reorderEnv(env, famid, rownames(mat.snp))
    	if (any(is.na(env))) {
        	tmpEnv <- rep(env, e = 3)
        	mat.snp <- mat.snp[!is.na(tmpEnv), ]
        	env <- env[!is.na(env)]
    	}
	type <- match.arg(model)
	stat <- colGxE(mat.snp, env, model=type, size=size)$stat
	ids.in <- which(rowSums(is.na(stat)) == 0)
	tmpsplit <- rep(1:ceiling(ncol(mat.snp)/size), e=size)[1:length(ids.in)]
	vec.ids <- split(ids.in, tmpsplit)
	env2 <- rep(env, e=3)==1
	matG <- matGxE <- matrix(NA, ncol(mat.snp), B)
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:length(vec.ids)){
		tmp <- gxe.null.stats(mat.snp[,vec.ids[[i]], drop=FALSE], env2, type, B=B) 
		matG[vec.ids[[i]],] <- tmp$G
		matGxE[vec.ids[[i]],] <- tmp$GxE
	}
	pval <- cbind(G=rowMeans(stat[,1] <= matG, na.rm=TRUE), GxE=rowMeans(stat[,2] <= matGxE, na.rm=TRUE))
	rownames(stat) <- rownames(pval) <- colnames(mat.snp)
	if(!addPerms)
		return(list(stat=stat, pval=pval))
	rownames(matG) <- rownames(matGxE) <- colnames(mat.snp)
	structure(list(stat=stat, pval=pval, matPermG=matG, matPermGxE=matGxE))
}
	

gxe.null.stats <- function(geno, env2, type, B=10000){
	matA0 <- getNullStats(geno[!env2,,drop=FALSE], type)
	matA1 <- getNullStats(geno[env2,,drop=FALSE], type)
	ids.het0 <- matA0[,3] == 0
	ids.het1 <- matA1[,3] == 0
	null.fun <- match.fun(paste("gxeBetaV", substring(type, 1, 3), sep=""))
	G <- GxE <- matrix(NA, ncol(geno), B)
	for(i in 1:ncol(geno)){
		out <- getGxENullStats(matA0[i,], matA1[i,], B, null.fun)
		G[i,] <- out$stat0
		GxE[i,] <- out$stat1
	}
	return(list(G=G, GxE=GxE))
}


getGxENullStats <- function(a0, a1, perm, null.fun){
	out0 <- null.fun(a0, perm)
	out1 <- null.fun(a1, perm)
	beta0 <- out0$beta
	beta1 <- out1$beta - out0$beta
	stat0 <- beta0 * beta0 / out0$v
	stat1 <- beta1 * beta1 / (out0$v + out1$v)
	return(list(stat0=stat0, stat1=stat1))
}
	
gxeBetaVadd <- function(vec, perm){
	beta <- v <- NA
	if(vec[3]==0){
		a2 <- rbinom(perm, vec[1], 0.5)
		if(vec[2]==0){
			beta <- log(a2/(vec[1]-a2))
			v <- vec[1] / ((vec[1]-a2) * a2)
		}
		else{
			a4 <- rbinom(perm, vec[2], 0.5)
			num <- a2 + a4
			denom <- vec[1] + vec[2]
			beta <- log(num/(denom-num))
			v <- denom/((denom-num)*num)
		}
	}
	else{
		a2 <- rbinom(perm, vec[1], 0.5)
		a4 <- rbinom(perm, vec[2], 0.5)
		matMulti <- rmultinom(perm, vec[3], c(0.25, 0.5, 0.25))
		num <- a2 + a4 + matMulti[2,] + 2*matMulti[3,]
		beta <- log(num/(vec[4]-num))
		v <- vec[4] / ((vec[4] - num) * num)
	}
	return(list(beta=beta, v=v))
}
  
gxeBetaVdom <- function(vec, perm){
	beta <- v <- NA
	if(vec[3]==0 & vec[1]>0){
		d2 <- rbinom(perm, vec[1], 0.5)
		d1 <- vec[1] - d2
		h <- (1/3*d1 - d2) / (2*d1)
		or <- sqrt(d2 / (3*d1) + h*h) - h
		beta <- log(or)
		v <- vec[1]*or/(or+1)^2
	}
	if(vec[3] > 0){
		d4 <- rbinom(perm, vec[3], 0.75)
		d3 <- vec[3] - d4
		d2 <- rbinom(perm, vec[1], 0.5)
		d1 <- vec[1] - d2
		h <- (1/3 * d1 - d2 + d3 - 1/3 * d4) / (2 * (d1 + d3))
		or <- sqrt((d2 + d4) / (3 * (d1 + d3)) + h * h) - h
		beta <- log(or)
		v <- vec[1] * or/(or+1)^2 + vec[3] * or / (3 * (or + 1/3)^2)
	}
	return(list(beta=beta, v=1/v))
}

gxeBetaVrec <- function(vec, perm){
	beta <- v <- NA
	if(vec[3]==0 & vec[2]>0){
		r2 <- rbinom(perm, vec[2], 0.5)
		r1 <- vec[2] - r2
		h <- (3*r1 - r2) / (2*r1)
		or <- sqrt(3*r2 / r1 + h*h) - h
		beta <- log(or)
		v <- vec[2]*or/(or + 1)^2
	}
	if(vec[3]>0){
		r4 <- rbinom(perm, vec[3], 0.25)
		r3 <- vec[3] - r4
		r2 <- rbinom(perm, vec[2], 0.5)
		r1 <- vec[2] - r2
		h <- (3*r1 - r2 + r3 - 3*r4) / (2 * (r1 + r3))
		or <- sqrt(3 * (r2+r4) / (r1+r3) + h*h) - h
		beta <- log(or)
		v <- vec[2] * or/(or + 1)^2 + 3*vec[3] * or/(or + 3)^2
	}
	return(list(beta=beta, v=1/v))
}
	

			