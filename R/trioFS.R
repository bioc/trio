trioFS <- function(x, ...) UseMethod("trioFS")

trioFS.default <- function(x, y, B=20, nleaves=5, replace=TRUE, sub.frac=0.632, control=lrControl(), fast=FALSE, 
		addMatImp=TRUE, addModels=TRUE, verbose=FALSE, rand=NA, ...){
  	if(nrow(x)%%4 != 0)
    		stop("x does not seem to contain trio data with three pseudo-controls for each affected children\n",
      			"since the number of rows is not a multiple of 4.", call.=FALSE)
	if(missing(y))
		y <- rep.int(c(3,0,0,0), nrow(x)/4)
	if(is.null(colnames(x))){
		colnames(x) <- paste("S", 1:ncol(x), sep="")
		warning("Since x has no column names, generic ones are added.", call.=FALSE)
	}
	bagg.out <- trioBagging(x, y, B=B, nleaves=nleaves, replace=replace, sub.frac=sub.frac, control=control,
		fast=fast, verbose=verbose, rand=rand)
	vim.out <- vim.trioFS(bagg.out, x, y, addInfo=TRUE, addMatImp=addMatImp)
	if(addModels){
		vim.out$logreg.model <- bagg.out$logreg.model
		vim.out$inbagg <- bagg.out$inbagg
	}
	vim.out
}

trioFS.formula <- function(formula, data, recdom=TRUE, ...){
	(require(logicFS, quietly=TRUE) && packageVersion("logicFS") >= "1.28.1") || 
		stop("Package logicFS >= 1.28.1 is required.")
	xy <- getXy(formula, data, recdom=recdom)
	trioFS(xy$x, xy$y, ...)
}

trioFS.trioPrepare <- function(x, ...){
	trioFS(x$bin[,-1], x$bin[,1], ...)
}	

trioBagging <- function(x, y, B=20, nleaves=5, replace=TRUE, sub.frac=0.632, control=NULL, fast=FALSE, 
		verbose=FALSE, rand=NA){
	(require(logicFS, quietly=TRUE) && packageVersion("logicFS") >= "1.28.1") || 
		stop("Package logicFS >= 1.28.1 is required.")
	if(B < 1)
		stop("B must be a positive integer.")
  	n.trios <- nrow(x)/4
	if(missing(y))
		y <- rep.int(c(3,0,0,0), n.trios)
	if(length(nleaves)!=1 || nleaves<1 || nleaves!=round(nleaves))
		stop("nleaves must be a	single positive integer.")
	choice <- ifelse(fast, "greedy", "sa")
	weights <- rep.int(1, n.trios)
	select <- checkTLRinput(x, y, choice, n.trios, nleaves=nleaves, weights=weights, control=control, rand=rand)
	if(!replace){
		if(sub.frac < 0.1 | sub.frac > 0.9)
			stop("sub.frac must be between 0.1 and 0.9.")
		n.sub <- ceiling(sub.frac * n.trios)
		sampling <- paste("Subsampling (", 100 * sub.frac, "% of the trios)",
			sep="")
	}
	else
		sampling <- "Bagging"
	if(verbose)
		cat("Starting Building Trio Logic Regression Models at ", date(), ".", "\n",
			"Depending on B and anneal.control this might take up to several hours.\n\n", sep="")
	list.trees <- list.bagg <- vector("list", B)
	if(!is.na(rand))
		set.seed(rand)
	for(i in 1:B){
		if(verbose)
			cat("Model ", i, ". ", sep="")
		ids.trio <- if(replace) sample(n.trios, n.trios, replace=TRUE)
			else sample(n.trios, n.sub)
		bagg <- rep(4*(ids.trio-1), e=4) + rep(1:4, length(ids.trio))
		tmp.out <- triologreg(x[bagg,], y[bagg], weights[ids.trio], select, nleaves=nleaves, 
			control=control, notBuiltStop=FALSE)
		if(!fast)
			list.trees[[i]] <- tmp.out$model
		else{
			ids.min <- which.min(tmp.out$allscores[,1])
			list.trees[[i]] <- tmp.out$alltrees[[ids.min]]
		}
		list.bagg[[i]] <- bagg
		if(verbose)
			cat("Done.\n", sep="")
	}
	if(verbose)
		cat("\n", "Finished Building Models at ", date(), ".", "\n\n", sep="")
	idsNull <- sapply(list.trees, is.null)
	if(any(idsNull)){
		n.null <- sum(idsNull)
		if(n.null>B/2 || B-n.null<5)
			stop("In less than ", ifelse(B-n.null<5, "five", "half of the"), " iterations, the trio logic regression",
				"model\n", "has been fitted properly. Please rerun the analysis with trioFS.\n",
				"If this error persists, please contact the maintainer of the trio package.", call.=FALSE)
		warning("In ", n.null, " of the iterations, the fitting of the trio logic regression model failed\n",
			"for an unknown reason so that only ", B-n.null, " models are considered in the detection\n",
			"of the SNP interactions in these models and the quantification of their importance.", call.=FALSE)
		list.trees <- list.trees[!idsNull]
		list.bagg <- list.bagg[!idsNull]
	}
	out <- list(logreg.model=list.trees, inbagg=list.bagg, sampling=sampling, nleaves=nleaves, type=0)
	class(out) <- "trioBagg"
	out
}

vim.trioFS <- function(object, x, y, neighbor=FALSE, addInfo=FALSE, addMatImp=TRUE){
	require(survival, quietly=TRUE) || stop("Package survival is required.")
	if(!is(object, "trioBagg"))
		stop("object must be an object of class trioBagg.")
	list.primes <- trioPImp(object)
	B <- length(list.primes)
	vec.primes <- unlist(list.primes)
	prop <- table(vec.primes)/B
	vec.primes <- unique(vec.primes)
	dat <- as.data.frame(x)
	colnames(dat) <- paste("X", 1:ncol(dat), sep="")
	mat.eval <- trioPImpEval(dat, vec.primes)
	if(ncol(mat.eval) < length(vec.primes)){
		ids <- which(!vec.primes %in% colnames(mat.eval))
		mono <- vec.primes[ids]
		vec.primes <- vec.primes[-ids]
		rmMonoPI <- function(lpi, mono) lpi[!lpi %in% mono]
		list.primes <- lapply(list.primes, rmMonoPI, mono=mono)
	}
	cl <- as.numeric(y > 0)
	inbagg <- object$inbagg
	mat.imp <- matrix(0, length(vec.primes), B)
	rownames(mat.imp) <- colnames(mat.eval)
	for(i in 1:B)
		mat.imp[,i] <- vimTrio(list.primes[[i]], mat.eval, inbagg[[i]], cl)
	vim <- rowMeans(mat.imp, na.rm=TRUE)
	prop <- prop[vec.primes]
	if(neighbor)
		warning("Neighborhood measures are currently not available. Therefore,\n",
			"neighbor = TRUE is ignored.", call.=FALSE)
	primes <- getSNPdummyNames(vec.primes, colnames(x))
	param <- if(addInfo) list(B=B, nleaves=object$nleaves, sampling=object$sampling) else NULL
	if(!addMatImp)
		mat.imp <- NULL
	measure <- "Trio"
	vim.out <- list(vim=vim, prop=prop, primes=primes, param=param, mat.imp=mat.imp)
	class(vim.out) <- "trioFS"
	vim.out
}


trioPImp <- function(log.out){
	require(mcbiopi, quietly=TRUE) || stop("Package mcbiopi is required.")
	lmodel <- log.out$logreg.model
	list.primes <- vector("list", length(lmodel))
	for(i in 1:length(lmodel))
		list.primes[[i]] <- getPImps(lmodel[[i]]$trees[[1]], 0)
	whichNull <- which(sapply(list.primes, is.null))
	if(length(whichNull)>0){
		list.primes <- list.primes[-whichNull]
		warning("Since ", length(whichNull), " of the models contain no variables, ",
			"they are removed.", call.=FALSE)
	}
	list.primes	
}

trioPImpEval <- function(dat, vec.primes){
	mat.eval <- with(dat, sapply(vec.primes, function(x) eval(parse(text=x))))
	cs <- colSums(mat.eval)
	if(any(cs %in% c(0, nrow(mat.eval)))){
		ids <- which(cs %in% c(0, nrow(mat.eval)))
		warning("For ", length(ids), "of the ", ncol(mat.eval), " interactions, ",
			"all cases and pseudo-controls show the same value.\n",
			"These interactions are therefore removed.", call.=FALSE)
		mat.eval <- mat.eval[,-ids]
	}
	mat.eval
}

predLL <- function(beta, x){
	mat.x <- matrix(x, ncol=4, byrow=TRUE)
	expbeta <- exp(mat.x * beta)
	vec.ll <- expbeta[,1] / rowSums(expbeta)
	sum(log(vec.ll))
}

vimTrio <- function(primes, mat.eval, inbagg, cl){
	n.primes <- length(primes)
	id.primes <- colnames(mat.eval) %in% primes
	oob <- which(!(1:nrow(mat.eval)) %in% inbagg)
	vec.improve <- numeric(ncol(mat.eval))
	names(vec.improve) <- colnames(mat.eval)
	if(n.primes == 1){
		vec.improve[id.primes] <- vimTrio.oneprime(mat.eval[,id.primes], cl,
			oob, inbagg)
		return(vec.improve)
	}
	mat.in <- 1 - diag(n.primes)
	mat.pred <- mat.eval[,id.primes] %*% mat.in
	mat.pred[mat.pred > 1] <- 1
	vec.ll <- numeric(n.primes)
	st <- rep(1:(length(inbagg) / 4), e=4)
	for(i in 1:n.primes){
		x <- mat.pred[,i]
		coef <- clogit(cl[inbagg] ~ x[inbagg] + strata(st))$coefficients
		vec.ll[i] <- predLL(coef, x[oob])
	}
	x <- (rowSums(mat.eval[,id.primes]) > 0) * 1
	coef <- clogit(cl[inbagg] ~ x[inbagg] + strata(st))$coefficients
	ll.full <- predLL(coef, x[oob])
	vec.improve[id.primes] <- vec.ll - ll.full
	-2 * vec.improve	
}

vimTrio.oneprime <- function(oneprime, cl, oob, inbagg){
	n.oobfam <- length(oob) / 4
	ll.null <- log(0.25) * n.oobfam
	st <- rep(1:(length(inbagg) / 4), e=4)
	coef <- clogit(cl[inbagg] ~ oneprime[inbagg] + strata(st))$coefficients
	ll.prime <- predLL(coef, oneprime[oob])
	-2 * (ll.null - ll.prime)
}

getSNPdummyNames<-function(vec.primes,col.names){
	n.col<-length(col.names)
	vec.primes<-gsub("X","XtendedName",vec.primes)
	coded<-paste("XtendedName",1:n.col,sep="")
	for(i in n.col:1)
		vec.primes<-gsub(coded[i],col.names[i],vec.primes)
	vec.primes
}

print.trioFS<-function(x,topX=5,show.prop=TRUE,coded=FALSE,digits=2,...){
	param<-x$param
	cat("Identification of Interactions Using Trio Logic Regression\n\n")
	cat("Number of Iterations: ", param$B,"\n")
	cat("Sampling Method:      ", param$sampling,"\n")
	cat("Number of Trees:      ", 1,"\n")
	cat("Max. Number of Leaves:", param$nleaves,"\n\n")
	vim<-sort(x$vim,decreasing=TRUE)
	topX<-min(topX,length(vim))
	names.vim<-if(coded) names(vim) else x$primes[order(x$vim,decreasing=TRUE,na.last=NA)]
	out<-data.frame(Importance = vim,
		Proportion = if(!is.null(x$prop)) x$prop[order(x$vim,decreasing=TRUE,na.last=NA)] else NA,
		Expression = names.vim)
	if(!is.null(x$name))
		names(out)[3]<-x$name
	rownames(out)<-1:nrow(out)
	if(!show.prop | is.null(x$prop))
		out<-out[,-2]
	out<-format(out[vim>=vim[topX],],digits=digits,nsmall=2)
	cat("The",nrow(out),"Most Important Interactions:\n\n")
	print(out)
}

plot.trioFS<-function(x,topX=15,show.prop=FALSE,coded=TRUE,cex=.9,pch=16,col=1,force.topX=FALSE,
		include0=TRUE,add.v0=TRUE,v0.col="grey50",main=NULL,...){
	if(!show.prop){
		vim<-x$vim
		if(is.null(main))
			main<-paste(x$measure,"Measure")
		xlab<-"Importance"
	}
	else{
		vim<-x$prop
		if(is.null(vim))
			stop("Data for the Ad Hoc Measure are not available.")
		if(is.null(main))
			main<-"Ad Hoc Measure"
		xlab<-"Proportion"
	}
	if(!coded)
		names(vim)<-x$primes
	vim<-sort(vim,decreasing=TRUE)
	topX<-min(topX,length(vim))
	vim<-if(force.topX) vim[1:topX] else vim[vim>=vim[topX]]
	rangex<-if(!show.prop) range(if(include0) 0,vim) else c(0,1)
	dotchart(rev(vim),main=main,xlab=xlab,pch=pch,color=col,cex=cex,xlim=rangex)
	if(!show.prop & add.v0 & rangex[1]<=0)
		abline(v=0,lty="dotted",col=v0.col)
}
	