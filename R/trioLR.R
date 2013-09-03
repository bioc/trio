trioLR <- function(x, ...) UseMethod("trioLR")

trioLR.formula <- function(formula, data, recdom=TRUE, ...){
	(require(logicFS, quietly=TRUE) && packageVersion("logicFS") >= "1.28.1") || 
		stop("Package logicFS >= 1.28.1 is required.")
	xy <- getXy(formula, data, recdom=recdom)
	trioLR(xy$x, xy$y, ...)
}

trioLR.trioPrepare <- function(x, ...){
	trioLR(x$bin[,-1], x$bin[,1], ...)
}

trioLR.default <- function(x, y, search=c("sa", "greedy", "mcmc"), nleaves=5, penalty=0, weights=NULL,
		control=lrControl(), rand=NA, ...){
	(require(logicFS, quietly=TRUE) && packageVersion("logicFS") >= "1.28.1") || 
		stop("Package logicFS >= 1.28.1 is required.")
	call <- match.call()
  	if(nrow(x)%%4 != 0)
    		stop("x does not seem to contain trio data with three pseudo-controls for each affected children\n",
      			"since the number of rows is not a multiple of 4.", call.=FALSE)
  	n.trios <- nrow(x)/4
	if(missing(y))
		y <- rep.int(c(3,0,0,0), n.trios)
	if(is.null(weights))
		weights <- rep.int(1, n.trios)
	choice <- match.arg(search)
	select <- checkTLRinput(x, y, choice, n.trios, nleaves=nleaves, penalty=penalty, weights=weights,
		control=control, rand=rand)
	out <- triologreg(x, y, weights, select, nleaves=nleaves, penalty=penalty, control=control, rand=rand)
	out$call <- call
	class(out) <- "trioLR"
	out
}

print.trioLR <- function(x, asDNF=FALSE, posBeta=FALSE, digits=3, ...){
	if(asDNF)
		require(mcbiopi, quietly=TRUE) || stop("Package mcbiopi is required.")
	cat("         Trio Logic Regression\n\n", sep="")
	if(x$choice==7){
		cat("Search Algorithm: MCMC\n\n", sep="")
		if(x$mc.control$output>0)
			cat("Visited models are stored in triolrlisting.tmp.\n")
		else
			cat("Visited models are not stored.\n")
	}
	if(x$choice %in% 1:2){
		cat("Search Algorithm: Simulated Annealing\n\n", sep="")
		if(x$choice==1)
			cat("A single model has been fitted:\n", sep="")
		if(x$choice==2)
			cat(x$nmodels, " models have been fitted:\n", sep="")
	}
	if(x$choice == 6)
		cat("Search Algorithm: Greedy\n\n", x$nmodels, " models have been fitted:\n", sep="")
	cn <- if(is.null(x$binnames)) paste("S", 1:ncol(x$binary), sep="") else x$binnames
	if(x$choice==1){
		printTrioTree(x$model$tree[[1]], cn, asDNF=asDNF, posBeta=posBeta, digits=digits)
		cat("Score: ", signif(x$model$score, digits), "\n")
	}
	if(x$choice %in% c(2,6)){
		score <- signif(x$allscores[,1], digits)
		nleaves <- x$allscores[,2]
		for(i in 1:x$nmodels){
			cat("\n", "Model with ", ifelse(nleaves[i]>0, "a maximum of ", ""), nleaves[i], 
				ifelse(nleaves[i]==1, " leaf", " leaves"), ":\n", sep="")
			printTrioTree(x$alltrees[[i]]$trees[[1]], cn, asDNF=asDNF, posBeta=posBeta, digits=digits)
			cat("Score: ", score[i], "\n", sep="")
		}
	}
}

printTrioTree <- function(ltree, cn, asDNF=FALSE, posBeta=FALSE, digits=3){
	mat <- as.matrix(ltree$trees)
	coef <- signif(ltree$coef, digits)
	n.row <- nrow(mat)
	vec <- character(nrow(mat))
	ids3 <- which(mat[,2] == 3)
	if(length(ids3)==0)
		cat("Model contains no variables.\n", sep="")
	if(length(ids3)==1)
		cat(coef, " * ", ifelse(mat[1,4]==1, "!", ""), cn[mat[1,3]], "\n", sep="")
	if(length(ids3)>1){
		if(asDNF){
			if(posBeta)
				coef <- abs(coef)
			else
				ltree$coef <- abs(ltree$coef)
			tmppi <- getPImps(ltree, 1)
			tmpdnf <- paste("(", getSNPdummyNames(tmppi, cn), ")", sep="", collapse=" | ")
			vec[1] <- if(length(tmppi)==1) tmpdnf else paste("(", tmpdnf, ")", sep="")
		}
		else{			
			for(i in ids3)
				vec[i] <- paste(ifelse(mat[i,4]==1, "!", ""), cn[mat[i,3]], sep="")
			nhalf <- floor(n.row/2)
			for(i in nhalf:1){
				if(mat[i,2] == 1)
					vec[i] <- paste("(", vec[2*i], " & ", vec[2*i+1], ")", sep="")
				if(mat[i,2] == 2)
					vec[i] <- paste("(", vec[2*i], " | ", vec[2*i+1], ")", sep="")
			}
		}
		cat(coef, " * ", vec[1], "\n", sep="")
	}
}


plot.trioLR <- function(x, whichTree=NA, freqType=1, useNames=FALSE, addStats=TRUE, digits=3, main=NULL, cexOper=1.5, 
		cexLeaf=1.5, sizeLeaf=7, cexPar=1.3, ...){
	if(useNames){
		if(is.null(x$binnames))
			stop("If useNames = TRUE, the column names of the genotype matrix must have been specified.")
		cn <- x$binnames
	}
	else
		cn <- 1:ncol(x$binary)
	if(x$choice==1)
		plotTrioTree(x$model$trees[[1]], cn, x$model$score, addStats=addStats, digits=digits, main=main, cexOper=cexOper,
			cexLeaf=cexLeaf, sizeLeaf=sizeLeaf, cexPar=cexPar)
	if(x$choice %in% c(2,6)){
		if(is.na(whichTree))
			stop("If multiple models have been fitted, whichTree must be specified.")
		if(!whichTree %in% 1:x$nmodels)
			stop("whichTree must be an integer between 1 and ", x$nmodels, ".")
		plotTrioTree(x$alltrees[[whichTree]]$trees[[1]], cn, x$allscores[whichTree, 1], addStats=addStats,
			digits=digits, main=main, cexOper=cexOper, cexLeaf=cexLeaf, sizeLeaf=sizeLeaf, cexPar=cexPar)
	}
	if(x$choice == 7){
		if(!freqType %in% 1:3)
			stop("For plotting the results of a MC Trio Logic Regression, freqType must be set\n",
				"either to 1 (for plotting the percentage of visited models containing the individual variables),\n",
				"or to 2 (for plotting the percentage of models containing the different pairs of variables),\n",
				"or to 3 (for the observed-to-expected ratio for the pairs being jointly in the models).")
		if(freqType==1){
			if(length(x$single) == 0)
				stop("No information available in x on the number of models containing the binary variables.")
			if(is.null(main))
				main <- "Marginal Frequency of Being in the Models"
			las <- ifelse(useNames, 3, 0)
			xlab <- ifelse(useNames, "", "Variable")
			barplot(as.vector(x$single)/sum(x$size[,2]) * 100, names=cn, main=main, xlab=xlab,
				ylab="Percentage", las=las)
		}
		else{
			if(length(x$double) == 0)
				stop("No information available in x on the number of models containing the pairs of variables.")
			if(is.null(main))
				main <- paste(ifelse(freqType==2, "Frequency", "Observed/Expected Ratio"), "of Being Jointly in the Models")
			tmp2 <- x$double
			tmp3 <- 1:length(tmp2[,1])
			tmp4 <- outer(tmp3, tmp3, "-")
			tmp2[tmp4 <= 0] <- NA
			if(freqType==2)
				image(tmp3, tmp3, tmp2, main=main, xlab="", ylab="",xaxt="n", yaxt="n")
			else{
				tmp <- x$single
				tmp5 <- outer(tmp/sum(tmp), tmp)
				image(tmp3, tmp3, tmp2/tmp5, main=main, xlab="", ylab="", xaxt="n", yaxt="n")
			}
			axis(1, at=1:length(tmp3), labels=cn, las=ifelse(useNames, 2, 0))
			axis(2, at=1:length(tmp3), labels=cn, las=ifelse(useNames, 1, 0))
		}
	}
}

plotTrioTree <- function(ltree, cn, score, addStats=TRUE, digits=3, main=NULL, cexOper=1.5, cexLeaf=1.5, 
		sizeLeaf=7, cexPar=1.3){
	dat <- ltree$trees
	coef <- signif(ltree$coef, digits)
	score <- signif(score, digits)
	if(is.null(main))
		main <- "Logic Tree"
	if(dat[1,5] == 0){
		cat("The trio logic regression model contains no variable.\n")
		plot(0:1, 0:1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", main=main)
		text(0.5, 0.5, "Empty Tree", cex=2)
		if(addStats)
			title(sub=paste("Score:", score, "\n"), cex.sub=cexPar)
	}
	else{
		level <- ceiling(log2(1 + max(which(dat[,5]==1))))
		max.p <- 2^level - 1
		dat <- dat[1:max.p,]
		conc <- dat[,2]
		cn <- c("", cn)
		knot <- cn[dat[,3]+1]
		neg <- dat[,4]
		num <- dat[,1]
		pick <- dat[,5]
		x <- y <- rep(0, max.p)
		left.b <- 2^(0:(level-1))
		right.b <- 2^(1:level)-1
		l.max <- max(num[left.b][pick[left.b]==1])
		r.max <- max(num[right.b][pick[right.b]==1])
		l.shift <- 1/(l.max+1)
		r.shift <- 1/(r.max+1)
		w.shift <- min(l.shift, r.shift)
		plot(c(0.85 * l.shift, 0.95*(1-r.shift)), c(-1, 1.05*(-level)), type="n", xlab="", ylab="", xaxt="n", yaxt="n", col=0, main=main)
		if(addStats)
			title(sub = paste("Parameter:", coef, "      Score:", score, "\n"), cex.sub=cexPar)
		for(k in level:1){
			y.val <- (-k)
			min.val <- 2^(k-1)
			max.val <- 2^k - 1
			diff.x <- 1/(min.val + 1)
			x.val <- seq(diff.x, 1-diff.x, diff.x)
			if(level==1)
				x.val <- x.val * (0.95 * (1-r.shift) - 0.85 * l.shift) + 0.85 * l.shift
			for(j in 1:length(x.val)){
				ids <- 2^(k-1) + j - 1
				if(pick[ids] == 1){
					x[ids] <- x.val[j]
					y[ids] <- y.val
					if(conc[ids] < 3)
						text(x.val[j], y.val, ifelse(conc[ids]==1, "AND", "OR"), cex=cexOper, font=2)
					else{
						points(x.val[j], y.val, pch=ifelse(neg[ids]==0, 0, 15), cex=sizeLeaf)
						text(x.val[j], y.val, knot[ids], cex=cexLeaf, col=ifelse(neg[ids]==0, 1, "white"))
					}
				}
			}
		}
		for(k in 1:max.p){		
			if((ceiling(log2(k)) < level) && (2*k<max.p) && (pick[k] + pick[2*k] == 2)){
				lines(x[c(k, 2*k)], c(y[k] - 0.2, y[2*k] + 0.2))
				lines(x[c(k, 2*k+1)], c(y[k]-0.2, y[2*k+1] + 0.2))
			}
		} 
		invisible()
	}
}				

