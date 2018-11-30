util.randomCaseAssign <- function (segLen, segCt){
  	inSegPos <- sample(seq(1, segLen, by=1), segCt, replace=TRUE)
  	addLen <- seq(0, (segCt-1)*segLen, by = segLen) 
  	casePos <-  as.numeric(inSegPos) + addLen
  	randomCol <- rep(0, length=segCt*segLen)
  	randomCol[casePos] <- 1 
  	randomCol
}

 
util.shuffleIdx <- function(segLen, segCt){
	sampleFun <- function(item, numIn){
   		ranIdx <- sample(1:numIn, size=1)
   		allIdx <- 1:numIn
   		leftIdx <- allIdx[allIdx != ranIdx]
   		shuIdx <- c(ranIdx, leftIdx)
  	}
	newIdx <- unlist(lapply(1:segCt, FUN=sampleFun, numIn = segLen))
  	if(segCt>1)
    		baseSeq <- rep(seq(0, segLen*segCt-1, by=segLen), each=segLen)
  	
	else
    		baseSeq <- 0
  	respMa <- baseSeq + newIdx
  	respMa
}


trio.permTest <- function(object, conditional=FALSE, n.perm=10, nleaves=NULL, control=NULL, rand=NA){
	# require(mcbiopi, quietly=TRUE) || stop("Package mcbiopi is required.")
    	if(n.perm <= 0)
		stop("n.perm must be at least 1.")
	if(!is(object, "trioLR"))
		stop("object must be an object of class trioLR.")
	if(!object$choice == 1)
		stop("In the original analysis leading to object, a single model must have been fitted using simulated annealing.")
	resp <- object$response
	nr <- length(resp)
	bina <- object$binary
	weights <- object$weight[resp==3]
	if(length(weights) != length(resp)/4)
		stop("The number of weights is incorrect.")
	if(is.null(nleaves))
		nleaves <- object$nleaves
	if(length(nleaves) != 1)
		stop("nleaves must be a single integer.")
	if(nleaves < 1)
		stop("nleaves must be at least 1.")
	if(is.null(control))
		control <- object$control
	else
		checkControlPars(control, nr/4)
    	setId <- rep(1:(nr/4), each=4)
    	dIdx <- seq(1, nr, by=4)
	if(!is.na(rand))
		set.seed(rand)
	if(!conditional){
      		logreg.pred <- util.randomCaseAssign(segLen=4, segCt=nr/4)
      		logreg.score <- NA
	}
	else{
		logreg.score <- object$model$score
  		logreg.pred <- mcbiopi::evalTree(object$model$trees[[1]]$trees, bina)
  	}
  	scoreVec <- rep(NA, times=n.perm+1)
  	scoreVec[1] <- logreg.score
	
	tmpfoo1 <- function(jj, predStatus){
        	setPred <- predStatus[jj:(jj+3)]
              	if(predStatus[jj]==1){
                	oneIdx <- which(setPred==1)
                	zeroIdx <- which(setPred==0)
              	}
		else{
                	oneIdx <- which(setPred==0)
                	zeroIdx <- which(setPred==1)
              	}
              	lens <- c(length(oneIdx), length(zeroIdx))
              	lens <- lens[which(lens!=0)]
    		re <- list(a=oneIdx, b=zeroIdx)
 	}
  	tmpfoo2 <- function(kk, pl){
          	ll <- pl[[(kk-1)/4+1]]
          	oneIdx <- ll$a
          	idx.new1<-oneIdx
            	if(length(oneIdx)!=0){
            		idx.shuffle1 <- util.shuffleIdx(segLen=length(oneIdx), segCt=1)
            		idx.new1 <- oneIdx[idx.shuffle1]
          	}
            	zeroIdx <- ll$b
          	c(kk + (idx.new1-1), kk + (zeroIdx-1) )
	}
	idxSetIdxShuffle <- NULL
	cat("\n", sep="")
  	for(i in 1:n.perm){
		cat("Permutation ", i, ".\n", sep="")
       		if(is.null(idxSetIdxShuffle)){
          		idxSetIdxShuffle <- lapply(dIdx, FUN=tmpfoo1, predStatus=logreg.pred)
        		#oneIdxes <- idxSetIdxShuffle[seq(1, length(dIdx)*2, by=2)]
          		#zeroIdxes <- idxSetIdxShuffle[seq(2, length(dIdx)*2, by=2)]
      		}
       		newIdx <- lapply(dIdx, FUN=tmpfoo2, pl=idxSetIdxShuffle)
   		binaNew <- bina[unlist(newIdx),]
   		logregRe <- NULL
       		logregRe <- triologreg(binaNew, resp, weights, 1, nleaves=nleaves, penalty=object$penalty, control=control)
  		scoreVec[1+i] <- logregRe$model$score
  	} 
	structure(list(origScore=scoreVec[1], permScore=scoreVec[-1]))
}
