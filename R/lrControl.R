lrControl <- function(start=0, end=0, iter=0, earlyout=0, update=0, treesize=8, opers=1, minmass=0,
		nburn=1000, hyperpars=0, output=4){
	list(start=start, end=end, iter=iter, earlyout=earlyout, update=update, treesize=treesize,
		opers=opers, minmass=minmass, nburn=nburn, hyperpars=hyperpars, output=output)
}

checkControlPars <- function(parlist, n.trios){
	if(!is.list(parlist))
		stop("control must be a list.")
	if(is.null(names(parlist)))
		stop("All elements in control must have a name.")
	if(any(sapply(parlist, length) != 1))
		stop("All elements in parlist must have length 1.")
	namesPar <- names(lrControl())
	namesList <- names(parlist)
	if(any(!namesList %in% namesPar))
		stop("Some of the parameters specified in control are not control parameters for a Trio Logic Regression.")
	if(any(!namesPar %in% namesList))
		stop("Some control parameters are missing in control.")
	if(parlist$start < parlist$end)
		stop("start must be larger than end.")
	if(parlist$iter!=0 && parlist$iter<1000)
		stop("iter must be at least 1000 (and should actually be much, much larger than 1000).")
	if(parlist$earlyout < 0)
		stop("earlyout must be a non-negative integer.")
	if(parlist$treesize<1)
		stop("treesize must be a positive integer.")
	if(!parlist$opers %in% (1:3))
		stop("opers must be either 1, 2, or 3.")
	if(parlist$minmass < 0)
		stop("minmass must be at least 0.")
	if(parlist$minmass > n.trios)
		stop("minmass must be at most equal to the number of trios, i.e. 25% of the number of (pseudo-)observations.")
	if(parlist$nburn < 0)
		stop("nburn must be at least 0.")
	if(abs(parlist$output) > 4)
		stop("output must be between -4 and 4.")
} 
	 

checkTLRinput <- function(x, y, choice, n.trios, nleaves=5, penalty=0, weights=NULL, control=NULL, rand=NA){
	if(!is.matrix(x))
    		stop("x must be a matrix.")
	if(ncol(x) < 1)
		stop("x must contain at least one variable.")
	if(ncol(x) > 1000)
		stop("x is not allowed to contain more than 1000 variables.")
  	if(any(is.na(x)))
    		stop("No missing values allowed in x.")
  	if(any(!x %in% c(0,1)))
    		stop("All values in x must be either 0 or 1.")
	if(n.trios>5000)
		stop("x is not allowed to contain data from more than 5000 trios.")
  	if(!is.numeric(y))
      		stop("y must be numeric.")
    	if(any(is.na(y)))
      		stop("No missing values are allowed in y.")
    	if(nrow(x) != length(y))
      		stop("The number of rows of x must be equal to the length of y.")
    	if(any(y != rep.int(c(3,0,0,0), n.trios)))
      		stop("y is not specified correctly.")
  	if(any(nleaves<1) || any(nleaves>32))
    		stop("nleaves must be an integer (or vector of integers) between 1 and 32.")
  	if(choice == "sa"){
		if(length(nleaves)>2)
			stop("nleaves is not allowed to consist of more than two values.")
    		select <- ifelse(length(nleaves)==1, 1, 2)
	}
  	else{
  		if(length(nleaves)>1)
    			stop("nleaves must be a single integer if either a greedy search or MC logic regression is performed.")
    		select <- ifelse(choice=="greedy", 6, 7)
	}
	if(select==2 && nleaves[2]<=nleaves[1])
		stop("When fitting several models, the second entry in nleaves must be larger than the first.")
  	if(penalty < 0)
    		stop("penalty must be a non-negative value.")
  	if(length(weights)!=n.trios)
    		stop("The number of weigths is not equal to the number of trios.")
  	if(any(weights<=0))
		stop("All weights must be larger than zero.")
	checkControlPars(control, n.trios)
	if(!is.na(rand)){
		if(rand < 1)
			stop("If specified, rand must be larger than zero.")
		if(rand != round(rand))
			stop("rand must be an integer.")
	}
	select
}
