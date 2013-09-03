splitBlocks <- function(blocks){
	if(!is(blocks, "LDblocks"))
		stop("blocks must be the output of findLDblocks.")
	ids.block <- blocks$blocks
	tab.block <- tabulate(ids.block)
	if(all(tab.block <= 7))
		return(ids.block)
	tmp.ids <- which(tab.block > 7)
	for(j in 1:length(tmp.ids))
		ids.block[ids.block==tmp.ids[j]] <- splitBlock(tmp.ids[j],
			tab.block[tmp.ids[j]])
	ids.block
}

splitBlock <- function(b, n){
	nb <- (n-1)%/%7 + 1
	if(nb==1)
		return(rep(b, n))
	minn <- n%/%nb
	rb <- n%%nb
	num <- sample(rep(c(minn,minn+1), c(nb-rb,rb)))
	rep(b+(1:nb)/10, num)
}
