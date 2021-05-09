qingEnv <- new.env()

augBkFrame <-
function(hapBkMap, key, probLeftover = .01){
	keyIndex = which(hapBkMap$keys==key)
	bk = hapBkMap$bks[[keyIndex]]
	
	## extract the original exp and prob, then restandard prob
	estBkExp = as.character(bk[,hapBkMap$expCol])
	estBkProbRe = bk[,hapBkMap$probCol]*(1-probLeftover)
	
	## create the exhaust haplotypes 
	bkSnpLen = bk[1,hapBkMap$hapLenCol]
	exhaustExp = exhaustHapExp(lociCt=bkSnpLen, snpCoding=c(1,2))$hapStr
	exhaustExpCt = length(exhaustExp)
	
	## augmate the additional expression
	rematch = match(estBkExp, exhaustExp)
	
	resiProb = probLeftover / (exhaustExpCt-length(estBkExp))
	resiProb = rep(resiProb, times=exhaustExpCt)
	resiProb[rematch]=estBkProbRe
	
	base = bk[1,]
	newBk = NULL
	for( i in 1:exhaustExpCt ){
		newBk = rbind(newBk, base)
	}
	## replace the expression and probability
	newBk[,hapBkMap$probCol]=resiProb
	newBk[,hapBkMap$expCol]=as.integer(exhaustExp)
	
	## replace the necessary part in hapBkMap
	hapBkMap$bks[[keyIndex]]=newBk
	
	hapBkMap$bkLens[keyIndex]=exhaustExpCt
	
	## other parts, like df, dfStr,  in hapBkMap is not updated because it is not necessary
	
	return(hapBkMap)
	
}

binaTree.1Level.changeSignal <-
function(bina.col, model.signal.cur ){

    bina.col.num = sapply(bina.col, FUN=function(item) {as.numeric(substr(item, 1, nchar(item)-1))})
    bina.col.ct = table(bina.col.num)
    bina.snp.single = paste(dimnames(bina.col.ct)[[1]][bina.col.ct==1], "b", sep="")

    bTree = binaTree.parser(str=model.signal.cur)
    leaves.all = unlist(bTree$elm.vlist)
    leaves.comb = unique(leaves.all[  is.element(unlist(bTree$elm.vlist), bina.snp.single)  ])
    if(length(leaves.comb)>=1){
      for (i in 1:length(leaves.comb)){
        tmp.str = leaves.comb[i]
        ##print(tmp.str)
        model.signal.cur = util.str.replace(str=model.signal.cur,
                           replaced=tmp.str,
                           new=paste(substr(tmp.str, 1, nchar(tmp.str)-1), "a", sep=""),
                           replace.all=TRUE)
        ##print(model.signal.cur)
      }
      
    }
    return(model.signal.cur)
  }

binaTree.addElm <-
function(binaTree,  elm.level, elm.var,  elm.bool){
  # objects for variable only elment for memory reason
  # 1)elm.vlist: a list of vector of variable names, one vector for one element
  # 2)elm.vMa: A matrix for variable only element
  #   matrix col 1:id: list id
  #          col 2:level: level
  #          col 3:bool: boolean term: not(-1) vs 1 ## CHANGED 8JAN09 before:[and(0); or(1); not(-1)]
  
  # two objects for state info
  # curElm: whether current element is a variable only(0) or node(1), or missing(NA)
  # curElmId: the index of current element, or missing(0)
  
  # when the element is a "not" element, the level has to be -1
  ifD = FALSE
 
  if (ifD) print("binaTree.addElm")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print(qp("elm.var=", elm.var))
  if (ifD) print(qp("elm.bool=", elm.bool))
  if(elm.bool==-1) elm.level=elm.level-1
  
  elm.vlist = binaTree$elm.vlist
  elm.vMa = binaTree$elm.vMa
  
  elm.vlist = c(elm.vlist, list(elm.var))
  
  if(binaTree$curElmId!=0){
    
    elm.lastIdx = max(elm.vMa[,1])
    elm.vMa = rbind(elm.vMa, c(elm.lastIdx+1, elm.level, elm.bool))
  }else{
    ## initialize the empty variables
    elm.lastIdx=0
    elm.vMa[1,]= c(elm.lastIdx+1, elm.level, elm.bool)
  }
  binaTree$elm.vlist = elm.vlist
  binaTree$elm.vMa = elm.vMa

  binaTree$curElm = 0;  #variable element
  binaTree$curElmId = elm.lastIdx+1 # id
  if (ifD) print("after addElm, binaTree")
  if (ifD) print(binaTree)
  return(binaTree)
}

binaTree.addNode <-
function(binaTree,  elm.level, elm.bool){
  ifD = FALSE
  if (ifD) print("binaTree.addNode")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print(qp("elm.bool=", elm.bool))
   
  # objects for node element. (contain node itself or variable only element)
  # 1)elm.nlist: a list of matrix, each vector for one element
  #    matrix col 1: idx for node or variable element
  #    matrix col 2: is a node(1) or variable element(0)
  # 2)elm.nMa: A matrix for node
  #    matrix col 1:id: list id
  #           col 2:level: levle
  #           col 3:bool: boolean term: and(0); or(1)
  #           col 4:open: open(1)/close(0)
  # two objects for state info
  # curElm: whether current element is a variable only(0) or node(1), or missing(NA)
  # curElmId: the index of current element
  
  ## called when the node is open or adding element
  elm.nlist = binaTree$elm.nlist
  elm.nMa = binaTree$elm.nMa

  ## if elm.nMa has the node open
  if (length(elm.nlist)>=1){
      fil = elm.nMa[,2]==elm.level & elm.nMa[,4]==1
      if(sum(fil, na.rm=TRUE)==1){
        # append to old node
        open.id = elm.nMa[,1][fil]
        
        nMa = elm.nlist[[fil]]
        nMa = rbind(nMa, c(binaTree$curElmId, binaTree$curElm))
        elm.nlist[fil]=list(nMa)
      }else{
        # new node
        elm.lastIdx = max(elm.nMa[,1])
        open.id = elm.lastIdx + 1
        elm.nMa = rbind(elm.nMa, c(open.id, elm.level, elm.bool, 1))

        nMa = matrix(c(binaTree$curElmId, binaTree$curElm), nrow=1, ncol=2)
        elm.nlist=c(elm.nlist, list(nMa))
      }
  }else{
      ## initialize the empty variables
      open.id = 1
      elm.nMa[1,]=c(open.id, elm.level, elm.bool, 1)
      nMa = matrix( c(binaTree$curElmId, binaTree$curElm), nrow=1, ncol=2)
      elm.nlist = c(elm.nlist, list(nMa))
  }
  
  binaTree$elm.nlist = elm.nlist
  binaTree$elm.nMa = elm.nMa
  
  binaTree$curElm = 1;  #node element
  binaTree$curElmId = open.id # id
  if (ifD) print("after addNode, binaTree:")
  if (ifD) print(binaTree)
  return(binaTree)
}

binaTree.apply <-
function(binaTree, data, colnames, re.final=TRUE){
  ifD = FALSE
  ori.rowCt = nrow(data)
  col.idx = 1:ncol(data)
  
  ## data.v: a matrix of boolean (0/1) for each row in binaTree$vMa, ordered by id
  ## data.n: a matrix of boolean (0/1) for each row in binaTree$nMa, ordered by id
  data.v=NULL
  data.n=NULL

  ## data.f: a vector of boolean (0/1) for the node at level 0
  ## or for the variable element if no node is defined.
  data.f=NULL

  data.f = rep(NA, ori.rowCt)
  
  ## process the variable only element first, regardless the level
  vCt = nrow(binaTree$elm.vMa)
  data.v = matrix(NA, ncol=vCt, nrow=ori.rowCt) 
  for( i in 1:vCt){
    cur.row = binaTree$elm.vMa[i,]
    ## col seq= id, level, bool
    ## find the var list
    cur.elm = binaTree$elm.vlist[[ i ]]
    cur.col = col.idx[is.element(colnames, cur.elm)]
    if (length(cur.col)!=length(cur.elm)) stop()

    cur.bool = cur.row[3]
    if(length(cur.col)>1){
     
      if(cur.bool==0){ ## and
        data.v[,i]=as.numeric(apply(data[,cur.col], 1, all))
      }
      if(cur.bool==1){ ## or
        data.v[,i]=as.numeric(apply(data[,cur.col], 1, any))
      }
    }else{
      if(cur.bool==-1){ ## not
        data.v[,i]=as.numeric(1-data[,cur.col])
      }else{
        data.v[,i]=data[,cur.col]
      }
    }
  } # for( i in 1:vCt){

  ## if have node element, process the node from deeper level to zero
  data.n=NULL
  nCt = length(binaTree$elm.nlist)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n=matrix(NA, ncol=nCt, nrow=ori.rowCt)

    new.order = order(nMa[,2], decreasing = TRUE)
    if(ifD) print(nMa)
    nMa = nMa[new.order, ,drop=FALSE]
    if(ifD) print(nMa)
    #nlist = binaTree$elm.nlist[new.order]
    nlist = binaTree$elm.nlist
    if(ifD) print(nlist)
    for( i in 1:nCt ){
      cur.row = nMa[i, ,drop=FALSE]
      if(ifD) print("cur node row after sorting")
      if(ifD) print(cur.row)
      ## col seq= id, level, bool, open
      
      ## find the var list
      cur.elm = nlist[[ cur.row[1] ]]
      if(ifD) print("cur node element")
      if(ifD) print(cur.elm)
      dd.velm=NULL
      dd.nelm=NULL
      # if variable element, get from data.v
      if(sum(cur.elm[,2]==0)>=1){
        dd.velm = data.v[, cur.elm[,1][cur.elm[,2]==0]]
        
        if (ifD) {
          if(is.null(dim(dd.velm))){
            print(dd.velm[1:5])
          }else{
            print(dd.velm[1:5,])
          }
        }
      }
      # if node element, get from data.n
      if(sum(cur.elm[,2]==1)>=1){
        dd.nelm = data.n[, cur.elm[,1][cur.elm[,2]==1]]
        
        if (ifD) {
          if(is.null(dim(dd.nelm))){
            print(dd.nelm[1:5])
          }else{
            print(dd.nelm[1:5,])
          }
        }
      }

      dd.cur = cbind(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]
     
      if(cur.bool==0){ ## and
        data.n[,cur.row[1]]=as.numeric(apply(dd.cur, 1, all))
      }
      if(cur.bool==1){ ## or
        data.n[,cur.row[1]]=as.numeric(apply(dd.cur, 1, any))
      }

      if( (cur.bool!=1) & (cur.bool!=0)) stop("Wrong boolean value.")
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

  if(re.final){
    if(nCt==0){
      # only one variable only element
      return(data.v[,1])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[,t.id]
      return(re)
    }

  }else{
    return(list(data.v=data.v, data.n=data.n))
  }
}

binaTree.closeNode <-
function(binaTree,  elm.level){
  ## if the current element is an element with not as boolean,
  ## it will have have row in the node map

  # objects for node element. (contain node itself or variable only element)
  # 1)elm.nlist: a list of vector, each vector for one element
  #    matrix col 1: idx for node or variable element
  #    matrix col 2: is a node(1) or variable element(0)
  # 2)elm.nMa: A matrix for node
  #    matrix col 1:id: list id
  #           col 2:level: levle
  #           col 3:bool: boolean term: and(0); or(1)
  #           col 4:open: open(1)/close(0)
  # two objects for state info
  # curElm: whether current element is a variable only(0) or node(1), or missing(NA)
  # curElmId: the index of current element

  ifD = FALSE
  if (ifD) print("binaTree.closeNode")
  if (ifD) print(qp("elm.level=", elm.level))
  if (ifD) print("binaTree::")
  if (ifD) print(binaTree)
  
  elm.nlist = binaTree$elm.nlist
  elm.nMa = binaTree$elm.nMa
  open.id = NA

  ## print(elm.nMa)
  ## see if any node is open (should be)
  if (nrow(elm.nMa)>=1){
      fil = (elm.nMa[,2]==elm.level) & (elm.nMa[,4]==1)
      if (ifD) print(fil)
      if(sum(fil, na.rm=TRUE)==1){
        # close the old node
        open.id = elm.nMa[,1][fil]
        elm.nMa[fil, 4]=0

        nMa = elm.nlist[[which(fil)]]
        if(ifD) print(nMa)
        nMa = rbind(nMa, c(binaTree$curElmId, binaTree$curElm))
        elm.nlist[which(fil)]=list(nMa)
      }else{
        stop("no node match the requirement")
      }
  }else{
      stop("no open node")
  }
  
  binaTree$elm.nlist = elm.nlist
  binaTree$elm.nMa = elm.nMa
  binaTree$curElmId = open.id#
  binaTree$curElm = 1 #node

  if(ifD) print("after close node, binaTree:")
  if(ifD) print(binaTree)
  
  return(binaTree)
}

binaTree.constr <-
function(){

  elm.vlist = list()
  elm.vMa = matrix(NA, ncol=3, nrow=1)
  colnames(elm.vMa) = c("id", "level", "bool")

  elm.nlist = list()
  elm.nMa = matrix(NA, ncol=4, nrow=1)
  colnames(elm.nMa) = c("id", "level", "bool", "open")

  binaTree = list(elm.vlist=elm.vlist, elm.vMa=elm.vMa, elm.nlist=elm.nlist, elm.nMa=elm.nMa,
               curElm=NA, curElmId=0)
  return(binaTree)
}

binaTree.merge <-
function(treeForSearch){
 fN = "binaTree.merge:"
 ifD = FALSE
 
 searchWholeTree = TRUE
 tt = 0
 search = FALSE
 # do an iterative search for ( or ) and , keep search until all the "or" have no parent with "and" operator
 while(searchWholeTree){
     tt = tt+1
     if(tt>100) stop(paste(fN, "iterative procdure exceed 100 times. Sometime is wrong!"))
     
     # if the tree is updated, use the updated one. Otherwise, search the next possible choice
     if(!search){
       if(ifD) {
         print("Updated the tree:")
         #print(binaTree.toStr(treeForSearch))
         print(treeForSearch)
       }
       tree = treeForSearch
       elm.nMa = tree$elm.nMa
  
       ## search for "and"
       filter = elm.nMa[,3]==0
       if(sum(filter)==0) {
  #          search = FALSE
  #          searchWholeTree = FALSE
           return(treeForSearch)
       }
       search = TRUE
     }

     if( sum(filter)>=1){
       
       newMa = elm.nMa[filter, , drop=FALSE]

       if(ifD) print(newMa)
       
       j = 1 # loop though node with "and" operator
       remove.list = NULL
       while ( j <= sum(filter) & search){
         if(ifD) print(paste("j=", j))
         curCh.idx = newMa[j, 1]
         curChRow = match(curCh.idx, elm.nMa[,1])
         curNode.level = newMa[j, 2]
         par.idx = binaTree.shuffle.findNodePar(tree, curCh.idx)

         if(is.null(par.idx)) stop (paste(fN, " No parent for a child node!"))

         if(par.idx != 0) {

            parRow = match(par.idx, elm.nMa[,1])
            if(ifD){
              
              print(paste("Matched parent row =", parRow, " parent id =", par.idx))
            }
            # find parent node is using "and"
            if( elm.nMa[parRow, 3] == 0){
              if(ifD) print("find a need for merge")
              search = FALSE
       
              ## find out the element of current node and the other node id
              curNodeMa = tree$elm.nlist[[curCh.idx]]
       
              parNodeMa = tree$elm.nlist[[par.idx]]

              ## need to be the same id as well a node
              filter.curkidRow = parNodeMa[,1]==curCh.idx & parNodeMa[,2]==1

              othNode.row = parNodeMa[!filter.curkidRow,]
              curNode.rows = curNodeMa
              tree$elm.nlist[[par.idx]] = rbind(othNode.row, curNode.rows)

              ## remove from the nMa the current node, will NOT remove the node from nlist for simplicity
              remove.list = c(remove.list, curChRow)

              if(ifD) {
                print("elm.nlist and nMa after updating the partent of old node")
                #print(tree$elm.nlist)
                #print(tree$elm.nMa)
              }
   
              treeForSearch = tree
              tree = NULL
            
            } ## if( elm.nMa[match.id, 3] == 0){
          } # if(par != 0) {
         j = j +1
       } ##  while ( j <= sum(filter) & search){

       if(length(remove.list)>0){
         elm.nMa = elm.nMa[-c(remove.list),,drop=FALSE]
         print(paste("removed:", paste(c(remove.list), collapse=";")))
         if(ifD) print(elm.nMa)
         treeForSearch$elm.nMa = elm.nMa
# 
#          remove.nlist  = elm.nMa[remove.list,1]
#          treeForSearch$elm.nlist = treeForSearch$elm.nlist [-c(remove.nlist)]
       }

       if(j > sum(filter) & search) searchWholeTree = FALSE
       if(ifD & j > sum(filter) & search ) print("stopped")
     } ## if( sum(filter)>=1){
     
   } ## while(searchWholeTree){

 ##return(NULL)
 return(treeForSearch)
}

binaTree.parser <-
function(str){

  ifD = FALSE

  ## the first character is "(", "not", "and", "or", or variable name
  b.lev=0
  b.open = FALSE
  cur.idx = 0
  cur.str = str

  binaTree = binaTree.constr()

  binaTree = binaTree.parser.proc(str, binaTree, curLevel=0, first.seg=TRUE)

  # check wether the level 0 is closed
  nMa = binaTree$elm.nMa

  if (length(binaTree$elm.nlist)>0){
    cur.level = nMa[, 2]
    if(sum(cur.level==0)>1) stop()
    zero.level = which(cur.level==0)
    if( nMa[zero.level, 4]==1){
      ## close the last node
      binaTree=binaTree.closeNode(binaTree,  elm.level=0)
    }
    
  }
  
  return(binaTree)

}

binaTree.parser.proc <-
function(str, binaTree, curLevel, first.seg=FALSE){
  ## need to grow the list
  ifD = FALSE
  
  str = util.str.rmSpacePadder(str, rm=c("b", "e"))
  
  if(ifD) print(qp("binaTree.parser.proc::str=", str, "|end"))
  if(ifD) print(qp("binaTree.parser.proc::curLevel=", curLevel, "|end"))
  #if(ifD) print(binaTree)
  
   ## patterns: parse out the seg separated by (/)
   ## "" hides the part has already been processed, ')' hides the part not necessary there
   ## 7 "("A and B and C')'
   ## 5  and A and B')'
   ## 8 ")") and A 
   ## 1 "((A and "(B and C)) and
   ## 2   A and (
   ## 4  and (
   
   b.level = util.findBrace(str)
   if(ifD) print(b.level)
   ## if no more brace, it is ending!!!
   t.elm=NULL
   if(is.null(b.level)){

     if(first.seg){
       ## #7:"("A and B and C')' only one elm
       ## can only be and/or/not
       ## add elm, elment only
       if(ifD) print("first.seg")
       t.elm = binaTree.parser.pVarElm(str, type="middle")
       if(is.null(t.elm)) stop()

       ## changed, so it would not allow two elments as one element variable
       if(t.elm$mixed==2){
         if(ifD) print("t.elm$mixed==2")
         # add the node
         
         ## add the element and add the node as well. Note that at the maximum, only two elements is allowed
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[1], 0)
         binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[2], 0)
         ## close the node
         binaTree = binaTree.closeNode(binaTree, curLevel)
       }else{       
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
       }
       
       return(binaTree)
     }else{
       ## #5: and A and C')', ending
       ## can be and/or,  must be mix node and
       ## ENDING
       ## add elm, add node, close node
       t.elm = binaTree.parser.pVarElm(str, type="start")
       if(is.null(t.elm)) stop("binaTree.parser.pVarElm(str, type='start')")
       if(t.elm$mixed==1){
         ## add the node
         binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
         ## add elm first
         binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
         ## close the node
         binaTree = binaTree.closeNode(binaTree, curLevel)
       }else{
         stop("binaTree.parser.pVarElm(str, type='start'):: not node" )
       }
       return(binaTree)
     } # if(first.seg){

   } # if(is.null(b.level)){
  
   ## if more brace, keep subsetting
   if(!is.null(b.level)){
     if (b.level$brace==1){
       ## #8 or #5 or #7 or having ")"
       if(first.seg){
         stop(" Can not have ')' before '('.")
       }
       
       if(b.level$idx==1){
           ## #8
           ## process node, then subset
           ## close the node
           
           binaTree = binaTree.closeNode(binaTree, curLevel)
           curLevel = curLevel-1
           
           ## subset
           if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx+1, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=FALSE)
           }
           return(binaTree)
       }else{ # if(b.level$idx==1){
         ## #7 or #5
         ## #7: "("A and B and C')', #5:  and A and B')'
         ## depends on the first is boolean or factor
         ## process node, then subset
         segStr = substr(str, 1, b.level$idx-1)
         t.elm = binaTree.parser.pVarElm(segStr, type="un")
         if(ifD) print(t.elm) 
         if(is.null(t.elm)) stop( "binaTree.parser.pVarElm(str, type='un')" )
         if(t.elm$mixed==1){
           ## add the node
           binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)           
           ## print(t.elm)
           ## add elm first
           binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
           
           binaTree = binaTree.closeNode(binaTree, curLevel)
           curLevel = curLevel-1
         }else{
           if(t.elm$mixed==2){
             if(ifD) print("second:: t.elm$mixed==2")
             # add the node
             
             ## add the element and add the node as well. Note that at the maximum, only two elements is allowed
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[1], 0)
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx[2], 0)
             ## close the node
             binaTree = binaTree.closeNode(binaTree, curLevel)
             curLevel = curLevel - 1
           }else{
               
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
             ##print(" binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool) "  )
             ##print(binaTree)
             curLevel = curLevel-1
             ## if the last is an variable only element, need to close the node
           }
         }

         ## subset
         if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx+1, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=FALSE)
         }else{
             if( curLevel==0){
               binaTree = binaTree.closeNode(binaTree, curLevel)
             }else{
               stop(qp("curLevel=", curLevel, ", should be 0"))
             }
         }
         return(binaTree)
       } # if(b.level$idx==1){
     } # if (b.level$brace==1){

     if(b.level$brace==0){
       ## having "("
       if(b.level$idx==1){
         ## #1: "((A and "(B and C)) and
         ## see how many levels it opens, then subset
           add.lev = util.char.1stIdx(str, "(", match=FALSE)
           curLevel = curLevel+add.lev
           curStr = substr(str, add.lev+1, nchar(str))
           binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=FALSE)
           return(binaTree)
       }else{ # if(b.level$idx==1){
         ## #4: and (, #2: A and (

         ## process the node, then subset
         segStr = substr(str, 1, b.level$idx-1)
         t.elm = binaTree.parser.pVarElm(segStr, type="end")
         ##print(paste("t.elm=", t.elm))
         if(is.null(t.elm)) stop( " binaTree.parser.pVarElm(str, type='end')" )
         if(t.elm$mixed==1){
           ## if the boolean term itself
           if(is.null(t.elm$idx)){
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
           }else{
             binaTree = binaTree.addElm(binaTree, curLevel, t.elm$idx, t.elm$bool)
             binaTree = binaTree.addNode(binaTree, curLevel, t.elm$bool)
           }

         }else{
           stop( " binaTree.parser.pVarElm(str, type='end'):: not node")        
         }           
         
         ## subset
         if(b.level$id < nchar(str)){
             curStr = substr(str, b.level$idx, nchar(str))
             binaTree = binaTree.parser.proc(curStr, binaTree, curLevel, first.seg=FALSE)
         }else{
           stop ("find ( without )")
         }
         return(binaTree)
       }  # if(b.level$idx==1){
     }# if(b.level$brace==0){
   } # if(!is.null(b.level)){

}

binaTree.parser.pVarElm <-
function(str, type="un", check.format=TRUE){
    ifD = FALSE
    fN = "binaTree.parser.pVarElm:"

    str = util.str.rmSpacePadder(str, rm=c("b", "e"))
    if(ifD) print(qp("binaTree.parser.pVarElm::str=", str, "|begin"))
    if(ifD) print(qp("binaTree.parser.pVarElm::type=", type, "|begin"))
    
    idx = NULL
    bool = NA
    mixed = 0
    
    if(type=="un"){
      # see if starting with 'and' or 'or'or 'not'
      idx.and = util.char.1stStrIdx(str, find="and ")
      idx.or = util.char.1stStrIdx(str, find="or ")
      idx.not = util.char.1stStrIdx(str, find="not ")
      if(idx.not==1){ ## start with not
          mixed = 0
          bool = -1
          idx = util.str.splitAndOrNot(str, "not", check.space=check.format)
          if(!is.null(idx)){
            if(check.format){
              if(length(idx$items)==1){
                return(list(idx=idx$items, bool=bool, mixed=mixed))
              }
            }else{
              if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
  
              return(list(idx=idx$items, bool=bool, mixed=mixed))
            }
          }
      }

      if(idx.and==1){ ## start with and
          mixed = 1
          bool = 0
          idx = util.str.splitAndOrNot(str, "and", check.space=check.format)

          if(is.null(idx)) return(NULL)

          if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
  
          return(list(idx=idx$items, bool=bool, mixed=mixed))
      }
      if(idx.or==1){ ## start with or
          mixed = 1
          bool = 1
          idx = util.str.splitAndOrNot(str, "or", check.space=check.format)

          if(is.null(idx)) return(NULL)
          if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
          return(list(idx=idx$items, bool=bool, mixed=mixed))
      }

      list.and = util.str.splitAndOrNot(str, " and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, " or ", check.space=check.format)

      ## operator is in the middle
      #print(list.and)
      #print(list.or)
      
      if (!is.null(list.and)){
        if (length(list.and$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
        return(list(idx=list.and$items, bool=0, mixed=2))
      }
      
      if (!is.null(list.or)){
        if (length(list.or$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
        return(list(idx=list.or$items, bool=1, mixed=2))
      }
    }
    
    if(type=="start"){
      list.and = util.str.splitAndOrNot(str, "and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, "or ", check.space=check.format)

      if (!is.null(list.and)){
        if (length(list.and$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.and$items, bool=0, mixed=1))
      }
      if (!is.null(list.or)){
        if (length(list.or$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.or$items, bool=1, mixed=1))  
      }
    }

    if(type=="end"){
      if(str=="and")
        return(list(idx=NULL, bool=0, mixed=1))
      
      if(str=="or")
        return(list(idx=NULL, bool=1, mixed=1))
      
      list.and = util.str.splitAndOrNot(str, "and", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, "or", check.space=check.format)

      if (!is.null(list.and)){
        if (length(list.and$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.and$items, bool=0, mixed=1))
      }
      
      if (!is.null(list.or)){
        #print(list.or)
        if (length(list.or$items)>1) stop(paste(fN, "do not allow more than two variables without a brace."))
        return(list(idx=list.or$items, bool=1, mixed=1))  
      }
    }    

    
    if(type=="middle"){
      idx.not = util.char.1stStrIdx(str, find="not ")
      if(idx.not==1){ ## start with not
          mixed = 0
          bool = -1
          idx = util.str.splitAndOrNot(str, "not", check.space=check.format)
          if(!is.null(idx)){
            if(check.format){
              if(length(idx$items)==1) return(list(idx=idx$items, bool=bool, mixed=mixed))
            }else{
              if (length(idx$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))
              return(list(idx=idx$items, bool=bool, mixed=mixed))
            }
          }
      }
      list.and = util.str.splitAndOrNot(str, " and ", check.space=check.format)
      list.or  = util.str.splitAndOrNot(str, " or ", check.space=check.format)

      #print(list.and)
      #print(list.or)
      if (!is.null(list.and)){
        if (length(list.and$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))      
        return(list(idx=list.and$items, bool=0, mixed=2))
      }
      if (!is.null(list.or)){
        if (length(list.or$items)>2) stop(paste(fN, "do not allow more than two variables in an element."))      
        return(list(idx=list.or$items, bool=1, mixed=2))
      }
      if(is.null(list.and) & is.null(list.or)){
        # could be just one element
         if(util.char.1stIdx(str, " ")==0){
           return(list(idx=str, bool=0, mixed=0))
         }
      }
                
    }    

    return(NULL) 
}

binaTree.patternForm <-
function(binaTree, elmMList, bkIdx.RuleOrder){
  ifD = FALSE
  fN = "binaTree.patternForm:"

  if(ifD) {
    print("binaTree")
    print(binaTree)
    print("elmMList")
    print(elmMList)
    print(bkIdx.RuleOrder)

  }


  # place to hold the matching hap pairs
  mgrid = NULL

  unique.bk = unique(bkIdx.RuleOrder)
  unique.bk = sort(unique.bk)
  
  # total number of HRCB
  bkCt = length(unique.bk)

  nMa = binaTree$elm.nMa
  nCt = nrow(binaTree$elm.nMa)

  ## in binaTree, all the op among node is "OR".
  ## but within a node, it is one element or a boolean term of multiple element
  ## with more than one element.

  mgrid = NULL
  ##process node with "and" op
  filter = nMa[,3]==0
  # find out the node with "and" as operator
  range = (1:nCt)[filter]
  ## need to parse the node with "and" as operator
  if(sum(filter)>=1){
    for ( i in range){
      # for each node 
      node.idx = nMa[i,1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      if(sum(nodeMa[,2])>=1) stop(paste(fN, "'and' node with node as element."))

      # nodeMa[,1] gives the element's id, which follow the elmBkIdx
      bkIdx.seqInRule = nodeMa[,1]

      bkm.ct = table(bkIdx.seqInRule)
      bkm.bench = NULL
      bkm.bench = as.integer(dimnames(bkm.ct)[1][[1]])

      ## if two elements in one bk need to be meet, need to update elmMList and bkIds
      newMList = NULL
      
      for(j in 1:length(bkm.bench)){
        ## for each bkIdx related to the element in the node
        tmp.ct = bkm.ct[j]
        elm.m.ids =  bkm.bench[j]
        if(tmp.ct>1){ # if more than one element from the same bk
          stop (paste(fN, "not implemented"))
#           elm.m.ids = nodeMa[ bkm.idx==bkm.cur.idx, 1]
#           #obtain all the elmMlist for this set of elm
#           tmp.allMatch = elmMList[elm.m.ids]
#           ## keep the first set of meeting one
#           tmp.finalSet = tmp.allMatch[[1]]
#           final.fit = rep(1, times=nrow(tmp.allMatch))
#           for(ii in 2:tmp.ct){
#             ## do an intersection operation, first assume all are in
#             ## then match two indexes
#             tmp.cur = tmp.allMatch[[ii]]
#             filter = is.element(tmp.finalSet[,1], tmp.cur[,1])
#             ## match the first idx
#             final.fit = final.fit & filter
#             ## match the second idx
#             filter = match(tmp.finalSet[,2], tmp.cur[,2])
#             final.fit = final.fit & filter            
#           } # for(ii in 2:tmp.ct){
#           if(sum(filter)<=0) stop (paste(fN, "no fitted hap pairs for the more than two matching elment in a block."))
# 
#           tmp.ma = tmp.finalSet[filter,, drop=FALSE ]
#           newMList = c(newMList, list(tmp.ma))
        }else{
          newMList = c(newMList, elmMList[elm.m.ids])
        } # if(tmp.ct>1){ # one element for one bk
        
      } # for(j in 1:length(bkm.bench))
      ## form the pattern for add node, add to the match.grid1 and 2
      bkMVec = bkIdx.RuleOrder[ bkm.bench]
      # match the bkm.cur.idx to a position in the whole pattern
      mat.idx = match(bkMVec, unique.bk)
      hapPair.matched=hapPair.match(mat.idx, bkCt, newMList)
      if(ifD) print(hapPair.matched)
      mgrid = rbind(mgrid, hapPair.matched)
    } # for ( i in range){ parse node with "and" as operator
  } # if(sum(filter)>=1){

  ## parse all the element in "or" node

  for (i in 1:nCt){
     curNode = nMa[i,]
     ## only for "or" node
     if(curNode[3]==1){
      node.idx = curNode[1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      ## only considering the element variables
      if(sum(nodeMa[,2])<nrow(nodeMa)){
        if(ifD) print("find the node linked with one element!!")
        elm.ids = nodeMa[nodeMa[,2]==0,1]
        #if(length(elm.ids)!=1) stop(paste(fN, "number of element in a node is not one"))

        # it is OK to have more than one element
        for (ij in 1:(length(elm.ids))){
          bkm.idx = bkIdx.RuleOrder [ elm.ids[ij] ]
          mat.idx = match(bkm.idx, unique.bk)
          if (ifD) print(bkm.idx)
          if (ifD) print(unique.bk)
          #print(elmMList[elm.ids[ij]])
          hapPair.matched=hapPair.match(mat.idx, bkCt, elmMList[elm.ids[ij]])
          mgrid = rbind(mgrid, hapPair.matched)
        }
        
        if(ifD) print(mgrid)
      }# else it is a all node node      
     } # if(curNode[3]==1){
  } # for i in (1:nCt){
  
  return( mgrid )
}

binaTree.patternFormOLD <-
function(binaTree, elmMList, bkIdx.RuleOrder){
  ifD = TRUE
  fN = "binaTree.patternForm:"

  if(ifD) {
    print("binaTree")
    print(binaTree)
    print("elmMList")
    print(elmMList)
    print(bkIdx.RuleOrder)

  }
  # place to hold the matching hap pairs
  mgrid = NULL

  unique.bk = unique(bkIdx.RuleOrder)
  unique.bk = sort(unique.bk)
  
  # total number of HRCB
  bkCt = length(unique.bk)

  nMa = binaTree$elm.nMa
  nCt = nrow(binaTree$elm.nMa)

  ## in binaTree, all the op among node is "OR".
  ## but within a node, it is one element or a boolean term of multiple element
  ## with more than one element.

  mgrid = NULL
  ##process node with "and" op
  filter = nMa[,3]==0
  # find out the node with "and" as operator
  range = (1:nCt)[filter]
  ## need to parse the node with "and" as operator
  if(sum(filter)>=1){
    for ( i in range){
      # for each node 
      node.idx = nMa[i,1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      if(sum(nodeMa[,2])>=1) stop(paste(fN, "'and' node with node as element."))

      # nodeMa[,1] gives the element's id, which follow the elmBkIdx
      bkIdx.seqInRule = nodeMa[,1]

      bkm.ct = table(bkIdx.seqInRule)
      bkm.bench = NULL
      bkm.bench = as.integer(dimnames(bkm.ct)[1][[1]])

      ## if two elements in one bk need to be meet, need to update elmMList and bkIds
      newMList = NULL
      
      for(j in 1:length(bkm.bench)){
        ## for each bkIdx related to the element in the node
        tmp.ct = bkm.ct[j]
        elm.m.ids =  bkm.bench[j]
        if(tmp.ct>1){ # if more than one element from the same bk
          stop (paste(fN, "not implemented"))
#           elm.m.ids = nodeMa[ bkm.idx==bkm.cur.idx, 1]
#           #obtain all the elmMlist for this set of elm
#           tmp.allMatch = elmMList[elm.m.ids]
#           ## keep the first set of meeting one
#           tmp.finalSet = tmp.allMatch[[1]]
#           final.fit = rep(1, times=nrow(tmp.allMatch))
#           for(ii in 2:tmp.ct){
#             ## do an intersection operation, first assume all are in
#             ## then match two indexes
#             tmp.cur = tmp.allMatch[[ii]]
#             filter = is.element(tmp.finalSet[,1], tmp.cur[,1])
#             ## match the first idx
#             final.fit = final.fit & filter
#             ## match the second idx
#             filter = match(tmp.finalSet[,2], tmp.cur[,2])
#             final.fit = final.fit & filter            
#           } # for(ii in 2:tmp.ct){
#           if(sum(filter)<=0) stop (paste(fN, "no fitted hap pairs for the more than two matching elment in a block."))
# 
#           tmp.ma = tmp.finalSet[filter,, drop=FALSE ]
#           newMList = c(newMList, list(tmp.ma))
        }else{
          newMList = c(newMList, elmMList[elm.m.ids])
        } # if(tmp.ct>1){ # one element for one bk
        
      } # for(j in 1:length(bkm.bench))
      ## form the pattern for add node, add to the match.grid1 and 2
      bkMVec = bkIdx.RuleOrder[ bkm.bench]
      # match the bkm.cur.idx to a position in the whole pattern
      mat.idx = match(bkMVec, unique.bk)
      hapPair.matched=hapPair.match(mat.idx, bkCt, newMList)
      if(ifD) print(hapPair.matched)
      mgrid = rbind(mgrid, hapPair.matched)
    } # for ( i in range){ parse node with "and" as operator
  } # if(sum(filter)>=1){

  ## parse all the element in "or" node

  for (i in 1:nCt){
     curNode = nMa[i,]
     ## only for "or" node
     if(curNode[3]==1){
      node.idx = curNode[1]
      nodeMa = binaTree$elm.nlist[[node.idx]]
      ## only considering the element variables
      if(sum(nodeMa[,2])<nrow(nodeMa)){
        if(ifD) print("find the node linked with one element!!")
        elm.ids = nodeMa[nodeMa[,2]==0,1]
        if(length(elm.ids)!=1) stop(paste(fN, "number of element in a node is not one"))

        bkm.idx = bkIdx.RuleOrder [ elm.ids ]
        mat.idx = match(bkm.idx, unique.bk)
        if (ifD) print(bkm.idx)
        if (ifD) print(unique.bk)
        hapPair.matched=hapPair.match(mat.idx, bkCt, elmMList[elm.ids])
        mgrid = rbind(mgrid, hapPair.matched)
      }# else it is a all node node      
     } # if(curNode[3]==1){
  } # for i in (1:nCt){
  
  return( mgrid )
}

binaTree.setApply <-
function(binaTree, setList){
  ifD = FALSE
  # the setList already contain the set for each variable

  data.n=NULL
  
  nCt = length(binaTree$elm.nlist)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n = rep(list(NULL), nrow(nMa))

    new.order = order(nMa[,2], decreasing = TRUE)
    if(ifD) print(nMa)
    nMa = nMa[new.order, ,drop=FALSE]
    if(ifD) print(nMa)
    #nlist = binaTree$elm.nlist[new.order]
    nlist = binaTree$elm.nlist
    if(ifD) print(nlist)
    for( i in 1:nCt ){
      cur.row = nMa[i, ,drop=FALSE]
      if(ifD) print("cur node row after sorting")
      if(ifD) print(cur.row)
      ## col seq= id, level, bool, open
      
      ## find the var list
      cur.elm = nlist[[ cur.row[1] ]]
      if(ifD) print("cur node element")
      if(ifD) print(cur.elm)
      dd.velm=NULL
      dd.nelm=NULL
      
      # if variable element, get from data.v
      ## only three possibilties: both are variable, both are node, one variable and one node
      if( sum(cur.elm[,2]==0)  >=1){
        # dd.velm = data.v[, cur.elm[,1][cur.elm[,2]==0]]
        dd.velm = setList[ c(cur.elm[,1][cur.elm[,2]==0]) ]
        
        if (ifD) {
          if(is.null(dim(dd.velm))){
            print(dd.velm[1:5])
          }else{
            print(dd.velm[1:5,])
          }
        }
      }
      # if node element, get from data.n
      if(sum(cur.elm[,2]==1)>=1){
        dd.nelm = data.n[ c(cur.elm[,1][cur.elm[,2]==1]) ]
        
        if (ifD) {
          if(is.null(dim(dd.nelm))){
            print(dd.nelm[1:5])
          }else{
            print(dd.nelm[1:5,])
          }
        }
      }

      dd.cur = c(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]
     
      if(cur.bool==0){ ## and -> intersection

        if(length( dd.cur[[1]] )==0) {
          warning("length( dd.cur[[1]] )==0")
          data.n[[ cur.row[1]   ]] = NULL
        }
        if(length( dd.cur[[2]] )==0) {
          warning("length( dd.cur[[2]] )==0")
          data.n[[ cur.row[1]   ]] = NULL
        }
        if(length(dd.cur[[1]])!=0 & length(dd.cur[[2]])!=0){
          filter = match(dd.cur[[1]], dd.cur[[2]], nomatch=0)
          filter = filter[filter>0]
          data.n[[  cur.row[1]  ]] = dd.cur[[2]][ filter ]
        }
      }
      if(cur.bool==1){ ## or -> union
        
        data.n [[  cur.row[1]  ]] = unique(c(dd.cur[[1]], dd.cur[[2]]))
        if( length(dd.cur[[1]])==0 & length(dd.cur[[2]])==0 ){
          data.n[[  cur.row[1]  ]] = NULL
        }
      }
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

    if(nCt==0){
      # only one variable only element
      return(setList[[1]])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[[t.id]]
      return(re)
    }

}

binaTree.shuffle.findElmPar <-
function(tree, id.elm){
  fN = "binaTree.shuffle.findElmPar:"
  elm.nMa = tree$elm.nMa

  ## loop through each node in nMa
  search = TRUE
  i = 1
  child.node=NULL
  while ( i <= length(tree$elm.nlist) & search){
    curMa = tree$elm.nlist[[i]]
    filter = curMa[ ,2]==0 & curMa[,1]==id.elm
    if(sum(filter)>=1){
      search = FALSE
      child.node = i
    }
    i = i+1
  }

  if(search){
    stop(paste(fN, " cannot find the node using the current element."))
  }

  return(child.node)

}

binaTree.shuffle.findNodePar <-
function(tree, id.node){
  fN = "binaTree.shuffle.findNodePar:"
  elm.nMa = tree$elm.nMa

  par.level = elm.nMa[elm.nMa[,1]==id.node,2]-1
  if(par.level<0) {
    #warning(paste(fN, " top level has no parent for id =", id.node))
    return(0)
  }

  par.posiId = elm.nMa[elm.nMa[,2]==par.level,1]

  if(length(par.posiId)<=0 ) stop( paste(fN, " no node belongs to the upper level"))

  i = 1

  posId = NULL
  ## multiple parents may have it as a child
  search =TRUE
  while(i<=length(par.posiId) & search){

    posParId = par.posiId[i]
    nodeMa = tree$elm.nlist[[ posParId ]]
    yesNode = nodeMa[ ,2]==1
    childId = nodeMa[yesNode,1]
    if(sum(is.element(id.node, childId))==1){
      posId = posParId
      search = FALSE
    }
    i = i +1;
  }

  return(posId)
 
}

binaTree.shuffleNode <-
function(treeForSearch){
 fN = "binaTree.shuffleNode:"
 ifD = FALSE
 
 searchWholeTree = TRUE
 tt = 0
 search = FALSE
 # do an iterative search for ( or ) and , keep search until all the "or" have no parent with "and" operator
 while(searchWholeTree){
     tt = tt+1
     if(tt>100) stop(paste(fN, "iterative procdure exceed 100 times. Sometime is wrong!"))
     
     # if the tree is updated, use the updated one. Otherwise, search the next possible choice
     if(!search){
       if(ifD) {
         print("Updated the tree:")
         print(treeForSearch)
       }
       tree = treeForSearch
       elm.nMa = tree$elm.nMa
  
       ## search for "or"
       filter = elm.nMa[,3]==1
       if(sum(filter)==0) {
  #          search = FALSE
  #          searchWholeTree = FALSE
           return(treeForSearch)
       }
       search = TRUE
     }

     if( sum(filter)>=1){
       
       newMa = elm.nMa[filter, , drop=FALSE]

       if(ifD) print(newMa)
       
       j = 1 # loop though node with "or" operator
       while ( j <= sum(filter) & search){
         if(ifD) print(paste("j=", j))
         curCh.idx = newMa[j, 1]
         curNode.level = newMa[j, 2]
         par.idx = binaTree.shuffle.findNodePar(tree, curCh.idx)

         if(is.null(par.idx)) stop (paste(fN, " No parent for a child node!"))

         if(par.idx != 0) {
           
            parRow = match(par.idx, elm.nMa[,1])
            if(ifD){
              
              print(paste("Matched parent row =", parRow, " parent id =", par.idx))
            }
            # find parent node is using "and"
            if( elm.nMa[parRow, 3] == 0){
              if(ifD) print("find a case")
              search = FALSE
       
              ## find out the element of current node and the other node id
              curNodeMa = tree$elm.nlist[[curCh.idx]]
       
              parNodeMa = tree$elm.nlist[[par.idx]]

              filter.curKid = parNodeMa[,1]==curCh.idx & parNodeMa[,2]==1
              
              twoKid = parNodeMa[,1]
              othNode.idx = twoKid[!filter.curKid]
              othNode.type = parNodeMa[!filter.curKid, 2]

              ## if the other node is a nodel
              ## change the other kid level one down
              if(othNode.type==1){
                othNodeRow = match(othNode.idx, elm.nMa[,1])
                elm.nMa[othNodeRow, 2] =  curNode.level+1

                ## create a duplicate of the othNode
                othNodeMaDup  = tree$elm.nlist[[othNode.idx]]
                othNodeDup.New = length(tree$elm.nlist)+1
                tree$elm.nlist = c(tree$elm.nlist, list(othNodeMaDup))
                elm.nMa = rbind(elm.nMa, c(othNodeDup.New, curNode.level+1, elm.nMa[othNode.idx, 3], 0))
                ma =  matrix(c(curNodeMa[2,1], othNodeDup.New, curNodeMa[2,2], othNode.type), ncol=2, nrow=2, byrow=FALSE)
              }else{
                ma =  matrix(c(curNodeMa[2,1], othNode.idx, curNodeMa[2,2], othNode.type), ncol=2, nrow=2, byrow=FALSE)
              }
              
              ## create the new node
              othNode.New = length(tree$elm.nlist)+1
              tree$elm.nlist = c(tree$elm.nlist, list(ma))
   
              if(ifD) {
                print("elm.nlist after adding the new node")
                #print(tree$elm.nlist)
   
              }
   
              ## added the new node to the nMa
              elm.nMa = rbind(elm.nMa, c(othNode.New, curNode.level, 0, 0))
              
       
              ## update the old node with "and"
              curNodeMa[2,] = c(othNode.idx, othNode.type)
              elm.nMa[curCh.idx, 3]=0
              tree$elm.nlist[[curCh.idx]]=curNodeMa
   
              if(ifD){
                print("elm.nlist after updating the current node")
                #print(tree$elm.nlist)
              }
              
       
              ## update the parent of old node
              parNodeMa =  matrix(c(curCh.idx, othNode.New, 1, 1), ncol=2, nrow=2, byrow=FALSE)
              elm.nMa[parRow, 3]=1
              tree$elm.nlist[[par.idx]]=parNodeMa
              tree$elm.nMa = elm.nMa
              
              if(ifD) {
                print("elm.nlist and nMa after updating the partent of old node")
                #print(tree$elm.nlist)
                #print(tree$elm.nMa)
              }
   
              treeForSearch = tree
              tree = NULL
            
            } ## if( elm.nMa[match.id, 3] == 0){
          } # if(par != 0) {
         j = j +1
       } ##  while ( j <= sum(filter) & search){
       if(j > sum(filter) & search) searchWholeTree = FALSE
       if(ifD & j > sum(filter) & search ) print("stopped")
     } ## if( sum(filter)>=1){
     
   } ## while(searchWholeTree){

 ##return(NULL)
 return(treeForSearch)

}

binaTree.toStr <-
function(binaTree){
  
  ifD = FALSE

  ## process the variable only element first, regardless the level
  vCt = nrow(binaTree$elm.vMa)
  data.v = rep("", times=vCt)
  
  for( i in 1:vCt){
    # print(i)
    cur.row = binaTree$elm.vMa[i,]
    ## col seq= id, level, bool
    ## find the var list
    cur.elm = binaTree$elm.vlist[[ i ]]

    cur.bool = cur.row[3]
    if(length(cur.elm)>1){
     
      if(cur.bool==0){ ## and
        data.v[i]= paste(cur.elm, collapse=" and ")
      }
      if(cur.bool==1){ ## or
        data.v[i]= paste(cur.elm, collapse=" or ")
      }
    }else{
      if(cur.bool==-1){ ## not
        data.v[i]= paste("(not ", cur.elm, ")", sep="")
      }else{
        data.v[i]= cur.elm
      }
    }
  } # for( i in 1:vCt){

  ## if have node element, process the node from deeper level to zero
  data.n=NULL
  nCt = nrow(binaTree$elm.nMa)
  if(nCt>=1){
    nMa = binaTree$elm.nMa
    
    data.n=rep("", times=nCt)

    new.order = order(nMa[,2], decreasing = TRUE)
    if(ifD) print(nMa)
    nMa = nMa[new.order, ,drop=FALSE]
    if(ifD) print(nMa)
    #nlist = binaTree$elm.nlist[new.order]
    nlist = binaTree$elm.nlist
    if(ifD) print(nlist)
    for( i in 1:nCt ){
      cur.row = nMa[i, ,drop=FALSE]
      if(ifD) print("cur node row after sorting")
      if(ifD) print(cur.row)
      ## col seq= id, level, bool, open
      
      ## find the var list
      cur.elm = nlist[[ cur.row[1] ]]
      if(ifD) print("cur node element")
      if(ifD) print(cur.elm)
      dd.velm=NULL
      dd.nelm=NULL
      # if variable element, get from data.v
      if(sum(cur.elm[,2]==0)>=1){
        dd.velm = data.v[cur.elm[,1][cur.elm[,2]==0]]
      }
      # if node element, get from data.n
      if(sum(cur.elm[,2]==1)>=1){
        dd.nelm = data.n[cur.elm[,1][cur.elm[,2]==1]]
      }

      
      dd.cur = cbind(dd.velm, dd.nelm)
      if(ifD) print(cur.row)
      cur.bool = cur.row[3]

      if(i!=nCt){
        if(cur.bool==0){ ## and
          data.n[cur.row[1]]= paste("(" , paste(dd.cur, collapse= " and "), ")", sep="")
        }
        if(cur.bool==1){ ## or
          data.n[cur.row[1]]= paste("(" , paste(dd.cur, collapse= " or "), ")", sep="")
        }
      }else{ # for the top node, no need for outside brace
        if(cur.bool==0){ ## and
          data.n[cur.row[1]]= paste(dd.cur, collapse= " and ")
        }
        if(cur.bool==1){ ## or
          data.n[cur.row[1]]= paste(dd.cur, collapse= " or ")
        }  
      }
    } ## for( i in 1:nCt ){
    
  } ##  if(nCt>=1){

    if(nCt==0){
      # only one variable only element
      return(data.v[,1])
    }else{
      # return the last element (must be a node)
      t.id = binaTree$curElmId
      re = data.n[t.id]
      if(ifD) print(data.n)
      return(re)
    }

}

bindHapBkGenoMaps <-
function (hapBkMap=NULL, genoMap) 
{
	if(is.null(hapBkMap)){
		# when all the marker are single marker
		snpCtNum = nrow(genoMap$df)/3
		genomeMarkerInfo = matrix(NA, ncol=5, nrow=snpCtNum)
		genomeMarkerInfo[,1]=rep(1, snpCtNum)
		genomeMarkerInfo[,2]=1:snpCtNum
		genomeMarkerInfo[,3]=1:snpCtNum
		genomeMarkerInfo[,4]=1:snpCtNum
		genomeMarkerInfo[,5]=rep(1, snpCtNum)

		genoIndex = which(genomeMarkerInfo[, 5] == 1)
		
		genoMapDf = cbind(genoMap$df, hapLens = rep(1, nrow(genoMap$df)))
		genoOnlyMap = dfToHapBkMap(genoMapDf, keyCol = NULL, chCol = genoMap$chCol, 
				blockCol = genoMap$genomeSeqCol, expCol = genoMap$expCol, 
				probCol = genoMap$probCol, hapLenCol = genoMap$genomeSeqCol, 
				beginCol = NULL, endCol = NULL, snpBase = genoMap$snpBase, 
				re.bf = TRUE, re.javaGUI = TRUE)		
		
		re = list(hapBkOnlyMap=NULL, genoOnlyMap = genoOnlyMap, genomeMarkerInfo = genomeMarkerInfo, genoIndex = genoIndex)
		return(re)
		
	}

	if (is.null(hapBkMap$beginCol)) {
		stop("hapBkMap$beginCol is null, methods not implemented.")
	}
	if (is.null(hapBkMap$dfStr)) {
		stop("hapBkMap$dfStr is null, methods not implemented.")
	}
	hapOnlyMap = hapBkMap
	genoMapDf = cbind(genoMap$df, hapLens = rep(1, nrow(genoMap$df)))
	genoOnlyMap = dfToHapBkMap(genoMapDf, keyCol = NULL, chCol = genoMap$chCol, 
			blockCol = genoMap$genomeSeqCol, expCol = genoMap$expCol, 
			probCol = genoMap$probCol, hapLenCol = genoMap$genomeSeqCol, 
			beginCol = NULL, endCol = NULL, snpBase = genoMap$snpBase, 
			re.bf = TRUE, re.javaGUI = TRUE)
	qu = NULL
	qu.type = NULL
	markers_gb = NULL
	markers_ge = NULL
	markersIndex_genome = 0
	oldCh = "ooo"
	newCh = "ooo"
	hapIdx = 0
	chkey = unlist(lapply(genoOnlyMap$keys, FUN = function(item) {
						re = util.str.seqCutter(item, delims = "-")[1]
						re
					}))
	chkey.hap = unlist(lapply(hapOnlyMap$keys, FUN = function(item) {
						re = util.str.seqCutter(item, delims = "-")[1]
						re
					}))
	chkeyUni = unique(chkey)
	chkeyUni.hap = unique(chkey.hap)
	chkeyDiff = setdiff(chkeyUni, chkeyUni.hap)
	if (length(chkeyDiff) == 0) {
		chkeyDiff = NULL
	}
	genomeMarkerInfo = NULL
	
	for (ikey in chkeyUni) {
		#print(ikey)
		snp.allCt = sum(chkey == ikey)
		if (is.element(ikey, chkeyDiff)) {
			tmp.endCt = sum(chkey == ikey)
			tmp.ma = NULL
			if (tmp.endCt > 0) {
				for (tt2 in 1:tmp.endCt) {
					tmp.ma = rbind(tmp.ma,
							c(rep(markersIndex_genome + tt2, 3), 1))
				}
				tmp.ma = matrix(tmp.ma, ncol = 4)
				hap.ma = data.frame(ch = I(rep(ikey, times = nrow(tmp.ma))), 
						markers_gb = tmp.ma[, 1], markers_ge = tmp.ma[, 1], 
						qu = tmp.ma[, 1], qu.type = tmp.ma[, 4])
			}
			else {
				stop("This is wrong.")
			}
			hap.ma = hap.ma[order(hap.ma[, 2]), ]
			genomeMarkerInfo = rbind(genomeMarkerInfo, hap.ma)
			markersIndex_genome = markersIndex_genome + snp.allCt
		}
		else {
			ibks.f <- (chkey.hap == ikey)
			
			tmp.ma = NULL
			## need to know whether there are some singleton appearing at the beginning
			if( hapOnlyMap$markers_b[ibks.f][1] > 1){
				for (tt in 1:(hapOnlyMap$markers_b[ibks.f][1]-1) ) {
					tmp.ma = rbind(tmp.ma, c(rep(markersIndex_genome+tt, 3), 1))
				}
			}
			
			markers_gb = markersIndex_genome + hapOnlyMap$markers_b[ibks.f]
			markers_ge = markersIndex_genome + hapOnlyMap$markers_e[ibks.f]
			qu = (hapIdx + 1):(hapIdx + sum(ibks.f))
			qu.type = rep(0, times = sum(ibks.f))
			hap.ma = data.frame(ch = I(rep(ikey, times = sum(ibks.f))), 
					markers_gb, markers_ge, qu, qu.type)
			hapIdx = hapIdx + sum(ibks.f)
			
			if (sum(ibks.f) > 1) {
				tmp.diff = markers_gb[-1] -
						markers_ge[-(sum(ibks.f))] - 1
				tmp.pos1 = which(tmp.diff > 0)
				if (length(tmp.pos1)>0) {
					tmp.pos = hap.ma[tmp.pos1, c(1, 3), drop = FALSE]
					for (tt in 1:(length(tmp.pos1))) {
						for (tt2 in 1:(tmp.diff[tmp.pos1][tt])) {
							tmp.ma = rbind(tmp.ma, c(rep(tmp.pos[tt, 
															2] + tt2, 3), 1))
						}
					}
				}
			}
			tmp.endCt = snp.allCt - hapOnlyMap$markers_e[ibks.f][sum(ibks.f)]
			if (tmp.endCt > 0) {
				for (tt2 in 1:tmp.endCt) {
					tmp.ma = rbind(tmp.ma, c(
									rep(markers_ge[sum(ibks.f)] + 
													tt2, 3), 1))
				}
			}
			if (!is.null(tmp.ma)) {
				tmp.ma = matrix(tmp.ma, ncol = 4)
				single.ma = data.frame(ch = I(rep(ikey, times = nrow(tmp.ma))), 
						markers_gb = tmp.ma[, 1], markers_ge = tmp.ma[, 1], 
						qu = tmp.ma[, 1], qu.type = tmp.ma[, 4])
				hap.ma = rbind(hap.ma, single.ma)
				hap.ma = hap.ma[order(hap.ma[, 2]), ]
			}
			markersIndex_genome = markersIndex_genome + snp.allCt
			genomeMarkerInfo = rbind(genomeMarkerInfo, hap.ma)
		}
	}
	genomeMarkerInfo = genomeMarkerInfo[order(genomeMarkerInfo[, 2]), ]
	re = list(hapBkOnlyMap = hapOnlyMap, genoOnlyMap = genoOnlyMap, 
			genomeMarkerInfo = genomeMarkerInfo)
	hapIndex = which(genomeMarkerInfo[, 5] == 0)
	genoIndex = which(genomeMarkerInfo[, 5] == 1)
	re = c(re, list(hapIndex = hapIndex, genoIndex = genoIndex))
	return(re)
}

bkMap.constr <-
function(data, keyCol, hapLenCol=NULL, expCol, probCol, alleleCode=1:2, ...){
	fN = "bkMap.constr"
	
	if( is.character(data)){
		data = read.csv(data, ...)
		#print(qp(fN, ": Info on the ", data, " file."))
		#print(str(data))
	}
	
	colNum = ncol(data)
	#print(keyVal)
	keyVal = unique(data[,keyCol]) # not reordered, following the original order
	#print(keyVal)
	
	bks = NULL
	snpLen = 0
	bkLens = NULL
	keys = NULL
	bkEndingIdx = NULL
	if(is.null(hapLenCol)){
		# if the input dataset does NOT contain info. on haplotype block length
		for( i in 1:length(keyVal)){
			if(i==1){
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")          
				tmp.prob = tmp.prob/sum(tmp.prob)          
				r[,probCol]=tmp.prob
				
				bkLens = dim(r)[1]
				snpCtVec = rep(snpCt1, bkLens)
				
				r = cbind(r, snpCtVec)
				
				bks = list(r)          
				snpLen = snpLen + snpCt1
				keys = as.character(r[1,keyCol])
				snpCt = snpCt1
			}else{
				
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				r[,probCol]=tmp.prob
				
				bkLens = c(bkLens, dim(r)[1])
				snpCtVec = rep(snpCt1, dim(r)[1])
				
				r = cbind(r, snpCtVec)
				
				bks = c(bks, list(r))
				snpLen = snpLen + snpCt1
				keys = c(keys, as.character(r[1,keyCol]))
				snpCt = c(snpCt, snpCt1)
			}
			
			hapLenCol = dim(r)[2]
		}
	}else{
		# if the input dataset contain info. on haplotype block length
		for( i in 1:length(keyVal)){
			if(i==1){
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				if (min(r[,hapLenCol]==snpCt1)==0) stop("At least one of the block provides inconsistant haplotype size.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				
				r[,probCol]=tmp.prob
				
				bks = list(r)
				bkLens = dim(r)[1]
				snpLen = snpLen + snpCt1
				keys = as.character(r[1,keyCol])
				snpCt = snpCt1
			}else{
				r = data[data[,keyCol]==keyVal[i], ]
				
				## normalize the prob and check the haplen
				tmp.expLen  = unlist(lapply(r[,expCol], FUN=function(i){nchar(as.character(i))}  ))
				snpCt1 = tmp.expLen[1]
				if (min(tmp.expLen==snpCt1)==0) stop("At least one of the block has different lengths for haplotypes.")
				if (min(r[,hapLenCol]==snpCt1)==0) stop("At least one of the block provides inconsistant haplotype size.")
				
				tmp.prob = r[,probCol]
				if (max(tmp.prob<0)==1) stop("At least one of the block has negative haplotype frequencies.")
				if (sum(tmp.prob)==0) stop("At least one of the block has haplotype frequencies sum up to 0.")
				tmp.prob = tmp.prob/sum(tmp.prob)
				
				r[,probCol]=tmp.prob
				
				bks = c(bks, list(r))
				bkLens = c(bkLens, dim(r)[1])
				snpLen = snpLen + snpCt1
				keys = c(keys, as.character(r[1,keyCol]))
				snpCt = c(snpCt, snpCt1)
			}
		}
	}
	bkMap = list(bks = bks, bkLens = bkLens, keys = keys, snpLen = snpLen, expCol = expCol, hapLenCol=hapLenCol, probCol = probCol, snpCt=snpCt, alleleCode=alleleCode )
	return(bkMap)
	
}

bkMap.ESp.apply1Rule <-
function(bkMap, signalRule){
	ifD = FALSE
	fN = "bkMap.ESp.apply1Rule"
	if(ifD) print(paste(fN, "begin::"))
	
	signal = signalRule$signal
	
	total.bkct = length( bkMap$snpCt )
	
	# construct the beginning/ending index of snp in original map
	idx.be = matrix(NA, nrow=total.bkct, ncol=2)
	idx.be[,2]=cumsum(bkMap$snpCt)
	idx.be[,1]=c(1, cumsum(bkMap$snpCt)+1)[1:total.bkct]
	
	
	causalBkUniqOrderedIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=TRUE)
	## find out # of haplotype in each HRCB
	HRCB.hapCt = bkMap$bkLens [causalBkUniqOrderedIdx]
	HRCB.hapCtProd = cumprod(HRCB.hapCt)[length(causalBkUniqOrderedIdx)]
	
	causalBkIdx.ruleOrder = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=FALSE)
	
	#if(length(causalBkUniqOrderedIdx)==1){
	
	if(length(causalBkIdx.ruleOrder)==1){
		
	}else{
		## need to update the rule a little
		if(ifD) print( binaTree.toStr(signal))
		signal.n = binaTree.shuffleNode(signal)    
		if(ifD) print( binaTree.toStr(signal.n))
		signal.f = binaTree.merge(signal.n)
		if(ifD) print( binaTree.toStr(signal.f))
		
	}
	
	## have to find the matching with the HRCB only bkMap
	#causalBkIdx.ruleOrder = match(causalBkIdx.ruleOrder, causalBkIdx)
	
	elmMList = list()
	# for each block, need to find the hap pairs that are associated with higher risk
	
	for ( j in 1:length(causalBkIdx.ruleOrder)){
		bkIdx = causalBkIdx.ruleOrder[ j ]
		if(ifD) print(paste("bkIdx inside the bkMap =", bkIdx ))
		
		##  20Dec08Change!!!: straighten coding: 1, 2, 3 for 1-digitCoding, and 1, 2 for 2-digit
		
		genoMap =  hap2GenoBlock(hapExp = as.character(bkMap$bks[[bkIdx]][, bkMap$expCol]),
				hapProb = bkMap$bks[[bkIdx]][, bkMap$probCol], snpCt = bkMap$snpCt[bkIdx],
				alleleCode = bkMap$alleleCode )
		if(ifD) print(genoMap)
		item = idx.be[bkIdx,]
		## use recycling rules
		re = paste("v", rep(item[1]:item[2], each=2), sep="")
		newData.col = paste(re,  c("a", "b"), sep="")
		
		## singleSNPRule follow the original order in the genome
		## HARD CODE!!!HARD CODE: genoMa from hap2GenoBlock keep digits  1, 2 for 2-digit, need to change to 0/1 binary coding
		## 20Dec08Change!!! to use 2- for 0/1 binay coding, so the first alleleCode (changed to 1) correspond to the minor allele
		## and the second  (changed to 2) for major allele. Just need to coordinate with signal Rule.
		
		varIdx = match (signalRule$signal$elm.vlist[[j]], newData.col)
		treePred = (2-genoMap$genoMa)[,varIdx]
		## if the variable is not proceed by "not" operator
		## CHANGED on 8JAN09, from ==0 to !=-1. -1 means not
		if (signalRule$signal$elm.vMa[j,3]!=-1){
			nct = sum(treePred==1)
			if(ifD) print("not proceed by not")
			pairHRCB = genoMap$tb[  treePred==1, c(1,2), drop=FALSE ]
		}else{
			## if the variable is proceed by "not" operator
			nct = sum(treePred==0)
			if(ifD) print("proceed by not")
			pairHRCB = genoMap$tb[  treePred==0, c(1,2), drop=FALSE ]
		}
		if (ifD){
			print(paste("matching hap pairs in block ", j))
			print(pairHRCB)
		}
		elmMList=c(elmMList, list(pairHRCB))
	} # for ( j in 1:length(causalBkIdx.new)){
	#print("This is elmMList")
	#print(elmMList)
	#print(causalBkIdx.ruleOrder)
	#print("qing mark")
	
	## simplied if only one HRCB
	
	if (length(causalBkIdx.ruleOrder)==1){
		
		superDipIdx = apply(pairHRCB, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		
		superDipIdx = sort(superDipIdx)
		
		## FIX LATER!!!FIX : can further reduce the workload
		HRCBGrp = t(sapply(superDipIdx, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=TRUE))
		
		# print(str(HRCBGrp))
		
		HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=FALSE]
		
		HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
		
		#print(HRCBGrp.A)
		#print(HRCBGrp.BIdx)
		return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
		
	}
	
	hapGrids = binaTree.patternForm(binaTree=signal.f, elmMList=elmMList, bkIdx.RuleOrder=causalBkIdx.ruleOrder)
	if(ifD) print(str(hapGrids))
	##write.csv(hapGrids, file="hapGrids.csv") 
	
	## check whether no matching BK
	apply(hapGrids, 2, FUN=function(coo){
				qing.check = mean(is.na(coo))
				if(qing.check==1) stop(qp(fN, ": no matching pattern for one block"))
			})
	
	
	## need to find super-hap index for the matching patterns
	uniBkCt = ncol(hapGrids) / 2
	
	HRCB.hapCt2 = HRCB.hapCt
	
	if(ifD) print(hapGrids)
	supIdx = NULL
	
	for(n in 1:nrow(hapGrids)){
		
		if(ifD) {
			print(paste("proc row for pattern: row=", nrow(hapGrids)))
			#print(n)
			#print(hapGrids[n,])
		}
		if( n>240) {
			#print(hapGrids[n,])
		}
		hap1 = hapGrids[n, 1:uniBkCt] 
		m1 = filterHaps2SHaps(hap1, HRCB.hapCt2) 
		
		hap2 = hapGrids[n, (uniBkCt+1):(2*uniBkCt)]
		m2 = filterHaps2SHaps(hap2, HRCB.hapCt2)  
		
		hapPairs = qExpandTable(listOfFactor =list(m1, m2), removedRowIdx=NULL, re.row=FALSE)
		hapPairs = t(apply(hapPairs, 1, FUN=range))
		superDipIdx = apply(hapPairs, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		superDipIdx = unique(unlist(superDipIdx))
		
		supIdx = c(supIdx, superDipIdx)
		supIdx = unique(supIdx)
		
		if(length(supIdx)>10^10) stop( qp(fN,": Too many matching super-haplotypes, exceeding 10^10."))
		#if(ifD) print(supIdx)
	}
	
	supDipAll = NULL
	if( HRCB.hapCtProd < 250 ){
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="radix")]
	}else{
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="quick", na.last=NA)]
	}
	
	
	#print(paste(fN, "QingMark2::sorted superDip:"))
	#print(length(supDipAll))
	## parse out the four stratum
	if (length(supDipAll)==0){
		warning(paste(fN, ":no matching allele on the HRCB"))
	}
	
	#print(supDipAll)
	## FIX LATER!!!FIX : can further reduce the workload
	HRCBGrp = t(sapply(supDipAll, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=TRUE))
	
	
	#print(HRCBGrp)
	# print(str(HRCBGrp))
	
	HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=FALSE]
	
	HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
	#save(HRCBGrp.A, file="HRCBGrp.A.RData")
	#save(HRCBGrp.BIdx, file="HRCBGrp.BIdx.RData")
	
	return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
	
}

bkMap.ESp.apply1Rule.stepBy <-
function(bkMap, signalRule){
	ifD = FALSE
	fN = "bkMap.ESp.apply1Rule"
	if(ifD) print(paste(fN, "begin::"))
	
	signal = signalRule$signal
	
	total.bkct = length( bkMap$snpCt )
	
	# construct the beginning/ending index of snp in original map
	idx.be = matrix(NA, nrow=total.bkct, ncol=2)
	idx.be[,2]=cumsum(bkMap$snpCt)
	idx.be[,1]=c(1, cumsum(bkMap$snpCt)+1)[1:total.bkct]
	
	
	causalBkUniqOrderedIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=TRUE)
	## find out # of haplotype in each HRCB
	HRCB.hapCt = bkMap$bkLens [causalBkUniqOrderedIdx]
	HRCB.hapCtProd = cumprod(HRCB.hapCt)[length(causalBkUniqOrderedIdx)]
	
	causalBkIdx.ruleOrder = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=FALSE)
	
	#if(length(causalBkUniqOrderedIdx)==1){
	
	if(length(causalBkIdx.ruleOrder)==1){
		
	}else{
		## need to update the rule a little
		if(ifD) print( binaTree.toStr(signal))
		signal.n = binaTree.shuffleNode(signal)    
		if(ifD) print( binaTree.toStr(signal.n))
		signal.f = binaTree.merge(signal.n)
		if(ifD) print( binaTree.toStr(signal.f))
		
	}
	
	## have to find the matching with the HRCB only bkMap
	#causalBkIdx.ruleOrder = match(causalBkIdx.ruleOrder, causalBkIdx)
	
	elmMList = list()
	# for each block, need to find the hap pairs that are associated with higher risk
	
	for ( j in 1:length(causalBkIdx.ruleOrder)){
		bkIdx = causalBkIdx.ruleOrder[ j ]
		if(ifD) print(paste("bkIdx inside the bkMap =", bkIdx ))
		
		##  20Dec08Change!!!: straighten coding: 1, 2, 3 for 1-digitCoding, and 1, 2 for 2-digit
		
		genoMap =  hap2GenoBlock(hapExp = as.character(bkMap$bks[[bkIdx]][, bkMap$expCol]),
				hapProb = bkMap$bks[[bkIdx]][, bkMap$probCol], snpCt = bkMap$snpCt[bkIdx],
				alleleCode = bkMap$alleleCode )
		if(ifD) print(genoMap)
		item = idx.be[bkIdx,]
		## use recycling rules
		re = paste("v", rep(item[1]:item[2], each=2), sep="")
		newData.col = paste(re,  c("a", "b"), sep="")
		
		## singleSNPRule follow the original order in the genome
		## HARD CODE!!!HARD CODE: genoMa from hap2GenoBlock keep digits  1, 2 for 2-digit, need to change to 0/1 binary coding
		## 20Dec08Change!!! to use 2- for 0/1 binay coding, so the first alleleCode (changed to 1) correspond to the minor allele
		## and the second  (changed to 2) for major allele. Just need to coordinate with signal Rule.
		
		varIdx = match (signalRule$signal$elm.vlist[[j]], newData.col)
		treePred = (2-genoMap$genoMa)[,varIdx]
		## if the variable is not proceed by "not" operator
		## CHANGED on 8JAN09, from ==0 to !=-1. -1 means not
		if (signalRule$signal$elm.vMa[j,3]!=-1){
			nct = sum(treePred==1)
			if(ifD) print("not proceed by not")
			pairHRCB = genoMap$tb[  treePred==1, c(1,2), drop=FALSE ]
		}else{
			## if the variable is proceed by "not" operator
			nct = sum(treePred==0)
			if(ifD) print("proceed by not")
			pairHRCB = genoMap$tb[  treePred==0, c(1,2), drop=FALSE ]
		}
		
		#print(paste("matching hap pairs in block ", bkIdx))
		#print(pairHRCB)
		
		
		elmMList=c(elmMList, list(pairHRCB))
	} # for ( j in 1:length(causalBkIdx.new)){
	#print("This is elmMList")
	#print(elmMList)
	#print(causalBkIdx.ruleOrder)
	#print("qing mark")
	
	## simplied if only one HRCB
	
	if (length(causalBkIdx.ruleOrder)==1){
		
		superDipIdx = apply(pairHRCB, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		
		superDipIdx = sort(superDipIdx)
		
		## FIX LATER!!!FIX : can further reduce the workload
		HRCBGrp = t(sapply(superDipIdx, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=TRUE))
		
		# print(str(HRCBGrp))
		
		HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=FALSE]
		
		HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
		
		#print(HRCBGrp.A)
		#print(HRCBGrp.BIdx)
		return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
		
	}
	
	hapGrids = binaTree.patternForm(binaTree=signal.f, elmMList=elmMList, bkIdx.RuleOrder=causalBkIdx.ruleOrder)
	if(ifD) print(str(hapGrids))
	##write.csv(hapGrids, file="hapGrids.csv") 
	
	## check whether no matching BK
	apply(hapGrids, 2, FUN=function(coo){
				qing.check = mean(is.na(coo))
				if(qing.check==1) stop(qp(fN, ": no matching pattern for one block"))
			})
	
	
	## need to find super-hap index for the matching patterns
	uniBkCt = ncol(hapGrids) / 2
	
	HRCB.hapCt2 = HRCB.hapCt
	
	if(ifD) print(hapGrids)
	supIdx = NULL
	
	for(n in 1:nrow(hapGrids)){
		
		if(ifD) {
			print(paste("proc row for pattern: row=", nrow(hapGrids)))
			#print(n)
			#print(hapGrids[n,])
		}
		if( n>240) {
			#print(hapGrids[n,])
		}
		hap1 = hapGrids[n, 1:uniBkCt] 
		m1 = filterHaps2SHaps(hap1, HRCB.hapCt2) 
		
		hap2 = hapGrids[n, (uniBkCt+1):(2*uniBkCt)]
		m2 = filterHaps2SHaps(hap2, HRCB.hapCt2)  
		
		hapPairs = qExpandTable(listOfFactor =list(m1, m2), removedRowIdx=NULL, re.row=FALSE)
		hapPairs = t(apply(hapPairs, 1, FUN=range))
		superDipIdx = apply(hapPairs, 1, FUN=util.it.triMatch, len=HRCB.hapCtProd)
		superDipIdx = unique(unlist(superDipIdx))
		
		supIdx = c(supIdx, superDipIdx)
		supIdx = unique(supIdx)
		
		if(length(supIdx)>10^10) stop( qp(fN,": Too many matching super-haplotypes, exceeding 10^10."))
		#if(ifD) print(supIdx)
	}
	
	supDipAll = NULL
	if( HRCB.hapCtProd < 250 ){
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="radix")]
	}else{
		supDipAll = supIdx[sort.list(as.integer(supIdx), method="quick", na.last=NA)]
	}
	
	
	#print(paste(fN, "QingMark2::sorted superDip:"))
	#print(length(supDipAll))
	## parse out the four stratum
	if (length(supDipAll)==0){
		warning(paste(fN, ":no matching allele on the HRCB"))
	}
	
	#print(supDipAll)
	## FIX LATER!!!FIX : can further reduce the workload
	HRCBGrp = t(sapply(supDipAll, FUN=util.it.triMatch2, len=HRCB.hapCtProd, re.homo=TRUE))
	
	
	#print(HRCBGrp)
	# print(str(HRCBGrp))
	
	HRCBGrp.A = HRCBGrp[HRCBGrp[,3]==1, c(1,2), drop=FALSE]
	
	HRCBGrp.BIdx = HRCBGrp[HRCBGrp[,3]==0,1]
	#save(HRCBGrp.A, file="HRCBGrp.A.RData")
	#save(HRCBGrp.BIdx, file="HRCBGrp.BIdx.RData")
	
	return(list(A=HRCBGrp.A, B=HRCBGrp.BIdx))
	
}

bkMap.findHRCBIdx <-
function(bkMap, snpIdx, re.keys=TRUE, unique=TRUE, complement=FALSE){
  if(!is.integer(snpIdx)) stop()
  # find out the bk bin the snp fall into by comparing the idx with the ending idx of block
  endIdx = cumsum(bkMap$snpCt)
  idx.bk.belong = qing.cut(val=snpIdx, cutPt=endIdx, cutPt.ordered = TRUE, right.include=TRUE)

  if ((max(idx.bk.belong<0)==1) | (max(idx.bk.belong>length(endIdx))==1)) stop("One of the SNP index in sigStr is out of range.")

  if(unique){
    idx.bk.belong = sort(unique(idx.bk.belong)) 
  }
  if(complement){
    total = 1:(length(bkMap$keys))
    idx.bk.belong.unique =  sort(unique(idx.bk.belong)) 
    idx.bk.belong = total[ -idx.bk.belong.unique ]
  }
  
  if(re.keys){
    idx.bk.belong = bkMap$keys[idx.bk.belong]
  }
  return(idx.bk.belong)
}

bkMap.genoFreq <-
function(bkMap, keys){

  # only one key
  idx = match(keys, bkMap$keys)

  genoFreq = NULL

  for( i in idx){
    if( is.na(i) ) stop()
    if(i>1){
      snpBIdx = sum(bkMap$snpCt[ 1:(i-1) ])
    }else{
      snpBIdx = 0
    }
    bks = bkMap$bks[[ i ]]
    genoFreq.byBk = hap.2geno( as.character(bks[, bkMap$expCol]), bks[,bkMap$probCol], snpBIdx)
    genoFreq = rbind(genoFreq, genoFreq.byBk)
  }
  return(genoFreq)

}

bkMap.HRCB.Esp1Rule.Base <-
function(bkMap, rule, baseName=NULL, dig1Code=0:3){
  ifD = FALSE
  fn = "bkMap.HRCB.Esp1Rule.Base::"
  if(ifD) print(paste(fn, "begin"))

  if (length( rule$slist)>1) stop("This function work only with one SignalRule list")

  
  signalRule = rule$slist[[1]]

  HRCBIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=TRUE)
  superHRCBMap = bkMap.superHRCB(bkMap, uniBkIndexes=HRCBIdx, re.probOnly = TRUE)
  
  # obtain super-hap and their probabilities
  supHapProb = superHRCBMap
  supLen = length(supHapProb)

    
  # if only one term in the rule, it doesn't matter whether the coef is positive or negative
  # if more than one term, need to choose the negative as the reference
  ## 20Dec08Change!!!
  HRCBGrps = bkMap.ESp.apply1Rule(bkMap, signalRule)

  if(ifD){
    print( HRCBGrps )
    #return(NULL)
  }

  ## construct the objects for four sampling stratum
  
  if (nrow(HRCBGrps$A)==0){
    ## if no heter pairs associated with high risk
    #print("check")
    #print(HRCBGrps$A)
  }
  
  HRCBStra.AD = HRCBSpGrp.cons(supHapProb, HRCBGrps$A, type="A")
  #print(HRCBStra.AD)

  HRCBStra.A = HRCBStra.AD$grpA

  HRCBStra.D = HRCBStra.AD$grpD

  #print(supHapProb)
  #print(HRCBGrps$B)
  HRCBStra.B = HRCBSpGrp.cons(supHapProb, HRCBGrps$B, type="B")
  #print(HRCBStra.B)
  
  if(ifD) {
    print(HRCBStra.A)
    print(HRCBStra.D)
    print(HRCBStra.B)
  }
  if(length( HRCBGrps$B )>0){
    tmp.Cidx = (1:supLen)[-HRCBGrps$B]
  }else{
    tmp.Cidx = (1:supLen)
  }
  HRCBStra.C = HRCBSpGrp.cons(supHapProb, tmp.Cidx, type="C")

  if(ifD){
    # check the sum of prob
    print(paste(fn, " check the sum of prob:"))
    print(sum(HRCBStra.A$idProb[,2],
        HRCBStra.B$idProb[,2],
        HRCBStra.C$idProb[,2],
        HRCBStra.D$idProb[,2]))
     print(HRCBStra.A)
     print(HRCBStra.B)
     print(HRCBStra.C)
     print(HRCBStra.D)
  }

  ## restandadize the prob
  straCt = c(HRCBStra.A$ct, HRCBStra.B$ct, HRCBStra.C$ct, HRCBStra.D$ct)
  straCumCt = cumsum(straCt)
  sttProb = c(HRCBStra.A$idProb[,2],
        HRCBStra.B$idProb[,2],
        HRCBStra.C$idProb[,2],
        HRCBStra.D$idProb[,2])
  sttProb = sttProb/sum(sttProb)

  #print(paste(fn, " risk allele freq = ", sum( sttProb [1: straCumCt[2] ])))

  
  if(HRCBStra.A$ct>0) HRCBStra.A$idProb[,2] = sttProb[1:straCumCt[1]]
  if(HRCBStra.B$ct>0) HRCBStra.B$idProb[,2] = sttProb[(straCumCt[1]+1): straCumCt[2] ]
  if(HRCBStra.C$ct>0) HRCBStra.C$idProb[,2] = sttProb[(straCumCt[2]+1): straCumCt[3] ]
  HRCBStra.D$idProb[,2] = sttProb[(straCumCt[3]+1): straCumCt[4] ]

  HRCBStra = list(HRCBStra.A = HRCBStra.A, HRCBStra.B = HRCBStra.B,
                  HRCBStra.C = HRCBStra.C, HRCBStra.D = HRCBStra.D)

  if(!is.null(baseName))  save(HRCBStra, file=paste(baseName, ".RData", sep=""))

  return(HRCBStra)

}

bkMap.HRCB.Esp1Rule.genoSeq <-
function(bkMap, rule, re.probOnly = TRUE){

  ifD = FALSE
  signalRule = rule$slist[[1]]

  HRCBIdx = bkMap.findHRCBIdx(bkMap, as.integer(signalRule$var.idx), re.keys=FALSE, unique=TRUE)

  causalInfo = bkMap.updateByRule(bkMap, rule)
  
  # after process signalRule, 
  keys.new = bkMap$keys[c(causalInfo$causalBkIdx)]
  bkMapS = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=FALSE)

  superHRCBMap = bkMap.superHRCB(bkMap, uniBkIndexes=HRCBIdx, re.probOnly=re.probOnly)

  if(re.probOnly){
    supHapProb = superHRCBMap
    #print(paste("# of prob: ", length( superHRCBMap)))

    return(list(bkMapS=bkMapS, supHapProb=supHapProb))
  }else{
   #print("Re additional stuff")
   #print(superHRCBMap)
   supHapProb = superHRCBMap$prob
   superHapExp = apply(superHRCBMap$linkBkExp, 1, paste, collapse="")

   if(ifD) print( cbind(superHapExp, superHRCBMap$prob, cumsum(superHRCBMap$prob)))
  
   # obtain super-hap and their probabilities

   #print(paste("# of prob: ", length( superHRCBMap$prob)))

   return(list(bkMapS=bkMapS, supHapExp=superHapExp, supHapProb=supHapProb))

  }

}

bkMap.HRCB.famMap <-
function(bkMapS, rule, newColName,  ifS="simuDirectInfo", baseName=NULL){
     ifD = FALSE

     ct.shap = cumprod(bkMapS$bkLens)[length(bkMapS$bkLens)]
     ct.sdip = .5*ct.shap*(1+ct.shap)
     ct.mrow = .5*ct.sdip*(1+ct.sdip)

     if(ct.mrow > 5*10^6) warning(qp("Total number of super-haplotype exceed 5*10^6. It is likely to exceed the memory limit."))

     if(ct.mrow > 5*10^7) stop(qp("Total number of super-haplotype exceed 5*10^7. It exceeds the memory limit."))

     linkedHap = bkMap.superHRCB(bkMapS)
     hapExp = apply(linkedHap$linkBkExp, 1, paste, collapse="")
     prob = linkedHap$prob/sum(linkedHap$prob)
   
     hapCt = length(prob)
     
     ab = hap2GenoBlock(hapExp, prob, snpCt = bkMapS$snpLen, alleleCode = bkMapS$alleleCode)

     risk = HRCB.applyRule(ab$tb$genoExp, rule, newColName)

     ## ruleMet only return the applied result for the last rule in the list
     ## associate the hap pair index with prob is used for the kids
     if(length(rule$slist)==1){
       tb = cbind(ab$tb, ruleMet=risk$treePred, risk.Prob=risk$riskProb)
     }else{
       tb = data.frame(ab$tb, risk.Prob=risk$riskProb)
     }

     if(ifD) print(str(tb))
     if(ifD) print(tb)

     #print("test")

     ## generat kids hap ids pair
     rowCt = nrow(tb)

     if(ifD) print(rowCt)
     matRowCt = (rowCt^2-rowCt)/2+rowCt
     matingRowIdxCorn = util.it.smallLargeIdx(rowCt, keep.same=FALSE)
     if(ifD) print(str(matingRowIdxCorn))
     matingRowIdxDiag = matrix(rep(1:rowCt, each=2), ncol=2, byrow=TRUE)
     if(ifD) print(str(matingRowIdxDiag))

     matingRowIdx = rbind(matingRowIdxCorn, matingRowIdxDiag)
     if(ifD) print(matingRowIdx[1:20,])

     tmpProb = tb$prob
     matingPCor = matrix(tmpProb[matingRowIdxCorn], ncol=2, byrow=FALSE)
     matingPCorn = 2*matingPCor[,1]*matingPCor[,2]
     matingPDiag = tmpProb^2
     matingP = c(matingPCorn, matingPDiag)
     if(ifD) print(matingP[1:20])
     
     ## HARD CODE!!! NEED to be blocked
##****  subs = 1:10
##****  matingRowIdx = matingRowIdx[subs,]
##****  matingP = matingP[subs]

     tmpTb = matrix(unlist(tb[,c(1,2)]), ncol=2, byrow=FALSE)
     if(ifD) print(str(tmpTb))
    
     fa.hapIdx = tmpTb[matingRowIdx[,1], ]
     ma.hapIdx = tmpTb[matingRowIdx[,2], ]

     if(ifD) {
       print("Sample parents hap idxes:")
       print(fa.hapIdx[1:10,])
       print(ma.hapIdx[1:10,])
     }

     fam.map = matrix(NA, ncol=12, nrow=matRowCt )
     fam.map[,c(1,2)]=fa.hapIdx
     fam.map[,c(3,4)]=ma.hapIdx
     fam.map[,5]=pmin( fa.hapIdx[,1], ma.hapIdx[,1])
     fam.map[,6]=pmax( fa.hapIdx[,1], ma.hapIdx[,1])
     fam.map[,7]=pmin( fa.hapIdx[,1], ma.hapIdx[,2])
     fam.map[,8]=pmax( fa.hapIdx[,1], ma.hapIdx[,2])
     fam.map[,9]=pmin( fa.hapIdx[,2], ma.hapIdx[,1])
     fam.map[,10]=pmax( fa.hapIdx[,2], ma.hapIdx[,1])
     fam.map[,11]=pmin( fa.hapIdx[,2], ma.hapIdx[,2])
     fam.map[,12]=pmax( fa.hapIdx[,2], ma.hapIdx[,2])

     if(ifD) print(fam.map[1:10,])

     ## obtain the disease prob given the hap idx
     kids.hapIdx = matrix(t(fam.map[,5:12]), nrow=2, byrow=FALSE, ncol=4*matRowCt)
     if(ifD) print(kids.hapIdx[,1:10])
     kids.hapIdx = t(kids.hapIdx)
     if(ifD) print(kids.hapIdx[1:10,])
     if(ifD) print(str(kids.hapIdx))

### instead of using the inefficient matching, we will use the virtual mapping function 
#     kids.matchedRow = unlist(apply(kids.hapIdx[1:20,], 1, util.vec.matchVecIdx, vec=t( tmpTb ),  vecLen=rowCt*2, #benchLen=2))

     kids.matchedRow = unlist(apply(kids.hapIdx, 1, util.it.triMatch, len=hapCt))

     if(ifD) print(kids.matchedRow[1:20])
     # kids.p = matrix(  rep(.25*matingP, each=4), ncol=4, byrow=TRUE)
     # print( paste("Check:: kids p sum (before standardization)=", sum(kids.p)))
     kids.p = matingP/sum(matingP)
  
     ## check
     #if(ifD)  print( paste("Check:: kids p sum=", sum(kids.p)))
     #fam.map = data.frame(fa.hapIdx, ma.hapIdx, ch1.h, ch2.h, ch3.h, ch4.h, rowProb = matingP, kids.p)
     #if(ifD) print(fam.map[1:10,])
  
     kids.risk = matrix( risk$riskProb [kids.matchedRow], ncol=4, byrow=TRUE)
     kids.matchedRow = matrix(kids.matchedRow, ncol=4, byrow=TRUE)
     if(ifD) {
       print("kids risk")
       print(str(kids.risk))
     }
     kids.risk[,1] = (.25*kids.p)* kids.risk[,1]
     kids.risk[,2] = (.25*kids.p)* kids.risk[,2]
     kids.risk[,3] = (.25*kids.p)* kids.risk[,3]
     kids.risk[,4] = (.25*kids.p)* kids.risk[,4]
  
     ## check
     #print( paste("Check:: kids risk sum=", sum(kids.risk)))
     
     fam.map = cbind(fam.map, matingP, kids.risk)
     colnames(fam.map) = c("f.hap1", "f.hap2", "m.hap1", "m.hap2",
                         "c1.hap1", "c1.hap2", "c2.hap1", "c2.hap2",
                         "c3.hap1", "c3.hap2", "c4.hap1", "c4.hap2",
                         "rowProb", 
                         "c1Risk", "c2Risk", "c3Risk", "c4Risk")
  
     if(ifD) print(round(fam.map[1:10,],5))

     if(!is.null(ifS)) {
       
       #write.csv(tb, file=paste(ifS, "hap2geno.csv", sep=""))
       #write.csv(fam.map,  file=paste(ifS, "hapMating.csv", sep=""))
     }
  
     matingTbInfo = list(matRowCt=matRowCt, kids.risk=kids.risk,
                 matingRowIdx = matingRowIdx,
                 hap2genoMap = ab$tb,
                 snpLen = bkMapS$snpLen,
                 kids.matchedRow = kids.matchedRow)
     #print(str(matingTbInfo))

     if(!is.null(baseName))  save( matingTbInfo, file=paste(baseName, ".RData", sep=""))

     return( matingTbInfo )
}

bkMap.HRCB.LRCB.split <-
function (bkMap, rule){

  causalInfo = bkMap.updateByRule(bkMap=bkMap, rule=rule)
  col.shuffled=causalInfo$new.colname

  nocausalInfo = bkMap.updateByRule(bkMap=bkMap, rule=rule,  complement=TRUE )
  col.shuffled.nocausal =nocausalInfo$new.colname

  
  # after process signalRule, 
  keys.new = bkMap$keys[c(causalInfo$causalBkIdx)]
  
  # one map for associated blocks
  bkMapS  = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=FALSE)


  bkMapNS = bkMap.shuffle(bkMap, bkCt=NULL, keys.new=keys.new, exclude=TRUE)
        
  return(bkMaps = list(bkMapS=bkMapS,  bkMapNS=bkMapNS, newColName=col.shuffled, noCausalColName=col.shuffled.nocausal))
}

bkMap.LRCB.spTrio <-
function(bkMap, caseNo, ifS = NULL, reControl=FALSE){

   ifD = FALSE
   # first sample parents
   f.samples = bkMap.spHap(bkMap, caseNo)
   f.samples2 = bkMap.spHap(bkMap, caseNo)

   m.samples = bkMap.spHap(bkMap, caseNo)
   m.samples2 = bkMap.spHap(bkMap, caseNo)

   #dim(par)
   fa = covDipStr2CodedGeno(f.samples$subjects, f.samples2$subjects, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
   ma = covDipStr2CodedGeno(m.samples$subjects, m.samples2$subjects, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))

   ## shuffle the 
   par = cbind(as.vector(f.samples$subjects), as.vector(f.samples2$subjects),
               as.vector(m.samples$subjects), as.vector(m.samples2$subjects))

   parIn =  cbind(as.vector(f.samples$subjectsIn), as.vector(f.samples2$subjectsIn),
               as.vector(m.samples$subjectsIn), as.vector(m.samples2$subjectsIn))

   #print(paste("!reControl: par dim:", paste(dim(par), collapse=" by ", sep="")))

   
   ## shuffle the random kid
   ## get index, !!!take the same happair combination for all blocks and all trio
   baseIdx =  matrix(c(1,3,2,3,1,4,2,4), nrow=2, byrow=FALSE)

   randomRowIdx = sample(1:4, size=1, replace=FALSE)
   newchild = baseIdx[,randomRowIdx]

   chExp = par[, newchild]
   chIn = parIn[, newchild]

   #print(paste("!reControl: chExp dim:", paste(dim(chExp), collapse=" by ", sep="")))

   #print(dim(chExp))
   #print(dim(chIn))
   #print(dim(newchild))

   child1 = matrix(chExp[,1], nrow=caseNo, byrow=FALSE)
   child2 = matrix(chExp[,2], nrow=caseNo, byrow=FALSE)

   affChild = covDipStr2CodedGeno(child1, child2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))


   trioSetts = rbind(fa, ma, affChild)
   #print(paste("!reControl: trioSetts dim:", paste(dim(trioSetts), collapse=" by ", sep="")))

   if(!reControl){     
     rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=FALSE))
     geno.FMCMa = trioSetts[rowNewSeq,]
   }else{

     othRowIdx = (1:4)[ -randomRowIdx ]
     
     exHapIdx = baseIdx[,othRowIdx]

     tt = FALSE
     if(tt) print(paste("othRowIdx=", paste(othRowIdx, collapse=";")))
     if(tt) print(paste("expHapIdx=", paste(exHapIdx, collapse=";")))

     for( i in 1:3){
          childExpIdx = exHapIdx[,i]
          #print(childExpIdx)
          othChildExp = par[,childExpIdx]
          #print(paste("!reControl: othChildExp dim:", paste(dim(othChildExp), collapse=" by ", sep="")))
          childOth1 = matrix(othChildExp[,1], nrow=caseNo, byrow=FALSE)
          childOth2 = matrix(othChildExp[,2], nrow=caseNo, byrow=FALSE)

          affChild =
            covDipStr2CodedGeno(childOth1, childOth2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
          trioSetts = rbind(trioSetts, affChild)
     }
     #print(paste("!reControl: trioSetts dim:", paste(dim(trioSetts), collapse=" by ", sep="")))
     
     rowNewSeq = t(matrix(1:(caseNo*6), ncol=6, nrow=caseNo, byrow=FALSE))
     #print(trioSetts)
     geno.FMCMa = trioSetts[rowNewSeq,]
     

   }

     ## rearrange the data, so the three subject in a family goes together.
   if(!is.null(ifS)){
       ## retain trio Index
       childIn1 = matrix(chIn[,1], nrow=caseNo, byrow=FALSE)
       childIn2 = matrix(chIn[,2], nrow=caseNo, byrow=FALSE)
       faIn1 = matrix(parIn[,1],  nrow=caseNo, byrow=FALSE)
       faIn2 = matrix(parIn[,2],  nrow=caseNo, byrow=FALSE)
       maIn1 = matrix(parIn[,3],  nrow=caseNo, byrow=FALSE)
       maIn2 = matrix(parIn[,4],  nrow=caseNo, byrow=FALSE)
  
       childIn = util.matrix.col.shuffle2(childIn1, childIn2)
       faIn = util.matrix.col.shuffle2(faIn1, faIn2)
       maIn = util.matrix.col.shuffle2(maIn1, maIn2)
  
       allIn = rbind(faIn, maIn, childIn)

       rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=FALSE))
       allIn = allIn[rowNewSeq,]

       write.table(allIn,  file=paste(ifS, "supHap.csv", sep=""), col.names=FALSE, row.names=FALSE, sep=",")         
   }
   return(geno.FMCMa)
   
#      ## not implemented
#      ## shuffle the random kid
#      child = apply(par, 1, FUN=function(row){
#        idx = matrix(c(1,3,2,3,1,4,2,4), nrow=2, byrow=FALSE)
#        randomChildIdx = sample(1:4, size=4, replace=FALSE)
#        idx = idx[,randomChildIdx]
#   
#        newchild = row[idx]
#      })
#      child1.1 = matrix(child[1,], nrow=caseNo, byrow=FALSE)
#      child1.2 = matrix(child[2,], nrow=caseNo, byrow=FALSE)
# 
#      child2.1 = matrix(child[3,], nrow=caseNo, byrow=FALSE)
#      child2.2 = matrix(child[4,], nrow=caseNo, byrow=FALSE)
# 
#      child3.1 = matrix(child[5,], nrow=caseNo, byrow=FALSE)
#      child3.2 = matrix(child[6,], nrow=caseNo, byrow=FALSE)
# 
#      child4.1 = matrix(child[7,], nrow=caseNo, byrow=FALSE)
#      child4.2 = matrix(child[8,], nrow=caseNo, byrow=FALSE)
# 
#      affChild1 = covDipStr2CodedGeno(child1.1, child1.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild2 = covDipStr2CodedGeno(child2.1, child2.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild3 = covDipStr2CodedGeno(child3.1, child3.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#      affChild4 = covDipStr2CodedGeno(child4.1, child4.2, subjectCt=caseNo, snpLen=bkMap$snpLen, snpCoding=0:3, snpBase=c(0, bkMap$alleleCode))
#   
#      allChild = rbind(affChild1, affChild2, affChild3, affChild4)
# 
#      # child is SNPct*4 by caseNo matrix
#      newSeq = matrix(1:(4*caseNo), ncol=caseNo, byrow=TRUE)
#      newChild = allChild[newSeq,]
# 
#      trioSetts = rbind(fa, ma, affChild1)
#      rowNewSeq = t(matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=FALSE))
# 
#      geno.FMCMa = trioSetts[rowNewSeq,]
#      return(trio=list(geno.FMCMa=geno.FMCMa, cc = newChild))

}

bkMap.shuffle <-
function(bkMap, bkCt=1, keys.new=NULL, exclude=FALSE){
  keys.ori = bkMap$keys
  if(is.null(keys.new)){
    keys.new =  resample(keys.ori, size=bkCt, replace=FALSE)
  }

  # order are index in the keys.ori
  if(!exclude){
    keys.neworder = match(keys.new, keys.ori)
    leftkeys.neworder = (1:length(keys.ori)) [-keys.neworder]

    ## CHANGE. do not think keep the left over blocks is ever be used.
    #if(!del.left){
    #  keys.neworder = c(keys.neworder, leftkeys.neworder)
    #}else{
      keys.neworder = keys.neworder
    #}
  }else{
    # if want to exclude the keys.new
    keys.neworder = match(keys.new, keys.ori)
    keptkeys.neworder = (1:length(keys.ori)) [-keys.neworder]

    #if(!del.left){
    #  keys.neworder = c(keptkeys.neworder, keys.neworder)
    #}else{
      keys.neworder = keptkeys.neworder
    #}
  }

  ##print(keys.neworder)
  newBkMap = bkMap
  newBkMap$bks = bkMap$bks[keys.neworder]
  newBkMap$bkLens = bkMap$bkLens[keys.neworder]
  newBkMap$keys = keys.ori[keys.neworder]
  newBkMap$snpCt = bkMap$snpCt[keys.neworder]
  newBkMap$snpLen = sum(newBkMap$snpCt)
  newBkMap$alleleCode = bkMap$alleleCode

  if(exclude){
    newBkMap$selBkCt = length(keys.ori) - length(keys.new)
  }else{
    newBkMap$selBkCt = length(keys.new)
  }
  return(newBkMap)
  
}

bkMap.spGenoSeq <-
function(bkMap, subjectCt){

   # first set of hap
   samples = bkMap.spHap(bkMap, subjectCt)
   samples2 = bkMap.spHap(bkMap, subjectCt)
   if(F) print(cbind(samples[[1]], samples2[[1]]))

   # change hap to genotype
   samplesBina = hapBk2AlleleSeq(subjects=samples$subjects, subjectCt, snpLen=bkMap$snpLen)
   samplesBina2 = hapBk2AlleleSeq(subjects=samples2$subjects, subjectCt, snpLen=bkMap$snpLen)
  
   ## construct the genotypic data
   binaNew1 = pmin(samplesBina, samplesBina2)
   binaNew2 = pmax(samplesBina, samplesBina2)
   bina = util.matrix.col.shuffle2(binaNew1, binaNew2)

   return(bina)
}

bkMap.spHap <-
function(bkMap, subjectCt){

  len = length(bkMap$bkLens)

  re = NULL
  people = NULL
  people.in = NULL
  for( i.bk in 1:len ){
    curIn = sample(bkMap$bkLens[i.bk], size = subjectCt, replace = TRUE, prob = bkMap$bks[[i.bk]][,bkMap$probCol])

    curExp = lapply(curIn, FUN=function(item, bkMap, i.bk){
                              as.character(bkMap$bks[[i.bk]][item,bkMap$expCol])
                            }, bkMap = bkMap, i.bk=i.bk)
    people = c(people, list(unlist(curExp)))
    people.in = c(people.in, list(curIn))
  }
  subjects = matrix(unlist(people), ncol=length(people), byrow=FALSE)
  subjects.in = matrix(unlist(people.in), ncol=length(people.in), byrow=FALSE)

  re = list(subjects=subjects, subjectsIn=subjects.in)
  return(re)
}

bkMap.superHRCB <-
function(bkMap, uniBkIndexes=c(1:length(bkMap$keys)), re.probOnly = FALSE){

  ifD = FALSE
  ## find every block expression for every block
  linked = lapply(uniBkIndexes, FUN=function(index, bkMap){
                                            bkMap$bks[[index]]}, bkMap=bkMap)

  hapLen = bkMap$snpCt
  
  len = length(uniBkIndexes)
  if(ifD) print(paste("block len=", len, sep=""))
  
  gridProb = lapply(linked, FUN=function(bk, probCol){
                                 re = bk[, probCol]}, probCol=bkMap$probCol)


    ## only calculate the probability
     gridProb = lapply(linked, FUN=function(bk, probCol){
                                   re = bk[, probCol]}, probCol=bkMap$probCol)
    if(length(uniBkIndexes)==1){
       prob = as.matrix(gridProb[[1]])
    }else{
      prob = qExpandTable(listOfFactor = gridProb)
      prob = apply(prob, 1, prod)
    }
    prob = prob/sum(prob)


  if(re.probOnly){
    return(prob)
  }
  
  gridBase = lapply(linked, FUN=function(bk, expCol){
                                 re = as.character(bk[, expCol])}, expCol=bkMap$expCol)
  gridBaseIndex = lapply(bkMap$bkLens[uniBkIndexes], FUN=function(bkLen){
                                 re = 1:bkLen})
  
  if(length(uniBkIndexes)==1){
     mas = as.matrix(gridBase[[1]])
     mashIn = as.matrix(gridBaseIndex[[1]])
  }else{
    mas = qExpandTable(listOfFactor =  gridBase )
    mashIn = qExpandTable(listOfFactor =  gridBaseIndex )

  }

  if(ifD) {
    print(paste("mas index dim=", ""))
    print(str(mas))
    print(mas[1:2,])
    print(mashIn[1:2,])
    print(prob[1:2])
  }
  #updateBkMap = specifyBks(oriData, uniBkIndexes, bkMap, keyCol, hapLenCol)
  #print(mas)

  re = list(linkBkExp=mas, linkBkIn = mashIn, prob=prob)

  ## print(mashIn)
  return(re)  
}

bkMap.updateByRule <-
function(bkMap, rule, complement=FALSE, is.hapGenoMap = FALSE){

  # shuffle the genotype blocks, but keep the columns
  # so need to shuffle the column names

  if(is.hapGenoMap) stop("Method not implemented for hapGenoMap listing yet!")

  # based on the snp idx, map them with block idx
  if(!is.hapGenoMap){
    causalBkIdx.new = bkMap.findHRCBIdx(bkMap, as.integer(rule$snpIdx), re.keys=FALSE, unique=TRUE, complement=complement)
    total.bkct = length(bkMap$keys)
  
    # construct the beginning/ending index of snp in original map
    idx.be = matrix(NA, nrow=total.bkct, ncol=2)
    idx.be[,2]=cumsum(bkMap$snpCt)
    idx.be[,1]=c(1, cumsum(bkMap$snpCt)+1)[1:total.bkct]
  
    id.be.shuffled = idx.be[ c(causalBkIdx.new,  (1:total.bkct)[!is.element(causalBkIdx.new, 1:total.bkct)]),, drop=FALSE ]
  
    newData.col = apply(id.be.shuffled, 1, FUN= function(item){
      ## use recycling rules
      re = paste("v", rep(item[1]:item[2], each=2), sep="")
      re = paste(re,  c("a", "b"), sep="")
      re})

  }else{
    ## not implemented. Because later we added the argu, complement. 
    ## causalBkIdx.new = hapBkGenoMap.findHRCBIdx(allmap=bkMap, snpIdx=as.integer(rule$snpIdx), re.key=FALSE, unique=TRUE)

    ## snp.idxSeq = hapBkGenoMap.findHRCBSnpIdx(allmap=bkMap, snpIdx=as.integer(rule$snpIdx), re.1digit=TRUE)
    ## newData.col = t( paste( "v", rep(snp.idxSeq, each=2), c("a", "b"), sep=""))

  }

  info = list(causalBkIdx = causalBkIdx.new,
            new.colname = unlist(newData.col))
  return(info)
  

}

bkMap.updateHapFreq <-
function(bkMap, bkIndex, expression=NULL, freq){
  bkFrame = bkMap$bks[[bkIndex]]

  if( length(freq)!= nrow(bkFrame)) stop(paste("updateBlockBinaFreq:: freq length doesn't match the ", bkIndex, " block config.", sep=""))

  newfreq = freq/(sum(freq))

  if(!is.null(expression)) {
    if( length(expression)!= nrow(bkFrame)) stop(paste("updateBlockBinaFreq:: expression length doesn't match the ", bkIndex, " block config.", sep=""))

    bkFrame[, bkMap$expCol ] = expression
  }

  bkFrame[, bkMap$probCol ] = newfreq

  bkMap$bks[[bkIndex]] = bkFrame

  return(bkMap)
}

calHapIdx2SHap <-
function(hapIdxes, hapCts){

  # return one index
  bkCt = length(hapCts)

  hapCts= c(1, hapCts[-bkCt])
  cumRows = cumprod(hapCts)
  enuStart = cumRows*(hapIdxes-1)
  idx = cumsum(enuStart)[bkCt] + 1

  return(idx)

}

calHapIdx2SHapSet <-
function( bkIdx, hapIdx, hapCts){

  ifD = FALSE
  
  ## find the number of row for each stratum
  cumRows = c(1, hapCts[1:bkIdx])
  cumRows = cumprod(cumRows)

  it = cumRows[bkIdx]
  matchSet = 1:it
  if(ifD) print(paste("matchSet=", paste(matchSet, collapse=";")))
  set = it*(hapIdx-1)+matchSet
  if(ifD) print(paste("set=", paste(set, collapse=";")))
  
  strataOffset = seq.int(from=0, to=cumprod(hapCts)[length(hapCts)]-1, by = it*hapCts[bkIdx])
  if(ifD) print(paste("strataOffset=", paste(strataOffset, collapse=";")))
  set = rep(strataOffset, each=it) + set;
  return(set)

}

checkMendelianError <-
function(codedSNPTrio, snpCoding=c(0,1,2,3)){

  # is the child homo with 1?
  if(codedSNPTrio[3]==snpCoding[2]){
    onePHaveNone = FALSE
    if(codedSNPTrio[1]==snpCoding[3]) onePHaveNone = TRUE
    othPHaveNone = FALSE
    if(codedSNPTrio[2]==snpCoding[3]) othPHaveNone = TRUE
    if(onePHaveNone | othPHaveNone)
      stop("Medelian error for homozygous (1) child with at least one parent homozygous (2)")
  # is the child homo with 2?
  }else if(codedSNPTrio[3]==snpCoding[3]){
    onePHaveNone = FALSE
    if(codedSNPTrio[1]==snpCoding[2]) onePHaveNone = TRUE
    othPHaveNone = FALSE
    if(codedSNPTrio[2]==snpCoding[2]) othPHaveNone = TRUE
    if(onePHaveNone | othPHaveNone)
      stop("Medelian error for homozygous (2) child with at least one parent homozygous (1)")
  # is the child hetero
  }else if(codedSNPTrio[3]==snpCoding[4]){

    if(codedSNPTrio[1]==snpCoding[2]){
      if(codedSNPTrio[2]==snpCoding[2]){
        stop("Medelian error for heterozygous child with two parents homozygous (1)") 
      }
    }else if(codedSNPTrio[1]==snpCoding[3]){
      if(codedSNPTrio[2]==snpCoding[3]){
        stop("Medelian error for heterozygous child with two parents homozygous (2)") 
      }
    }
  # is the child missing
  }

  return(NULL)
}

covDipBinaMa2CodedGeno <-
function(dip1, dip2, subjectCt, snpCoding=c(0,1,2,3), snpBase=c(0,1,2)){

  dipSum = dip1+dip2

  ##  print(dipSum)
  a1 = sum(snpBase[2:3] * c(2,0)) # 11-> to 1
  a2 = sum(snpBase[2:3] * c(1,1)) # 12-> to 3
  a3 = sum(snpBase[2:3] * c(0,2)) # 22-> to 2

  snp1d.f = factor(dipSum, levels=c(a1, a3, a2), labels = snpCoding[2:4])
  snp1d.f = as.character(snp1d.f)

  if(is.null(dim(dipSum))){
    ncol = length(dipSum)
    nrow = 1
  }else{
    ncol = dim(dipSum)[2]
    nrow = dim(dipSum)[1]
  }

  codedGeno = matrix(as.integer(snp1d.f), nrow=subjectCt, ncol=ncol)

  
##   codedGeno = matrix(NA, nrow=subjectCt, ncol=ncol)
##   for( i in 1:dim(dipSum)[1]){
##     for( j in 1:dim(dipSum)[2]){
##        if (dipSum[i,j] == 2*snpBase[1]) codedGeno[i,j] = snpCoding[1]
##        if (dipSum[i,j] == 2*snpBase[2]) codedGeno[i,j] = snpCoding[2]
##        if (dipSum[i,j] == 2*snpBase[3]) codedGeno[i,j] = snpCoding[3]
##        if (dipSum[i,j] == snpBase[2]+snpBase[3]) codedGeno[i,j] = snpCoding[4]
##     }
##   }
  return(codedGeno)
}

covDipStr2CodedGeno <-
function(dipStr1, dipStr2, subjectCt, snpLen, snpCoding=c(0,1,2,3), snpBase=c(0,1,2)){
  dip1 =  hapBk2AlleleSeq(dipStr1, subjectCt, snpLen, markdownOne = FALSE)
  dip2 =  hapBk2AlleleSeq(dipStr2, subjectCt, snpLen, markdownOne = FALSE)
  re = covDipBinaMa2CodedGeno(dip1, dip2, subjectCt, snpCoding=snpCoding, snpBase=snpBase)
  return(re)
}

dfToGenoMap <-
function(df, dataHeaders=NULL, genotype=c("11", "12", "22"), snpBase=1){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = df
	}else{
		## assume that column names of the data is passed from outside
		data = df
		colnames(data) = dataHeaders
	}
	
	m = match(c("prekey", "seq", "freq", "genotype"), colnames(data), 0)
	
	## assume except for markers_b and marekers_e, other variable should be presented
	if(min(m[2:3])==0) stop("One or more required variable(s) missing")
	
	
	if( !is.element("prekey", colnames(data)) ){
		## assume the seq of homo, hetero, homo
		data$prekey = I(rep("prekey", times=nrow(data)))
	}
	if( !is.element("genotype", colnames(data)) ){
		## assume the seq of homo, hetero, homo
		data$genotype = I(rep(genotype, times=nrow(data)/3))
	}
	mheader=c("prekey", "seq", "freq", "genotype")
	
	m = match(c("prekey", "seq", "freq", "genotype"), colnames(data), 0)
	## cast the data into the right format
	for ( i in 1:4 ){
		if(!is.na(m[i])){
			if (i==1){
				if(class(data[,m[i]])=="factor") data[,i]=as.character(data[,m[i]])
			}
			if (i==4){
				if(class(data[,m[i]])=="factor") data[,i]=as.numeric(as.character(data[,m[i]]))
			}
			if( i==2 | i==3){
				if(class(data[,m[i]])=="factor") data[,i]=as.numeric(as.character(data[,m[i]]))
			}
		}
	}
	
	chLastSNPIndex = util.sql.groupby(data[,m], groupCols=1,  varCol=2, type=c("max"))
	
	genoMap = list(df = data[,m], chLastSNPIndex=chLastSNPIndex, chCol=1, genomeSeqCol=2, probCol=3, expCol=4, snpBase=snpBase)
	return(genoMap)
	
}

dfToHapBkMap <-
function(data, keyCol=NULL, chCol, blockCol, expCol, probCol, hapLenCol, beginCol=NULL, endCol=NULL, snpBase=1, re.bf = TRUE, re.javaGUI = TRUE){
	
	if(is.null(keyCol)){
		data = cbind(data, key=paste(data[,chCol], data[,blockCol], sep="-"))
		keyCol = ncol(data)
	}
	
	allKeys = unique(data[,keyCol])
	keyCt = length(allKeys)
	
	bks = NULL
	snpLen = 0
	bkLens = NULL
	bkSnpLens = NULL
	markers_b = NULL
	markers_e = NULL
	
	df = data[,c( keyCol, chCol, blockCol, expCol, probCol, hapLenCol, beginCol, endCol)]
	
	if(!is.null(beginCol)){
		for( i in 1:keyCt){
			r = df[data[,keyCol]==allKeys[i], ]
			bks = c(bks, list(r))
			bkLens = c(bkLens, dim(r)[1])
			bkSnpLens = c(bkSnpLens, r[1, 6])
			snpLen = snpLen + r[1,6]
			markers_b = c(markers_b, r[1, 7])
			markers_e = c(markers_e, r[1, 8])
		}
	}else{
		for( i in 1:keyCt){
			r = df[data[,keyCol]==allKeys[i], ]
			bks = c(bks, list(r))
			bkLens = c(bkLens, dim(r)[1])
			bkSnpLens = c(bkSnpLens, r[1, 6])
			snpLen = snpLen + r[1, 6]
		}
	}
	
	
	dfStr = cbind(as.character(data[,keyCol]),
			as.character(data[,expCol]),
			as.character(data[,probCol]))
	
	dfMaster = cbind(as.character(allKeys),
			as.character(bkLens),
			as.character(bkSnpLens),
			
			as.character(markers_b),
			as.character(markers_e))
	
	hapBkMap = NULL
	if(re.bf){
		hapBkMap = c(hapBkMap, list(df=df))
	}
	
	if(re.javaGUI){
		hapBkMap = c(hapBkMap, list(dfStr=dfStr,       dfStrDim    = c(nrow(dfStr), ncol(dfStr)),
						dfMaster=dfMaster, dfMasterDim = c(nrow(dfMaster), ncol(dfMaster)),
						bkCumIndex=cumsum(bkLens)))
	}
	
	hapBkMap = c(hapBkMap, list(
					bks = bks,
					bkLens = bkLens,
					keys = as.character(allKeys),
					bkSnpLens = bkSnpLens, 
					snpBase =snpBase, snpLen = snpLen, keyCol=1, chCol=2, blockCol=3,
					expCol = 4, probCol = 5, hapLenCol=6))
	
	if(!is.null(beginCol)){
		hapBkMap = c(hapBkMap, list(markers_b = markers_b, markers_e = markers_e,
						beginCol=7, endCol=8))
	}
	
	return(hapBkMap)
}

ESp.impu1Par <-
function(othParPairs, childPairs,  semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, reType=FALSE, job=1 ){
  ifD = FALSE
  fStr = "[ESp.impu1Par]"
  if(ifD) {
    print(fStr)
    print(semiMapFrame)
    print(paste("snpLen=", snpLen))
  }
  
  if(is.null(dim(othParPairs))) othParPairs = matrix(othParPairs, ncol=2)
  if(is.null(dim(childPairs))) childPairs = matrix(childPairs, ncol=2)
  
  allParIdx= unique(as.vector(othParPairs))

  if(ifD){
    print(othParPairs)
    print(childPairs)

  }
  
  
  ## generate all combination of par and child pairs.
  prob.ps = getHapProb2(selIdx=allParIdx,
    semiMapFrame,
    resiProbCol=resiProbCol,
    augIdxCol=augIdxCol,
    probCol= probCol,
    snpLen, restandard=TRUE)
  if(ifD) {
    print(prob.ps)
  }
  
  n.par = length(othParPairs)/2
  n.ch = length(childPairs)/2

  sp.ma = NULL
  for( i in 1:n.par){
    par = othParPairs[i,]
    
    for( j in 1:n.ch){
      ch = childPairs[j,]
      
      ## find whether the pair match on at lease one hap
      if(sum(is.na(match(par, ch)))==2) {
        if(ifD) print( paste("Ignor par=", paste(par, collapse=";"), " and ch=", paste(ch, collapse=";") ))

      }else{
          sp.prob=NULL
          type = 0
          
          ##1) B/E/C vs. A/D
          if (par[1]!=par[2]){
            ##2) B/E/C, B/E vs. C
            if(ch[1]!=ch[2]){
              
              ch = sort(ch)
              par= sort(par)
              ##3) B/E, B vs. E
              if (sum( ch==par )==2){
                if(ifD) print("B")
                # B
                sp.prob = ESp.impu1Par.B(hap=ch, prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
                type=2

                ## additional col, no use
                add.col = rep(NA, 5)
                
              }else{
                if(ifD) print("E")
                # E: need to keep the common one at the front
                match.t = is.na(match(ch, par))
                if(match.t[1]){ ch = c(ch[2], ch[1])}
                match.t = is.na(match(par, ch))
                if(match.t[1]){ par = c(par[2], par[1]) }             
                
                sp.prob = ESp.impu1Par.E(hap=ch[1], prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
                type=5
                add.col = rep(par[2], 2)
                if(ifD) print(sp.prob)
                
              } ##B vs E::if (sum( ch==par )==2){
              
            }else{ ##2) B/E/C, B/E vs. C
              # C: need to keep the common one at front
              if(ifD) print("C")
              match.t = is.na(match(par, ch))
              if(match.t[1]){par = c(par[2], par[1]) }

              #print(par)
              #print( prob.ps[match(par, allParIdx)] )

              sp.prob = ESp.impu1Par.C(hap=par, prob.p=prob.ps[match(par, allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              # print(sp.prob)
              type=3
              add.col = rep(par[2], 3)

              
            } ## B/E vs C::if(ch[1]!=ch[2]){
            
          }else{ ## B/E/C vs D/A if (par[1]!=par[2]){
            ##2) D/A, D vs. A
            if(ch[1]!=ch[2]){
              # D need to keep the common one at front
              if(ifD) print("D")
              match.t = is.na(match(ch, par))
              if(match.t[1]){
                  ch = c(ch[2], ch[1])              
              }
              sp.prob = ESp.impu1Par.D(hap=ch[2], prob.p=prob.ps[match(par[1], allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              type=4
              
              add.col = rep(NA, 2)
            }else{
              # A
              if(ifD) print("A")
              sp.prob = ESp.impu1Par.A(hap=ch[1], prob.p=prob.ps[match(par[1], allParIdx)], semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
              type=1
              add.col = rep(NA, 2)
              
            } ## D vs. A::if(ch[1]!=ch[2]){
    
          } ## B/E/C vs D/A if (par[1]!=par[2]){
          ch.colA = matrix(rep(ch, length(sp.prob)), ncol=2, byrow=TRUE)
          
          rows = cbind(rep(type, length(sp.prob)), 1:length(sp.prob), sp.prob, add.col, ch.colA)
          if(ifD) print(rows)
          
          sp.ma = rbind(sp.ma, rows)
          # print(sp.ma)

        } # if(sum(is.na(match(par, ch)))<2) {
    } # for( j in 1:n.ch){
  } # for( i in 1:n.par){

  colnames(sp.ma)=c("type", "straSeq", "prob", "add", "ch1", "ch2")


  hap6idx = matrix(NA, nrow=job, ncol=7)
  
  for(ss in 1:job){
    # now sample it
    row.sp = sample(1:nrow(sp.ma), size=1, prob=sp.ma[,3])
  
    type.sp = sp.ma[row.sp,1]
    straSeq = sp.ma[row.sp,2]
    
    if(ifD) {
      print("Sampling")
      print(sp.ma)
    }
    if(type.sp==1){
      sp6 = ESp.impu1Par.A.sp(hap=sp.ma[row.sp, 5],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==2){
      sp6 = ESp.impu1Par.B.sp(hap=sp.ma[row.sp, 5:6],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==3){
      #print("type.sp==3")
      sp6 = ESp.impu1Par.C.sp(hap=sp.ma[row.sp, 5], hap.p=sp.ma[row.sp, 4], straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
      #print(sp6)
    }
    if(type.sp==4){
      sp6 = ESp.impu1Par.D.sp(hap=sp.ma[row.sp, 5:6],  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }
    if(type.sp==5){
      sp6 = ESp.impu1Par.E.sp(hap=sp.ma[row.sp, 5:6], hap.p=sp.ma[row.sp, 4] ,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen )
    }

    hap6idx[ss, ]=c(sp6, type.sp)
  }

  
  if(reType){
    return(hap6idx )
  }else{
    return(hap6idx[,1:6,drop=FALSE])
  }
}

ESp.impu1Par.A <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	prob.c = getHapProb2(selIdx=hap, semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=FALSE)
	stra.1 = (prob.c^2)*(prob.p^2)
	
	stra.2 = (prob.c)*(1-prob.c)*(prob.p^2)
	
	return(c(stra.1, stra.2))
	
}

ESp.impu1Par.A.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ii|ii
		re = rep(hap, 6)
	}
	
	if(straSeq==2){
		## i, !=i|ii|ii
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)
		re = c(sort(c(hap, sp)), rep(sort(hap), 4))
	}
	
	return(re)
	
}

ESp.impu1Par.B <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	
	prob.c1 = getHapProb2(selIdx=hap[1], semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=FALSE)
	prob.c2 = getHapProb2(selIdx=hap[2], semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=FALSE)
	
	stra.1 = (prob.c1^2)*(prob.p[1]*prob.p[2])
	stra.2 = (prob.c1)*(1-prob.c1-prob.c2)*(prob.p[1]*prob.p[2])
	stra.3 = 2*(prob.c1*prob.c2)*(prob.p[1]*prob.p[2])
	stra.4 = (prob.c2)*(1-prob.c1-prob.c2)*(prob.p[1]*prob.p[2])
	stra.5 = (prob.c2^2)*(prob.p[1]*prob.p[2])
	
	return(c(stra.1, stra.2, stra.3, stra.4, stra.5))
	
}

ESp.impu1Par.B.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ik|ik
		re = c(rep(hap[1], 2), sort(hap), sort(hap))
	}
	
	if(straSeq==2){
		## i, !=ik|ik|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)
		re = c(sort(c(hap[1], sp)), sort(hap), sort(hap))
	}
	if(straSeq==3){
		## ik|ik|ik
		re = rep(sort(hap), 3)
	}
	
	if(straSeq==4){
		## k, !=ik|ik|ik
		sp =  sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=hap)    
		re = c(sort(c(hap[2], sp)),  sort(hap), sort(hap))   
	}
	if(straSeq==5){
		## k,k|ik|ik
		re = c(rep(hap[2], 2),  sort(hap), sort(hap))
	}
	
	return(re)
	
}

ESp.impu1Par.C <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	prob.c1 = getHapProb2(selIdx=hap[1], semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol= probCol,
			snpLen, restandard=FALSE)
	prob.c2 = getHapProb2(selIdx=hap[2], semiMapFrame,
			resiProbCol= resiProbCol,
			augIdxCol= augIdxCol,
			probCol=  probCol,
			snpLen, restandard=FALSE)
	
	stra.1 = (prob.c1^2)*(prob.p[1]*prob.p[2])
	stra.2 = (prob.c1)*(1-prob.c1-prob.c2)*(prob.p[1]*prob.p[2])
	stra.3 = (prob.c1*prob.c2)*(prob.p[1]*prob.p[2])
	
	return(c(stra.1, stra.2, stra.3))
	
}

ESp.impu1Par.C.sp <-
function(hap, hap.p, straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	#print(paste("straSeq=", straSeq))
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## i,i|ij|ii
		re = c(hap, hap,   sort(c(hap, hap.p)), hap, hap)
	}
	
	if(straSeq==2){
		## i, !=ij|ij|ii
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap, hap.p))
		#print(paste("sp=", sp))
		re = c(sort(c(hap, sp)), sort(c(hap, hap.p)), hap, hap)
	}
	if(straSeq==3){
		## ij|ij|ii
		re = c(sort(c(hap, hap.p)), sort(c(hap, hap.p)), hap, hap)
	}
	
	return(re)
	
}

ESp.impu1Par.D <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	prob.c = getHapProb2(selIdx=hap, semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol= probCol,
			snpLen, restandard=FALSE)
	
	
	stra.1 = (prob.c^2)*(prob.p^2)
	stra.2 = (prob.c)*(1-prob.c)*(prob.p^2)
	
	return(c(stra.1, stra.2))
	
}

ESp.impu1Par.D.sp <-
function(hap,  straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## kk|ii|ik
		re = c(hap[2], hap[2], hap[1], hap[1], sort(hap))
	}
	
	if(straSeq==2){
		## k, !=k|ii|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap[2]))
		re = c(sort(c(hap[2], sp)),     hap[1], hap[1], sort(hap))
	}
	
	return(re)
	
}

ESp.impu1Par.E <-
function(hap, prob.p, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	
	prob.c = getHapProb2(selIdx=hap, semiMapFrame,
			resiProbCol=resiProbCol,
			augIdxCol=augIdxCol,
			probCol=probCol,
			snpLen, restandard=FALSE)
	
	
	stra.1 = (prob.c^2)*(prob.p[1]*prob.p[2])
	stra.2 = (prob.c)*(1-prob.c)*(prob.p[1]*prob.p[2])
	
	return(c(stra.1, stra.2))
	
}

ESp.impu1Par.E.sp <-
function(hap, hap.p, straSeq, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen ){
	# print(straSeq)
	re = NULL
	## give out 6 item vectors for the hap pairs for imputed par, given par and child
	if(straSeq==1){
		## k,k|ij|ik
		re = c(hap[2], hap[2], sort(c(hap[1], hap.p)), sort(hap))
	}
	
	if(straSeq==2){
		## k, !=k|ij|ik
		sp = sampleHapSemiAugMap2(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=c(hap[2]))
		re = c(sort(c(hap[2], sp)),     sort(c(hap[1], hap.p)), sort(hap))
	}
	
	
	return(re)
	
}

ESp.impuParent.heterKid <-
function(hapPair, hapProb, test=FALSE){
	ifD = FALSE
	fN = "ESp.impuParent.heterKid"
	hapCt = length(hapProb)
	
	## matrix to store the set info
	# first 2 col are hap idx, third is the probability
	setB = matrix(NA, ncol=3, nrow=5)
	setB[2,c(1,2)]=hapPair
	setB[3,c(1,2)]=rep(hapPair[1], 2)
	setB[5,c(1,2)]=rep(hapPair[2], 2)
	
	setB[2,3] = 2*hapProb[hapPair[1]]*hapProb[hapPair[2]] 
	setB[3,3] = (hapProb[hapPair[1]])^2
	setB[5,3] = (hapProb[hapPair[2]])^2
	setB[1,3] = 2*hapProb[hapPair[1]]*(1-sum(hapProb[hapPair]))
	setB[4,3] = 2*hapProb[hapPair[2]]*(1-sum(hapProb[hapPair]))
	
	
	sampleTb = matrix(NA, nrow=9, ncol=5)
	## first 2 col ae idx matched with setB
	sampleTb[,1]=c(1,1,1,2,2,2,3,3,3)
	sampleTb[,2]=c(4,2,5,4,2,5,4,2,5)
	## mating ratio
	sampleTb[,3]=c(2,2,2,2,1,2,2,2,2)
	## kids ratio
	sampleTb[,4]=c(.25, .25, .5, .25, .5, .5, .5, .5, 1)
	## last col is the sample prob
	
	## get the sample prob
	sampleTb[,5]= (setB[ sampleTb[,1], 3]) * (setB[ sampleTb[,2], 3]) * (sampleTb[,3]) * (sampleTb[,4])
	if(test){
		return(sampleTb)
	}
	if(ifD){
		print(paste(fN, "::sampling probability"))
		print(sampleTb)
	}
	
	sampleTb[,5]=sampleTb[,5]/sum(sampleTb[,5])
	chooseRow = sample(1:9, size=1, prob=sampleTb[,5])
	
	hap4idx = rep(NA, times=4)
	if(is.element(chooseRow, c(5,6,8,9))){
		#no need for resampling, the hapPair should come in as ordered
		hap4idx = c( setB[ sampleTb[chooseRow,1], c(1,2)], setB[ sampleTb[chooseRow,2], c(1,2)] )
		return(hap4idx)
	}else{
		hapProb[ hapPair ]=0
		sample.idx =  sample(1:hapCt, size=2, replace=TRUE, prob=hapProb)
		# need to resampling A:A2
		if (chooseRow==1){
			hap4idx[ c(1,3) ] = hapPair
			hap4idx[ c(2,4) ] = sample.idx
			hap4idx = c(range(hap4idx[c(1,2)]), range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# A:B
		if (chooseRow==2){
			hap4idx[ 1 ] = hapPair[1]
			hap4idx[ 2 ] = sample.idx[1]
			hap4idx[ c(3,4) ] = setB[2, c(1,2)]
			hap4idx = c(range(hap4idx[c(1,2)]), hap4idx[c(3,4)])
			return(hap4idx)
		}
		# A:C2
		if (chooseRow==3){
			hap4idx[ 1 ] = hapPair[1]
			hap4idx[ 2 ] = sample.idx[1]
			hap4idx[ c(3,4) ] = setB[5, c(1,2)]
			hap4idx = c(range(hap4idx[c(1,2)]), hap4idx[c(3,4)])
			return(hap4idx)
		}
		# B:A2
		if (chooseRow==4){
			hap4idx[ 3 ] = hapPair[2]
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx[ c(1,2) ] = setB[2, c(1,2)]
			hap4idx = c(hap4idx[c(1,2)], range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# C:A2
		if (chooseRow==7){
			hap4idx[ 3 ] = hapPair[2]
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx[ c(1,2) ] = setB[3, c(1,2)]
			hap4idx = c(hap4idx[c(1,2)], range(hap4idx[c(3,4)]))
			return(hap4idx)
		}   
		
		
	}
	
	
}

ESp.impuParent.homoKid <-
function(hapPair, hapProb, test=FALSE){
	
	hapCt = length(hapProb)
	
	## matrix to store the set info, no need anymore
	# first 2 col are hap idx, third is the probability
	
	sampleTb = matrix(NA, nrow=3, ncol=5)
	## first 2 col ae idx matched with setB
	sampleTb[,1]=c(1, 1, 2)
	sampleTb[,2]=c(1, 2, 2)
	## mating ratio, #not used anymore
	
	## kids ratio
	sampleTb[,4]=c(.25, .5, 1)
	## last col is the sample prob
	
	## get the sample prob
	
	## current approach: get the mating prob
	tp = hapProb[hapPair]
	matingProb = c(4*tp^2*(1-tp)^2, 4*tp^3*(1-tp), tp^4  )
	
	sampleTb[,5]= matingProb*sampleTb[,4]
	if(test){
		return(sampleTb)
	}
	sampleTb[,5]=sampleTb[,5]/sum(sampleTb[,5])
	chooseRow = sample(1:3, size=1, prob=sampleTb[,5])
	
	hap4idx = rep(NA, times=4)
	if(chooseRow==3){
		#no need for resampling
		hap4idx = rep(hapPair, times=4)
		return(hap4idx)
	}else{
		
		hapProb[ hapPair ]=0
		sample.idx =  sample(1:hapCt, size=2, replace=TRUE, prob=hapProb)
		# need to resampling A:A2
		if (chooseRow==1){
			hap4idx[ c(1,3) ] = rep(hapPair, 2)
			hap4idx[ c(2,4) ] = sample.idx
			hap4idx = c(range(hap4idx[c(1,2)]), range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
		# A:B
		if (chooseRow==2){
			hap4idx[ c(1,2,3) ] = rep(hapPair, 3)
			hap4idx[ 4 ] = sample.idx[1]
			hap4idx = c(hap4idx[c(1,2)],  range(hap4idx[c(3,4)]))
			return(hap4idx)
		}
	}
	
}

ESp.imputBlock <-
function(appVarNames,  trioBlock, snpLen=ncol(trioBlock),  bkIdx, job=1, snpCoding, snpBase, reType=FALSE, logF=NULL,  hapBkOnlyMap.vars){

  ## TODO!!! making missed only disappear
  fStr ="[ESp.imputBlock:]"
  ifD = FALSE
  # inside the function, assume snpCoding as c( 0, 1, 2, 3 )for NA, homo, homo, heter and snpBase as c(0, 1, 2) for NA, allele1, allele2

  
  if(ifD) print( paste(fStr, " processing block index:", bkIdx))

  if( min( c(snpCoding==c(0,1,2,3),  snpBase ==c(0,1,2))) <1 )  
     stop (paste("\nData configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""),
                                                       "] snpBase=[", paste(snpBase, collapse=";", sep=""), "]", sep=""))
  allhapKeys = get(appVarNames$freqMap, envir = qingEnv)$hapIndex
  semiMapFrame = get(appVarNames$freqMap, envir = qingEnv)$hapBkOnlyMap$bks[[ match(bkIdx, allhapKeys)]]

  ## the third row is the child
  child = trioBlock[3,]
  father = trioBlock[1,]
  mother  = trioBlock[2,]

  compMissing.trio = as.logical(apply(trioBlock, 1, sum)==0)

  cHapBkInfoMap = hapGenoBlockProc(child,  snpCoding=snpCoding)
  reqIn =  cHapBkInfoMap$homoIn
  reqDig = cHapBkInfoMap$homoDigit
  

  ## if no parents is completely missing
  if( !(compMissing.trio[1] | compMissing.trio[2])   ){
      
      fHapBkInfoMap = hapGenoBlockProc(father,  snpCoding=snpCoding)
      mHapBkInfoMap = hapGenoBlockProc(mother,  snpCoding=snpCoding)
    
    
      fHapFiltered = procSemiAugMap(appVarNames, fHapBkInfoMap, snpLen)
      mHapFiltered = procSemiAugMap(appVarNames, mHapBkInfoMap, snpLen)
      
      ## obtain completed parent filtered dip, knock off the one doesn't match with homozygous index
      fDipTb = fHapBkIdx2DipTb(appVarNames, fHapBkInfoMap, idxList=fHapFiltered, snpCoding,
                                      reqIn=reqIn, reqDigits = reqDig, expression=fHapBkInfoMap$ori, snpLen)
     
      mDipTb = fHapBkIdx2DipTb(appVarNames, mHapBkInfoMap, idxList=mHapFiltered, snpCoding,
                                      reqIn=reqIn, reqDigits = reqDig, expression=mHapBkInfoMap$ori, snpLen)


    
    if (compMissing.trio[3]){

      if(ifD) print( "no parents is compMiss, child is compMiss, sample parents.")
      if(ifD) print( "Type 2: F(d)M(d)C(CM)")
      ## no parents is compMiss, child is compMiss, sample parents.

      ## sample the non missing parent
      fProb = exDipProbSemiAugMap(fDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
      ## sample the non missing parent
      mProb = exDipProbSemiAugMap(mDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(2, job)
      for( ss in 1:job){
         ## sample fa then sample mo...
         chooseRow = sample(length(fProb), size=1, prob=fProb)
         tmpFather = fDipTb[chooseRow,]
   
         chooseRow = sample(length(mProb), size=1, prob=mProb)
         tmpMother = mDipTb[chooseRow,]
   
         tmpChild = matrix( c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                              tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=FALSE)
         t.choice = sample(4, size=1)
   
         trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice])]))
   
       }
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12, drop=FALSE])
      }
      
    }else{
      ## no parents is compMiss, child is NOT compMiss, get all and exclude the fDipTb and mDipTb for those doesn't fit cDipTb

      if(ifD) print( "no parents is compMiss, child is NOT compMiss, get all and exclude the fDipTb and mDipTb for those doesn't fit cDipTb")
      if(ifD) print( "Type 6: F(d)M(d)C(d)")
      ## if none is completed missing
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)

      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                      reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

##       print("###")
##       print(cDipTb)
##       print(fDipTb)
##       print(mDipTb)
      
      ## need to exclude the fDipTb and mDipTb for those doesn't fit cDipTb 
      tmpDip = exParentDip(cDipTb, fDipTb, fHapBkInfoMap,  snpLen)
      cDipTb = tmpDip$childTb
      fDipTb = tmpDip$parTb
#      print(tmpDip)
    
      tmpDip = exParentDip(cDipTb, mDipTb, mHapBkInfoMap, snpLen)
      cDipTb = tmpDip$childTb
      mDipTb = tmpDip$parTb
#      print(tmpDip)
      
      if(!is.null(logF)){
           logl(logF, paste("Possible dip for father with some data: row=", length(fDipTb)/2))
           logl(logF, paste(util.matrix.cat(fDipTb, 1:2, sep="."), collapse="; ", sep=""))
      
           logl(logF, paste("Possible dip for mother with some data: row=", length(mDipTb)/2))
           logl(logF, paste(util.matrix.cat(mDipTb, 1:2, sep="."), collapse="; ", sep=""))
    
           logl(logF, paste("Possible dip for child with some data: row=", length(cDipTb)/2))
           logl(logF, paste(util.matrix.cat(cDipTb, 1:2, sep="."), collapse="; ", sep=""))
      }
    
      ## need to refilter/standandize the pair probability
      ## obtain the diplotype map and probability
      
      prob1 = exDipProbSemiAugMap(fDipTb, semiMapFrame,
                         hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
      
      prob2 = exDipProbSemiAugMap(mDipTb, semiMapFrame,
                         hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
    
      if(ifD) {
        print("Possible diplotype (hap pair) indexes for father (fDipTb); and probability")
        print(fDipTb)
        print(prob1)
        print("Possible diplotype (hap pair) indexes for mother (mDipTb); and probability")
        print(mDipTb)
        print(prob2)
        print("Possible diplotype (hap pair) indexes for child  (cDipTb); ")
        print(cDipTb)
      }

      ## otherwise, build the mating table based on children geno      
      mating6hap = genMatingTBCondOnChild3(appVarNames,  child= cDipTb,
                             par1=fDipTb, par2=mDipTb, prob1, prob2, logF=logF, job=job)

      othChildtt = apply(mating6hap, 1, FUN=find.PsudoControlHap)
      
      trioHapIdx =  cbind(mating6hap[ ,1:4, drop=FALSE],  t(othChildtt))

      #if (ifD) print("@@@@@ return obj@@@@@")
      if(reType){
        #if(ifD) print(  cbind(trioHapIdx, rep(6, job)) )
        return(cbind(trioHapIdx, rep(6, job)))
      }else{
        #if(ifD) print(trioHapIdx)
        return(trioHapIdx)
      }
      
    }
  }

  ## if both parents are completely missing
  if( compMissing.trio[1] & compMissing.trio[2] ){
    if (compMissing.trio[3]){

      if(ifD) print("both parents are compMiss, child is compMiss, sample parents")
      if(ifD) print("Type 1: F(CM)M(CM)C(CM)")

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(1, job)
      for(ss in 1:job){
        ## both parents are compMiss, child is compMiss, sample parents.
        tmpFather = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
        tmpMother = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)
  
        tmpChild = matrix( c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                             tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=FALSE)
        t.choice = sample(1:4, size=1)
  
        trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice]) ]))
  
      }

      
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[, 1:12, drop=FALSE])
      }
      
    }else{

      if(ifD) print("both parents are compMiss, child is NOT compMiss, sample kids first, then imput parent (ESp method)")
      if(ifD) print("Type 4: F(CM)M(CM)C(d)")
      ## both parents are compMiss, child is NOT compMiss, sample kids first, then imput parent (ESp method)

      supHapProb = getHapProb.semiMapFrame( semiMapFrame, hapBkOnlyMap.vars$resiProbCol,
                           hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      ## for child not completely missing, need to restandardize the other parents' hap freq
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)
      ## need one function to FIX!!!FIX it 
      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                  reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

      ## sample child
      cProb = exDipProbSemiAugMap(cDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol,
                                  hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(4, job)
      for(ss in 1:job){
        chooseRow = sample(length(cProb), size=1, prob=cProb)
        hapPair = cDipTb[chooseRow,]
  
        if (hapPair[1]==hapPair[2]){
          matingTbl = ESp.impuParent.homoKid(hapPair[1], supHapProb)
        }else{
          matingTbl = ESp.impuParent.heterKid(hapPair, supHapProb)
        }
  
        trioHapIdxtt = c(matingTbl, hapPair)
        trioHapIdx[ss,1:12] =  c(matingTbl,  find.PsudoControlHap(trioHap6=trioHapIdxtt))
      }

      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12,drop=FALSE])
      }
    }
  }
  

  ## if only one parent is completely missing
  if( sum(compMissing.trio[1:2])==1 ){
    if(ifD) print( "only one parent is completely missing, need to sample the non-comp-missing parents from a restricted list")
    
    ## need to sample the non-comp-missing parents from a restricted list
    nonMparent = trioBlock[!compMissing.trio, ,drop=FALSE][1, ]
    pHapBkInfoMap = hapGenoBlockProc(nonMparent,  snpCoding=snpCoding)
    pHapFiltered = procSemiAugMap(appVarNames, pHapBkInfoMap, snpLen)

#     print(pHapBkInfoMap)
#     print(pHapFiltered)
    ## obtain completed parent filtered dip, knock off the one doesn't match with homozygous index
    pDipTb = fHapBkIdx2DipTb(appVarNames, pHapBkInfoMap,
                                     idxList=pHapFiltered, snpCoding = snpCoding,
                                     reqIn=reqIn, reqDigits = reqDig, expression=pHapBkInfoMap$ori, snpLen)

    if (compMissing.trio[3]){
      ## one parents is compMiss, child is compMiss, sample parents.
      ## just use the popu hap freq for the missing parent
      if(ifD) print( "one parents is compMiss, child is compMiss, sample parents ")
      if(ifD) print("Type 3: F(CM)M(d)C(CM)")
      missingParent = sampleDipSemiAugMap(semiMapFrame, hapBkOnlyMap.vars$resiProbCol, hapBkOnlyMap.vars$augIdxCol, hapBkOnlyMap.vars$probCol, snpLen)

      ## sample the non missing parent
      pProb = exDipProbSemiAugMap(pDipTb, semiMapFrame=semiMapFrame,
                                  hapBkOnlyMap.vars$resiProbCol,
                                  hapBkOnlyMap.vars$augIdxCol,
                                  hapBkOnlyMap.vars$probCol, snpLen)

      trioHapIdx = matrix(NA, nrow=job, ncol=13)
      trioHapIdx[,13]=rep(3, job)
      for(ss in 1:job){
        
        chooseRow = sample(length(pProb), size=1, prob=pProb)
        nomissingParent = pDipTb[chooseRow,]
  
        if(compMissing.trio[1]){
          tmpFather = missingParent
          tmpMother = nomissingParent
        }else{
          tmpFather = nomissingParent
          tmpMother = missingParent
        }
  
        tmpChild = matrix(  c(tmpFather[1], tmpMother[1], tmpFather[1], tmpMother[2],
                              tmpFather[2], tmpMother[1], tmpFather[2], tmpMother[2]), nrow=2, byrow=FALSE)
        
        t.choice = sample(4, size=1)
  
        trioHapIdx[ss,1:12] =  c(tmpFather, tmpMother, as.vector(tmpChild[, c(t.choice, (1:4)[-t.choice]) ]))
      }
      
      if(reType){
        return(trioHapIdx)
      }else{
        return(trioHapIdx[,1:12,drop=FALSE])
      }
      
    }else{

      if(ifD) print( "one parents is compMiss, child is NOT compMiss, use ESp method")
      if(ifD) print("Type 5: F(CM)M(d)C(d)")
      ## one parents is compMiss, child is NOT compMiss, use ESp method.

      ## for child not completely missing, need to restandardize the other parents' hap freq
      cHapFiltered = procSemiAugMap(appVarNames, cHapBkInfoMap, snpLen)
      ## need one function to FIX!!!FIX it 
      cDipTb = fHapBkIdx2DipTb(appVarNames, cHapBkInfoMap, idxList=cHapFiltered, snpCoding,
                                  reqIn=NULL, reqDigits = NULL, expression=cHapBkInfoMap$ori, snpLen)

      #print(cHapFiltered)
      #print(cDipTb) 
      trioHapIdx.f = ESp.impu1Par(
        othParPairs=pDipTb, childPairs=cDipTb,  semiMapFrame,
        resiProbCol=hapBkOnlyMap.vars$resiProbCol,
        augIdxCol=hapBkOnlyMap.vars$augIdxCol,
        probCol=hapBkOnlyMap.vars$probCol,
        snpLen, reType=reType, job=job)

      if(ifD) print(paste("Outcome hap:", paste(trioHapIdx.f, collapse=";")))
      
      ## need to know switch parents if the mom is missing
      if( compMissing.trio[1] ){
        othChildtt = apply(trioHapIdx.f, 1, FUN=find.PsudoControlHap)      
        trioHapIdx =  cbind(trioHapIdx.f[, 1:4, drop=FALSE],  t(othChildtt))
      }else{
        othChildtt = apply(trioHapIdx.f, 1, FUN=find.PsudoControlHap)
        trioHapIdx =  cbind(trioHapIdx.f[, c(3,4, 1,2), drop=FALSE],  t(othChildtt))
      }

      if(reType){
        return(cbind(trioHapIdx, 5+trioHapIdx.f[,7]/10))
      }else{
        return(trioHapIdx)
      }
 
    }
  }
  
}

exchangeDigit <-
function(ma, cols=NULL, dig1Code=c(0, 1, 3, 2), dig2Code = c(0, 1, 2), action=c("1to2", "2to1")){
	ifD = FALSE
	
	if(is.null(cols)){
		cols = c(1, ncol(ma))
	}
	
	fun.error = FALSE
	if(length(cols)!=2) fun.error=TRUE
	if(cols[2]>ncol(ma)) fun.error=TRUE
	
	outputColCt = cols[2]-cols[1]
	if( outputColCt <=0) fun.error=TRUE
	
	if(length(action)>=2) {
		
		if(length(action)>1) warning(paste("Two or more actions is request (", paste(action, collapse="; "),
							"). Only the first requested is performed."), sep="")
		action=action[1]
		
		if ((action!="1to2") & (action!="2to1"))
			stop(paste("Requested action, (",  action, "), is not implemented.", sep="" ))
	}
	ma.change = ma[, cols[1]:cols[2], drop=FALSE]
	## internally, do not use NA to represent missing,
	code.0 = 0
	code.012 = dig2Code
	code.0123 = dig1Code
	
	if( max( is.na(dig2Code[2:3]))==1) stop("NA is not allow to represent the non-missing allele!")
	if( max( is.na(dig1Code[2:4]))==1) stop("NA is not allow to represent the non-missing genotype!")
	
	if(action=="1to2"){
		if( is.na(dig1Code[1]) ){
			code.123 = dig1Code[2:4]
			## if the min is the same as the default code for NA(=0), then use the min -1 
			if( max(code.123==0)==1) code.0 = min(code.123)-1
			code.0123 = c(code.0, code.123)
			ma.change = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=code.0123)
		}
		outputColCt = outputColCt/2
		tmpBridge = dig2Code[c(1, 2, 2, 3)]
		tmpBridge2 = dig2Code[c(1, 2, 3, 3)]
		if(ifD) print(paste("tmpBridge:", paste(tmpBridge, collapse=";")))
		if(ifD) print(paste("tmpBridge2:", paste(tmpBridge2, collapse=";")))   
	}
	if(action=="2to1"){
		if( is.na(dig2Code[1]) ){
			#need to replace the given NA with 0 or others
			code.12 = dig2Code[2:3]
			if( max(code.12==0)==1) code.0 = min(code.12)-1
			code.012 = c(code.0, code.12)
		}
		ma.change = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = dig2Code, replaceBy=code.012)
		outputColCt = outputColCt*2
		sumBase = matrix(c(2, 0, 0,
						0, 2, 0,
						0, 1, 1,
						0, 0, 2), byrow=FALSE, nrow=3)
		sumCode = as.vector(matrix(code.012, ncol=3)%*%sumBase)
		if(ifD) print(paste("sumCode:", paste(sumCode, collapse=";")))
	}
	
	
	if (fun.error){
		stop("Argu, cols, refer to the ranges of the index for the columns in argu, ma.
						The current values are not valid.")
	}
	
	
	
	if(cols[1]>1) {
		ma.kept1 = ma[, 1:(cols[1]-1)]
	}else{
		ma.kept1 = NULL
	}
	
	ma.kept2=NULL
	if(cols[2]<ncol(ma)) ma.kept2 = ma[, (cols[2]+1):ncol(ma)]
	
	ma.ex=NULL
	if(action=="1to2"){
		
		
		ma.2dig1 = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = code.0123, replaceBy=tmpBridge)
		ma.2dig2 = apply(ma.change, 1:2, FUN= util.vec.replace, orignal = code.0123, replaceBy=tmpBridge2)
		
		## make sure the small digit is at the front
		if(dig2Code[2]<dig2Code[3]){
			ma.ex = util.matrix.col.shuffle2(ma.2dig1, ma.2dig2)
		}else{
			ma.ex = util.matrix.col.shuffle2(ma.2dig2, ma.2dig1)
		}
		
	}
	
	if(action=="2to1"){
		
		seq1 = seq.int(from=1, to=ncol(ma.change), by=2)
		
		ma.sum = matrix(ma.change[, seq1]+ma.change[,seq1+1], ncol=length(seq1), byrow=FALSE)
		
		## convert to Qing's 1-digit coding system. 3 is for hetero
		
		ma.ex = apply(ma.sum, 1:2, FUN= util.vec.replace, orignal = sumCode, replaceBy=dig1Code)
		
	}
	if(!is.null( ma.kept1)){
		all = cbind(ma.kept1, ma.ex)
	}else{
		all = ma.ex
	}
	if(!is.null( ma.kept2)){
		all = cbind(all, ma.kept2)
	}
	return(all)
	
}

exDipProbSemiAugMap <-
function(dipMap, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = FALSE
	
	dipProb = as.vector(dipMap)
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	mapProb = semiMapFrame[ , probCol]
	
	commonProbIdx = length(mapProb) + 1
	
	curIdx = match(dipProb, mappedAugIdx, nomatch=commonProbIdx )
	fittedFilter = curIdx==commonProbIdx
	
	## find out how many un/specified haplotype will fit the observed data
	fittedInMap= unique(curIdx[!fittedFilter])
	fittedOutMap = unique(dipProb[fittedFilter])
	
	if(ifD) print("fittedInMap")
	if(ifD) print(fittedInMap)
	if(ifD) print("fittedOutMap")
	if(ifD) print(fittedOutMap)
	
	## restandarize the prob
	if(length(fittedOutMap)>0){
		mapProbRe=rep(0, commonProbIdx)
		mapProbRe[commonProbIdx]=leftOverProb/length(fittedOutMap)
		mapProbRe[fittedInMap] = mapProb[fittedInMap]/sum(mapProb[fittedInMap])*(1-leftOverProb)
		
	}else{
		mapProbRe =rep(0, commonProbIdx-1)
		mapProbRe[fittedInMap] = mapProb[fittedInMap]/sum(mapProb[fittedInMap])
		
	}
	
	
	if(ifD) print(mapProbRe)
	
	dipProb = mapProbRe[curIdx]
	
	if(ifD) print(dipProb)
	
	## remember the dipMap is passed in as a vector filted by column
	idxSeq = 1:(length(curIdx)/2)
	reDipProb = dipProb[idxSeq] * dipProb[idxSeq+ (length(curIdx)/2) ]
	
	return(reDipProb)
	
}

exhaustHapExp <-
function(lociCt=3, re.str=TRUE,  snpCoding=c(1,2)){

  ma = matrix(snpCoding, ncol=1)

  ## grow the matrix on both sides
  if(lociCt==1) return(list(hap=ma, hapStr=as.character(ma)))
  for( i in 1:(lociCt-1) ){
    growing = util.matrix.clone(ma, 2)    
    addedBit = c(rep(snpCoding[1], times=2^i), rep(snpCoding[2], times=2^i))
    ma = cbind(growing, " "=addedBit)
  }
  if(re.str) {
    hapStr = util.matrix.cat(ma, 1:lociCt, sep="")
    return(list(hap=ma, hapStr=hapStr))
  }else{
    return(list(hap=ma))
  }
  
}

exIdxFromHap <-
function(lociCt){
  ## know the internal sequency of haplotype
  digitMap1 = matrix(NA, ncol=lociCt, nrow=2^(lociCt-1) );
  digitMap2 = matrix(NA, ncol=lociCt, nrow=2^(lociCt-1) );

  for( n in 1:lociCt ){
     col = seq.int(from=1, to=2^lociCt, by = 2^n)
     col2 = col+2^(n-1)

     colNext = col
     colNext2 = col2

     for( j in 1:(2^(n-1))){
       if( j > 1){
         newCol = col + j - 1
         colNext = rbind(colNext, newCol)

         newCol2 = col2 + j - 1
         colNext2 = rbind(colNext2, newCol2)
         ## print(colNext)
       }
     }

     digitMap1[,n]=as.vector(colNext)
     digitMap2[,n]=as.vector(colNext2)
   }

  
  return (list(digitMap1=digitMap1, digitMap2 = digitMap2))
}

exParentDip <-
function(childTb, parTb, pHapBkInfoMap, snpLen){
	## if a haplotype doesn't include in the parents or child, exclude
	
	if(is.null(dim(childTb)) | is.null(dim(parTb))  | length(pHapBkInfoMap$missingIn)==snpLen){
		return(list(childTb=matrix(childTb, ncol=2), parTb=matrix(parTb, ncol=2)))
	}
	
	child = unique(childTb)
	parent = unique(parTb)
	
	childLeft = child[is.element(child, parent)]
	parentLeft = parent[is.element(parent, child)]
	
	if(length(childLeft)==0 | length(parentLeft)==0 ) {
		## print(paste("ChildTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
		stop(paste("\n ChildTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
	}
	
	
	childIdx = apply(childTb, 1, FUN = function(rowItem, bench){
				re = as.logical(max(is.element(rowItem, bench)))
			}, bench = childLeft)
	
	
	parIdx = apply(parTb, 1, FUN = function(rowItem, bench){
				re = as.logical(max(is.element(rowItem, bench)))
			}, bench = parentLeft)
	
	if(sum(childIdx)==0 | sum(parIdx)==0 ) {
		## print(paste("childTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
		stop(paste("\n childTb:\n", childTb, "parTb:\n", parTb, "\n Cannot find the pair matched"))
	}
	return(list(childTb=childTb[childIdx,, drop=FALSE], parTb=parTb[parIdx,,drop=FALSE]))
	
}

fHapBkIdx2DipTb <-
function(appVarNames, filteredBkInfo, idxList, snpCoding, reqIn=NULL, reqDigits = NULL, expression, snpLen=nchar(expression)){
	ifD = FALSE
	fStr = "[fHapBkIdx2DipTb]:"
	if(ifD) {
		print(qp(fStr, "start"))
		print(filteredBkInfo)
	}
	if( min( snpCoding==c(0,1,2,3)) <1 )  
		stop (paste("\nData configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""), "]", sep=""))
	
	
	hetoIn = filteredBkInfo$hetoIn
	
	if(!is.null(hetoIn)){
		## need to build up the dip pair following the heto digit requirment\
		bkCt = length(idxList)
		
		if(bkCt<2)
			stop(paste("\nCannot find more than two haplotype to build heto digit for data:(", expression, ").", sep=""))
		hetoSeq = hetoIn
	}else{
		## only need to process the digit requirement posed by child
		hetoSeq = NULL
		bkCt = length(idxList)
	}
	
	if(ifD & (!is.null(hetoSeq))) print(paste("hetoSeq:", paste(hetoSeq, collapse=";", sep="")))
	
	## if no heto digit and no requirment posed by child, just return the all possible dip pair
	if(is.null(reqIn) & is.null(hetoIn) ){
		dipMap = util.it.upTriCombIdx(idxList, diag=TRUE, re.ordered=TRUE)
		#print(dipMap)
		return(dipMap)
	}
	
	
	## if no heto digit but some requirement posed by child
	if(is.null(hetoIn) & (!is.null(reqIn)) ){
		dipMap = NULL
		for ( i in 1:bkCt ){
			curBkIdx = idxList[i]
			## because no heto digits restriction, then dip built by the same hap should be considered.
			j = i
			
			while ( j <= bkCt ){
				## compare the heto digits
				othBkIdx = idxList[j]
				
				## check the required digits
				meetReq = TRUE
				tmpReq = 1
				while( tmpReq <= length(reqIn)){
					curDigitReq = reqDigits[tmpReq]
					tmpSeq = seq.int(from=1, to=2^(reqIn[tmpReq]-1), by=1)
					oneMeetReq = isDigitAtLociIdx(hapIdx=curBkIdx, digit=curDigitReq,
							lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
					othMeetReq = isDigitAtLociIdx(hapIdx=othBkIdx, digit=curDigitReq,
							lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
					
					meetReq = meetReq & (oneMeetReq | othMeetReq)
					tmpReq = tmpReq + 1
					if(!meetReq) tmpReq = length(reqIn)+10
					
				}
				## keep the matched pair
				if(meetReq)  dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
				j = j + 1
			} # while ( j <= bkCt ){
		} # for ( i in 1:bkCt ){
		if(is.null(dipMap))
			stop(paste("\nCannot find the haplotype pairs for parent meet digit requirement, exp=(", expression, ").", sep=""))
		#print(dipMap)
		return(dipMap)   
	}
	
	## if( !is.null(hetoIn) & [ is.null(reqIn) | !is.null(reqIn) ] )
	
	idx4hapDigit = NULL
	## check out the global variables
	tryCatch({
				
				tmpGetObj = NULL
				tmpGetObj = get(appVarNames$digit, envir=qingEnv)
				
				idx4hapDigit$digitMap1 = tmpGetObj$digitMap1[1:(2^(snpLen-1)), 1:snpLen]
				idx4hapDigit$digitMap2 = tmpGetObj$digitMap2[1:(2^(snpLen-1)), 1:snpLen]
				
				rm(tmpGetObj)
				##gc()
			}, error=function(e){
				errTrace = paste(e, collapse=";", sep="") 
				stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$digit, " does not exisit."))
			})
	
	
	dipMap = NULL
	## check for all required heto
	
	## all the idxList meet the homo digit requirement, but for heto digits
	## we need to match up the pair for all the heto digits
	#print(idxList)
	#print(idx4hapDigit)
	for ( i in 1:bkCt ){
		curBkIdx = idxList[i]
		## because heto digits restriction, then dip built by the same hap should be excluded.
		j = i + 1 
		while( j <= bkCt ){
			## compare the heto digits
			othBkIdx = idxList[j]
			
			if(ifD) print(paste("i=", i, "; j=", j, sep=""))
			
			## suppose cur bk contribute a digit 1
			curLociMatch1 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap1[, hetoSeq, drop=FALSE], val=curBkIdx)
			## suppose oth bk contribute a digit 2
			othLociMatch2 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap2[, hetoSeq, drop=FALSE], val=othBkIdx)
			
			## suppose oth bk contribute a digit 1
			othLociMatch1 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap1[, hetoSeq, drop=FALSE], val=othBkIdx)
			## suppose cur bk contribute a digit 2
			curLociMatch2 = util.matrix.colIdx4Match(ma=idx4hapDigit$digitMap2[, hetoSeq, drop=FALSE], val=curBkIdx)               
			
			#print(curLociMatch1)
			#print(othLociMatch2)
			
			#print(othLociMatch1)
			#print(curLociMatch2)
			
			## change code, may cause trouble      
			##        if(length(curLociMatch1)==0) curLociMatch1 = 100
			##        if(length(curLociMatch2)==0) curLociMatch2 = 101
			##        if(length(othLociMatch1)==0) othLociMatch1 = 102
			##        if(length(othLociMatch2)==0) othLociMatch2 = 103
			##        tmp = sum(c( suppressWarnings(curLociMatch1==othLociMatch2),
			##                    suppressWarnings(curLociMatch2==othLociMatch1)))
			
			
			## check same length
			tmp=0
			
			if(length(curLociMatch1)==length(othLociMatch2) ){
				if(length(curLociMatch1)!=0 ){
					tmp=sum(curLociMatch1==othLociMatch2)
				}
			}
			
			if(length(curLociMatch2)==length(othLociMatch1) ){
				if(length(curLociMatch2)!=0 ){
					tmp=tmp+sum(curLociMatch2==othLociMatch1)
				}
			}
			
			
			if(tmp!=0){
				
				if(tmp == length(hetoSeq)){
					## find the matched pair
					## check the required digit
					if(is.null(reqIn)){
						## keep the matched pair
						dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
						if(ifD) print(paste("curBkIdx=", curBkIdx, ": othBkIdx=", othBkIdx, sep=""))
					}else{
						## check the required digits
						meetReq = TRUE
						tmpReq = 1
						while( tmpReq <= length(reqIn)){
							curDigitReq = reqDigits[tmpReq]
							tmpSeq = seq.int(from=1, to=2^(reqIn[tmpReq]-1), by=1)
							oneMeetReq = isDigitAtLociIdx(hapIdx=curBkIdx, digit=curDigitReq,
									lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
							othMeetReq = isDigitAtLociIdx(hapIdx=othBkIdx, digit=curDigitReq,
									lociCt=snpLen, lociIdx=reqIn[tmpReq], intVec = tmpSeq)
							
							meetReq = meetReq & (oneMeetReq | othMeetReq)
							tmpReq = tmpReq + 1
							if(meetReq) tmpReq = length(reqIn)+10
							
						}
						## keep the matched pair
						if(meetReq)  {
							dipMap = rbind(dipMap, range(curBkIdx, othBkIdx))
							if(ifD) print(paste("curBkIdx=", curBkIdx, "; othBkIdx=", othBkIdx, sep=""))
						}
					} # if(is.null(reqIn)){
				} # if(tmp == length(hetoSeq)){
			} #  if(tmp!=0){
			
			j = j + 1
			
		} # while( j <= bkCt ){
	} # for ( i in 1:bkCt ){
	rm(idx4hapDigit)
	##gc()
	if(is.null(dipMap))
		stop(paste("\n Cannot find the haplotype pairs for parent meet heto/digit requirement, exp=(", expression, ").", sep=""))
	
	#print(dipMap)
	
	return(dipMap)
	
}

filterHapIdx2SuperHapSet <-
function( bkIdxes, hapIdxes, hapCts){

  it = lapply(hapCts, 1, FUN=seq, from=1)
  # first build the whole list
  idxComb = qExpandTable(listOfFactor = it, removedRowIdx=NULL, re.row=FALSE )
  
  # eliminate the impossible rows.
  left = 1: (cumprod(hapCts)[length(hapCts)]  )
  for ( i in 1:length(bkIdxes)){
    matchedVal = hapIdxes[[i]]
    set = left[ is.element(idxComb [left ,bkIdxes[i]], matchedVal) ]
    #print(set)
    left = set
  }
  return(set)

}

filterHaps2SHaps <-
function(hapForBk, hapCts){
  ifD = FALSE
  fN = "filterHaps2SHaps:"
  if(ifD) {
    print(paste(fN, "begin"))
    print(hapCts)
  }
  bkCt = length(hapCts)

  if(bkCt==1) {
    if (is.na(hapForBk)) {
      return(1:hapCts)
    }else{
      return (hapForBk)
    }
  }

  segLen = cumprod(hapCts)
  
  filter = is.na(hapForBk)

  cut.idx = (1:bkCt)[filter]
  endWithNA = TRUE

  # if the whole pattern doesn't end with NA
  if(!filter[bkCt]){
    endWithNA = FALSE
    cut.idx.end = c(cut.idx, bkCt)
    cut.idx.sta = c(1, cut.idx+1)
  }else{
    # if the whole pattern ends with NA
    cut.idx.end = cut.idx
    cut.idx.sta = c(1, cut.idx[-length(cut.idx)]+1) 
  }

  if(ifD) print(cbind(cut.idx.sta, cut.idx.end))
  segCt = length(cut.idx.sta)
  j = 1
  grp.idx = NULL
  
  while( j <=segCt){
    if(ifD) print(paste("j=", j))
    seg = hapForBk[cut.idx.sta[j]:cut.idx.end[j]]
    if(ifD) print(seg)
    offset.q = 0
    segSeq = 0
    ## calculate the offset
    if( (cut.idx.end[j]-cut.idx.sta[j])!=0){
      ## if the seg has other than NA, we have offset
      if(j==1 & j==segCt & !endWithNA){
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]) ]
        nums = seg
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
        return(offset.q)
      }
      if(j==1 & j==segCt & endWithNA){
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = seg[ - length(seg)]
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
        n.prec = segLen[cut.idx.end[j]-1]
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
        return(segSeq+offset.q)
      }      
      if(j==1){
        ## if it is the first segment and with length >= 2, must have an offset
        hapCc = hapCts[ (cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = seg[ - length(seg)]
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)
      }
      if(j>1 & j< segCt){
        ## if it is the middle segment and with length >= 2, must have an offset
        hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]-1) ]
        nums = c(1, seg[ - length(seg)])
        hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
        ## the out of box offset added one already, but we want to add another offset on it 
        offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
      }
      if(j==segCt){
        if(endWithNA){
          ## if it is the last segment and with length >= 2 and have a sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]-1) ]
          nums = c(1, seg[ - length(seg)])
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
        }else{
          ## if it is the last segment and with length >= 2 and have no sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]) ]
          nums = c(1, seg)
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
        } ## if(endWithNA){
      } ## if(j==segCt){
      if(ifD) print(offset.q)

      if(j!=segCt | endWithNA){
        ## it has a sequence as well
        if(cut.idx.end[j]==1){
          n.prec = 1
          
        }else{
          n.prec = segLen[cut.idx.end[j]-1]
        }        
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
      }
        
    }else{ ###  if( (cut.idx.end[j]-cut.idx.sta[j])!=0){
      ## if only one element in the segment

      if (j==segCt & !endWithNA){
          ## if it is the last segment and with length >= 2 and have no sequence cutoff
          hapCc = hapCts[(cut.idx.sta[j]):(cut.idx.end[j]) ]
          nums = c(1, seg)
          hapCc = c(segLen[(cut.idx.sta[j]-1) ], hapCc)
          ## the out of box offset added one already, but we want to add another offset on it 
          offset.q = calHapIdx2SHap(nums, hapCts=hapCc)-1
          segSeq = 0
      }else{
        ## it has a sequence only
        if(cut.idx.end[j]==1){
          n.prec = 1
          ## for the first element to be a NA
          offset.q = 1
        }else{
          n.prec = segLen[cut.idx.end[j]-1]
        }        
        segSeq = seq.int(from=0, to=segLen[cut.idx.end[j]]-1, by = n.prec)
      } ## if (j==segCt & !endWithNA){

    } ###  if( (cut.idx.end[j]-cut.idx.sta[j])!=0){

    if(ifD) print("segSeq")
    if(ifD) print(segSeq)
    if(ifD) print(offset.q)

    if(j==1){
      ## first seg
      grp.idx = segSeq + offset.q
    }else{
      oth.idx = segSeq + offset.q
      all.idx = qExpandTable(listOfFactor =list(grp.idx, oth.idx), removedRowIdx=NULL, re.row=FALSE)
      if(ifD) print(all.idx)
      grp.idx = all.idx[,1]+all.idx[,2]
    }
    if(ifD) print(grp.idx)
    j = j+1
  }
  return(grp.idx)

}

filterHaps2SHaps.check <-
function(hapForBk, hapCts){
  ifD = FALSE
  fN = "filterHaps2SHaps.check::"
  hapList = lapply(hapCts, FUN=function(i) {1:i})
  aa = qExpandTable(listOfFactor =hapList, removedRowIdx=NULL, re.row=FALSE)

  filter = !is.na(hapForBk)
  comp = rep(TRUE, times=nrow(aa))

  
  for( i in 1:length(hapForBk)){
    if(!is.na(hapForBk[i])){
      ff = aa[,i]==hapForBk[i]
      comp =  comp & ff
    }
  }
  row.idx = (1:nrow(aa))[comp]
  if(ifD) print(row.idx)
  
  my.idx = filterHaps2SHaps(hapForBk, hapCts)

  dis = sum(is.na(match(my.idx, row.idx)))
  len.m = length(my.idx)-length(row.idx)

  if (sum(dis, len.m)==0) return (NULL)
  print("compare:my.idx with row.idx")
  print(my.idx)
  print(row.idx)
  stop(paste(fN, "ERROR: not matched index."))
  return(NULL)
}

find.PsudoControlHap <-
function(trioHap6){
  
  child.sort = sort(trioHap6[5:6])
  
  othChild = matrix(trioHap6[1:4][c(1, 3, 1, 4, 2, 3, 2, 4)], ncol=2, byrow=TRUE)
  othChild = apply(othChild, 1, sort)
  
  child.m = NULL
  for( i in 1:4){
    oth = othChild[,i]
    if(sum(child.sort==oth)==2){
      child.m = c(child.m, i)
    }
  }
  if(length(child.m)<1) stop("No match for the child")

  t.choice = child.m[sample(length(child.m), size=1)]
 
  re = othChild[, c(t.choice, (1:4)[-t.choice]) ]

  return(re)
}

findLastPreviousDate <-
function(timeVec, probVec, cutoff){
  if(max(timeVec)< cutoff) return(NA)
  filter = which(timeVec<= cutoff)
  re = probVec[filter[length(filter)]]
  return(re)
}

findMissing <-
function(df, is.1digit=FALSE, snpStartLeftIndex, snpEndRightIndex, dig1Code=c(0, 1, 2, 3), dig2Code=c(0, 1, 2) ){
  ## use c(0,1,3,2) as the 1-digit coding, if not, exchange it
  snpEndLeftIndex = ncol(df)-snpEndRightIndex+1

  if(!is.1digit){
     snpNum = (snpEndLeftIndex - snpStartLeftIndex +1) /2
     snps = df[, snpStartLeftIndex:snpEndLeftIndex]

     snp1digit = exchangeDigit(ma=snps, cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
 
   }else{
     snpNum = snpEndLeftIndex - snpStartLeftIndex +1
     snp1digit = df[, snpStartLeftIndex:snpEndLeftIndex]
     if(min(dig1Code==c(0,1,3,2))==0){
       snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,3,2))
     }
   }
   nrow = nrow(df)
   missingPos = which(snp1digit==0)

   if(length(missingPos)==0) return(NULL)

   missingCol = (missingPos - (missingPos %% nrow))/nrow + 1

   miss.cord = matrix(NA, nrow=length(missingPos), ncol=2)

   miss.cord[,1]= (missingPos %% nrow)
   ## adjust the one at the last row
   miss.cord[  miss.cord[,1]==0, 1]= nrow
  
   miss.cord[,2]= (missingPos - miss.cord[,1])/nrow + 1

   colnames(miss.cord)=c("row", "snp")
   return(miss.cord)

}

freq.allele <-
function(ct, order = c("majorHM", "heter", "minorHM")){

  if ( sum(is.na(match(order, c("majorHM", "heter", "minorHM")))) != 0) stop ("Wrong order of count.")
  freq = (2*ct[1]+ct[2])/( 2*sum(ct) )
  freqs = c(freq, 1-freq)
  #print(ct)
  if ( freq>1 ) warning("Allele frequency is greater than 1")
  return(freqs)

}

freq.build <-
function(hap=NULL, geno){
	
	## make it can take .csv file
	if(!is.null(hap)){
		if(is.character(hap)){
			hap = read.csv(hap, header=TRUE, sep=",", as.is=TRUE)
		}
	}
	
	## make it can take .csv file
	if(is.character(geno)){
		geno = read.csv(geno, header=TRUE, sep=",", as.is=TRUE)
	}
	
	## check input data
	test = match(c("prekey", "seq", "freq"), colnames(geno))
	if (sum(is.na(test))>=1)
		stop("Input argument, geno, does not have the required columns, i.e., prekey, seq, and freq.") 
	
	if(is.null(hap)){
		warning("NULL value is provided for input argument, hap.")
	}
	
	freq = c(genoMap.info = list(geno[, match(c("prekey", "seq", "freq"), colnames(geno))]),
			hapMap.info = list(hap))
	
	return(freq)
}

freqbuild.haponly <-
function(hap, alleleCode = 1:2){
	## only hap is provided.
	## need to find out the hap with 1 loci, and get the geno freq from the 1 loci
	ifD = FALSE
	
	## check the maximum length for imputation
	bkMap = bkMap.constr(data=hap, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=alleleCode )
	if(sum(bkMap$snpCt >=8 )>0)  stop("At least one of the block size exceeds the maximum value of 7 loci.") 
	
	prekey = rep("ch", nrow(hap))
	cumidx = cumsum(bkMap$snpCt)
	markers.e = cumidx
	
	markers.b = c(1, (cumidx+1)[-length(cumidx)])
	others = lapply(  1:(length(bkMap$keys)), FUN=function(i, rep1, rep2, rep3, reptimes){
				outData1 = rep(rep1[i], times=reptimes[i])
				outData2 = rep(rep2[i], times=reptimes[i])
				outData3 = rep(rep3[i], times=reptimes[i])
				ooo = cbind(outData1, outData2, outData3)
			},
			rep1=bkMap$snpCt, rep2=markers.b, rep3=markers.e, reptimes=bkMap$bkLens)
	
	others.ma =NULL
	for( dd in others){
		others.ma = rbind(others.ma, dd)
	}
	others = matrix(others, nrow=nrow(hap), ncol=3, byrow=TRUE)
	
	hapFrame = cbind(prekey, hap, others.ma)
	colnames(hapFrame)=c("prekey", "block", "hap", "freq", "hapLen", "markers_b", "markers_e")
	
	## need to remove the singletones
	hapFrame = hapFrame[ (hapFrame[,5]!=1), ]
	
	
	genoFrame = data.frame(prekey=rep("ch", bkMap$snpLen*3),
			seq=rep(1:bkMap$snpLen, each=3),
			freq=rep(NA, bkMap$snpLen*3))
	
	single = which(bkMap$snpCt == 1)
	if(length(single)>0){
		
		## find singletons.
		for( ss in single ){
			freq.info = bkMap$bks[[ss]]
			allele.freq=freq.info[,3][match(freq.info[,2], alleleCode)]
			geno1 = allele.freq[1]^2
			geno2 = 2*allele.freq[1]*allele.freq[2]
			geno3 = allele.freq[2]^2
			genoFrame[ ((ss-1)*3+1):(ss*3) , 3]=c(geno1, geno2, geno3)
		}
		
		
	}
	if(nrow(hapFrame)>0){
		freq.haponly = freq.build(hap=hapFrame, geno=genoFrame)
	}else{
		freq.haponly = freq.build(hap=NULL, geno=genoFrame)
	}
	if(ifD) {
		print(hapFrame)
		print(dim(genoFrame))
		print(genoFrame[1:6,])
	}
	return(freq.haponly)
}

freqmap.reconstruct <-
function(data, cols=NULL, loci.ct, is.1digit=FALSE, dig1Code=c(0, 1, 2, 3), dig2Code = c(0, 1, 2), key.prefix="", start.base=1, ...){
   ifD = FALSE

   ## is NA is error used to represent the missing, change it to 0, or the max(allele code)+1
   
   if (is.1digit){
     if(is.null(cols)) cols=c(1, ncol(data))
     data = exchangeDigit(ma=data, cols=cols, dig1Code=dig1Code,
       dig2Code = dig2Code, action=c("1to2")) [,(cols[1]: ((cols[2]-cols[1]+1)*2+cols[1]-1))]
     miss.code = dig1Code[1]
     #str(data)
   }else{
     data = data[, (cols[1]:cols[2])]
     miss.code = dig2Code[1]
   }
   
   ## need to construct the haplotype freq file
   ## with hap key, haplotype, haplotype prob 
   ct = length(loci.ct)
   geno.code = as.vector(matrix(dig2Code[2:3], ncol=2) %*% matrix(c(2,0, 1,1, 0,2), ncol=3, byrow=FALSE))

   #if(sum(loci.ct>7)>=1) stop("Cannot process haplotype block with 8 or more loci.")
   if(sum(loci.ct)!=(ncol(data)/2)) {
     if(ifD) print(sum(loci.ct))
     if(ifD) print(ncol(data)/2)
     stop("Dismatching number of loci with the number of column in data.")

   }

   hapMap.info= NULL
   genoMap.info=NULL
   
   startId = start.base
   
   ## process the file if only one block exist
   if(ct==1){

     if(loci.ct>=2){
       hap.re = haplo.stats::haplo.em(geno=data, ...)
       
       hap.prob = hap.re$hap.prob
    
       hap.exp =  hap.re$haplotype
    
       hap.expVec = as.integer(apply(hap.exp, 1, paste, collapse=""))
       #m = match(c("ch", "block", "hap", "freq", "hapLen","markers_b", "markers_e"), 
       df = data.frame( prekey= rep(key.prefix, length(hap.prob)),
                   block =rep(1, length(hap.prob)),
                   hap=as.vector(hap.expVec),
                   freq = hap.prob,
                   hapLen = rep(loci.ct, length(hap.prob)),
                   markers_b=rep(startId, length(hap.prob)),
                   markers_e=rep(startId+loci.ct-1, length(hap.prob))
                  )
       hapMap.info=df
     }

     dd = data[,1]+data[,2]
       
     tb = table(factor(unlist(dd), levels=geno.code))
     tb.freq = tb/sum(tb)

     #ch	seq	freq	genotype
     df = data.frame(prekey= I(rep(key.prefix, 3)),
                   seq =rep(startId, 3),
                   freq = as.vector(tb.freq)
                  )
     genoMap.info=df
     
     
     return(re=list(hapMap.info=hapMap.info, genoMap.info=genoMap.info))
     
   }

    ## process the file if only more than one block exist
   idx.end = cumsum(loci.ct)*2
   idx.start = c(1, (idx.end+1)[-length(loci.ct)])
   cut.rg = cbind(idx.start, idx.end)


   all.df.key2 = NULL
   all.df.oth2 = NULL

   snp.b = startId
   snp.e = startId-1
   tmp1 = 1
   
   allsnp = TRUE
   for( i in 1:ct){
      keys =NULL
      oth=NULL
      snp.b = snp.e+1
      snp.e = snp.b+loci.ct[i]-1
      if(loci.ct[i]!=1){
		  allsnp = FALSE
          ## for Hap
          hap.re = haplo.stats::haplo.em(geno=data[, ((cut.rg[i,1]):(cut.rg[i,2]))], ...)

          # print(hap.re)
          hap.prob = hap.re$hap.prob
       
          hap.exp  = hap.re$haplotype
       
          hap.expVec = as.integer(apply(hap.exp, 1, paste, collapse=""))

          keys = rep(key.prefix, length(hap.prob))
          
          oth =cbind(block =  rep(i, length(hap.prob)),
                   hap=hap.expVec,
                   freq = hap.prob,
                   hapLen = rep(loci.ct[i], length(hap.prob)),
                   markers_b=rep(snp.b, length(hap.prob)),
                   markers_e=rep(snp.e,length(hap.prob))
                  )

          all.df.key2 = c(all.df.key2, keys)
          all.df.oth2 = rbind(prekey=all.df.oth2, oth)
          #print(all.df.oth2[1:3,])
          
      } #if(loci.ct[i]!=1){
    } #for( i in 1:ct){

	if(!allsnp){
	    rownames(all.df.key2)=NULL
	    rownames(all.df.oth2)=NULL
	    hapMap.info = data.frame(prekey=I(all.df.key2), all.df.oth2)
	}else{
		hapMap.info = NULL
	}
 
    ## gen one genomap for every SNP
    dd.t = data[, (min(cut.rg)):(max(cut.rg))]
    # convert into two columns
    dd.ma = matrix(unlist(dd.t), ncol=2, byrow=FALSE)
    dd = dd.ma[,1]+dd.ma[,2]

    ## remove NA first
    #print(table(dd))
    ## convert into n columns each column is for one SNP
    dd = matrix(dd, nrow=nrow(data), byrow=FALSE)
    tb = apply(dd, 2, FUN= function(col, geno.code){
      tt = table(factor(unlist(col), levels=geno.code))
      tt.freq = tt/sum(tt)
      tt.freq
    }, geno.code=geno.code)

    #print(tb[1:3,])

    #ch	seq	freq	genotype
 
    genoMap.info = data.frame(prekey= I(rep(key.prefix,  3*ncol(dd))),
        seq =rep(startId : (startId + cumsum(loci.ct)[length(loci.ct)] -1), each=3),
        freq = as.vector(tb) )

   	return(re=list(hapMap.info=hapMap.info, genoMap.info=genoMap.info))
   
}

genGenoProb <-
function(){
   genoProb = matrix(c(1,1, 1, 0, 0,
                       2,1, 0, 0, 1,
                       3,1,.5, 0,.5,
                       1,2, 0, 0, 1,
                       2,2, 0, 1, 0,
                       3,2, 0,.5,.5,
                       1,3,.5, 0,.5,
                       2,3, 0,.5,.5,
                       3,3,.25,.25, .5), ncol=5, byrow=TRUE)
   ## test
   ## apply(genoProb[,3:5], 1, sum)
   ## apply(genoProb[,3:5], 2, sum)
   return (genoProb)
}

genMatingTBCondOnChild3 <-
function(appVarNames, child, par1, par2, prob1, prob2, logF = NULL, job=1){

  fStr ="[genMatingTBCondOnChild3]:"
  ifD = FALSE
  maxTbl = NULL
  if(ifD) print(par1)
  if(ifD) print(par2)
  tryCatch({
    maxTbl = get(appVarNames$tbl, envir=qingEnv)
    maxRow = nrow(maxTbl)
  }, error=function(e){
    errTrace = paste(e, collapse=";", sep="")
    stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$freqMap, " does not exisit."))
  })
  
  tryCatch({
    maxMateTbl = get(appVarNames$mateTbl, envir=qingEnv)
  }, error=function(e){
    errTrace = paste(e, collapse=";", sep="")
    stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$freqMap, " does not exisit."))
  })

  childPairNo = nrow(child)

  if(is.null(childPairNo)) {
    child = matrix(child, ncol=2)
    childPairNo = 1
  }

  if(is.null(dim(par1))){
    par1 = matrix(par1, ncol=2)
  }
  if(is.null(dim(par2))){
    par2 = matrix(par2, ncol=2)
  }
  #print("OOOOOLLLDDDDD")
  
  #print(par1)
  #print(par2)
  
  ## get the most probable ones, order the pair by their prob
  par1 = par1[order(prob1,decreasing=TRUE),,drop=FALSE]
  par2 = par2[order(prob2,decreasing=TRUE),,drop=FALSE]

  prob1 = sort(prob1,decreasing=TRUE)
  prob2 = sort(prob2,decreasing=TRUE)
  
  #print(par1)
  #print(par2)
  #print(prob1)
  #print(prob2)
  
  ## would assume that all pairs with different indexes will come at the beginning of the list
  ## and pairs with same indexes will come at the end of the list

  par1short = nrow(par1)<=nrow(par2)

  if(par1short){
    baseNo = nrow(par1)
    othNo = nrow(par2)
    basePair = par1
    othPair = par2
    baseCol = 1:2
    othCol = 3:4
    mateCol = 1
    
  }else{
    baseNo = nrow(par2)
    othNo = nrow(par1)
    basePair = par2
    othPair = par1
    baseCol = 3:4
    othCol = 1:2
    mateCol = 2
  }

  probVec = rep(NA, times=maxRow)
  counter = 0
  ## iterate choice of one parent and the child
  for ( i in 1:baseNo){
    if(ifD) print(paste("i=", i, " out of ", baseNo))
    
    childMeet = apply(child, 1, FUN = function(rowItem, bench){
        re = as.logical(max(is.element(rowItem, bench)))
    }, bench = basePair[i,])

    childSelRowIdx = NULL
    if( sum(childMeet) >=1 ){
        childSelRowIdx = (1:childPairNo)[childMeet]
        for(j in childSelRowIdx){
             ## if two hap in child is the same
             sameHapChild = child[j,1]==child[j,2]
             if(sameHapChild) {
                par2Meet = apply(othPair, 1, FUN = function(rowItem, bench){
                   re = as.logical(max(is.element(rowItem, bench)))
                }, bench = child[j,1])               

             }else{

                childMatch = is.element(child[j,], basePair[i,])
                if(sum(childMatch)==2){
                  othHapReq = child[j,]
                }else if(sum(childMatch)==1){
                  othHapReq = child[j,][!childMatch]
                }else{
                  stop ("programming error")
                }
                if(ifD) print(childMatch)

                par2Meet = apply(othPair, 1, FUN = function(rowItem, bench){
                   re = as.logical(max(is.element(rowItem, bench)))
                }, bench = othHapReq)                  
             }

             if(ifD) print(par2Meet)
             if( sum(par2Meet) >=1 ){
               par2SelRowIdx = (1:othNo)[par2Meet]
               
               if(ifD)  print(paste("par2SelRowIdx=", par2SelRowIdx))
               addedRow = length(par2SelRowIdx)
               if(ifD) print(paste("addedRow=", addedRow))
               if( addedRow != 0){
                  if( (counter + addedRow) > maxRow ) {
                       #rm(maxTbl)
                       #gc()
                       # stop (paste("\n", fStr, "List length (", counter + addedRow , ") exceed the maximum (", maxRow, ") for i =", i, se
                       simpleResult = matTBCondOnChild3.finalProc(maxTbl=maxTbl, counter=counter, maxMateTbl=maxMateTbl,
                               probVec=probVec, prob1=prob1, prob2=prob2, job=job, logF=logF  )
                       return(simpleResult)

                  }
                  tmpPrb = rep(.25, times=addedRow)
                  othSameHap = NULL
                  othSameHap = othPair[par2SelRowIdx,1]==othPair[par2SelRowIdx,2]
  
                  if(ifD)  print(sum(othSameHap))
                  if(sum(othSameHap)>0) {
                    # if other parent are homo, double =.5
                    tmpPrb[othSameHap]=tmpPrb[othSameHap]*2
                  }
                  # if base parent are homo, double for the second time =1
                  if(basePair[i,1]==basePair[i,2]) tmpPrb = tmpPrb*2

                  ## depend on the child is heto, same pairs for parents and child has different prob than different pairs for parents
                  ## parent must be same heto, homo kid will have .25, heter kid will have .5 
                  if( (sum(othSameHap)<addedRow) && (basePair[i,1]!=basePair[i,2]) && ( child[j,1]!=child[j,2])  ){
                    ## 1st rule, there are some hapPairs not the same for the other parent
                    ## 2nd rule, hapPairs for base parent are not the same
                    ## 3nd rule, hapPairs for child are not the same 
                    tmpMaRow = (1:addedRow)[!othSameHap]
                    othPairMa = othPair[par2SelRowIdx[tmpMaRow],,drop=FALSE]

                    ## compare the hapPairs for the other parents with the hapPairs for the base parents.
                    rowTmpMeet = apply(othPairMa, 1, FUN = function(rowItem, bench){
                      re = as.logical(sum(is.element(rowItem, bench))==2)
                    }, bench = basePair[i,])
                    if(sum(rowTmpMeet)>0) tmpPrb[tmpMaRow[rowTmpMeet]]=.5
                  }
                  
                  probVec[(counter+1):(counter+addedRow)]= tmpPrb
                  
                  ##print(matrix(.Internal(rep(par1[i,], addedRow)), ncol=2, byrow=TRUE))
                  ##print( matrix(par2[par2Meet, ], ncol=2))
                  maxTbl[(counter+1):(counter+addedRow), baseCol] = matrix(rep(basePair[i,], addedRow), ncol=2, byrow=TRUE)
                  maxTbl[(counter+1):(counter+addedRow), othCol] = othPair[par2SelRowIdx, , drop=FALSE]
                  ## newly added
                  maxTbl[(counter+1):(counter+addedRow), 5:6] = matrix(rep(child[j,], addedRow), ncol=2, byrow=TRUE)

                  if(ifD) print(maxTbl[(counter+1):(counter+addedRow),])
                  maxMateTbl[(counter+1):(counter+addedRow), mateCol] = rep(i, addedRow)
                  maxMateTbl[(counter+1):(counter+addedRow), 3-mateCol] = par2SelRowIdx
         
                  counter = counter + addedRow
                  if(ifD)print(paste("counter=", counter))
               }
           } ## if( sum(par2Meet) >=1 ){

        } ## for(j in 1:childSelRowIdx){
    } ## if( sum(childMeet) >=1 ){
  }
  simpleResult = matTBCondOnChild3.finalProc(maxTbl=maxTbl, counter=counter, maxMateTbl=maxMateTbl,
                               probVec=probVec, prob1=prob1, prob2=prob2, job=job, logF=logF  )
  return(simpleResult)
}

geno.2dStr2BinaMa <-
function(expStr, subjectCt=length(expStr), snpLen){

  samplesBina = hapBk2AlleleSeq(subjects=expStr, subjectCt=subjectCt, snpLen=snpLen*2)
  bina1= samplesBina[, seq.int(from=1,to=snpLen*2, by=2), drop=FALSE]
  bina2 = samplesBina[, seq.int(from=2,to=snpLen*2, by=2), drop=FALSE]
   binaNew1 = pmin(bina1, bina2)
   binaNew2 = pmax(bina1, bina2)
   bina = util.matrix.col.shuffle2(binaNew1, binaNew2)

   if(subjectCt==1) bina=matrix(bina, nrow=1)
   return(bina)  

}

geno2hapFreq <-
function(prekey, block, genoFreq, snpBase=1, snpSeq){
	df = data.frame( key=I(rep(paste(prekey, block, sep="-"), 2)),
			hapF.ch = rep(prekey, 2),
			hapF.block = rep(block, 2),
			hapF.hap = c(snpBase, snpBase+1),
			freq = c(genoFreq[1]+genoFreq[2]/2, 1- (genoFreq[1]+genoFreq[2]/2) ),
			hapF.hapLen = rep(1, 2),
			hapF.markers_b = rep(snpSeq, 2),
			hapF.markers_e = rep(snpSeq, 2))
	
	return(df)
}

getBackParentGeno <-
function(trioDf, famCol=1, memCol=2, snpIdx=NULL, prefix=NULL, re.child=FALSE){

    if(is.null(snpIdx)){
      snpIdx = (1:ncol(trioDf))[-c(famCol, memCol)]
    }
      ## HARD CODE!!!HARD CODE: know the last three digits are for covariate values
     rowCt = nrow(trioDf)
     trioSNP=trioDf[, c(famCol, memCol, snpIdx)]

     if(re.child) {
       selMa = seq.int(from=3, to=rowCt, by = 3)

     }else{
       ## we know the first two rows in every three rows are parents
       selSeq = seq.int(from=1, to=rowCt, by = 3)
       selSeq2 = selSeq + 1
       selMa = matrix(c(selSeq, selSeq2), ncol=2, byrow = FALSE)
       selMa = as.vector(t(selMa))
     }

     #print(length(selMa))
     trioSNP.id = paste(trioSNP[selMa, 1], trioSNP[selMa, 2], sep="-")
     unique.par = unique(trioSNP.id)

     unique.par.idx = match(unique.par, trioSNP.id)

     #print(length(unique.par))
     
     if(is.null(prefix)){
       if(re.child) {
         re = trioSNP[unique.par.idx,]
       }else{
         re = trioSNP[unique.par.idx,]
       }
       return(re)
     }
     
     if(re.child) {
       write.table(trioSNP[unique.par.idx,], file=paste(prefix, "_child.txt", sep=""), sep=" ",
                   append = FALSE, row.names = FALSE, col.names = FALSE)
     }else{
       write.table(trioSNP[unique.par.idx,], file=paste(prefix, "_parent.txt", sep=""), sep=" ",
                   append = FALSE, row.names = FALSE, col.names = FALSE)
     }
     return(NULL)
}

getHapProb.semiMapFrame <-
function( semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = FALSE
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb)
	
	lef = leftOverProb/(2^snpLen-commonProbIdx)
	allProb = rep(lef, length=2^snpLen)
	
	allProb[mappedAugIdx]=mapProb
	
	return(allProb)
}

getHapProb2 <-
function(selIdx, semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, restandard=FALSE){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = FALSE
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb)
	
	findMatchMajor = match(selIdx, mappedAugIdx)
	selMajor = findMatchMajor[!is.na(findMatchMajor)]
	if(ifD) print(semiMapFrame)
	if(ifD) print(selMajor)
	
	selMajor.ct = sum(!is.na(findMatchMajor))
	selMinor.ct = sum(is.na(findMatchMajor))
	
	if(commonProbIdx == (2^snpLen)){
		## if the map is complete
		
		if(selMinor.ct>0){
			stop (paste("non existing selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
		}else if(selMajor.ct>0){
			newProb = mapProb[ selMajor ]
		}else {
			stop (paste("not valide selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
		}
	}else{
		## if the map is incomplete
		
		augProb = c(mapProb, leftOverProb/(2^snpLen - commonProbIdx))
		idxMatched = selMajor
		if(selMinor.ct>0){
			validIdx = selIdx[is.na(findMatchMajor)]
			tt = (validIdx>=1 & validIdx <=2^snpLen)
			if(sum(tt)!=length(validIdx)){
				stop (paste("non existing selIdx=[", paste(selIdx, collapse=";", sep=""), "]", sep=""))
			}
			
			idxMatched = findMatchMajor
			idxMatched[is.na(findMatchMajor)] = commonProbIdx+1
			
		}
		newProb = augProb[idxMatched]
	}
	
	if(restandard){
		newProb = newProb/sum(newProb)
	}
	
	return(newProb)
}

getP.logodds <-
function(logodds){
	return (exp(logodds)/(1+exp(logodds)))
}

grp.CI <-
function (maUp, maLow, position, barLen, col=NULL,...) 
{
    vLen = length(position)
    lapply(1:vLen, FUN = function(i, up, low, xPos, barLen, col,...) {
        segments((xPos - barLen/2)[i], low[i], (xPos + barLen/2)[i], 
            low[i], col,...)
        segments(xPos[i], low[i], xPos[i], up[i],col, ...)
        segments((xPos - barLen/2)[i], up[i], (xPos + barLen/2)[i], 
            up[i],col, ...)
    }, up = maUp, low = maLow, xPos = position, barLen = barLen, col)
    return(NULL)
}

grp.kmStep <-
function(timeVec, probVec, lty=1, lwd=1,  col ="black", plot.dot=FALSE){

  len = length(timeVec)

  for( i in 1:(len-1)){
    segments(timeVec[i], probVec[i], timeVec[i+1], probVec[i], col=col, lty=lty, lwd=lwd)
    segments(timeVec[i+1], probVec[i], timeVec[i+1], probVec[i+1], col=col, lty=lty, lwd=lwd)
    if(plot.dot) points(timeVec[i+1], probVec[i], pch=1, col="black")
  }

  return(NULL)
}

grp.palette <-
function(name, width = 500, height =500){
     jpeg(filename = name, width = width, height = height,
          pointsize = 10, quality = 100, bg = "white", res = NA)
     
     
     recL = 1
     colMa = matrix(colors(), ncol = 25, nrow = 27, byrow = TRUE)
     
     my.xlim = c(0, 25)
     my.ylim = c(0, 27)
     plot(my.xlim, my.ylim, type = "n", xlab="", ylab="", axes=FALSE)
     
     
     for(row in 1:27){
       for(col in 1:25){
         rect( (col-1)* recL + .1, (row-1)* recL + .1, col*recL - .1, row*recL - .1, col = colMa[row, col] )
         #print(paste("x=", (col-1)* recL + .1, ", y=", (row-1)* recL + .1))
         #print(paste("col=", col, "row=", row))
         
       }
     
     }
     axis(2, seq(.5,27.5, by=2), seq(0,670, by=50), cex=.5)
     axis(1)
     grid(nx = NULL, ny = NULL, col = colors()[31], lty = "dotted",
          lwd = NULL, equilogs = TRUE)
     
     dev.off()
     return(NULL)
}

hap.2geno <-
function(hapExp, hapProb, snpBIdx=NA){

    if(is.na(snpBIdx)) snpBIdx = 0
    bina1 = hapBk2AlleleSeq(subjects=hapExp,
                            subjectCt=length(hapExp),
                            snpLen=nchar(hapExp[1]), markdownOne = FALSE)
    
    geno1dFreq = lapply(1:ncol(bina1), FUN=function(col, prob, snpMa){
      snpF = factor(snpMa[,col],levels=1:2)
      
      gFreq = tapply( prob, snpF, sum)
      #print(gFreq)
      gFreq[ is.na(gFreq) ] = 0
      gFreq1 = c(gFreq, 2*gFreq[1])
      gFreq2 = c(gFreq, gFreq[2])
      re = gFreq1*gFreq2
      #print(re)
    }, prob=hapProb  , snpMa = bina1)
    
    genoFreq = util.list.2matrix(list=geno1dFreq, byRow = TRUE)
    colnames(genoFreq)=c("11","22","12")
    rownames(genoFreq)=paste("snp", snpBIdx+1:ncol(bina1), sep="")
    return(genoFreq)
}

hap2GenoBlock <-
function(hapExp = c("11", "12", "21", "22"), hapProb = rep(1/length(hapExp), length(hapExp)), snpCt = nchar(hapExp[1]), alleleCode = c(1, 2)){

     ## first standardize??
     
     ct = length(hapExp)
#      hapIdxCorner=NULL
#      ## 2 hap combination to form the diplotype, identify each variation by its index
#      for ( i in 1:ct ){
#        for ( j in 1:ct){
#          if(i<j) hapIdxCorner = rbind(hapIdxCorner, c(hap1=i, hap2=j))
#        }
#      }
#      
# 
     ifD = FALSE
     if(ifD) print(paste("hapLen:", ct))
     hapIdxCorner = util.it.smallLargeIdx(ct, keep.same=FALSE)
     hapIdxDiag = matrix(rep(1:ct, times=2), ncol=2, byrow=FALSE)
  
     probCorner = matrix(as.vector(hapProb)[as.vector(hapIdxCorner)], ncol=2, byrow=FALSE)
     probDiag = matrix(as.vector(hapProb)[as.vector(hapIdxDiag)], ncol=2, byrow=FALSE)
     joinCorner = 2*probCorner[,1]*probCorner[,2]
     joinDiag = probDiag[,1]*probDiag[,2]
     
     expCorner = matrix(hapExp[hapIdxCorner], ncol=2, byrow=FALSE)
     expDiag = matrix(hapExp[hapIdxDiag], ncol=2, byrow=FALSE)
     
     hapIdx = rbind(hapIdxCorner, hapIdxDiag)
     colnames(hapIdx) = c("id1", "id2")
     hapExp = rbind(expCorner, expDiag)
     colnames(hapExp) = c("hExp1", "hExp2")

     if(ifD) print(hapExp)
     # get the one-digit coding
     # when add to digit, get c(2,3,4)-1=c(1,2,3) the common SNP coding
     # HARD CODE!!!HARD CODE
     bina1 = hapBk2AlleleSeq(hapExp[,1], subjectCt=dim(hapExp)[1], snpLen=snpCt, markdownOne = FALSE)
     bina2 = hapBk2AlleleSeq(hapExp[,2], subjectCt=dim(hapExp)[1], snpLen=snpCt, markdownOne = FALSE)
      if(ifD) print(bina1)
     
     ##  20Dec08Change!!!: straighten coding: output must be of digit 1, 2, 3 for 1-digitCoding, and 1, 2 for 2-digit
     a1 = sum(alleleCode * c(2,0)) # 11-> to 1
     a2 = sum(alleleCode * c(1,1)) # 12-> to 3
     a3 = sum(alleleCode * c(0,2)) # 22-> to 2

     snp1d = bina1+bina2
     snp1d.f = factor(snp1d, levels=c(a1, a3, a2))

     if(ifD) {
       print(c(a1, a2, a3))
       print(";;")
       print(snp1d)
       print(snp1d.f)
       print(is.na(snp1d.f))
     }
     
     if( max(is.na(snp1d.f))==1) stop("Input allelCode does not match with other input.")
     
     snp1d = as.integer(snp1d.f)
     snp1d = matrix(snp1d, ncol=ncol(bina1), nrow=nrow(bina1))

     # get the two-digit coding
     bina1.f = factor(bina1, levels=alleleCode)
     if( max(is.na(bina1.f))==1) stop("Input allelCode does not match with other input.")
     bina1.f = as.integer(bina1.f)

     bina2.f = factor(bina2, levels=alleleCode)
     if( max(is.na(bina2.f))==1) stop("Input allelCode does not match with other input.")
     bina2.f = as.integer(bina2.f)

     bina1 = matrix(bina1.f, ncol=ncol(bina1), nrow=nrow(bina1))
     bina2 = matrix(bina2.f, ncol=ncol(bina1), nrow=nrow(bina1))
     
     binaNew1 = pmin(bina1, bina2)
     binaNew2 = pmax(bina1, bina2)

     if(ifD) print(binaNew1)
     
     genoTypes.m = util.matrix.col.shuffle2(binaNew1, binaNew2)
     genoExp = util.matrix.cat(genoTypes.m, 1:(snpCt*2), sep="")
      
     prob = matrix(c(joinCorner, joinDiag), ncol=1)

     tb = data.frame(hapIdx, hapExp, prob, genoExp)
     
     return(list(tb=tb, genoMa=genoTypes.m, geno1d = snp1d))
}

hap2genotype.bk <-
function(bkDframe, chCol, blockCol, keyCol, expCol, probCol, hapLenCol){
  hapExp = bkDframe[,expCol]
  prob = bkDframe[,probCol]
  dipExp =expand.grid(hapExp, hapExp)
  probExp = expand.grid(prob, prob)
  varCt =  bkDframe[1, hapLenCol]
  genotype = hap2genotype.m(as.character(dipExp[,1]), as.character(dipExp[,2]), nrow(dipExp), varCt)
  probExp = probExp[,1]*probExp[,2]
  geno = util.matrix.cat( genotype, 1:(2*varCt))
  collaps = tapply(probExp, geno, sum)

  nrow = length(collaps)
  
  re = data.frame(ch=rep(bkDframe[1,chCol], nrow), block=rep(bkDframe[1,blockCol], nrow),
                  exp=names(collaps), varCt=rep(2*varCt, nrow), prob=collaps,
                  key=rep(bkDframe[1,keyCol], nrow))
 ## re = list(probList=collaps, varCt = varCt)
  return(re)
}

hap2genotype.m <-
function(hapMa1, hapMa2, subjectCt, snpLen){
  bina1 = hapBk2AlleleSeq(hapMa1, subjectCt, snpLen, markdownOne = FALSE)
  bina2 = hapBk2AlleleSeq(hapMa2, subjectCt, snpLen, markdownOne = FALSE)
  
  binaNew1 = pmin(bina1, bina2)
  binaNew2 = pmax(bina1, bina2)
  re = util.matrix.col.shuffle2(binaNew1, binaNew2)

  return(re)
}

hapBk2AlleleSeq <-
function(subjects, subjectCt, snpLen, markdownOne=TRUE){

  if(is.null(dim(subjects))){
     subjects.snp = lapply(subjects, FUN=util.str.2CharArray,  len=snpLen)    
  }else{     
     bkCt = dim(subjects)[2]
     
     subjects.hap = util.matrix.cat(subjects, 1:bkCt, sep="")
     subjects.snp = lapply(subjects.hap, FUN=util.str.2CharArray,  len=snpLen)
   }
  
  if(markdownOne){
    subjects.snp = matrix(as.numeric(unlist(subjects.snp))-1, ncol = snpLen, byrow = TRUE)
  }else{
    subjects.snp = matrix(as.numeric(unlist(subjects.snp)), ncol = snpLen, byrow = TRUE)
  }
  return(subjects.snp)
  
}

hapBkGenoMap2HapMap <-
function(hapBkGenoMap, reSimuMap=TRUE){
	
	singleton = hapBkGenoMap$genomeMarkerInfo[hapBkGenoMap$genomeMarkerInfo$qu.type==1, 2 ]
	singletonGeno = hapBkGenoMap$genoOnlyMap$bks[ singleton ]
	newFrame = NULL
	for( i in 1: (length(singletonGeno))){
		g = singletonGeno[[i]]
		newF = geno2hapFreq(prekey=g[1,hapBkGenoMap$genoOnlyMap$chCol] ,
				block =g[1,hapBkGenoMap$genoOnlyMap$blockCol] ,
				genoFreq=g[,hapBkGenoMap$genoOnlyMap$probCol] ,
				snpBase=hapBkGenoMap$hapBkOnlyMap$snpBase,
				snpSeq=singleton[i])
		newFrame=c(newFrame, list(newF))
	}
	
	singletonSeqId = which(hapBkGenoMap$genomeMarkerInfo$qu.type==1)
	allSeqId = 1: (nrow(hapBkGenoMap$genomeMarkerInfo))
	leftSeqId = allSeqId[ -singletonSeqId ]
	
	## bind with old
	allBkMap = c(hapBkGenoMap$hapBkOnlyMap$bks, newFrame)
	## reshuffle the bks
	newBkMap = allBkMap[ match(allSeqId, c(leftSeqId, singletonSeqId)) ]
	
	
	## get new DF
	newDf = matrix(NA, nrow=length(singletonGeno)*2+nrow(hapBkGenoMap$hapBkOnlyMap$df), ncol=7)
	newDf = data.frame(key=I(rep("-", nrow(newDf))), newDf)
	
	rowIdx = 0
	for( j in 1:(length(newBkMap))){
		hapDf = newBkMap[[j]]
		newDf[ (rowIdx+1):(rowIdx+nrow(hapDf)), ]= hapDf
		rowIdx = rowIdx+nrow(hapDf)
	}
	
	#str(newDf)
	if(reSimuMap){
		simuDf = newDf[ , c(1, 4, 5)]
		simuMap = bkMap.constr(data=simuDf, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3)
		return(list( newBkMap=newBkMap, newDf=newDf, simuMap=simuMap))
	}
	
	return(list(  newBkMap=newBkMap, newDf=newDf ))
}

hapGenoBlockProc <-
function(hapGenoSub,  snpCoding=c(0, 1, 2, 3)){
	ifD = FALSE
	if(ifD) print(hapGenoSub)
	
	snpCt = length(hapGenoSub)
	
	homoIn = NULL
	homoDigit = NULL
	missingIn = NULL
	hetoIn = NULL
	for ( i in 1: snpCt ){
		if(ifD) print(hapGenoSub[i])
		matchIndex = which(as.numeric(hapGenoSub[i])==snpCoding)
		if(matchIndex==1){
			missingIn = c(missingIn, i)
		}else if(matchIndex==2 | matchIndex==3){
			homoIn = c(homoIn, i)
			homoDigit = c(homoDigit, snpCoding[matchIndex])
		}else if(matchIndex==4){
			hetoIn = c(hetoIn, i)
		}  
	}
	return(list(homoIn=homoIn, homoDigit = homoDigit, missingIn=missingIn, hetoIn=hetoIn, ori=hapGenoSub))
}

hapPair.match <-
function(listNewOrder, bkMax, pairList){
  ifD=FALSE
  fN = "hapPair.match:"
  if (ifD)   print(paste(fN, " start"))
  if (ifD)   print(listNewOrder)
  if (ifD)   print(bkMax)
  if (ifD)   print(pairList)
    
  bkCt = length(listNewOrder)

  map.row = 2^(bkCt-1)
  if(bkCt ==1) {
    row.ct = nrow(pairList[[1]])
    mgrid1 = matrix(NA, nrow=row.ct, ncol=bkMax)
    mgrid2 = matrix(NA, nrow=row.ct, ncol=bkMax)

    mgrid1[, listNewOrder]= pairList[[1]][,1]
    mgrid2[, listNewOrder]= pairList[[1]][,2]
    #print(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
    return(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
  }else{
    row.ct = cumprod( unlist(lapply(pairList, FUN=nrow)) )[bkCt]
  }

  # map of index for which hap to pull: each column for each index in bkIds,
  idx.map = exhaustHapExp(lociCt=bkCt)$hap
  idx.map = idx.map[1:(nrow(idx.map)/2),,drop=FALSE]

  if(ifD) print(paste("map.row=", map.row, " row.ct=", row.ct))
  
  mgrid1 = matrix(NA, nrow=row.ct*map.row, ncol=bkMax)
  mgrid2 = matrix(NA, nrow=row.ct*map.row, ncol=bkMax)
  
  for(i in 1:map.row){
    if(ifD) print(i)
    # for each row in the map, grap the hap list and do an expand
    grab.hapSet = idx.map[i,]
    ex.list = lapply(1:bkCt, FUN=function(ii, maps, pList){
      ma = pList[[ii]] 
      re = ma[,maps[ii]]
      return(re)
    }, maps=grab.hapSet, pList = pairList)
    ex.grid = qExpandTable(listOfFactor =ex.list, removedRowIdx=NULL, re.row=FALSE)
    if(ifD) print(ex.grid)
    
    row.seq = ((i-1)*row.ct+1) : (i*row.ct) 
    mgrid1[ row.seq,listNewOrder ] = ex.grid

    grab.hapSet = 3-idx.map[i,]
    ex.list = lapply(1:bkCt, FUN=function(ii, maps, pList){
      ma = pList[[ii]] 
      re = ma[,maps[ii]]
      return(re)
    }, maps=grab.hapSet, pList = pairList)
    if(ifD) print(ex.list)
    ex.grid = qExpandTable(listOfFactor =ex.list, removedRowIdx=NULL, re.row=FALSE)
    if(ifD) print(ex.grid)
    
    mgrid2[ row.seq,listNewOrder ] = ex.grid

  }
  #print(cbind(mgrid1=mgrid1, mgrid2=mgrid2))
  return(cbind(mgrid1=mgrid1, mgrid2=mgrid2))

}

HRCB.applyRule <-
function(subjectsExpMa, rule, col.shuffled){
  ifD = FALSE
  len = length(rule)
  if(is.null(dim(subjectsExpMa))){
    nsub = length(subjectsExpMa)
  }else{
    nsub = nrow(subjectsExpMa)
  }
  varCt = length(col.shuffled)

  if(ifD) print("subjectsExpMa:")
  if(ifD) print(subjectsExpMa[1:10])
  subjectListMa = hapBk2AlleleSeq(subjectsExpMa, nsub, varCt, markdownOne=TRUE)
  print(str(subjectListMa))

  subjectList = util.matrix.2list(subjectListMa)
  
  linear = unlist(rep(rule$inter, nsub))
  slist = rule$slist
  
  for ( i in 1:length(slist)){
    curRule = slist[[i]]
    treePred = binaTree.apply(binaTree=curRule$signal, data=subjectListMa, colnames=col.shuffled, re.final=TRUE)
    #print(treePred)
    linear = linear + treePred * unlist(rep(curRule$coef, nsub))
  }

  if(ifD){
    aaaa = cbind(subjectsExpMa, subjectListMa, treePred )
    write.csv(aaaa, file="HRCB.applyRule.diagnostic.csv")

  }
  #print(linear)
  riskProb =  getP.logodds(linear)
  attributes(riskProb)=NULL
  #print(riskProb)
  return(list(riskProb=riskProb, treePred=treePred))
}

HRCB.Esp1Rule.sampleKid <-
function(rule, caseNo, spStrata, supHapProb){
	ifD =FALSE
	fn = "HRCB.Esp1Rule.sampleKid::"
	
	tmpObj = NULL
	if( is.character(spStrata)){
		
		a = load(paste(spStrata, ".RData", sep=""))
		tmpObj = get(a[1])
		
	}else{
		tmpObj = spStrata
	}
	
	
	HRCBStra.A = tmpObj$HRCBStra.A
	HRCBStra.B = tmpObj$HRCBStra.B
	HRCBStra.C = tmpObj$HRCBStra.C
	HRCBStra.D = tmpObj$HRCBStra.D
	
	signalRule = rule$slist[[1]]
	inter = rule$inter
	coef = signalRule$coef
	
	HR = exp(inter+coef)/(1+exp(inter+coef))
	LR = exp(inter)/(1+exp(inter))
	
	
	if(HRCBStra.A$ct>0) {
		risk.A = HRCBStra.A$idProb[,2]*unlist(HR)
	}else{
		risk.A = 0
	}
	if(HRCBStra.B$ct>0) {
		risk.B = HRCBStra.B$idProb[,2]*unlist(HR)
	}else{
		risk.B = 0
	}
	if(HRCBStra.C$ct>0) {
		risk.C = HRCBStra.C$idProb[,2]*unlist(LR)
	}else{
		risk.C = 0
	}
	risk.D = HRCBStra.D$idProb[,2]*unlist(LR)
	
	risk.sumHR = lapply(list(risk.A, risk.B), FUN=sum)
	risk.sumLR = lapply(list(risk.C, risk.D), FUN=sum)
	
	
	# check the sum of prob
	#print(paste(fn, " population risk =", sum(unlist(c(risk.sumHR, risk.sumLR)))))
	
	
	stra.sampled = rmultinom(1, size=caseNo, prob=c(risk.sumHR, risk.sumLR) )
	stra.cumsam = cumsum(stra.sampled)  
	if(ifD){
		print(paste("sampled group:", paste(stra.sampled, collapse="; ")))
	}
	
	
	simMa = matrix(NA, nrow=caseNo, ncol=2)
	## simu for group A
	if(stra.sampled[1]>0){
		#print("a")
		simMa[1:stra.cumsam[1],] = HRCBSpGrp.sp(grp=HRCBStra.A, size=stra.sampled[1], hapProb=supHapProb, riskProb=risk.A)
		#print("a")
	}
	## simu for group B
	if(stra.sampled[2]>0){
		#print("b")
		simMa[(stra.cumsam[1]+1):stra.cumsam[2] ,] = HRCBSpGrp.sp(grp=HRCBStra.B, size=stra.sampled[2], hapProb=supHapProb, riskProb=risk.B)
		
	}
	## simu for group C
	if(stra.sampled[3]>0){
		#print("c")
		simMa[(stra.cumsam[2]+1):stra.cumsam[3], ] =
				HRCBSpGrp.sp(grp=HRCBStra.C, size=stra.sampled[3], hapProb=supHapProb)
		
	}
	## simu for group D
	if(stra.sampled[4]>0){
		#print("d")
		simMa[(stra.cumsam[3]+1):stra.cumsam[4], ] =
				HRCBSpGrp.sp(grp=HRCBStra.D, size=stra.sampled[4], hapProb=supHapProb, grpA =HRCBStra.A )
		
	}
	
	return(simMa)
}

HRCB.Esp1Rule.spTrioOnBase <-
function(bkMap=NULL, preObj=NULL, spStrata, rule, caseNo, ifS="simuInfo",  reControl=FALSE){
  FN = "HRCB.Esp1Rule.spTrioOnBase"
  if(is.null(preObj)){
    if(is.null(bkMap)) stop(" No object is passed as bkMap nor as preObj. The function need at least one to work.")
    warning( "No object is passed as preObj. The function will need to generate the preObj.")
    
    preObj = bkMap.HRCB.Esp1Rule.genoSeq(bkMap, rule, re.probOnly = FALSE)
  }

  #print(FN)
  spStrataObj  = get(spStrata, envir = qingEnv)
  #print(str(spStrataObj))

  
  bkMapS=preObj$bkMapS
  supHapExp=preObj$supHapExp
  supHapProb=preObj$supHapProb

  if(is.null(preObj$supHapExp)) stop("Object, preObj, has the variable $supHapExp as NULL!")
  
  kids = HRCB.Esp1Rule.sampleKid(rule=rule, caseNo=caseNo, spStrata=spStrataObj, supHapProb=supHapProb)
    
#  print("return kids")
#  return(kids)
  parIdx = matrix(NA, nrow = caseNo, ncol=4)
  for( eachChild in 1:caseNo){
     hapPair = kids[eachChild,]
     if (hapPair[1]==hapPair[2]){
       parIdx[eachChild,] = ESp.impuParent.homoKid(hapPair[1], supHapProb)
     }else{
       ## testing  hapPair=c(1,2)
       parIdx[eachChild,] = ESp.impuParent.heterKid(hapPair, supHapProb)
     }
   }

    ## shuffle the trio, so risk groups are mixed
    tmpRandomShu = sample(1:caseNo, size=caseNo, replace=FALSE)

   if(!reControl){
     ## need to combine family and run.
     fam.hapIdx.unsort = rbind(parIdx[,1:2], parIdx[,3:4], kids)
     fam.hapIdx = fam.hapIdx.unsort[t( matrix(1:(caseNo*3), ncol=3, nrow=caseNo, byrow=FALSE)), ]
  
     ## find out the exact string expression
     hap.str1 = supHapExp[fam.hapIdx[,1]]
     hap.str2 = supHapExp[fam.hapIdx[,2]]
  
     ## convert string into digits
     geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*3, snpLen=bkMapS$snpLen, snpCoding=0:3, snpBase=c(0, bkMapS$alleleCode))

     ## shuffle the trio, so risk groups are mixed
     allRowIdx = matrix(1:(3*caseNo), nrow=3, byrow=FALSE) 
     trioShuffledIdx = allRowIdx[,tmpRandomShu]
     geno.FMCMa=geno.FMCMa[trioShuffledIdx,]

     if(!is.null(ifS)) {
       ## check!!!
       supDipIdx = matrix(unlist(apply(fam.hapIdx.unsort, 1, FUN=util.it.triMatch, len=length(supHapProb))),
                    ncol=3, nrow=caseNo, byrow=FALSE)
  
       sample.idx = data.frame(parIdx, kids, supDipIdx)[tmpRandomShu,]
       colnames(sample.idx) = c("hapIdx_f1", "hapIdx_f2", "hapIdx_m1", "hapIdx_m2",
                                "hapIdx_c1", "hapIdx_c2", "dipIdx_f", "dipIdx_m", "dipIdx_c")
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }
   }else{ # if(!reControl){
     ## need to find out which one is for the affected child
     othKids =  matrix(NA, nrow = caseNo, ncol=6)
     for( eachChild in 1:caseNo){
       find.idx1 = which( kids[eachChild, 1] == parIdx[eachChild, 1:4])
       if(length(find.idx1)==0) stop("Affected child hap indexes do not match with parents'.")
       find.idx2 = which( kids[eachChild, 2] == parIdx[eachChild, 1:4])
       if(length(find.idx2)==0) stop("Affected child hap indexes do not match with parents'.")
   
       if(max(find.idx1)<=2){
         # idx1 can only come from F, then idx2 must come from M
         fix.idx1 = find.idx1[1]
         fix.idx2 = find.idx2[ find.idx2>=3 ][1]
         comp.idx1 = (1:2)[-fix.idx1]
         comp.idx2 = (3:4)[-(fix.idx2-2)]
       }else{ ## if(max(find.idx1)<=2){
         if( min(find.idx1)>=3){
           # idx1 can only come from M, then idx2 must come from F
           fix.idx1 = find.idx1[1]
           fix.idx2 = find.idx2[ find.idx2<=2 ][1]
           comp.idx1 = (3:4)[-(fix.idx1-2)]
           comp.idx2 = (1:2)[-fix.idx2]
         }else{ ### if( min(find.idx1)>=3){
           # idx1 can come from either F or M, then need to see where idx2 come from
           if(max(find.idx2)<=2){
             # idx2 must come from F, then idx1 come from M
             fix.idx1 = find.idx1[ find.idx1>=3 ][1]
             fix.idx2 = find.idx2[1]
             comp.idx1 = (3:4)[-(fix.idx1-2)]
             comp.idx2 = (1:2)[-fix.idx2] 
           }else{ #### if(max(find.idx2)<=2){
#             if(min(find.idx2)>=3){
#               # idx2 must come from M, then idx1 come from F
#               fix.idx1 = find.idx1[ find.idx1<=2 ][1]
#               fix.idx2 = find.idx2[1]
#               comp.idx1 = (1:2)[-fix.idx1]
#               comp.idx2 = (3:4)[-fix.idx2]
#             }else{
               # set idx1 come from F, then idx2 come from M
               fix.idx1 = find.idx1[ find.idx1<=2 ][1]
               fix.idx2 = find.idx2[ find.idx2>=3 ][1]
               comp.idx1 = (1:2)[-fix.idx1]
               comp.idx2 = (3:4)[-(fix.idx2-2)]               
#             }
           } #### if(max(find.idx2)<=2){
         } ### if( min(find.idx1)>=3){
       } ## if(max(find.idx1)<=2){
       # fill in the other three
       
       if(  max( is.na(c(fix.idx1, comp.idx2, comp.idx1, comp.idx2, comp.idx1, fix.idx2)))==1 ) print("##############")
       othKids[eachChild, ] =  parIdx[eachChild,1:4] [c(fix.idx1, comp.idx2, comp.idx1, comp.idx2, comp.idx1, fix.idx2)]
       #print( c(parIdx[eachChild,1:4], kids[eachChild,], othKids[eachChild,]) )
       
     } # for( eachChild in 1:caseNo){

     ## need to combine family and run.
     fam.hapIdx.unsort = rbind(parIdx[,1:2], parIdx[,3:4], kids, othKids[,1:2], othKids[,3:4], othKids[,5:6])
     fam.hapIdx = fam.hapIdx.unsort[t( matrix(1:(caseNo*6), ncol=6, nrow=caseNo, byrow=FALSE)), ]
  
     ## find out the exact string expression
     hap.str1 = supHapExp[fam.hapIdx[,1]]
     hap.str2 = supHapExp[fam.hapIdx[,2]]
  
     ## convert string into digits
     geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*6, snpLen=bkMapS$snpLen, snpCoding=0:3, snpBase=c(0, bkMapS$alleleCode))

     ## shuffle the trio, so risk groups are mixed
     allRowIdx = matrix(1:(6*caseNo), nrow=6, byrow=FALSE) 
     trioShuffledIdx = allRowIdx[,tmpRandomShu]
     geno.FMCMa=geno.FMCMa[trioShuffledIdx,]

     if(!is.null(ifS)) {
       ## check!!!
       supDipIdx = matrix(unlist(apply(fam.hapIdx.unsort, 1, FUN=util.it.triMatch, len=length(supHapProb))),
                    ncol=6, nrow=caseNo, byrow=FALSE)
  
       sample.idx = data.frame(parIdx, kids, supDipIdx)[tmpRandomShu,]
       colnames(sample.idx) = c("hapIdx_f1", "hapIdx_f2", "hapIdx_m1", "hapIdx_m2",
                                "hapIdx_c1", "hapIdx_c2", "dipIdx_f", "dipIdx_m", "dipIdx_c", "dipIdx_cA", "dipIdx_cB" ,"dipIdx_cC" )
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }
   } # if(!reControl){
  
     
   return(geno.FMCMa)

}

HRCB.famMap.spTrio <-
function(caseNo, matingTbName, ifS="simuDirectInfo", reControl =FALSE ){
     ifD = FALSE

     FN = "HRCB.famMap.spTrio" 
     #print(FN)
     matingTbInfo  = get(matingTbName, envir = qingEnv)
     #print(str(matingTbInfo))

     
     matRowCt = matingTbInfo$matRowCt
     #print(qp("matRowCt=", matRowCt*4))
     #print(str(matingTbInfo))
     
     kids.risk = matingTbInfo$kids.risk
     matingRowIdx = matingTbInfo$matingRowIdx
     genoMap = matingTbInfo$hap2genoMap
     snpLen = matingTbInfo$snpLen
     kids.matchedRow = matingTbInfo$kids.matchedRow

     ## sample row
     ##caseNo = 5
     kid.idx.sample = sample( 1:(matRowCt*4), size = caseNo, prob=kids.risk, replace=TRUE)

     ## find the matrix idx for the kids
     row.sample = kid.idx.sample%%matRowCt
     row.sample[row.sample==0]=matRowCt
     case.idx.matchedRow = kids.matchedRow[kid.idx.sample]
     
     if(reControl){
       ## need to calculate the control id
       ## recreat the child id for sampled rows
       controlIdx = unlist(lapply(1:caseNo, FUN=function(i, totalRow, rowIdx, cIdx){
             childIDX = rep(rowIdx[i], 4)+ (0:3)*totalRow
             leftControl = childIDX[ childIDX!=cIdx[i] ]
           }, totalRow = matRowCt, rowIdx=row.sample, cIdx=kid.idx.sample))

       tt.newIdxSeq = cbind( case.idx.matchedRow, matrix( kids.matchedRow[controlIdx], ncol=3, byrow=TRUE))
       genoOth = genoMap[,6][ tt.newIdxSeq  ]
       geno.CC = genoOth[t( matrix(1:(caseNo*4), ncol=4, byrow=FALSE) )   ]
       geno.CCMa = t(sapply(geno.CC, FUN=geno.2dStr2BinaMa, subjectCt=caseNo*4, snpLen=snpLen))
       
     }


     ## generate parents, kids, pesudo controls
     fa.idx.matchedRow = matingRowIdx[,1][row.sample]
     ma.idx.matchedRow = matingRowIdx[,2][row.sample]
     

     if(!is.null(ifS)) {
       sample.idx = data.frame(kid.idx.sample, row.sample, fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow)
       colnames(sample.idx) = c("idx_kids", "rowIdx_matingTb", "fRowIdx_genoTb",
                                 "mRowIdx_genoTb", "cRowIdx_genoTb")
       write.csv(sample.idx,  file=paste(ifS, "supHap.csv", sep=""))
     }

     ## convert string into digits
     hap.str1 = genoMap[,3][ c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow) ]
     hap.str2 = genoMap[,4][ c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow) ]

     geno = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=caseNo*3, snpLen=snpLen, snpCoding=0:3, snpBase=c(0, 1, 2))

     #geno = genoMap[,6][c(fa.idx.matchedRow, ma.idx.matchedRow, case.idx.matchedRow)]

     geno.FMCMa = geno[ t( matrix(1:(caseNo*3), ncol=3, byrow=FALSE) ),   ]
     if(ifD) print(matrix(geno.FMCMa, ncol=1)[1:10,])
     ## convert string into digits
     # geno.FMCMa = t(sapply(geno.FMC, FUN=geno.2dStr2BinaMa, subjectCt=caseNo*3, snpLen=snpLen))
   
     if(!reControl){
        return(geno.FMCMa)
     }

     
    stop ("NOT implemented yet")
    return(list(cc=geno.CCMa, geno.FMCMa=geno.FMCMa))

}

HRCBSpGrp.cons <-
function(hapProb, ids, type="A"){
	ifD = FALSE
	## contain the first index and the haplotype probabilities for hap1, hap2 and joint
	idProb = NULL
	id2 = NULL
	
	if (type == "A"){
		## for type A, the ids are two column indexes
		## need to find the list of second indexes for the same first index
		tmp.row = nrow(ids)
		if (tmp.row==0) {
			grpA = list(type="A",  ct=0)
			
		}else{
			
			idsAll = matrix(NA, nrow=tmp.row*2, ncol=2)
			idsAll[1:tmp.row,]=ids
			idsAll[(1:tmp.row)+tmp.row, 1]= ids[,2]
			idsAll[(1:tmp.row)+tmp.row, 2]= ids[,1]
			
			new.order = order(idsAll[,1])
			ids = NULL
			ids = idsAll[new.order,]
			
			uniqueIdx1 = unique(ids[,1])
			if(ifD) {
				print("unique idx:")
				print(uniqueIdx1)
			}
			beginIdx = match(uniqueIdx1, table=ids[,1])
			
			## rely on the natural order of ids[,2]
			id2 = ids[,2]
			ids = NULL
			
			senCt = length(uniqueIdx1)
			idProb = matrix(NA, nrow=senCt, ncol=3)
			idProb[,1]=uniqueIdx1
			idProb[,3]=beginIdx
			
			tmp = matrix(NA, nrow=senCt, ncol=2)
			tmp[,1]=hapProb[ uniqueIdx1 ]
			
			
			## to avoid processing too many if, leave the last one out for now
			
			tmp[1:(senCt-1), 2] = sapply(1:(senCt-1), FUN=function(i, idxcut, allmatch, prob){
						
						prob.idxR = (idxcut[i]):(idxcut[i+1]-1)
						prob.sum = sum(prob[  allmatch[ prob.idxR ]  ])
						return(prob.sum)
					}, idxcut=idProb[,3], allmatch=id2, prob=hapProb)
			
			## append the last one
			tmp[senCt, 2]= sum(hapProb[ id2[  (idProb[senCt,3]:length(id2))  ] ])
			
			idProb[,2]= tmp[,1]*tmp[,2]
			
			idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
			
			if(max(idProb[,2]<0)==1) stop("Sampling probability for strata A is wrong.")
			
			grpA = list(type="A", idProb=idProb, id2=id2, ct=senCt)
		}
		
		## also get D
		idProb=NULL
		len = length(hapProb)
		
		idProb = matrix(NA, nrow = len, ncol=2)
		idProb[,1] = 1:(len)
		
		tmp2 = rep(0, times=len)
		
		## if the first id is found in group A, we need to take out the probabilities
		if (grpA$ct > 0)  tmp2[ grpA$idProb[,1] ] = tmp[,2]
		
		idProb[,2] = hapProb*(1-hapProb-tmp2)
		tmp=NULL
		tmp2=NULL
		
		idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
		if(max(idProb[,2]<0)==1) stop("Sampling probability for strata D is wrong.")
		
		grpD = list(type="D", idProb=idProb, ct = len)
		
		return(list(grpA=grpA, grpD=grpD))
	}
	
	if (type == "B" | type=="C"){
		if(length(ids)==0) return(list(type=type, ct = 0))
		
		idProb = matrix(NA, nrow=length(ids) , ncol=2)
		# for type B or C, the ids are one column for homogenuous pair
		idProb[,1] = ids
		idProb[,2] = hapProb[ids]^2
		
		idProb[  (idProb[,2] > ((-1)*10^(-15))) &  (idProb[,2] < 0), 2]=0
		
		if(max(idProb[,2]<0)==1) stop("Sampling probability for strata B/C is wrong.")
		
		
		return(list(type=type, idProb=idProb, ct=length(ids)))
	}
	
}

HRCBSpGrp.sp <-
function(grp, size, hapProb, riskProb=NULL, grpA=NULL){
	ifD = FALSE
	samMa = matrix(NA, nrow=size, ncol=2)
	#print(grp)
	#print(hapProb)
	
	
	straCt = nrow(grp$idProb)
	stra.sampled=NULL
	if(grp$type=="A" | grp$type=="B"){
		if (is.null(riskProb))   riskProb=grp$idProb[,2]
		stra.sampled = rmultinom(1, size=size, prob= riskProb )
	}else{
		stra.sampled = rmultinom(1, size=size, prob= grp$idProb[,2])
	}
	if(ifD) print(grp$type)
	if(ifD) print("strata ct")
	if(ifD) print(stra.sampled)
	
	#print("sampled")
	tmp.ma = cbind(1:straCt, stra.sampled)
	tmp.ma = tmp.ma[tmp.ma[,2]>0, ,drop=FALSE]
	#print(tmp.ma)
	
	samMa[,1] =  unlist(apply(tmp.ma, 1, FUN = function(row, ma){
						rep(ma[row[1]], times=row[2])}, ma=grp$idProb[,1]))
	
	#print(samMa)
	
	if (grp$type=="A"){
		#id2 = grp$id2
		samMa[,2] = unlist(apply(tmp.ma, 1, FUN=function(row, id2, hapProb, idxct ){
							#print(hapProb)
							if(ifD) print(row)
							if(row[1]!=length(idxct)){
								tmp.range = idxct[row[1]]:(idxct[row[1]+1]-1)
							}else{
								tmp.range = idxct[row[1]]:length(id2)
							}
							tmp.sam = id2 [tmp.range]
							
							id2.sel = resample(tmp.sam, size=row[2], prob = hapProb[tmp.sam], replace=TRUE)
							#print(id2.sel)
							if(length(id2.sel)==1) id2.sel = rep(id2.sel, times=row[2])
							#print(id2.sel)
							return(id2.sel)
							
						},  id2 = grp$id2, hapProb=hapProb, idxct = grp$idProb[,3]))
		if(ifD) print(samMa)
		samMa = cbind(pmin(samMa[,1], samMa[,2]), pmax(samMa[,1], samMa[,2]))
	}
	if (grp$type=="B" | grp$type=="C"){
		samMa[,2] = samMa[,1]
	}
	
	if (grp$type=="D"){
		## sample the second idx by rejection
		if(is.null(grpA)) stop("HRCBSpGrp.sam::error")
		samMa[,2] = unlist(
				apply(tmp.ma, 1, FUN = function(row, grpAid, idxct, id2){
							maByIdx = match(row[1], grpAid)
							if( !is.na(maByIdx)){
								if(ifD) print(row)
								if(maByIdx!=length(idxct)){
									tmp.range = idxct[maByIdx]:(idxct[maByIdx+1]-1)
								}else{
									tmp.range = idxct[maByIdx]:length(id2)
								}
								## ask the program to trace up
								tmp.sam = id2 [tmp.range]
								removeList = c(tmp.sam, row[1])
								
								# sample for the entire list, except the one already in A
								idx2 = sam.reject(allIdx=grp$idProb[,1], idxProb = hapProb, rejectIdx = removeList, size=row[2])
							}else{
								# sample for the others, exept the one that is the same as the column 1           
								idx2 = sam.reject(allIdx=grp$idProb[,1], idxProb = hapProb, rejectIdx = row[1], size=row[2])
							}
							return(idx2)
						}, grpAid = grpA$idProb[,1], idxct=grpA$idProb[,3], id2=grpA$id2))
		
		if(ifD) print(samMa)
		samMa = cbind(pmin(samMa[,1], samMa[,2]), pmax(samMa[,1], samMa[,2]))
		if(ifD) print(samMa)
	}  
	
	return(samMa)
	
}

impuBk.scheduler <-
function(raw, idx, job=1, toolname=NULL, freqMaps=NULL, dir="", is.1digit=TRUE, dig1Code, dig2Code, reType=FALSE, reHap=NULL, logF=NULL, logErr=""){

  fStr ="[impuBk.scheduler:]"
  ifD = FALSE

  if(!is.null(dir)){
    if(dir==""){
      print(paste("Imputation information files(s) will be saved under current working directory:",  getwd()))
    }else{
      print(paste("Create directory:", dir, ", where imputation information file(s) are saved.", sep=""))
      dir.create(path=dir, showWarnings = TRUE)
      if (!is.null(logF)){
        logF = file.path(dir, qp(logF, ".txt"))
      }
      if (!is.null(reHap)){
        reHap = file.path(dir, reHap)
      }
      logErr = file.path(dir, logErr)
    }
  }else{
    if(job==1) {
#      print(
#          paste("Value for argument dir is NULL. Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))
#       paste("Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))
      reHap = NULL
      logF = NULL
    }else{
      if (!is.null(logF)){
        logF = qp(logF, ".txt")
      }
#      print(
#          paste("Value for argument dir is NULL. However, users ask to do multiple imputation. ",
#                "Imputation information files(s) will be saved under current working directory:",  getwd()))
    }
  }

  # inside the function, assume 0 1 2 3 for NA, homo, homo, heter and 0 1 2 for NA, allele1, allele2

  ## first check the existance of required parameters
  if(is.null(toolname)){
    if(is.null(freqMaps)){
      stop("Values for both arguments, toolboxName and freqMaps, are missing.")
    }else{
      toolname = toolbox.load(freqMaps=freqMaps)
    }
  }

  ## find the start and end position for the blocks of interest

  all.genomeMarkerInfo = get(toolname$freqMap, envir = qingEnv)$genomeMarkerInfo

  bd.start = all.genomeMarkerInfo[idx[1], 2, drop=TRUE]
  bd.end   = all.genomeMarkerInfo[idx[ length(idx) ], 3, drop=TRUE]

  snpOffset =bd.start-1

  ## change digit to 1 digit coding using Qing's coding scheme
  if(!is.1digit){
      ## excluding the parents cols
      genos  = raw[, 2+( ((bd.start-1)*2+1):(bd.end*2) ),drop=FALSE  ]
      snpNum = ncol(genos)/2
      snp1digit = exchangeDigit(ma=genos,
        cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
  }else{
      ## excluding the parents cols
      genos  = raw[, 2+(bd.start:bd.end), drop=FALSE]
      snpNum = ncol(genos)
      snp1digit = genos
      #dig1Code.inside = dig1Code[c(1, 2, 4, 3)]
      
      ## change to inside Qing's coding scheme
      if(min(dig1Code==c(0,1,3,2))==0){
        snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,3,2))
      }
  }

  trioCt = nrow(snp1digit)/3
  maxRow = trioCt*job

  ## 18 columns: bk, trio, x1, x2, y1, y2, hap_f, hap_m, hap_c1-4
  if(reType){
    imputBkRecord = matrix(NA, ncol = 19, nrow=maxRow)
    impuDummy = imputBkRecord 
  }else{
    imputBkRecord = matrix(NA, ncol = 18, nrow=maxRow)
    impuDummy = imputBkRecord 
  }
  imputBkRecord.ct = 0
  errorTrap = NULL
  
  tryCatchEnv = new.env(parent=qingEnv)
  assign("trapID", 0, envir=tryCatchEnv)
  assign("errorTrap", errorTrap, envir=tryCatchEnv)

  if(!is.null(get(toolname$freqMap, envir = qingEnv)$hapBkOnlyMap)){
	  all.hapIndex = get(toolname$freqMap, envir = qingEnv)$hapIndex
	  
	  hapBkOnlyMap.vars=list()
	  tmp.hapMap =  get(toolname$freqMap, envir = qingEnv)$hapBkOnlyMap
	  hapBkOnlyMap.vars$resiProbCol= tmp.hapMap$resiProbCol
	  hapBkOnlyMap.vars$augIdxCol= tmp.hapMap$augIdxCol
	  hapBkOnlyMap.vars$probCol= tmp.hapMap$probCol

  }else{
	  all.hapIndex = c(-1)
	  hapBkOnlyMap.vars=list()
  }
  snpCoding = 0:3
  snpBase = 0:2
  genoProb = genGenoProb()

    for( unit in idx){
      if (is.element(unit, all.hapIndex)){
        ## haplotype, for each block, search every trio for missingness
        ## find block boundary
        bk.bd = unlist(all.genomeMarkerInfo[unit, c(2,3)])
        bk.geno = snp1digit[, (bk.bd[1]:bk.bd[2])]

        snpCt = bk.bd[2]-bk.bd[1]+1

        exhaustHap = get(toolname$exp, envir = qingEnv)
        exhaustHap = exhaustHap[1:2^snpCt, 1:snpCt]

        bk.genoRowComp = rowSums(bk.geno!=0) == snpCt 
        bk.genoTrioComp = matrix(bk.genoRowComp, ncol=3, byrow=TRUE)
        bk.missTrioIdx = (1:trioCt) [ rowSums(bk.genoTrioComp) < 3 ]
    
        for (famId in bk.missTrioIdx){
          x1 = bk.bd[1]
          x2 = bk.bd[2]
          y1 = (famId-1)*3+1
          y2 = famId*3
          trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
          if(ifD) print( paste(fStr, " processing fam index:", famId))
          if(ifD) print(trioBlock)
          
          tryCatch({
            replace = ESp.imputBlock(appVarNames=toolname, 
                           trioBlock=trioBlock, snpLen=snpCt, bkIdx = unit, job=job,
                           snpCoding = snpCoding, snpBase=snpBase, reType = reType,  logF=logF, 
                           hapBkOnlyMap.vars=hapBkOnlyMap.vars
                            )

            imputBkRecord.ct = imputBkRecord.ct + 1
            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) , 1:6 ] = matrix(
                           rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=TRUE)
            
            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) ,
                          7:ncol( imputBkRecord ) ] = replace
            
            fam.hapIdx1 = replace[1, c(1,3,5)]
            fam.hapIdx2 = replace[1, c(2,4,6)]
            hap.str1 = exhaustHap[fam.hapIdx1, ]
            hap.str2 = exhaustHap[fam.hapIdx2, ]

            ## convert string into digits
            geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=3, snpLen=snpCt)
 
            if(ifD) print(geno.FMCMa)
            if(ifD) print(snp1digit[y1:y2, x1:x2-  snpOffset])

            if( sum(abs(geno.FMCMa - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

            snp1digit[y1:y2, x1:x2-  snpOffset] = geno.FMCMa

            if(!is.null(logF)){
              logl(logF, paste("For trio #", y2/3, ", Choose hap idx=[",
                               paste(replace[1:6], collapse=".", sep=""),
                               "]"))
              logl(logF, paste("Choose hap exp=[",
                    paste( apply(geno.FMCMa, 1, FUN=paste, collapse=""), collapse=".", sep=""), "]",
                               sep=""))
              logl(logF, paste("missing bk data:",
                    paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
      
            }
            if(ifD){
              print("-------------------")
              print(  c(y1, y2, x1, x2)  )              
            }

#            if( (imputBkRecord.ct!=0) & (imputBkRecord.ct %%cutpt==0) & (!is.null(reHap)) ){  
#                  ## if there is record and upto certain number, need to be outputed
#                  ## evaluated after each trio
#                  write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=FALSE],
#                              file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
#                              append = TRUE, row.names = FALSE, col.names=FALSE)
#                  imputBkRecord.ct = 0
#                  imputBkRecord = impuDummy
#             }                       
           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   b = get("trapID", envir=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, envir=tryCatchEnv)
               #trapID <<- trapID+1
	
               ## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
			   #errorTrap <<- rbind(errorTrap, errorInfo)
	           errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
			   a = get("errorTrap", envir=tryCatchEnv)
			   a = rbind(a, errorInfo)
			   assign("errorTrap", a, envir=tryCatchEnv)
			   
               
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
                  
                  #logl(logF, errTrace)
                }         
               }, warn =function(w) {
                 print("ImpuBlock::Warnings")
                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
              
        } ## for (famId in bk.missTrioIdx){
         if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){  
                  ## if there is record and upto certain number, need to be outputed
                  ## evaluated when the block is over
                 write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=FALSE],
                             file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                             append = TRUE, row.names = FALSE, col.names=FALSE)
         }
         imputBkRecord.ct = 0
         imputBkRecord = impuDummy
      }else{
        # print("A genotype")  ## genotype
        bk.bd = unlist(all.genomeMarkerInfo[unit, 2])
        tmp.Map =  get(toolname$freqMap, envir = qingEnv)$genoOnlyMap
        popuProb =  tmp.Map$bks[[bk.bd]][,tmp.Map$probCol]

        snpCt = 1

        bk.geno = snp1digit[,bk.bd -  snpOffset]        
        bk.genoRowComp = as.integer(bk.geno!=0) 
        bk.genoTrioComp = matrix(bk.genoRowComp, ncol=3, byrow=TRUE)
        bk.missTrioIdx = (1:trioCt) [ rowSums(bk.genoTrioComp) < 3 ]

        for (famId in bk.missTrioIdx){
           x1 = bk.bd
           x2 = bk.bd
           y1 = (famId-1)*3+1
           y2 = famId*3
           trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
           if(ifD) print( paste(fStr, " processing fam index:", famId))
           if(ifD) print(trioBlock)
           
           tryCatch({
 
             replace = imputGeno( trioBlock, job=job, genoProb=genoProb,
                    popuProb = popuProb, data.order="FMC", snpCoding=snpCoding)

             #print(replace)
             imputBkRecord.ct = imputBkRecord.ct + 1
             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) , 1:6 ]  =  matrix(
                            rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=TRUE)

             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job),
                           c(7, 9, 11, 13, 15, 17) ] = replace
             
             geno.FMCMa = replace[1, 1:3]
             if( sum(abs(geno.FMCMa - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

             snp1digit[y1:y2, x1:x2 -  snpOffset] = replace[1,1:3]

#              if( (imputBkRecord.ct!=0) & (imputBkRecord.ct %%cutpt==0) & (!is.null(reHap)) ){  
#                   ## if there is record and upto certain number, need to be outputed
#                   ## evaluated after each trio                
#                 write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=FALSE],
#                             file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
#                             append = TRUE, row.names = FALSE, col.names=FALSE)
#                 imputBkRecord.ct = 0
#                 imputBkRecord = impuDummy
#              }
             
           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   
			   b = get("trapID", envir=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, envir=tryCatchEnv)
			   #trapID <<- trapID+1
			   
	           ## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
			   #errorTrap <<- rbind(errorTrap, errorInfo)
			   errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
    		   a = get("errorTrap", envir=tryCatchEnv)
			   a = rbind(a, errorInfo)
			   assign("errorTrap", a, envir=tryCatchEnv)			   
     
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste(trioBlock, collapse="."),
                               sep=""))
                }
               }, warn =function(w) {
                 # print("ImpuBlock::Warnings")                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
         } ## for (famId in bk.missTrioIdx){
         if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){  
                ## if there is record, write out after each singletonn SNP 
                write.table(imputBkRecord[1:(imputBkRecord.ct*job) ,,drop=FALSE],
                            file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                            append = TRUE, row.names = FALSE, col.names=FALSE)
         }
         imputBkRecord.ct = 0
         imputBkRecord = impuDummy
     
      } ## if (is.element(unit, all.hapIndex)){
    } ## for( unit in idx){

	errorTrap=get("errorTrap", envir=tryCatchEnv)
    if(!is.null(errorTrap)){
      write.table(errorTrap, file=paste(logErr, "_errorTrap.csv", sep=""), sep=",",
                 append = FALSE, row.names = FALSE, col.names = TRUE)
     
    }
    
    return (snp1digit)
}

impuBkTDT.scheduler <-
function(raw=data, idx, job=1, toolname=NULL, freqMaps=NULL, dir=NULL, is.1digit=TRUE, dig1Code, dig2Code, reType=FALSE, reHap=NULL, logF=NULL, logErr=""){

  fStr ="[impuBkTDT.scheduler:]"
  ifD = FALSE
  #logF = "test"
  if(ifD) print(fStr)

  if(!is.null(dir)){
    if(dir==""){
      print(paste("Imputation information files(s) will be saved under current working directory:",  getwd()))

    }else{
      print(paste("Create directory:", dir, ", where imputation information file(s) are saved.", sep=""))
      dir.create(path=dir, showWarnings = TRUE)

      if (!is.null(logF)){
        logF = file.path(dir, qp(logF, ".txt"))
      }

      if (!is.null(reHap)){
        reHap = file.path(dir, reHap)
      }

      logErr = file.path(dir, logErr)

    }
  }else{

    if(job==1) {
#      print(
#          paste("Value for argument dir is NULL. Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))

#       paste("Function will not save any information about imputation. In case of error, log file will be saved under current working directory:", getwd()))
      reHap = NULL
      logF = NULL
    }else{
      if (!is.null(logF)){
        logF = qp(logF, ".txt")
      }
      print(
          paste("Value for argument dir is NULL. However, users ask to do multiple imputation. ",
                "Imputation information files(s) will be saved under current working directory:",  getwd()))
    }

  }


  # inside the function, assume 0 1 2 3 for NA, homo, homo, heter and 0 1 2 for NA, allele1, allele2

  ## first check the existance of required parameters
  if(is.null(toolname)){
    if(is.null(freqMaps)){
      stop("Values for both arguments, toolboxName and freqMaps, are missing.")
    }else{
      toolname = toolbox.load(freqMaps=freqMaps)
    }
  }
  
  ## find the start and end position for the blocks of interest

  all.genomeMarkerInfo = get(toolname$freqMap, envir = qingEnv)$genomeMarkerInfo

  bd.start = all.genomeMarkerInfo[idx[1], 2, drop=TRUE]
  bd.end   = all.genomeMarkerInfo[idx[ length(idx) ], 3, drop=TRUE]

  snpOffset =bd.start-1
  #print("Get here")
  ## change digit to 1 digit coding using Qing's coding scheme
  if(!is.1digit){
      ## excluding the parents cols
      genos  = raw[, 2+( ((bd.start-1)*2+1):(bd.end*2) ), drop=FALSE  ]
      snpNum = ncol(genos)/2
      snp1digit = exchangeDigit(ma=genos,
        cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
  }else{
      ## excluding the parents cols
      genos  = raw[, 2+(bd.start:bd.end), drop=FALSE]
      snpNum = ncol(genos)
      snp1digit = genos
      ## change to inside Qing's coding scheme
      if(min(dig1Code==c(0,1,3,2))==0){
        snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,3,2))
      }
  }

  #print("Get snp1digit")
  trioCt = nrow(snp1digit)/3
  maxRow = trioCt*job

  snp1digitTDT = matrix(NA, nrow=trioCt*6, ncol=(bd.end-bd.start+1))
  
  ## 18 columns: bk, trio, x1, x2, y1, y2, hap_f, hap_m, hap_c1-4
 
  if(reType){
    imputBkRecord = matrix(NA, ncol = 19, nrow=maxRow)
    impuDummy = imputBkRecord
  }else{
    imputBkRecord = matrix(NA, ncol = 18, nrow=maxRow)
    impuDummy = imputBkRecord
  }
  imputBkRecord.ct = 0
  errorTrap = NULL
  
  tryCatchEnv = new.env(parent=qingEnv)
  assign("trapID", 0, envir=tryCatchEnv)
  assign("errorTrap", errorTrap, envir=tryCatchEnv)
  
  #print("Get imputBkRecord.ct")
  #print(dim(imputBkRecord))
  
	if(!is.null(get(toolname$freqMap, envir = qingEnv)$hapBkOnlyMap)){
		all.hapIndex = get(toolname$freqMap, envir = qingEnv)$hapIndex
		
		hapBkOnlyMap.vars=list()
		tmp.hapMap =  get(toolname$freqMap, envir = qingEnv)$hapBkOnlyMap
		hapBkOnlyMap.vars$resiProbCol= tmp.hapMap$resiProbCol
		hapBkOnlyMap.vars$augIdxCol= tmp.hapMap$augIdxCol
		hapBkOnlyMap.vars$probCol= tmp.hapMap$probCol
		
	}else{
		all.hapIndex = c(-1)
		hapBkOnlyMap.vars=list()
	}


  snpCoding = 0:3
  snpBase = 0:2
  genoProb = genGenoProb()

  #print("Get loaded item")
  #print(idx)
  for( unit in idx){
      if (ifD) print( paste(fStr, " processing block index:", unit))
      if (is.element(unit, all.hapIndex)){
        ## haplotype, for each block, search every trio for missingness
        ## find block boundary
        bk.bd = unlist(all.genomeMarkerInfo[unit, c(2,3)])

        snpCt = bk.bd[2]-bk.bd[1]+1

        exhaustHap = get(toolname$exp, envir = qingEnv)
        exhaustHap = exhaustHap[1:(2^snpCt), 1:snpCt]

        bk.missTrioIdx = (1:trioCt)
    
        for (famId in bk.missTrioIdx){
          x1 = bk.bd[1]
          x2 = bk.bd[2]
          y1 = (famId-1)*3+1
          y2 = famId*3
          trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
          if(ifD) print( paste(fStr, " processing fam index:", famId))
          if(ifD) print(trioBlock)
          
          tryCatch({

            replace = ESp.imputBlock(appVarNames=toolname, 
                           trioBlock=trioBlock, snpLen=snpCt, bkIdx = unit, job=job,
                           snpCoding = snpCoding, snpBase=snpBase, reType = reType,  logF=logF, 
                           hapBkOnlyMap.vars=hapBkOnlyMap.vars
                            )
            
            imputBkRecord.ct = imputBkRecord.ct + 1
            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job), 1:6 ] = matrix(
                           rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=TRUE)

#             print(raw[, 1:4])
#             print(y1)
#             print(raw[y1, 1])
#             print(      matrix(rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=TRUE)[,1:4]
#                   )

            imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job), 7:ncol( imputBkRecord ) ] = replace
            fam.hapIdx1 = replace[1, c(1, 3, 5, 7, 9, 11 )]
            fam.hapIdx2 = replace[1, c(2, 4, 6, 8,10, 12 )]
            hap.str1 = exhaustHap[fam.hapIdx1, ]
            hap.str2 = exhaustHap[fam.hapIdx2, ]

            ## convert string into digits
            geno.FMCMa = covDipStr2CodedGeno(hap.str1, hap.str2, subjectCt=6, snpLen=snpCt)
 
            if(ifD) print(geno.FMCMa)
            if(ifD) print(snp1digit[y1:y2, x1:x2-  snpOffset])

            if( sum(abs(geno.FMCMa[1:3,] - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

            snp1digitTDT[ ((famId-1)*6+1) :( famId*6  ), x1:x2-  snpOffset] =
              geno.FMCMa[ c(1,2,3, sample(1:3, size=3, replace=FALSE)+3), ]
            
            if(!is.null(logF)){
      
              logl(logF, paste("For trio #", y2/3, ", Choose hap idx=[",
                               paste(replace[1:6], collapse=".", sep=""),
                               "]"))
              logl(logF, paste("Choose hap exp=[",
                               paste( apply(geno.FMCMa, 1, FUN=paste, collapse=""), collapse=".", sep=""), "]",
                               sep=""))
              logl(logF, paste("missing bk data:",
                               paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
      
            }
            if(ifD){
              print("-------------------")
              print(  c(y1, y2, x1, x2)  )              
            }

           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   
			   b = get("trapID", envir=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, envir=tryCatchEnv)			   
			   #trapID <<- trapID+1
		
				## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
				#errorTrap <<- rbind(errorTrap, errorInfo)	
	            errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
				a = get("errorTrap", envir=tryCatchEnv)
				a = rbind(a, errorInfo)
				assign("errorTrap", a, envir=tryCatchEnv)
	
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste( apply(trioBlock, 1, FUN=paste, collapse=""), collapse="."),
                               sep=""))
                  
                  #logl(logF, errTrace)
                }         
               }, warn =function(w) {
                 print("ImpuBlock::Warnings")
                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
              
        } ## for (famId in bk.missTrioIdx){

        if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){
          #print("Write")
          write.table(imputBkRecord[1:(imputBkRecord.ct*job) , , drop=FALSE],
                      file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                      append = TRUE, row.names = FALSE, col.names=FALSE)
        }
        imputBkRecord.ct  = 0
        imputBkRecord = impuDummy

      }else{
        ## genotype
        # print("A genotype")

        bk.bd = unlist(all.genomeMarkerInfo[unit, 2])

        tmp.Map =  get(toolname$freqMap, envir = qingEnv)$genoOnlyMap
        popuProb =  tmp.Map$bks[[bk.bd]][,tmp.Map$probCol]

        snpCt = 1
        bk.missTrioIdx = (1:trioCt) 

        for (famId in bk.missTrioIdx){
           x1 = bk.bd
           x2 = bk.bd
           y1 = (famId-1)*3+1
           y2 = famId*3
           trioBlock = snp1digit[y1:y2, x1:x2 -  snpOffset]
           if(ifD) print( paste(fStr, " processing fam index:", famId))
           if(ifD) print(trioBlock)
           
           tryCatch({
 
             replace = imputGeno( trioBlock, job=job, genoProb=genoProb,
                    popuProb = popuProb, data.order="FMC", snpCoding=snpCoding)

             #print(replace)
             imputBkRecord.ct = imputBkRecord.ct + 1
             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) ,
                           1:6 ] = matrix(rep(c(unit, raw[y1, 1], y1, y2, x1, x2), times=job), ncol=6, byrow=TRUE)

             imputBkRecord[((imputBkRecord.ct-1)*job+1):(imputBkRecord.ct*job) ,
                           c(7, 9, 11, 13, 15, 17) ] = replace
             geno.FMCMa = replace[1, 1:6]

             if( sum(abs(geno.FMCMa[1:3] - trioBlock)[ trioBlock!=0 ])!=0 ) stop("No matching")

             snp1digitTDT[ ((famId-1)*6+1):( famId*6  ), x1:x2-  snpOffset] =
               geno.FMCMa [ c(1,2,3, sample(1:3, size=3, replace=FALSE)+3)]
                     
           },  error = function(e) {
               ##print(qTraceback())
               traceback()
			   b = get("trapID", envir=tryCatchEnv)
			   b = b +1 
			   assign("trapID", b, envir=tryCatchEnv)			   
			   #trapID <<- trapID+1
			
               ## HARD CODE!!!HARD CODE: trio id is assumed to the be first one
			   #errorTrap <<- rbind(errorTrap, errorInfo)
			   errorInfo = c(trapID=b, bkIdx=unit, pedgree=raw[y1,1], case=raw[(y1+2),2], c(y1, y2, x1, x2))
			   a = get("errorTrap", envir=tryCatchEnv)
			   a = rbind(a, errorInfo)
			   assign("errorTrap", a, envir=tryCatchEnv)
			   
                if(!is.null(logF)){
                  logl(logF, paste("\nError trap id=(", b, ") and details for errors:", sep=""))
                  logl(logF, paste("Error block index {idx=", unit, "}------------", sep=""))
                  logl(logF, paste("Error trap famId=(", famId, ") and details for errors:", sep=""))
                  
                  logl(logF, paste("missing bk data:",
                               paste(trioBlock, collapse="."),
                               sep=""))
                  
                  #logl(logF, errTrace)
                }         
               }, warn =function(w) {
                 # print("ImpuBlock::Warnings")
                 
            } ) ## tryCatch

           if(!is.null(logF)){
             logl(logF, paste("end imputing block index {idx=", unit, "}------------\n", sep=""))
           }
         } ## for (famId in bk.missTrioIdx){

         if( (imputBkRecord.ct!=0) & (!is.null(reHap)) ){
          #print("Write")
          write.table(imputBkRecord[1:(imputBkRecord.ct*job) , , drop=FALSE],
                      file=paste(reHap, "_imputBkRecord.csv", sep=""), sep=",",
                      append = TRUE, row.names = FALSE, col.names=FALSE)
        }
        imputBkRecord.ct  = 0
        imputBkRecord = impuDummy
        
      } ## if (is.element(unit, all.hapIndex)){

#       print( qp("unit=", unit))
#       print( imputBkRecord[1:(imputBkRecord.ct*job), 1:4])
    } ## for( unit in idx){

	errorTrap=get("errorTrap", envir=tryCatchEnv)
    if(!is.null(errorTrap)){
      write.table(errorTrap, file=paste(logErr, "_errorTrap.csv", sep=""), sep=",",
                 append = FALSE, row.names = FALSE, col.names = TRUE)
     
    }
    return (snp1digitTDT)
}

imputGeno <-
function(trioBlock, job=1, genoProb, popuProb, data.order, snpCoding){
  fStr = "[imputGeno]:"
  ifD = FALSE
  if(ifD) print(qp(fStr, "start"))
  if( min( c(snpCoding==c(0,1,2,3),data.order=="FMC")) <1 )
     stop (paste(fStr, "Data configuration is not right:\n", "snpCoding=[", paste(snpCoding, collapse=";", sep=""),
                 "]", "trioBlock order=[", data.order, "]", sep=""))

  ##CAUTION, need to flip the popuProb, it is in the sequence for 11, 22, 12
  ## and re-standardize if freq==0
  tmpPop = popuProb/sum(popuProb)
  zeroPop = tmpPop==0
  #print(zeroPop)

  if(sum(zeroPop)>0){
    tmpPop = tmpPop/1.0001
    
    tmpPop[zeroPop]=.0001/sum(zeroPop)

#    if(ifD){
#      print(popuProb)
#      print(tmpPop)
#      print(sum(tmpPop))
#    }
    popuProb=tmpPop
  }
  
  popuProb=popuProb[c(1,3,2)]
  
## first process the block with complete genotypes
  if( sum(trioBlock==snpCoding[1])==0 ){
    seqpar = paste(trioBlock[1], trioBlock[2], sep="-")
    bench = paste(genoProb[,1], genoProb[,2], sep="-")
    matchPos = match(seqpar, bench)

    famRe = matrix(NA, nrow=job, ncol=6)

    
    for( i in 1:job){
      famRe[i, 1:2]= trioBlock[1:2]
      kid = unlist(lapply(1:3, FUN=function(item, repp){
        rep(item, repp[item])
      }, repp=genoProb[matchPos, 3:5]*4))
  
      newlyM = match(trioBlock[3], kid)
  
      allKid = kid[ c(newlyM[1], kid[-newlyM][sample(1:3, size=3, replace=FALSE)])  ]

      famRe[i, 3:6]=allKid      
    }

    return(famRe)
  }

  
  finalIdx = rep(T, nrow(genoProb))
  ## filter out the parents
  if(trioBlock[1]!=snpCoding[1]){
    matched = genoProb[,1]==trioBlock[1]
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)
  if(trioBlock[2]!=snpCoding[1]){
    matched = genoProb[,2]==trioBlock[2]
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)
  cProb = rep(1, nrow(genoProb))
  if(trioBlock[3]!=snpCoding[1]){
    ## hard code!!! hard code geno coding need to be 1,2,3 to correspond to the child geno
    ## indicated by genoProb
    matched = genoProb[, trioBlock[3]+2]!=0

    cProb = genoProb[, trioBlock[3]+2]
    
    finalIdx = finalIdx & matched
  }
  if(ifD) print(finalIdx)

  if(sum(finalIdx)==0) stop(paste(fStr,
          "Mendelian error! trio geno=[", paste(trioBlock, collapse=".", sep=""),
          "] for order =[", order, "]", sep=""))
#  print("pop")
#  print(popuProb)
  fProb = popuProb[genoProb[,1]]
  mProb = popuProb[genoProb[,2]]
  jProb = fProb*mProb
  jProb = jProb/sum(jProb)

  ## check if estimated genotype return no genotype for this one
  jointProb = sum(jProb[finalIdx])

  if(jointProb==0) {
    jProb[!finalIdx]=0
    jProb = jProb*cProb 
    jProb[finalIdx]=rep(1/sum(finalIdx), sum(finalIdx))
    print("should not happen")
    #print(trioBlock)
    stop("should not happen")
    
  }else{
    jProb[!finalIdx]=0
    jProb = jProb*cProb
    jProb = jProb/sum(jProb)
  }
  
  if(ifD) print(cbind(fProb, mProb, jProb, finalIdx, genoProb[,1:2]))

  allKid = t(apply(genoProb[, 3:5]*4, 1, FUN=function(row){
                  a = rep(c(1, 2, 3), times=row); a}))
  #print(allKid)
  #print(jProb)
  
  re6Geno = matrix(NA, ncol=6, nrow=job)

  for( ss in 1:job){
    
    chooseRow = sample(1:9, size=1, prob = jProb)
  
    re = genoProb[chooseRow, 1:2]
    if(trioBlock[3]!=snpCoding[1]){
      childGeno = trioBlock[3]
    }else{
      ## hard code!!! hard code, geno coding need to be 1,2,3 to correspond to the child geno
      ## indicated by genoProb
      childGeno = sample(1:3, size=1, prob = genoProb[chooseRow, 3:5])
  
    }

    oth = allKid[chooseRow,]
    #print(oth)
    othleft = oth [- which(oth==childGeno)[1] ]
    re6Geno[ss,] = c(re, childGeno, othleft)
  
  }
  #print(re6Geno)
  return(re6Geno)
}

isDigitAtLociIdx <-
function(hapIdx, digit=1, lociCt, lociIdx, intVec=seq.int(from=1, to=2^(lociIdx-1), by=1)){

  ifD = FALSE
  
  if(digit!=1 & digit!=2) stop(paste("Invalide digit imput: ", digit, sep=""))
  
  leftSide = (hapIdx - intVec)/2^lociIdx + 1

  if(ifD) print(leftSide)
  
  roundLeftSide = as.integer(leftSide)

  if(ifD) print(roundLeftSide)

  integerLeft = leftSide[roundLeftSide == leftSide]

  if(ifD) print(integerLeft)

  one = integerLeft >= 1
  oth = integerLeft <= 2^(lociCt-lociIdx)

  fittedCt = one & oth

  if(ifD) print(fittedCt)

  if(length(fittedCt)>1) stop("Error!. Should left with only one or none choice")

  if(length(fittedCt)==1){
    if(digit==2) fittedCt = !fittedCt
    return(fittedCt)
  }else{
    if(digit==2) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}

linkageFile.proc <-
function(data, snpIdxRange=NULL, key.prefix="", bk.sizes=NULL, action = c("outputTrio", "formTrio1digit", "missingReport", "Mendelian check", "freq estimate"), dig2Code=0:2, dig1Code=c(0,1,3,2), ... ){
	
	#dig2Code=0:2
	#dig1Code=c(NA,0,1,2)
	
	pedCol =1
	memCol=2
	affectCol = 6
	dadCol=3
	momCol=4
	txt.affect=2
	
	sep=""
	header=FALSE
	
	
	trioTwoDigit = snpPREFileMatchTrio(txtF=data, sep=sep, header=header,
			pedCol =pedCol, memCol=memCol, affectCol =affectCol,
			dadCol=dadCol, momCol=momCol, txt.affect=txt.affect,
			logF=NULL)
	#print(dim(trioTwoDigit))
	
	if(is.null(snpIdxRange)) snpIdxRange = c(7, ncol(trioTwoDigit))
	
	snpStartLeftIndex=snpIdxRange[1]
	snpEndRightIndex =ncol(trioTwoDigit) - snpIdxRange[2] + 1
	
	re = list()
	if(is.element("outputTrio", action)){
		re = c(re, trio2digit=list(trioTwoDigit))
	}
	
	genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
	snpNum = ncol(genos)/2
	
	
	## check linkage file code
	linkqing = unique(as.vector(unlist(genos)))
	tt.check = match(linkqing, dig2Code) 
	if(max(is.na(tt.check))==1)
		stop(paste("dig2Code is wrong.", "Existing codes in data:[",
						paste(linkqing, collapse=";", sep=""), "].")
		)
	
	snp1digit = exchangeDigit(ma=genos,
			cols=c(1,snpNum*2), dig1Code=dig1Code, dig2Code =dig2Code, action=c("2to1"))
	
	snp1digit.inside = exchangeDigit(ma=genos,
			cols=c(1,snpNum*2), dig1Code=c(0, 1, 3, 2), dig2Code =dig2Code, action=c("2to1"))
	
	trio1digit = cbind( trioTwoDigit[, c(pedCol, memCol)], snp1digit)
	if(is.element("formTrio1digit", action)){
		
		re = c(re, trio=list(trio1digit))
	}
	
	if(is.element("missingReport", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
			warning("Cannot report missing information. The argument, snpIdxRange, includes one of the special column for linkage file.")
		}else{
			missSNPPos = findMissing(df=trioTwoDigit, is.1digit=FALSE, snpStartLeftIndex=snpStartLeftIndex,
					snpEndRightIndex=snpEndRightIndex, dig1Code=NULL, dig2Code=dig2Code )
			re = c(re, missIdx = list(missSNPPos))
		}
	}
	
	if(is.element("Mendelian check", action)){
		snpTrio = matrix(snp1digit.inside, nrow=3, byrow=FALSE)
		
		MedErr = matrix(NA, ncol=4, nrow=ncol(snpTrio))
		colnames(MedErr)=c("y", "x", "trio", "SNP")
		
		tryCatchEnv = new.env(parent=qingEnv)
		assign("MedErr.ct", 0, envir=tryCatchEnv)	
		assign("MedErr", MedErr, envir=tryCatchEnv)
		
		trioCt = nrow(snp1digit)/3
		
		## CHANGED: 2009: report the index in the trio1digit
		tmpDigit = 1
		#print(str(snpTrio))
	
		for ( i in 1:ncol(snpTrio)){
			tryCatch({
						tt = checkMendelianError(codedSNPTrio=snpTrio[,i], snpCoding=c(0,1,2,3))
					}, error = function(e){
						#print(paste("trioCt=", trioCt))
						#print(paste("snpStartLeftIndex=", snpStartLeftIndex))
					    a = get("MedErr.ct", envir=tryCatchEnv)
					    a = a+1
					    assign("MedErr.ct", a, envir=tryCatchEnv)
					    #MedErr.ct <<- MedErr.ct +1
						#print(paste("MedErr.ct=", MedErr.ct, " i=", i))
						tttx = ceiling(i/trioCt)
						ttty = i%%trioCt
						if(ttty==0) ttty = trioCt
						## CHANGED: 2009: report the index in the trio1digit
						b = get("MedErr", envir=tryCatchEnv)
						b[a,] = c( (ttty-1)*3+1,  (tttx-1)*tmpDigit+3,  ttty,    tttx)
						assign("MedErr", b, envir=tryCatchEnv)
						#print( MedErr[MedErr.ct,,drop=FALSE] )
					})          
		}
		MedErr.ct = get("MedErr.ct", envir=tryCatchEnv)
		MedErr = get("MedErr", envir=tryCatchEnv)
		
		if(MedErr.ct==0) {
			MedErr=NULL
			#print("No Mendelian error.")
			re = c(trio=list(trio1digit), MedErr=list(MedErr[1:MedErr.ct,,drop=FALSE]))
		}else{
			#print("Found Mendelian error(s).")
			re = c(MedErr=list(MedErr[1:MedErr.ct,,drop=FALSE]), trio.err=list(trio1digit))
		}
	}
	
	if(is.element("freq estimate", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
			warning("Cannot provide frequencies estimation. The argument,  snpIdxRange, includes one of the special column for linkage file.")
		}else{
			
			tmp = ncol(trioTwoDigit)
			parTwoDigit = getBackParentGeno(trioDf=trioTwoDigit, famCol=pedCol, memCol=memCol,
					snpIdx=snpStartLeftIndex:(tmp-snpEndRightIndex+1), re.child=FALSE, prefix=NULL)
			#print("###")
			#print(dim(parTwoDigit))
			
			if( is.null(bk.sizes) ){
				warning("Cannot provide frequencies estimation. The argument, bk.sizes, is not provided for user-specify option.")
			}else{
				
				map=freqmap.reconstruct(data=parTwoDigit, cols=c(3, ncol(parTwoDigit)), loci.ct=bk.sizes, is.1digit=FALSE,
						dig1Code=NULL, dig2Code = dig2Code, key.prefix=key.prefix, start.base=1, ...)
				
				re = c(re, freq=list(map))
			}
			
		} ##if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
		
	}
	
	return(re)
	
}

listExtractor <-
function(nameList, varList, na.replace = NULL){
  ifD = FALSE

  varNames = names(varList)
  varLen = length(varList)
  exLen = length(nameList)

  m = match(nameList, varNames, 0)
  if(min(m)==0 & is.null(na.replace)){
    stop(paste("Variable(s) not found in the varList. NULL not allowed. Names in varList: ", paste(varNames, collapse=", "), sep=""))
  }

  re = NULL
  for(i in 1:length(m)){
    if(m[i]==0){
      re = c(re, na.replace)
    }else{
      re = c(re, varList[m[i]])
    }
  }
  return (re)
}

logBe <-
function (fileName, str=NULL){
       cat(paste("\n\n=============Begin of the Script::" , Sys.time(), "==================") , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       if(!is.null(str)) cat(str , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       return(NULL)
}

loge <-
function (fileName, str=NULL){
       if(!is.null(str)) cat(str , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       cat(paste("-------------End of the Script::" , Sys.time(), "--------------------") , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       return(NULL)
}

logErr <-
function (fileName, str){
       cat("!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!" , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       cat(str , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       cat("!!!!!!!!!!!! ERROR ------------------- ERROR !!!!!!!!!!!" , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       return(NULL)
}

logl <-
function (fileName, str){
       cat(str , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       return(NULL)
}

logs <-
function (fileName, str){
       cat(str , file = fileName, sep = " ", fill = FALSE, labels = NULL, append = TRUE)
       return(NULL)
}

logWarn <-
function (fileName, str){
       cat("************ WARNING ***************** WARNING ***********" , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       cat(str , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       cat("************ WARNING ----------------- WARNING ***********" , file = fileName, sep = " ", fill = TRUE, labels = NULL, append = TRUE)
       return(NULL)
}

ls.list <-
function(path=NULL, objname,  exam.ext = .3){
  
  myobj = NULL
 
  
  if( (!is.null(path)) & (!is.character(objname))){
    assign(objname, NULL, envir = qingEnv)
    load(path, envir = qingEnv)
    myobj = get(objname, envir = qingEnv)
  }else{
    myobj = objname
  }

  if( !is.list(myobj)){
    stop(paste("Object name (", objname, ") is not a list!", sep=""))
  }
  element = NULL
  l.len = length(myobj)
  names.all = names(myobj)
  for( i in 1:l.len){
    #print(i)
    #print(element)
    elm = myobj[[i]]
    # str(elm)
    if( !is.list(elm)){
      element = c(element, paste("Object", i, ") name=", names.all[i], "; value=", elm, ";", sep=""))
    }else{
      # check structure
      #str(elm[[1]])
      if (length(elm)==1){
        element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=1", sep=""))
      }else if(length(elm)>1){
        is.allSame = qStr.pattern(elm, exam.ext=exam.ext)
        ##print(is.allSame)
        if( is.null(is.allSame)){
          element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=", length(elm), ";", sep=""))
        }else{
          ## same objects 
          element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=", length(elm), ";", sep=""))
      
          element = c(element, paste("            structure::", is.allSame, ";", sep=""))

        }
      }else{
        element = c(element, paste("Object", i, ") List::name=", names.all[i], "; length=0", sep=""))
      }

    }
  }
  return(element)
}

matTBCondOnChild3.finalProc <-
function(maxTbl, counter, maxMateTbl, probVec, prob1, prob2, job, logF){
  #print("Run finalProc")
  ifD = FALSE

  cutoff = min(counter, 20)
  ## in the maxMateTbl, the first column is index within par1 and the second column is index within par2
  if (ifD) print(maxTbl[1:100,])
  if (ifD) print(cbind(maxMateTbl[1:cutoff, ], probVec[1:cutoff]))
  prob1n = rep(0, length(prob1))
  prob1n[unique(maxMateTbl[1:counter,1])]= prob1[unique(maxMateTbl[1:counter,1])]

  prob1n = prob1n/sum(prob1n)
  onePDipProb = prob1n[maxMateTbl[1:counter,1]]
  if(ifD) print(paste("onePDipProb=", paste(round(onePDipProb,3)[1:cutoff], collapse=";", sep="")))


  prob2n = rep(0, length(prob2))
  prob2n[unique(maxMateTbl[1:counter,2])]= prob2[unique(maxMateTbl[1:counter,2])]

  prob2n = prob2n/sum(prob2n)
  othPDipProb = prob2n[maxMateTbl[1:counter,2]]
  if(ifD) print(paste("othPDipProb=", paste(round(othPDipProb,3)[1:cutoff], collapse=";", sep="")))

  prob = onePDipProb * othPDipProb *  probVec[1:counter]
  
  prob = prob/sum(prob)

  if(ifD) print(paste("final prob=", paste(round(prob,3)[1:cutoff], collapse=";", sep="")))

  ## sample the row
  hap6idx = matrix(NA, nrow=job, ncol=6)
  if(ifD) print(hap6idx)
  for(ss in 1:job){
    chooseRow = sample(1:counter, size=1, prob=prob)
    #print(chooseRow)
    hap6idx[ss,] = maxTbl[chooseRow, , drop=FALSE]
    #print(hap6idx)
  }

#   if(!is.null(logF)){
#     if(cutoff==20){
#       bkMapStr = paste(fStr, "\nCutoffed possible dip mating table:\n",
#                      paste(wrComTbl( cbind(maxTbl[1:cutoff, ,drop=FALSE], prob[1:cutoff], probVec[1:cutoff]),
#                               colNames=c("f1", "f2", "m1", "m2", "c1", "c2", "prob", "probVec")), collapse="\n", sep=""),
#                      sep="")
#     }else{
#       bkMapStr = paste(fStr, "\nCompleted poss dip mating table:\n",
#                      paste(wrComTbl( cbind(maxTbl[1:cutoff, ,drop=FALSE], prob[1:cutoff], probVec[1:cutoff]),
#                               colNames=c("f1", "f2", "m1", "m2", "c1", "c2", "prob", "probVec")), collapse="\n", sep=""),
#                      sep="")
# 
#     }
#     logl(logF, bkMapStr)
#     logl(logF, paste(fStr, "choose row num=", chooseRow, sep=""))
#     #logl(logF, paste("Parents index=[", paste(as.vector(hap6idx), collapse=".", sep=""), "]", sep=""))
#   }

  if(ifD) {
    tttt = cbind(maxTbl[1:cutoff, ,drop=FALSE], prob[1:cutoff])
    print(tttt)
    print(hap6idx)
  }
  rm(maxTbl)
  rm(maxMateTbl)
  rm(probVec) 

  return( hap6idx )

}

moveDir <-
function(par=getwd(), fromRoot, toRoot, subDir=NULL){


  fListdir = list.files(path = file.path(par, fromRoot))
  dir.index = file.info(file.path(par, fromRoot, fListdir))$isdir
  dirs = fListdir[dir.index]
  for( j in dirs){
    newdir = file.path(par, toRoot, j)
    dir.create(path=newdir, showWarnings = TRUE)
  }
  ## addtional subdirectory
  if(!is.null(subDir))
    dir.create(path=file.path(par, toRoot, subDir[1], subDir[2]), showWarnings = TRUE)
  
  fList = list.files(path = file.path(par, fromRoot), all.files = TRUE,
           full.names = FALSE, recursive = TRUE)

  for( i in fList){
    f1 = file.path( par, fromRoot, i)
    f2 = file.path( par, toRoot, i)
    print(f1)
    print(f2)
    file.copy(from=f1, to=f2, overwrite=TRUE )

  }
}

procSemiAugMap <-
function(appVarNames, homoHetoInfo, snpLen){
	
	fStr="[procSemiAugMap]:"
	
	if(is.null(homoHetoInfo$homoIn)){
		## if not homo digit requirment, every hap is possible
		retainList = 1:(2^snpLen)
		return(retainList)
	}
	
	## if there is homo digit requirement
	idx4hapDigit = NULL
	## check out the global variables
	tryCatch({
				
				tmpGetObj = NULL
				tmpGetObj = get(appVarNames$digit, envir=qingEnv)
				idx4hapDigit$digitMap1 = tmpGetObj$digitMap1[1:(2^(snpLen-1)), 1:snpLen]
				idx4hapDigit$digitMap2 = tmpGetObj$digitMap2[1:(2^(snpLen-1)), 1:snpLen]
				
				rm(tmpGetObj)
				##gc()
			}, error=function(e){
				## traceback()
				## print(e)
				errTrace = paste(e, collapse=";", sep="")
				stop(paste("\n", fStr, errTrace, "\nApp-wise Global Variable ", appVarNames$digit, " does not exisit."))
			})
	
	
	
	## first process the homo digit
	
	retainList = NULL
	matchedIdx = NULL
	for( i in 1:length(homoHetoInfo$homoIn)){
		curDigit =  homoHetoInfo$homoDigit[i]
		tmpPos = homoHetoInfo$homoIn[i]
		if(curDigit == 1) {
			matchedIdx = idx4hapDigit$digitMap1[, tmpPos, drop=FALSE]
			if(is.null(retainList)){
				retainList = matchedIdx
			}else{
				retainList = retainList[is.element(retainList, matchedIdx)]
			}
			## print(retainList)
		}else{
			matchedIdx = idx4hapDigit$digitMap2[, tmpPos, drop=FALSE]
			if(is.null(retainList)){
				retainList = matchedIdx
			}else{
				retainList = retainList[is.element(retainList, matchedIdx)]
			}
			##  print(retainList)
		} ## if(curDigit == 1) {
	} ## for( i in homoHetoInfo$homoIn){
	
	rm(idx4hapDigit)
	##gc()
	if( is.null (retainList) )
		stop(paste("\nBlock infomation error. No population haplotype matches existing data:(", expression, ").", sep=""))
	
	return (retainList)
	
}

qExpandTable <-
function(listOfFactor =  list( 1:3, 10:13), removedRowIdx=NULL, re.row=FALSE ){
   tblDim = length(listOfFactor)
   tblDimSeq = unlist(lapply(listOfFactor, FUN=length))

   recyMa = matrix(listOfFactor[[1]], ncol=1, nrow=tblDimSeq[1])
   recyRow = tblDimSeq[1]
   ## grow the index list  
   for ( i in 2:tblDim ){
     growing = util.matrix.clone(recyMa, tblDimSeq[i])
     ## addecCol =  rep(listOfFactor[[i]], each=recyRow)
     recyMa = cbind(growing,  rep(listOfFactor[[i]], each=recyRow))
     recyRow = recyRow * tblDimSeq[i]  
   }

   rowIdx = NULL
   if ( (!is.null(removedRowIdx)) | re.row){
     lengthRev = c(0, cumprod(tblDimSeq)) [ tblDim:1 ]
     tmp = ( removedRowIdx[tblDim:1] - 1)* lengthRev
     rowIdx = sum(tmp)+removedRowIdx[1]
     ## print(lengthRev)
     ## print(tmp)
   }

   ## print(rowIdx)
   ## print(recyMa)
   if(re.row){
     return(list(ma=recyMa, rowNum = rowIdx))
   }else{
     if(!is.null(removedRowIdx)){
       recyMa = recyMa[-rowIdx, , drop=FALSE]
     }
     return(recyMa)
   }
}

qing.cut <-
function(val, cutPt, cutPt.ordered = TRUE, right.include=TRUE){

  ## !!! HARD CODE !!! -1 means the val >/>= the maximum value in the cutPt
  ## return the matched index
  if(!cutPt.ordered) cutPt = order(cutPt)

  cutPt.ct = length(cutPt)

  if(cutPt.ct<1) stop("Zero length cutPt.")
  if(cutPt.ct<1) stop("Zero length val.")
  
  val.ct = length(val)

  re = rep(NA, length=val.ct)
  
  numFalse = re
  cellMatchedIdx = re

  if(right.include){
    for( i in 1:val.ct){
      cVal = val[i]
      falseMat = cVal > cutPt
      numFalse[i] = sum(falseMat)
    }
  }else{
    for( i in 1:val.ct){
      cVal = val[i]
      falseMat = cVal >= cutPt
      numFalse[i] = sum(falseMat)
    }
  }
    
  cellMatchedIdx [ numFalse==cutPt.ct ] = -1
  cellMatchedIdx [ numFalse!=cutPt.ct ] = numFalse[ numFalse!=cutPt.ct ]+1  

  return(cellMatchedIdx)

}

qing.mulMatch <-
function(val, bench){
  
  matched = bench==val
  bench.seq = 1:length(bench)
  re = bench.seq[matched]
  if(length(re)>=1){
    return(re)
  }else{
    return(0)
  }
}

qp <-
function(..., sep=""){
  str = paste(..., sep=sep)
  return(str)
}

qStr.pattern <-
function(objList,  exam.ext = .3){
  l.len = length(objList)

  class.last = NULL
  len.last = NULL
  names.last = NULL
  i = 1
  search=TRUE
  while( i <= min((round(l.len*exam.ext, 0)+1), l.len) & search ){
    #print(i)
    elm = objList[[i]]
    names = names(elm)
    class = unlist(lapply(elm, FUN=class))
    len = unlist(lapply(elm, FUN=length))

    if(!is.null(names.last)){
      ## see if the names matched
      if(length(names.last)==length(names)){
        all.m = names.last==names
        if(sum(all.m)!=length(names)) search=FALSE
      }else{
        search=FALSE  
      }
      if(length(class.last)==length(class)){
        all.m = class.last==class
        if(sum(all.m)!=length(class)) search=FALSE
      }else{
        search=FALSE  
      }
      if(length(len.last)==length(len)){
        all.m = len.last==len
        if(sum(all.m)!=length(len)) search=FALSE  
      }else{
         search=FALSE  
      }      
      
    }
    names.last=names
    names=NULL
    class.last=class
    class=NULL
    len.last=len
    len=NULL    
    i = i+1
  }

  if(search){
    if (is.null(names.last[1])){
      reStr = paste(c("name", "class", "length"),
                  c("",  class.last[1],  len.last[1]), 
                  sep="=", collapse="; ")
    }else{
      reStr = paste(c("name", "class", "length"),
                  c(names.last[1],  class.last[1],  len.last[1]), 
                  sep="=", collapse="; ")
    }
  }else{
    reStr = NULL
  }
  return(reStr)
}

qstrsplit <-
function(str, delim, re.1st = FALSE){
  re = unlist(strsplit(str, delim))
  if (length(re)==1){
    # possible not find the delim
    if (nchar(re[1])==nchar(str)){
      # not find the delim
      return(NA)
    }else  if(nchar(re[1])==0){
      # str contain delim itself
      if(re.1st ){
        return(re)
      }else{
        return(NULL)
      }
    }else{
      # str ends with delim
      return(re)
    }
  }else{
    if(nchar(re[1])==0){
      # str starts with delim
      if(re.1st){
        return(re)
      }else{
        return(re[-1])      
      }
    }else{
      return(re)
    }
  }
}

qTraceback <-
function(x = NULL){
    if (is.null(x) && (exists(".Traceback", envir = baseenv()))) 
        x <- get(".Traceback", envir = baseenv())

    reStr = NULL
    if (is.null(x) || length(x) == 0) 
        cat(gettext("No traceback available"), "\n")
    else {
        n <- length(x)
        m0 <- getOption("deparse.max.lines")
        for (i in 1:n) {
            label <- paste(n - i + 1, ": ", sep = "")
            m <- length(x[[i]])
            if (m > 1) 
                label <- c(label, rep(substr("          ", 1, 
                  nchar(label, type = "w")), m - 1))
            if (is.numeric(m0) && m0 > 0 && m0 < m) {
                cat(paste(label[1:m0], x[[i]][1:m0], sep = ""), 
                  sep = "\n")
                cat(label[m0 + 1], " ...\n")
                reStr = c(reStr, paste(label[1:m0], x[[i]][1:m0], sep = "") )
                reStr = c(reStr, label[m0 + 1], " ...\n")
            }
            else {
              cat(paste(label, x[[i]], sep = ""), sep = "\n")
              reStr = c(reStr, paste(label, x[[i]], sep = ""))
            }
        }
    }
    return(reStr)

}

resample <-
function(x, size, ...)
  if(length(x) <= 1) { if(!missing(size) && size == 0) x[FALSE] else x
  } else sample(x, size, ...)

sam.Hap.old <-
function(exp, prob, subjectCt){
	
	if(length(exp)==1){
		reExp = rep(exp, subjectCt)
	}else{
		reExp = sample(exp, size=subjectCt, prob=prob, replace=TRUE)
	}
	
	return(reExp)
	
}

sam.reject <-
function(allIdx=NULL, idxProb=NULL, rejectIdx, size=1){
	
	if (is.null(allIdx)){
		## if allIdx==NULL, then allIdx = 1:length(idxProb)
		length = length
		allIdx = 1:length
	}else{
		## if idxProb==NUL, then prob are the same for all idx,
		length = length(allIdx)
		if(is.null(idxProb)) idxProb = rep(1/length, length)
	}
	
	if(length<=0) stop(paste("Invalid input value for allList with length=[", length, "]", sep=""))
	
	meet = match(rejectIdx, allIdx)
	idxProb[meet]=0
	idx = resample(allIdx, size=size, prob=idxProb, replace=TRUE)
	
	return(idx)
}

sampleDipSemiAugMap <-
function(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = FALSE
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	commonProbIdx = length(mapProb) 
	
	## need to take account the situation where all the hap is presented.
	if(commonProbIdx == 2^snpLen){
		chooseDip1 = sample(1:commonProbIdx, size=1, prob = mapProb)
		chooseDip2 = sample(1:commonProbIdx, size=1, prob = mapProb)
		
		chooseDip1 = mappedAugIdx[chooseDip1]
		chooseDip2 = mappedAugIdx[chooseDip2]
	}else{
		commonProbIdx = length(mapProb) + 1
		chooseDip1 = sample(1:commonProbIdx, size=1, prob = c(mapProb, leftOverProb))
		chooseDip2 = sample(1:commonProbIdx, size=1, prob = c(mapProb, leftOverProb))
		
		allIdx = NULL
		
		if(chooseDip1 == commonProbIdx){
			allIdx = 1 :(2^snpLen)
			chooseDip1 = sampleIdxOutsideList(allIdx, mappedAugIdx)
		}else{
			chooseDip1 = mappedAugIdx[chooseDip1]
		}
		
		if(chooseDip2 == commonProbIdx){
			if(is.null(allIdx)) allIdx = 1:(2^snpLen)
			chooseDip2 = sampleIdxOutsideList(allIdx, mappedAugIdx)
		}else{
			chooseDip2 = mappedAugIdx[chooseDip2]
		}
		
	}
	
	reDip = range(chooseDip1, chooseDip2)
	
	return(reDip)
	
}

sampleHapExclude <-
function(idxs, hapProb){
	newHapProb = hapProb
	newHapProb[idxs]=0
	idx = sample(1:length(newHapProb), prob=newHapProb, size=1)
	return(idx)
}

sampleHapSemiAugMap2 <-
function(semiMapFrame, resiProbCol, augIdxCol, probCol, snpLen, exHapIdx=NULL){
	## SemiAugMap only have the major haplotypes' expression and prob,
	## only a number for prob for each of the other haplotypes 
	## SemiAugMap also has a column showing the index of the haplotype in the fixed AugMap
	ifD = FALSE
	fStr = "[sampleHapSemiAugMap2]:"
	if(ifD) print(fStr)
	
	leftOverProb = semiMapFrame[1, resiProbCol]
	mapProb = semiMapFrame[ , probCol]
	mappedAugIdx = semiMapFrame[ , augIdxCol]
	
	exMajor = NULL
	exMinor = NULL
	exMajor.ct = 0
	exMinor.ct = 0
	commonProbIdx = length(mapProb)
	chooseHapIdx = NULL
	
	if(!is.null(exHapIdx)){
		
		matched = (  (exHapIdx >=1) & (exHapIdx <= (2^snpLen)))
		
		if(sum(matched)!=length(exHapIdx)) stop( paste("not valid exHapIdx=[", paste(exHapIdx, collapse=";", sep=""), "]", sep=""))
		
		majorMatch = is.element(exHapIdx, mappedAugIdx)
		exMajor = exHapIdx[majorMatch]
		exMinor = exHapIdx[!majorMatch]
		
		leftMajorRow = (1:commonProbIdx)[!is.element(mappedAugIdx, exMajor)]
		exMajor.ct = sum(majorMatch)
		exMinor.ct = sum(!majorMatch)
		
	}
	
	if(ifD){
		print(semiMapFrame)
		if(!is.null(exHapIdx)){
			print(paste("exMajor=[", paste(exMajor, collapse=";", sep=""), "]", sep=""))
			print(paste("leftMajorRow=[", paste(leftMajorRow, collapse=";", sep=""), "]", sep=""))
			print(paste("exMinor=[", paste(exMinor, collapse=";", sep=""), "]", sep=""))
			print(paste("majorMatch=[", paste(majorMatch, collapse=";", sep=""), "]", sep=""))
		}
	}
	
	
	## need to take account the situation where all the hap is presented.
	if(commonProbIdx == 2^snpLen){
		if(exMinor.ct>=1) stop(paste(fStr, " impossible to exclude minor when the map is completed."))
		## the only excluded are major
		if(exMajor.ct>0){
			## readjust the prob
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			
			if(ifD) print(paste("mapProb=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct), size=1, prob=mapProbN)
			chooseHapIdx = mappedAugIdxN[chooseRow]
			
		}else{
			## sample from completed map
			if(ifD) print(paste("mapProb=[", paste(round(mapProb, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx , "]", sep=""))
			
			chooseRow = sample(1:commonProbIdx, size=1, prob=mapProb)
			chooseHapIdx = mappedAugIdx[chooseRow]
		}
		
	}else{
		## not a complete map
		if(exMajor.ct==0 & exMinor.ct==0){
			## no exclusion
			
			mapProbN = c(mapProb,  leftOverProb)
			
			if(ifD) print("(exMajor.ct==0 & exMinor.ct==0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx+1), size=1, prob=mapProbN)
			
			if(chooseRow==commonProbIdx+1){
				allIdx = 1:(2^snpLen)
				chooseHapIdx = sampleIdxOutsideList(allIdx, mappedAugIdx)
			}else{
				chooseHapIdx = mappedAugIdx[chooseRow]
			}
			
		}else if(exMajor.ct==0 & exMinor.ct>0){
			## only exclude minor
			before.minor = 2^snpLen - commonProbIdx
			mapProbN = c(mapProb,  leftOverProb/before.minor*(before.minor-exMinor.ct))
			
			if(ifD) print("(exMajor.ct==0 & exMinor.ct>0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx+1), size=1, prob=mapProbN)
			
			
			if(chooseRow==commonProbIdx+1){
				allIdx = 1:(2^snpLen)
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx, exMinor))
			}else{
				chooseHapIdx = mappedAugIdx[chooseRow]
			}      
			
		}else if(exMajor.ct>0 & exMinor.ct==0){
			## only exclude major
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			mapProbN = c(mapProbN,  leftOverProb)
			
			if(ifD) print("(exMajor.ct>0 & exMinor.ct==0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct+1), size=1, prob=mapProbN)
			
			if(chooseRow==(commonProbIdx-exMajor.ct+1)){
				allIdx = 1:(2^snpLen)
				## tricky, can only sample minor, augIdx pool are not reduced
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx))
			}else{
				chooseHapIdx = mappedAugIdxN[chooseRow]
			}      
		}else if(exMajor.ct>0 & exMinor.ct>0){
			## exlude both major and minor
			mapProbN = mapProb[leftMajorRow]
			mappedAugIdxN = mappedAugIdx[leftMajorRow]
			
			before.minor = 2^snpLen - commonProbIdx
			
			mapProbN = c(mapProbN,  leftOverProb/before.minor*(before.minor-exMinor.ct))
			
			
			if(ifD) print("(exMajor.ct>0 & exMinor.ct>0):")
			if(ifD) print(paste("mapProbN=[", paste(round(mapProbN, 2), collapse=";", sep=""), "]", sep=""))
			if(ifD) print(paste("length=[", commonProbIdx-exMajor.ct+1 , "]", sep=""))
			
			chooseRow = sample(1:(commonProbIdx-exMajor.ct+1), size=1, prob=mapProbN)
			
			if(chooseRow==(commonProbIdx-exMajor.ct+1)){
				allIdx = 1:(2^snpLen)
				## tricky, can only sample minor, augIdx pool are not reduced
				chooseHapIdx = sampleIdxOutsideList(allIdx, c(mappedAugIdx, exMinor))
			}else{
				chooseHapIdx = mappedAugIdxN[chooseRow]
			}    
			
			
		}else{
			stop("Nothing here")
		}
		
	}
	
	return(chooseHapIdx)
	
}

sampleIdxOutsideList <-
function(allList, listVec, maxIt = 1000){
	
	length = length(allList)
	if(length<=0) stop(paste("Invalid input value for allList with length=[", length, "]", sep=""))
	
	keepSearch = TRUE
	idx = NULL
	count = 1 
	while(keepSearch & count<= maxIt){
		idx = sample(1:length, size=1)
		if(!is.element(allList[idx], listVec)){
			keepSearch = FALSE
			return(allList[idx])
		}   
		count = count+1
	}
	
	stop(paste("\nInefficient sampling. Iteration count exceeds maximum limit=[", maxIt, "].",
					"\nallList = [", paste(allList, collapse=";", sep=""),
					"] and listVect = [", paste(listVec, collapse=";", sep=""), "].",
					sep=""))
	return(NULL)
}

semiAugBkFrame <-
function(hapBkMap, key, probLeftOver = .01){
	
	ifD = FALSE
	
	keyIndex = which(hapBkMap$keys==key)
	bk = hapBkMap$bks[[keyIndex]]
	
	## extract the original exp and prob, then restandard prob
	estBkExp = as.character(bk[,hapBkMap$expCol])
	estBkProbRe = bk[,hapBkMap$probCol]/(sum( bk[,hapBkMap$probCol] ))*(1-probLeftOver)
	
	## create the exhaust haplotypes 
	bkSnpLen = bk[1,hapBkMap$hapLenCol]
	exhaustExp = exhaustHapExp(lociCt=bkSnpLen, snpCoding=c(1,2))$hapStr
	exhaustExpCt = length(exhaustExp)
	
	## match the expression in the short list to the exhaustive list
	rematch = match(estBkExp, exhaustExp)
	
	## replace the expression's probability
	bk[,hapBkMap$probCol]=estBkProbRe
	bk$augIdx = rematch
	bk$resiProb = rep(probLeftOver, times=nrow(bk))
	
	hapBkMap$bks[[keyIndex]] = bk
	
	hapBkMap$augIdxCol = ncol(bk)-1
	hapBkMap$resiProbCol = ncol(bk)
	
	#if(ifD) print(bkFrame)
	
	## other parts, like df, dfStr,  in hapBkMap is not updated because it is not necessary
	return(hapBkMap)
	
}

setTrioMissingSNP <-
function(trioDf, cord, snp1digit=FALSE, missingDigit = 0){
  ## cord has the beginning row number and col numbers for the trio with missing data

  cord = matrix(cord, ncol=4, byrow=FALSE)
  rowCt = nrow(cord)

  dataNew = trioDf
  if(!snp1digit){
    for( i in 1:rowCt){
      x.y = cord[i,1:2]
      dataNew[  x.y[1]: (x.y[1]+2),  x.y[2]: (x.y[2]+1) ] = missingDigit
      #print(paste("i=", i, " x.y=[", paste(x.y, collapse=";", sep=""), "]", sep=""))
      #print(paste(    x.y[1]: (x.y[1]+2), collapse=";", sep=""))
      #print(paste(    x.y[2]: (x.y[2]+1), collapse=";", sep=""))
    }
  }else{
    for( i in 1:rowCt){
      x.y = cord[i,3:4]
      #dataNew[  ((x.y[1]-1)*3+1): (x.y[1]*3),  x.y[2]+2 ] = 0
      #print(paste("i=", i, " x.y=[", paste(x.y, collapse=";", sep=""), "]", sep=""))
      #print(paste( ((x.y[1]-1)*3+1): (x.y[1]*3) , collapse=";", sep=""))
      #print(paste( x.y[2]+2, collapse=";", sep=""))

      dataNew[  ((x.y[1]-1)*3+1): (x.y[1]*3),  x.y[2]+2 ] = missingDigit
    }
  }
  return(dataNew)
}

signal.new <-
function(coef, type, signalStr){

  # var used internally
  dig2Code=c("00", "11", "12", "22")
  
  ifD = FALSE
  vars = NULL
  segReplace = NULL

  
  signal=binaTree.parser(str=signalStr)
  #signal =binaTree.parser.proc(str=signalStr, binaTree=binaTree, curLevel=0, first.seg=TRUE)
  vars = unique(unlist(signal$elm.vlist))
  #boolOp = signal$elm.vMa[,3]

  ori.vars = unlist(signal$elm.vlist)
  ori.match = match(ori.vars, vars)
  
  
  if(type=="geno2d"){
    ## translate into snp 2-digit coding
    pos.equalSign = sapply(vars, FUN=util.char.1stIdx, find="=")

    snpIdx = sapply(1:length(vars), FUN=function(i, name, pos) {
      ## variable name starts with g, followed by snp idx and "="
      a = substr(name[i], 2, pos[i]-1)
      a
    }, name=vars, pos=pos.equalSign)
    
    segReplace = sapply(1:length(vars), FUN=function(i, name, pos, dig2Code, snpId){
      ##  Dec08Change!!! Change the sigStr to two-digit string representation
      a = substr(name[i], pos[i]+1, pos[i]+2)
      ## change!!!change
      if(a== dig2Code[2]){
        ## less common homo is "11", (a) is dominant for the less common, (b) is recessive for the less common
        ## 11 is equivalent to recessive 
        re = paste("v", snpId[i], "b", sep="")  
      }
      if(a== dig2Code[4]){
        ## "22" is equivalent to not dominant
        re =  paste("(not v", snpId[i] ,"a)", sep="")
      }      
      if(a== dig2Code[3]){
        ## heto: "12" 
        re = paste("( (not v", snpId[i],"b) and v", snpId[i], "a )", sep="")    
      }
      re
    }, name=vars, pos=pos.equalSign, dig2Code=dig2Code, snpId = snpIdx)
  } ## if(type=="geno2d"){

  if(type=="D/R"){
    ## translate into snp 2-digit coding
    snpIdx = sapply(1:length(vars), FUN=function(i, name) {
      ## variable name starts with snp idx, followed by one character "D"/"R"
      a = substr(name[i], 1, nchar(name[i])-1)
      a
    }, name=vars)

    snpDR = sapply(1:length(vars), FUN=function(i, name) {
      ## variable name starts with snp idx, followed by one character "D"/"R"
      a = substr(name[i], nchar(name[i]), nchar(name[i]))
      a
    }, name=vars)


    ## need to take care of not
    #if (length(vars)!=length( boolOp )) stop("Error in signal.new: not unique markers in string.")
    
    segReplace = sapply(1:length(vars), FUN=function(i,  snpDR, snpId) {
      re = ""
      if(snpDR[i]== "D"){
        ## dominant:
        re = paste( "v", snpId[i] ,"a", sep="")
      }
      if(snpDR[i]== "R"){
        ## recessive:
        re = paste( "v", snpId[i] ,"b", sep="")
      }
      re
    },  snpDR=snpDR,  snpId = snpIdx)
  } ## if(type=="d/r"){
  
    ## replace the original str with new element within str
    changedStr = signalStr
    if(length(vars)>=1){
      for (i in 1:length(vars)){
        tmp.str = vars[i]
        ##print(tmp.str)
        changedStr = util.str.replace(str=changedStr,
                           replaced=vars[i],
                           new=segReplace[i],
                           replace.all=TRUE)
        ##print(model.signal.cur)
      }
    }

  if(ifD) print(changedStr)
  if( length(unlist(signal$elm.vlist))==1 & substr(changedStr,1,1)=="(")
    changedStr = substr(changedStr, 2, nchar(changedStr)-1)
  
  ## reconstruct the single
  #signal =binaTree.parser.proc(str=changedStr, binaTree=binaTree, curLevel=0, first.seg=TRUE)
  signal = binaTree.parser(str=changedStr)
  
  ## also return the snpIdx so reshuffle the column can be done later
  all.varSnp = snpIdx[ ori.match ]
  
  return(list(coef=coef, signal=signal, var.idx = all.varSnp, str=changedStr)) 

}

signalRule.2strata.build <-
function(sigStr="g9=11 and g13=11", sigType="geno2d", para=c(-5, 1)){

  ## Dec08Change!!!: allow D/R coding
  #if (sigType !="geno1d") stop( paste("sigType=", sigType, " is not implemented.", sep=""))
  if( !is.element( sigType, c("geno2d", "D/R") ) ) stop( paste("Wrong input: sigType=", sigType, ". It is not implemented.", sep=""))
  
  if( length(para)!=2) stop( paste("Wrong length of input para:", length(para), sep=""))
        
  rule = signalRule.contr()
  rule = signalRule.setInter(rule, para[1])
  sig1 = signal.new(para[2], type=sigType, signalStr = sigStr)
  rule = signalRule.addSignal(rule, sig1)

  return (rule)
}

signalRule.addSignal <-
function(rule, signal){

  rule$slist = c(rule$slist, list(signal))
  rule$snpIdx = c(rule$snpIdx, signal$var.idx)
  
  return (rule)
}

signalRule.contr <-
function(){
  rule = list()
  rule$inter = NA
  rule$slist = NULL
  rule$snpIdx = NULL
  return(rule)
}

signalRule.setInter <-
function(rule, inter){
  rule$inter = inter;
  return (rule)
}

simuHapMap.build <-
function(hapInfoFrame=NULL, genoInfoFrame){
	
	## make it can take .csv file
	if(!is.null(hapInfoFrame)){
		if(is.character(hapInfoFrame)){
			hapInfoFrame = read.csv(hapInfoFrame, header=TRUE, sep=",", as.is=TRUE)
		}
	}
	
	## make it can take .csv file
	if(is.character(genoInfoFrame)){
		genoInfoFrame = read.csv(genoInfoFrame, header=TRUE, sep=",", as.is=TRUE)
	}
	
	## check input data
	test = match(c("prekey", "seq", "freq"), colnames(genoInfoFrame))
	if (sum(is.na(test))>=1)
		stop("Input argument, genoInfoFrame, does not have the required columns, i.e., prekey, seq, and freq.") 
	
#  if(is.null(hapInfoFrame)){
#    warning("NULL value is provided for input argument, hapInfoFrame.")
#  }
	
	genoMap.info = genoInfoFrame[, match(c("prekey", "seq", "freq"), colnames(genoInfoFrame))]
	
	genoMap = dfToGenoMap(df = genoMap.info, dataHeaders = c("prekey", 
					"seq", "freq"), genotype = c("11", "12", "22"), snpBase = 1)
	
	
	if(is.null(hapInfoFrame)){
		#print("here")
		
		onlyGeno = genoInfoFrame
		## build the df first
		key = paste( onlyGeno[,1], onlyGeno[,2], sep="-")
		key = unique(key)
		
		allGenoBk = NULL
		for( i in 1:(length(key))){
			x.b = (i-1)*3+1
			x.e = i*3
			oneBk = geno2hapFreq(prekey=onlyGeno[x.b,1], block=onlyGeno[x.b,2],
					genoFreq = onlyGeno[x.b:x.e, 3], snpBase=1,
					snpSeq = i)
			allGenoBk = rbind(allGenoBk, oneBk[,c(1,4,5)])
		}
		colnames(allGenoBk)=c("key","haploytype","frequency")
		#print(allGenoBk)
		
		simuMap = bkMap.constr(data=allGenoBk, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3)
		simuSetup = list(NULL)
		simuSetup = c(simuSetup, newDf=list(allGenoBk), simuMap=list(simuMap)) 
	}else{
		bkFrame = hapInfoFrame
		bkFrame$key = paste(bkFrame[, 1], bkFrame[, 2], sep = "-h")
		
		hapBkMap = dfToHapBkMap(data = bkFrame, keyCol = ncol(bkFrame), 
				chCol = 1, blockCol = 2, expCol = 3, probCol = 4, 
				hapLenCol = 5, beginCol = 6, endCol = 7, snpBase = 1, 
				re.bf = TRUE, re.javaGUI = TRUE)
		hapBkGenoMap = bindHapBkGenoMaps(hapBkMap = hapBkMap, genoMap = genoMap)
		simuSetup = hapBkGenoMap2HapMap(hapBkGenoMap=hapBkGenoMap, reSimuMap=TRUE)
	}
	return(simuSetup)
}

snpPREFileMatchTrio <-
function(txtF, sep=" ", header = FALSE, pedCol=1, memCol=2, affectCol=6, dadCol=3, momCol=4, txt.affect=2, logF = NULL){

      if(is.character(txtF)){
        ## by default the text file has no header 
        df = read.table(file=txtF, header = header, sep=sep)
      }else{
        df = txtF
      }
      filter = df[, affectCol]==txt.affect 
      caseDf = df[filter,]

      ## it is ok for the parent to be affected
      ##controlDf = df[!filter,]
      controlDf = df

      ## don't sort the case, keep the original order
      #caseDf = caseDf[order(caseDf[,pedCol], caseDf[,memCol]),]
      #controlDf = controlDf[order(controlDf[,pedCol], controlDf[,memCol]),]

      numCase = nrow(caseDf)

      trioDf = NULL
      for( i in 1:numCase){
        caseRow = caseDf[i,]
        
        tryCatch({
             trioDf = rbind(trioDf, trio.match.single(caseRow, controlDf, pedCol, memCol, dadCol, momCol))
           },
             warning = function(warn){
                  if(!is.null(logF)){
                    #print(warn)
                    logWarn(logF, as.character(warn))
                  }else{
                    print("Throw warnings by case")
                    print(warn)
                }
        }) ## tryCatch({
      }
      return(trioDf)
}

toolbox.load <-
function(freqMaps){

	## HARD CODE!!!HARD CODE assuming the genotype is never used in the genoMap
	genoMap = dfToGenoMap(df=freqMaps$genoMap.info, dataHeaders=c("prekey", "seq", "freq"),  genotype=c("11", "12", "22"), snpBase=1)
	
    ## need to build the overall map
    if(!is.null(freqMaps$hapMap.info)){
	  if(max(freqMaps$hapMap.info$hapLen>7)==1) stop("Cannot process haplotype block with 8 or more loci.")  
      bkFrame = freqMaps$hapMap.info
      bkFrame$key = paste(bkFrame[,1], bkFrame[,2], sep="-") 

      hapBkMap = dfToHapBkMap(data=bkFrame,  keyCol=ncol(bkFrame),
                   chCol=1, blockCol=2,  expCol=3, probCol=4, hapLenCol=5, beginCol=6, endCol=7, snpBase=1, re.bf = TRUE, re.javaGUI = TRUE)
      
	  hapBkGenoMap = bindHapBkGenoMaps(hapBkMap=hapBkMap, genoMap=genoMap)
		   
      hapBkMap = hapBkGenoMap$hapBkOnlyMap
	  changedKey =  hapBkMap$keys
		   
	  newHapBkMap = hapBkMap
	  for( i in changedKey){
			   newHapBkMap =  semiAugBkFrame(newHapBkMap, key=i, probLeftOver = .01)
	  }
	  semiAugHapBkGenoMap = hapBkGenoMap
	  semiAugHapBkGenoMap$hapBkOnlyMap = newHapBkMap
    }else{
	  hapBkGenoMap = bindHapBkGenoMaps(hapBkMap=NULL, genoMap=genoMap)
	  semiAugHapBkGenoMap = hapBkGenoMap
	}
   
    ## HARD CODE!!!HARD CODE
    maxSNP=7
    idx4hapDigitAll = exIdxFromHap(maxSNP)
    exhaustHapExpAll = exhaustHapExp(maxSNP, re.str=FALSE)[[1]]
    
    maxRow = 4000
    maxTbl = matrix(NA, ncol =6, nrow = maxRow)
    maxMateTbl = matrix(NA, ncol =2, nrow = maxRow)
    
    varPrefix = paste(sample(letters[1:26], size=5), sep="", collapse="")
    varOriginalName = c("semiAugHapBkGenoMapSZ", "idx4hapDigitAll", "exhaustHapExpAll", "maxTbl", "maxMateTbl")
    appVarNames = paste(varPrefix, varOriginalName, sep="_")
    names(appVarNames) = c("freqMap", "digit", "exp", "tbl", "mateTbl")
    appVarNames = as.list(appVarNames)
    
    #print(paste("Application wise global environment var with names as", paste(appVarNames, sep="", collapse=";")))
    
    lapply(1:length(appVarNames), FUN=function(item, varNames, varList){
              assign(varNames[[item]], varList[[item]], envir=qingEnv); return(NULL)},
           varNames=appVarNames,
           varList = list(semiAugHapBkGenoMap, idx4hapDigitAll, exhaustHapExpAll, maxTbl, maxMateTbl))
    return(appVarNames)
}

trio.impu <-
function(triodd, freq,  impu.missingOnly=TRUE){
  dir=NULL

  dig1Code=c(NA,0,1,2)


  ## check number of loci match with map or not
  if( (ncol(triodd)-2) != nrow(freq$genoMap.info)/3 )
    stop("The number of SNPs in the frequency file does not equal the number of SNPs in trio dataset.")

  
  #if( subInF!="PPC" ) stop("Subject in the data should follow Father, Mother and child sequence.")
  
  if(is.character(triodd) ) {
    data = read.csv(triodd, header=FALSE)
  }else{
    data = triodd
  }
  code.Factor = apply( data[,c(-1,-2)], 2, FUN=is.factor)
  if(max(code.Factor)==1) stop("Genotypes cannot be factors.")

  code.Factor = apply( data[,c(-1,-2)], 2, FUN=is.numeric)
  if(min(code.Factor)==0) stop("Genotypes must be numerical.")
  
  ## check the coding for trio
  code.Check = unique(unlist(data[,c(-1,-2)]))
  
  if( sum(is.na(code.Check))==0 ){
    if (impu.missingOnly) {
      print("No missing genotype. Return original data.")
      return(triodd[,-c(1,2)])
    }
  }

  code.Left = code.Check[!is.na(code.Check)]

  code.Match = match(code.Left, dig1Code[2:4])
  if( sum(!is.na(code.Match)) != length(code.Match)){
    stop("Genotypes are not coded as 0, 1, nor 2")
  }

  toolboxNames = toolbox.load(freqMaps=freq)
  
  if(!is.null(  get(toolboxNames$freqMap, envir = qingEnv)$hapBkOnlyMap )){
	  bk.ssize =  get(toolboxNames$freqMap, envir = qingEnv)$hapBkOnlyMap$bkSnpLens
	  if (max(bk.ssize>=8)==1) stop("At least one haplotype block has 8 or more SNPs in the block. Method fails.")	  
  }

  # inside the function, assume 0 1 3 2 for NA, homo, hetero, homo, and 0 1 2 for NA, allele1, allele2
  dig1Default = c(0, 1, 3, 2)
  ## if the given code for 1-digit is not (0 1, 3, 2), change it.

  if (is.na(dig1Code[1])){
    ttt = data[,c(-1,-2)]
    ttt[is.na(ttt)]=max(dig1Code, na.rm=TRUE)+1

    data[,c(-1,-2)]=ttt
    dig1Code[1]=max(dig1Code, na.rm=TRUE)+1
  }

  
  if( sum(dig1Default == dig1Code)!=4 ){
    data.geno = apply(data[, c(-1, -2)], 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=dig1Default)
    data = cbind(data[, 1:2], data.geno)
  }
  ##TODO!!! confirm the number
  #print(toolboxNames)
  
  bkCt = nrow(get(toolboxNames$freqMap, envir = qingEnv)$genomeMarkerInfo)

  if(impu.missingOnly){

    ## get the hapPair without saving.
    imputed = impuBk.scheduler(raw=data, idx=1:bkCt, job=1,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=TRUE, dig1Code=dig1Default, dig2Code=0:2,
                reType=FALSE, reHap=NULL, logF=NULL, logErr="")
    if( sum(dig1Default == dig1Code)!=4 ){
      imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
    }
     return(imputed)
  }

  #print("get impuBkTDT.scheduler")
  ## get the hapPair without saving.
  imputed = impuBkTDT.scheduler(raw=data, idx=1:bkCt, job=1,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=TRUE, dig1Code=dig1Default, dig2Code=0:2,
                reType=FALSE, reHap=NULL, logF=NULL, logErr="")

  if( sum(dig1Default == dig1Code)!=4 ){
     imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
  }
  return(imputed)
 
}

trio.impuDev <-
function(trio1digit, freqMaps, bk.seq=NULL, dig1Code=0:3, subInF = "PPC", dir=NULL, job=1, prefix="", impu.missingOnly=TRUE){

  if( subInF!="PPC" ) stop("Subject in the data should follow Father, Mother and child sequence.")
  
  if(is.character(trio1digit) ) {
    data = read.csv(trio1digit, header=FALSE)
  }else{
    data = trio1digit
  }
  
  toolboxNames = toolbox.load(freqMaps=freqMaps)
  
  if(!is.null(  get(toolboxNames$freqMap, envir = qingEnv)$hapBkOnlyMap )){
	  bk.ssize =  get(toolboxNames$freqMap, envir = qingEnv)$hapBkOnlyMap$bkSnpLens
	  if (max(bk.ssize>=8)==1) stop("At least one haplotype block has 8 or more SNPs in the block. Method fails.")	  
  }  
  
  
  # inside the function, assume 0 1 3 2 for NA, homo, hetero, homo, and 0 1 2 for NA, allele1, allele2
  dig1Default = c(0, 1, 3, 2)
  ## if the given code for 1-digit is not (0 1, 3, 2), change it.

  if( sum(dig1Default == dig1Code)!=4 ){
    data.geno = apply(data[, c(-1, -2)], 1:2, FUN= util.vec.replace, orignal = dig1Code, replaceBy=dig1Default)
    data = cbind(data[, 1:2], data.geno)
    
  }
  ##TODO!!! confirm the number
  #print(toolboxNames)
  
  bkCt = nrow(get(toolboxNames$freqMap, envir = qingEnv)$genomeMarkerInfo)

  #print(data[1:3, 1:10])
  if(is.null(bk.seq)) bk.seq = 1:bkCt
  if(impu.missingOnly){

    ## get the hapPair without saving.
    imputed = impuBk.scheduler(raw=data, idx=bk.seq, job=job,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=TRUE, dig1Code=dig1Default, dig2Code=0:2,
                reType=FALSE,
                reHap=qp(prefix, "impuHap_bk", bk.seq[1]),
                logF=NULL,
                logErr=qp(prefix, "impuErr_bk", bk.seq[1])
      )
    if( sum(dig1Default == dig1Code)!=4 ){
      imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
    }
     return(imputed)
  }

  ## get the hapPair without saving.
  imputed = impuBkTDT.scheduler(raw=data, idx=bk.seq, job=job,
                toolname=toolboxNames,
                freqMaps=NULL, dir=dir,
                is.1digit=TRUE, dig1Code=dig1Default, dig2Code=0:2,
                reType=FALSE,
                reHap=qp(prefix, "impuHap_bk", bk.seq[1]),
                logF=NULL,
                logErr=qp(prefix, "impuErr_bk", bk.seq[1])
  )

  if( sum(dig1Default == dig1Code)!=4 ){
     imputed = apply(imputed, 1:2, FUN= util.vec.replace, orignal = dig1Default, replaceBy=dig1Code)
  }
  return(imputed)
 
}

trio.match.single <-
function(caseRow, controlDf, pedCol, memCol, dadCol, momCol){
   dadRow = which(controlDf[,pedCol]==unlist(caseRow[pedCol]) & controlDf[,memCol]==unlist(caseRow[dadCol]))
   momRow = which(controlDf[,pedCol]==unlist(caseRow[pedCol]) & controlDf[,memCol]==unlist(caseRow[momCol]))

   if( length(dadRow)!=1 ){
     warning(paste("\n The individual is affected, but we have no data on the father. Case information::\n",
                   paste(wrComTbl(caseRow[1:7], colNames=colnames(caseRow)[1:7]), collapse="\n"), sep=""))
     return(NULL)
   }
   if( length(momRow)!=1 ){
     warning(paste("\n The individual is affected, but we have no data on the mother. Case information::\n",
                   paste(wrComTbl(caseRow[1:7], colNames=colnames(caseRow)[1:7]), collapse="\n"), sep=""))
     return(NULL)
   }
  
   re = rbind(controlDf[dadRow,], controlDf[momRow,], caseRow)
   return(re)
 }

trio.simu.direct <-
function (bkMap, rule, caseNo, datasetCt=1, infoS="simuDirInfo", ddir=NULL, baseObj.saveFN=NULL, baseObj.name=NULL,  reControl=FALSE, dig1Code=0:3 ){

  info = bkMap.HRCB.LRCB.split(bkMap, rule)

  bkMapNS = info$bkMapNS

  varPrefix = paste(sample(letters[1:26], size=10), sep="", collapse="")
  varOriginalName = "stepstone"
  finalUse = paste(varPrefix, varOriginalName, sep="_")
  #print(qp("finalUse=", finalUse))

  if(is.null(baseObj.saveFN)){
    ## if the spStrata is not saved earlier, need to generate
    if(!is.null(baseObj.name)){
      ## then we can choose to save it
     
      print(qp("No stepstone object is saved. Choose to generate and save the object as:", baseObj.name))
      
      matingTbInfo = bkMap.HRCB.famMap(info$bkMapS, rule, newColName=info$newColName,  ifS=infoS, baseName=baseObj.name)

      #finalUse =  baseObj.name
      assign(finalUse, matingTbInfo, envir=qingEnv)
      
    }else{
      ## or not save it, just used. Not recommend.
      print(qp("No stepstone object is saved. Choose to generate but not save the object."))
      matingTbInfo = bkMap.HRCB.famMap(info$bkMapS, rule, newColName=info$newColName,  ifS=infoS, baseName=NULL)

      #finalUse = matingTblInfo
      assign(finalUse, matingTbInfo, envir=qingEnv)
    }
    
  }else{
    ## if the spStrata is saved earlier, just load the saved one
   
    if( is.character(baseObj.saveFN)){
      print(qp("Stepstone object is saved before as:", baseObj.saveFN, ". Do not need to generate it."))
      tmp = load(paste(baseObj.saveFN, ".RData", sep=""))
      assign(finalUse, get(tmp[1]), envir=qingEnv)
      
    }else{
      print(qp("Stepstone object is given."))
      assign(finalUse, baseObj.saveFN, envir=qingEnv)
    }
  }

  
  simuTrio = rep(list(NA), length=datasetCt)

  if(is.null(ddir) &  (datasetCt!=1)){
    print( "Request to generate multipe datasets, but do not give path for save. Will return as a list.")
  }

  for( i in 1:datasetCt){
      if(is.null(infoS)){
        trioData1 = HRCB.famMap.spTrio(caseNo=caseNo, matingTbName=finalUse,
          ifS =NULL , reControl=reControl)
      }else{
        trioData1 = HRCB.famMap.spTrio(caseNo=caseNo, matingTbName=finalUse,
          ifS = qp(infoS, "_", i) , reControl=reControl)
      }
      
      trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, reControl=reControl)
      
      simuTrioData = trioMerge(trioData1, trioData2, colName1=info$newColName, colName2=info$noCausalColName)

      if( is.null(dim(simuTrioData))){
        snp1recode =  util.vec.replace(simuTrioData, orignal = c(0,1,3,2), replaceBy= dig1Code)
      }else{
        snp1recode =  apply(simuTrioData, 1:2, FUN= util.vec.replace, orignal = c(0,1,3,2), replaceBy=dig1Code)
      }

      if(is.null(ddir)){
        simuTrio[[i]]=snp1recode
      }else{
         if (ddir==""){
           save(snp1recode, file=qp("trioDirSimu_", i, ".RData"))
         }else{
           save(snp1recode, file=file.path(ddir, qp("trioDirSimu_", i, ".RData")))
         } 
      }
    }

  if(is.null(ddir)){
    return(simuTrio)
  }else{
    return(NULL)
  }
}

trio.simu.proposed <-
function(bkMap, rule, caseNo, datasetCt=1, infoS="simuPropInfo", exInfoS="exSimuPropInfo", ddir=NULL, startIdx=NULL, spStrata.saveFN = NULL,  spStrata.name=NULL, reControl =FALSE, dig1Code=0:3 ){

  if (is.null(ddir)) {
    infoS=NULL
    exInfoS = NULL
  }
  ifD = FALSE
  
  info = bkMap.HRCB.LRCB.split(bkMap, rule)

  bkMapNS = info$bkMapNS

  varPrefix = paste(sample(letters[1:26], size=10), sep="", collapse="")
  varOriginalName = "stepstone"
  finalUse = paste(varPrefix, varOriginalName, sep="_")
  if(ifD) print(qp("finalUse=", finalUse))

  if(is.null(spStrata.saveFN)){
    ## if the spStrata is not saved earlier, need to generate
    if(!is.null(spStrata.name)){
      ## then we can choose to save it
      spStrata = bkMap.HRCB.Esp1Rule.Base(bkMap, rule, baseName=spStrata.name)
      if(ifD) print(qp("No stepstone object is previously saved. Generate and save the object as:", spStrata.name, ".RData"))

      #finalUse = spStrata.name
      assign(finalUse, spStrata, envir=qingEnv)
      
    }else{
      ## or not save it, just used. Not recommend.
      spStrata = bkMap.HRCB.Esp1Rule.Base(bkMap, rule, baseName=NULL)

      #finalUse = spStrata
      assign(finalUse, spStrata, envir=qingEnv)
      #print(qp("No stepstone object is previously saved. Generate but not save the object."))
    }
    
  }else{
    ## if the spStrata is saved earlier, just load the saved one
    
    if( is.character(spStrata.saveFN)){
      if(ifD) print(qp("Stepstone object is previously saved as:", spStrata.saveFN, ".RData. Do not need to generate it."))

      tryCatch({tmp = load(paste(spStrata.saveFN, ".RData", sep=""))
               assign(finalUse, get(tmp[1]), envir=qingEnv)  },
             error = function(e){
               stop(paste("Cannot open step-stone file '", spStrata.saveFN, ".RData'.", sep=""))
             })

    }else{
      #print(qp("Stepstone object is given."))
      #finalUse = spStrata.saveFN
      assign(finalUse, spStrata.saveFN, envir=qingEnv)
    }
  }
  
  preObj = bkMap.HRCB.Esp1Rule.genoSeq(bkMap, rule, re.probOnly = FALSE)

  
  simuTrio = rep(list(NA), length=datasetCt)

  if(is.null(ddir) &  (datasetCt!=1))
    if(ifD) print( "Request to generate multiple datasets, and will return as a list.")

  for( i in 1:datasetCt){

        if(!is.null(infoS)){
          trioData1 = HRCB.Esp1Rule.spTrioOnBase(bkMap=NULL, preObj=preObj, spStrata=finalUse,
                                       rule=rule, caseNo,
                                       ifS = qp(infoS, "_", i) , reControl=reControl)
        }else{
          trioData1 = HRCB.Esp1Rule.spTrioOnBase(bkMap=NULL, preObj=preObj, spStrata=finalUse,
                                       rule=rule, caseNo,
                                       ifS = NULL , reControl=reControl)
        }
    #tryCatch({
        if(!is.null(exInfoS)){
          trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, ifS = qp(exInfoS, "_", i), reControl=reControl)
        }else{
          trioData2 = bkMap.LRCB.spTrio(bkMapNS, caseNo=caseNo, ifS = NULL, reControl=reControl)
        }
     #}, warning=function(w){print(w)})
        
        simuTrioData = trioMerge(trioData1, trioData2, colName1=info$newColName, colName2=info$noCausalColName)

        ## get back to common coding scheme, 3 for heter
        if( is.null(dim(simuTrioData))){
          snp1recode =  util.vec.replace(simuTrioData, orignal = c(0,1,3,2), replaceBy= dig1Code)
        }else{
          snp1recode =  apply(simuTrioData, 1:2, FUN= util.vec.replace, orignal = c(0,1,3,2), replaceBy=dig1Code)
        }

        ## adding other stuff
        if(ifD) print("A")
        if(ifD) print(dim(snp1recode))
        famid = rep(1:caseNo, each=3)
        pid = rep(1:3, times=caseNo)
        trio.snp = cbind(famid, pid, snp1recode)
        colnames(trio.snp) = c("famid", "pid", paste("snp", 1:(ncol(snp1recode)), sep=""))
        
        if(is.null(ddir)){
          simuTrio[[i]]=trio.snp
        }else{

         if (ddir==""){
           if (is.null(startIdx)){
             save(trio.snp, file=qp("trioSimu_", i, ".RData"))
           }else{
             save(trio.snp, file=qp("trioSimu_", startIdx+i-1, ".RData"))
           }
         }else{
           if (is.null(startIdx)){
             save(trio.snp, file=file.path(ddir, qp("trioSimu_", i, ".RData")))
           }else{
             save(trio.snp, file=file.path(ddir, qp("trioSimu_", startIdx+i-1, ".RData")))
           }
         } 
       }

  }

  if(is.null(ddir)){
    return(simuTrio)
  }else{
    return(NULL)
  }

}

trio.simuDev <-
function(bkMap=NULL,  sigStr="g9=11 and g13=11", sigType, para=c(-1, .5),  caseNo=10, datasetCt=1,  dig1Code=0:3, ddF=NULL, startIdx=NULL, stepstone.saveFN=NULL, stepstone.name=NULL, spSupHap.namePrefix=NULL, verbose=TRUE, reControl=FALSE){

  print(paste("Try to simulated trio data: # datasets=", datasetCt, "; # trio =", caseNo, sep=""))


  if(!is.null(ddF)){
    if(ddF==""){
      print(paste("Simulated data file(s) and other result file(s) are saved under current working directory."))
      infoS.s = spSupHap.namePrefix
      exInfoS.s = paste(spSupHap.namePrefix, "LRCB", sep="")
    }else{
      print(paste("Create directory:", ddF, ", where simulated data file(s) and other result file(s) are saved.", sep=""))
      dir.create(path=ddF, showWarnings = TRUE)
      infoS.s = file.path(ddF, spSupHap.namePrefix)
      exInfoS.s = file.path(ddF, paste(spSupHap.namePrefix, "LRCB", sep=""))
    }

    if(is.null(spSupHap.namePrefix))  infoS.s = NULL
  }else{
    print(paste("Value for argument ddF is NULL. Function will return data as a list, ignoring input for argument, spSupHap.namePrefix"))
    infoS.s = NULL
    exInfoS.s = NULL
  }

  if(is.null(bkMap)){
    print(paste("Value for argument bkMap is NULL. Function will use package default object, simuBkMap."))
    data(simuBkMap, envir = environment())
    bkdata = get("simuBkMap")
    bkMap = bkMap.constr(data=bkdata, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=1:2)
  }

  #print("Info on data object for haplotype block frequencies:")
  #print(str(bkMap, max.level=1))

  ##  Dec08Change!!! allow D/R coding
  rule =signalRule.2strata.build(sigStr=sigStr, sigType=sigType, para=para)

  #print("Info on data object for risk factor, signalRule:")
  #print(rule$slist[[1]]$str)
  if(is.null(spSupHap.namePrefix)){
    infoS.s = NULL
    exInfoS.s = NULL
  }else{
    if(!is.null(ddF)){
      if(ddF==""){
        infoS.s = spSupHap.namePrefix
        exInfoS.s = paste(spSupHap.namePrefix, "LRCB", sep="")
      }else{
        infoS.s = file.path(ddF, spSupHap.namePrefix)
        exInfoS.s = file.path(ddF, paste(spSupHap.namePrefix, "LRCB", sep=""))
      }
    }else{
      infoS.s = NULL
      exInfoS.s = NULL
    }
  }

  ptm=proc.time()


  trioData = trio.simu.proposed(bkMap=bkMap,
                    rule=rule,
                    caseNo=caseNo,
                    datasetCt=datasetCt,
                    infoS=infoS.s,
                    exInfoS=exInfoS.s,
                    ddir=ddF,
                    startIdx=startIdx,
                    spStrata.saveFN = stepstone.saveFN,
                    spStrata.name = stepstone.name,
                    reControl =reControl,  dig1Code=dig1Code 
                    )
  
  if(verbose) print(paste("Time used to generate ", datasetCt,  " dataset(s).", sep=""))
  if(verbose) print(proc.time() - ptm)
  gc.e = gc()
 
  if(verbose) print("Info on memory usage:")
  if(verbose) print(gc.e)

  print(paste("Finished generating ", datasetCt,  " dataset(s).", sep=""))

  return(trioData)
}

trio.simuOLD <-
function(bkMap=NULL,  interaction="9R and 13R", alpha, beta, n=10, rep=1,  stepstone.saveFN=NULL, stepstone.name=NULL, verbose=TRUE){

  sigType="D/R"
  para = c(alpha, beta)
  dig1Code=c(4,0,1,2)
  ddF = NULL
  spSupHap.namePrefix=NULL
  alleleCode=1:2

  print(paste("Try to simulated trio data: # datasets=", rep, "; # trio =", n, sep=""))
  
  if(!is.null(ddF)){
    if(ddF==""){
      print(paste("Simulated data file(s) and other result file(s) are saved under current working directory:",  getwd()))
      infoS.s = spSupHap.namePrefix
    }else{
      print(paste("Create directory:", ddF, ", where simulated data file(s) and other result file(s) are saved.", sep=""))
      dir.create(path=ddF, showWarnings = TRUE)
      infoS.s = file.path(ddF, spSupHap.namePrefix)
    }

    if(is.null(spSupHap.namePrefix))  infoS.s = NULL
  }else{
    #print(paste("Value for argument ddF is NULL. Function will return data as a list, ignoring input for argument, spSupHap.namePrefix"))
    infoS.s = NULL
  }

  if(is.null(bkMap)){
    print(paste("Value for argument bkMap is NULL. Function will use package default object, simuBkMap."))
    data(simuBkMap, envir = environment())
    bkdata = get("simuBkMap")
    bkMap = bkMap.constr(data=bkdata, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=alleleCode)
  }

  #print("Info on data object for haplotype block frequencies:")
  #print(str(bkMap, max.level=1))
   ##  Dec08Change!!! allow D/R coding
  rule =signalRule.2strata.build(sigStr=interaction, sigType=sigType, para=para)

  #print("Info on data object for risk factor, signalRule:")
  #print(rule$slist[[1]]$str)

  ptm=proc.time()

  trioData = trio.simu.proposed(bkMap=bkMap,
                    rule=rule,
                    caseNo=n,
                    datasetCt=rep,
                    infoS=infoS.s,
                    exInfoS=NULL,
                    ddir=ddF,
                    spStrata.saveFN = stepstone.saveFN,
                    spStrata.name = stepstone.name,
                    reControl =FALSE,  dig1Code=dig1Code 
                    )
  
  if(verbose) print(paste("Time used to generate ", rep,  " dataset(s).", sep=""))
  if(verbose) print(proc.time() - ptm)
  gc.e = gc()
 
  if(verbose) print("Info on memory usage:")
  if(verbose) print(gc.e)

  print(paste("Finished generating ", rep,  " dataset(s).", sep=""))

  return(trioData)
}

trioFile.proc <-
function(data, key.prefix="", bk.sizes=NULL, dig2Code=0:2, dig1Code=c(0,1,3,2),
		action = c("missingReport", "Mendelian check", "freq estimate"), ... ){
	
	
	sep=""
	header =FALSE
	pedCol=1
	memCol=2
	snpIdxRange = c(3, ncol(data))
	is.1digit=TRUE
	
#   if(!is.null(txtF)){
#     if(is.character(txtF)){
#          ## by default the text file has no header 
#          #trioTwoDigit = read.csv(file=data, header = header, sep=sep)
#          stop("The argument, data, cannot be a string.")
#     }else{
#          trioTwoDigit = data
#     }
#   }
	
	trioTwoDigit = data
	
	snpStartLeftIndex=snpIdxRange[1]
	snpEndRightIndex =ncol(trioTwoDigit) - snpIdxRange[2] + 1
	
	
	if(!is.1digit){
		genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
		## check input, if is.1digit=TRUE, dig2Code must be c(0,1,2), if is.1digit=FALSE, dig1Code must be c(NA, 0,1,2,)
		code.Factor = apply( genos, 2, FUN=is.factor)
		if(max(code.Factor)==1) stop("Genotypes cannot be factors.")
		
		code.Factor = apply( genos, 2, FUN=is.numeric)
		if(min(code.Factor)==0) stop("Genotypes must be numerical.")
		
		code.Check = sort(unique(unlist(genos)))
		
		if( min(code.Check==dig2Code)==0 ){
			stop("Trio data is in 2-digit coding, but alleles are not represented by 0, 1, and 2.")
		}
		
		
		snpNum = ncol(genos)/2
		snp1digit = exchangeDigit(ma=genos,
				cols=c(1,snpNum*2), dig1Code=dig1Code, dig2Code =dig2Code, action=c("2to1"))
		
	}else{
		genos  = trioTwoDigit[, snpStartLeftIndex:(ncol(trioTwoDigit)-snpEndRightIndex+1)]
		
		## check input, if is.1digit=TRUE, dig2Code must be c(0,1,2), if is.1digit=FALSE, dig1Code must be c(NA, 0,1,2,)
		code.Factor = apply( genos, 2, FUN=is.factor)
		if(max(code.Factor)==1) stop("Genotypes cannot be factors.")
		
		code.Factor = apply( genos, 2, FUN=is.numeric)
		if(min(code.Factor)==0) stop("Genotypes must be numerical.")
		
		## check the coding for trio
		code.Check = unique(unlist(genos))
		
		code.Left = code.Check[!is.na(code.Check)]
		
		code.Match = match(code.Left, dig1Code[2:4])
		if( sum(!is.na(code.Match)) != length(code.Match)){
			stop("Non-missing genotypes are not coded as 0, 1, and 2")
		}
		
		snpNum = ncol(genos)
		snp1digit = genos
		## change to inside Qing's coding scheme, now assume the input code is c(NA, 0,1,2), no change, and no checking
		
#      if(min(dig1Code==c(0,1,2,3))==0){
#        snp1digit = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code, replaceBy=c(0,1,2,3))
		
	}
	
	snp1digit.inside = apply(snp1digit, 1:2,  FUN= util.vec.replace, orignal = dig1Code,
			replaceBy=c(0,1,3,2))
	
	re = list()
	
#   if(is.element("formTrio1digit", action)){
#     trio1digit = cbind( trioTwoDigit[, c(pedCol, memCol)], snp1digit)
#     re = c(re, trio1digit=list(trio1digit))
#   }
	
	if(is.element("missingReport", action)){
		## check necessary parameters
		if( is.element(snpStartLeftIndex, c(pedCol, memCol))){
			warning("Cannot report missing information. The argument, snpStartLeftIndex, is one of the special column for the trio file.")
		}else{
			
			if(is.1digit){
				#print(here)
				#print(snp1digit.inside)
				#print(snpStartLeftIndex)
				#print(snpEndRightIddex)
				missSNPPos = findMissing(df=snp1digit.inside, is.1digit=TRUE, snpStartLeftIndex=snpStartLeftIndex-2,
						snpEndRightIndex=snpEndRightIndex, dig1Code=c(0,1,3,2), dig2Code=dig2Code)
			}else{
				missSNPPos = findMissing(df=trioTwoDigit, is.1digit=FALSE, snpStartLeftIndex=snpStartLeftIndex,
						snpEndRightIndex=snpEndRightIndex, dig1Code=NULL, dig2Code=dig2Code )
			}
			
			if(is.null(missSNPPos)){
				
				re = c(re, missIdx = list(NULL))
			}else{
				re = c(re, missIdx = list(missSNPPos))
			}
		}
	}
	
	if(is.element("Mendelian check", action)){
		tmpDigit=2
		if(is.1digit) tmpDigit = 1
		snpTrio = matrix(snp1digit.inside, nrow=3, byrow=FALSE)
		
		MedErr = matrix(NA, ncol=4, nrow=ncol(snpTrio))
		colnames(MedErr)=c("y", "x", "trio", "SNP")
		tryCatchEnv = new.env(parent=qingEnv)
		assign("MedErr.ct", 0, envir=tryCatchEnv)	
		assign("MedErr", MedErr, envir=tryCatchEnv)

		trioCt = nrow(snp1digit)/3
		
		#print(str(snpTrio))
		for ( i in 1:ncol(snpTrio)){
			tryCatch({
						tt = checkMendelianError(codedSNPTrio=snpTrio[,i], snpCoding=c(0,1,2,3))
					}, error = function(e){
						a = get("MedErr.ct", envir=tryCatchEnv)
						a = a+1
						assign("MedErr.ct", a, envir=tryCatchEnv)
						tttx = ceiling(i/trioCt)
						ttty = i%%trioCt
						if(ttty==0) ttty = trioCt
						b = get("MedErr", envir=tryCatchEnv)
						b[a,] = c( (ttty-1)*3+1,  (tttx-1)*tmpDigit+1+snpStartLeftIndex-1,  ttty,    tttx)
						assign("MedErr", b, envir=tryCatchEnv)
						#print( MedErr[MedErr.ct,,drop=FALSE] )
					})          
		}
		MedErr.ct = get("MedErr.ct", envir=tryCatchEnv)
		MedErr = get("MedErr", envir=tryCatchEnv)
				
		if(MedErr.ct==0) {
			MedErr=NULL
			#print("No Mendelian error.")
		}else{
			#print("Found Mendelian error(s).")
		}
		re = c(re, MedErr=list(MedErr[1:MedErr.ct,,drop=FALSE]))
		
	}
	
	if(is.element("freq estimate", action)){
		## check necessary parameters
		if(is.null(bk.sizes)) warning("Cannot provide frequencies estimation. The argument, bk.sizes, is missing.")
		
		if( is.element(snpStartLeftIndex, c(pedCol, memCol))){
			warning("Cannot provide frequencies estimation. The argument, snpStartLeftIndex, is one of the special column for linkage file.")
		}else{
			
			
			tmp = ncol(trioTwoDigit)
			parTwoDigit = getBackParentGeno(trioDf=trioTwoDigit, famCol=pedCol, memCol=memCol,
					snpIdx=snpStartLeftIndex:(tmp-snpEndRightIndex+1), re.child=FALSE, prefix=NULL)
			#print(dim(parTwoDigit))
			
			if( is.null(bk.sizes) ){
				warning("Cannot provide frequencies estimation. The argument, bk.sizes, is not provided for user-specify option.")
			}else{
				
				map=freqmap.reconstruct(data=parTwoDigit, cols=c(3, ncol(parTwoDigit)), loci.ct=bk.sizes, is.1digit=is.1digit,
						dig1Code=dig1Code, dig2Code = dig2Code, key.prefix=key.prefix, start.base=1, ...)
				
				re = c(re, freq=list(map))
			}
			
		} ##if( is.element(snpStartLeftIndex, c(pedCol, memCol, affectCol, dadCol, momCol))){
		
	}
	
	return(re)
	
}

trioMerge <-
function(trioData1, trioData2, colName1, colName2,snpCoding=0:3){

  ifD = FALSE

  # print("trioMerge")
  # print(colName1)
  # print(colName2)
  # 
  # print(str(trioData1))
  # print(str(trioData2))
  # 
  # print(trioData1[1:6, 1:5])
  # print(trioData2[1:6, 1:5])

  if( sum(snpCoding == (0:3))!=4)
    stop ("Function use 1-digit coding for genotypes as the following, integer 1 to 3 to
                    represent the three genotypes: less common homozygous, common homozygous, and heterzygous.
                     And zero for missing value")

  trioData = cbind(trioData1, trioData2)

  ## note this is 1-digit coding.
  colName1.first = colName1[ seq.int(from=1, to=length(colName1), by=2) ]
  colName2.first = colName2[ seq.int(from=1, to=length(colName2), by=2) ]

  
  trioCol = c(colName1.first, colName2.first)

  id.all = 1:(length(trioCol)/2)
  
  ori.col = paste("v", id.all, "a", sep="")
  #print(ori.col)

  shuffle.order = match(ori.col, trioCol)

  #print( cbind(trioCol, ori.col, shuffle.order))

  trioData.shuffled = trioData[, shuffle.order]

  #print(shuffle.order)
  #print(colnames(trioData.shuffled))

  if(ifD) print(str(trioData.shuffled))
  
  return(trioData.shuffled)

}

txtToGenoMap <-
function(txtF, delim=",", dataHeaders=NULL, genotype=c("11", "12", "22"), snpBase=1){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = read.csv(file=txtF, header = TRUE, sep=delim, as.is=TRUE)
	}else{
		## assume that column names of the data is passed from outside
		data = read.csv(file=txtF, header = FALSE, sep=delim, as.is=TRUE)
		colnames(data) = dataHeaders
	}
	
	genoMap = dfToGenoMap(df=data, dataHeaders=dataHeaders,  genotype=genotype, snpBase=snpBase)
	
	return(genoMap)
	
}

txtToHapBkMap <-
function(txtF, delim=",", dataHeaders=NULL, sorted=FALSE, ...){
	##txtF = "tblBlock.csv"
	##dataHeaders = NULL
	
	if(is.null(dataHeaders)){
		## assume that txtF has the first row as column names
		data = read.csv(file=txtF, header = TRUE, sep=delim, na="missing")
	}else{
		## assume that column names of the data is passed from outside
		data = read.csv(file=txtF, header = FALSE, sep=delim, na="missing")
		colnames(data) = dataHeaders
	}
	
	m = match(c("ch", "block", "hap", "freq", "hapLen","markers_b", "markers_e"), 
			colnames(data), 0)
	
	## assume except for markers_b and marekers_e, other variable should be presented
	if(min(m[1:5])==0) stop("One or more required variable(s) missing")
	
	## cast the data into the right format
	for ( i in m ){
		if(class(data[,i])=="factor") data[,i]=as.numeric(as.character(data[,i]))
	}
	
	## create key as a combination of chromosome and block
	key = paste(data[,m[1]], data[,m[2]], sep="-")
	
	
	df = cbind(key, data[,m])
	
	if(!sorted){
		df = df[ order(df$ch, df$block),]
	}
	
	hapBkMap = NULL
	if(m[6]==0){
		hapBkMap = dfToHapBkMap(df, keyCol=1, chCol=2, blockCol=3,
				expCol=4, probCol=5, hapLenCol=6, ...)
	}else{
		hapBkMap = dfToHapBkMap(df, keyCol=1, chCol=2, blockCol=3,
				expCol=4, probCol=5, hapLenCol=6,
				beginCol=7, endCol=8, ...)
	}
	return(hapBkMap)             
}

util.array.rmEmptyStr <-
function(vec){
  lens= nchar(vec)
  re = vec[lens>=1]
  return(re)
}

util.array3d.2matrix <-
function(arr, dimLevel=dimnames(arr)[[3]], showDim3 = FALSE, re.num=TRUE, appAsRow = FALSE){
  re = NULL   
  for( i in 1:length(dimLevel)){
    ma = arr[,,i]
    if(showDim3){
      if(re.num){
        ma = cbind(ma, rep(as.numeric(dimLevel[i]), nrow(ma)))
      }else{
        mat = data.frame(ma)
        dimnames(mat)=dimnames(ma)
        ma = cbind(mat, rep(dimLevel[i], nrow(ma)))
      }
    }
    if(appAsRow){
      re = rbind(re, ma)
    }else{
      re = cbind(re, ma)
    }
  }
  return(re)
}

util.char.1stIdx <-
function(str, find, match=TRUE, is.array=FALSE){
   
   if(!is.array){
     charArray = util.str.2CharArray(str, len=nchar(str))
   }else{
     charArray=str
   }
   if(match){
     pos = which(charArray==find)
   }else{
     pos = which(charArray!=find)

     #print(pos)
   }
   if( length(pos)>=1){
     if(match){
       return(pos[1])
     }else{
       return(pos[1]-1)
     }
   }else{
     return (0)
   }
}

util.char.1stStrIdx <-
function(str, find){

   result = qstrsplit(str, find, re.1st=TRUE)
   if(length(result)==1){
     # NA==not found it
     if(is.na(result)) return(0)
     # NULL==contain find only
     if(is.null(result)) return(1)
   }

   return(nchar(result[1])+1)
}

util.dataframe.findColNo <-
function(data, colNameSearched){
  colNms = colnames(data)
  matchIdx = match(colNameSearched, colNms)
  return(matchIdx)
}

util.dataframe.merge <-
function(list, vertical=TRUE){

  len = length(list)
  data = NULL
  for( i in 1:len){
    cur = list[[i]]
    if(vertical){
      data = rbind(data, cur)
    }else{
      data = cbind(data, cur)
    }
  }
  return(data)
}

util.dataframe.round <-
function(vecData, keepZero=FALSE, na.replace=NULL, digits=getOption("digits"), ...){

  re = sapply(vecData, FUN=function(x, na.rep, digits, ...){
                                num = x
                                if(is.numeric(num)){

                                    if(is.na(num) & (!is.null(na.rep))){
                                      num = na.replace
                                    }else if(is.na(num) & is.null(na.rep)){
                                      num = NA
                                    }else{
                                      if(keepZero){
                                        if(digits==0){
                                          num = format(num, nsmall=0)
                                        }else{
                                          num = format(num, digits=digits, nsmall=digits)
                                        }
                                      }else{
                                        num = round(num, digits=digits, ...)
                                      }
                                    }
                                                                 
                                }else{
                                  if(num=="NA" & is.null(na.rep)){
                                    num = as.character(num)
                                  }else if(num=="NA" & (!is.null(na.rep))){
                                    num = na.rep
                                  }else{
                                    num = num
                                  }
                                }
                                  
                              }, na.rep = na.replace, digits=digits, ...)
  return(re)

}

util.findBrace <-
function(str){

  #left ==0, right=1
  brace=-1
  idx=0
  
  arr = util.str.2CharArray(str, len=nchar(str))
  l.find =util.char.1stIdx(arr, "(", match=TRUE, is.array=TRUE)
  r.find =util.char.1stIdx(arr, ")", match=TRUE, is.array=TRUE)

  if(l.find==0){
    if(r.find!=0)
      return (list(brace=1, idx=r.find))
    else
      return(NULL)
  }
  if(r.find==0){
    if(l.find!=0)
      return (list(brace=0, idx=l.find))
    else
      return(NULL)
  }
  if(l.find<r.find) return (list(brace=0, idx=l.find))
  if(l.find>r.find) return (list(brace=1, idx=r.find))
  
}

util.it.smallLargeIdx <-
function(len, keep.same=FALSE){

      if(keep.same){
               init = 1:len
                      idx.ma = matrix(NA, ncol=2 , nrow=(len^2-len)/2+len )
                      idx.init = unlist(lapply(init, FUN=function(i, ct){
                                 rep(i, times=ct-i+1)}, ct=len))
                      idx.next = unlist(lapply(init, FUN=function(i, ct){
                                 i:ct
                               }, ct=len))
                      idx.ma[,1]=idx.init
                      idx.ma[,2]=idx.next

             }else{
                      init=1:(len-1)
                      #print("util.it.smallLargeIdx2")
                      #print(init)
                      #print(len)
                             idx.ma = matrix(NA, ncol=2 , nrow=(len^2-len)/2 )
                             idx.init = unlist(lapply(init, FUN=function(i, ct){
                                        rep(i, times=ct-i)}, ct=len))
                             idx.next = unlist(lapply(init, FUN=function(i, ct){
                                        (i+1):ct
                                      }, ct=len))
                             idx.ma[,1]=idx.init
                             idx.ma[,2]=idx.next
                    }

          return(idx.ma)
    }

util.it.smallLargeIdxOOLD <-
function(idx, keep.same=FALSE){
    ct = length(idx)
    idxPair=NULL
    if(keep.same){
     for ( i in 1:ct ){
       for ( j in 1:ct){
         #if(i<j) IdxPair = rbind(hapIdxCorner, c(hap1=i, hap2=j))
         if(i<=j) idxPair = rbind(idxPair, c(idx1=i, idx2=j))
       }
     }
    return(idxPair)
   }else{
     for ( i in 1:ct ){
       for ( j in 1:ct){
         #if(i<j) IdxPair = rbind(hapIdxCorner, c(hap1=i, hap2=j))
         if(i<j) idxPair = rbind(idxPair, c(idx1=i, idx2=j))
       }
     }
     return(idxPair)
   }
}

util.it.triMatch <-
function(idxPair, len){
  ## rely on the structure of tri idx pair
  ## rely on the order to the pair
  if(idxPair[1]-idxPair[2]==0){
     ## two idx are the same
     endIdx = (len^2-len)/2
              rowIdx = endIdx + idxPair[1]
   }else{
              # two idx are not the same, old approach
              #idxLen = c(0, (len-1):1)
              #idxLen.upper = idxLen[1:idxPair[1]]
              #idx.added = idxPair[2]-idxPair[1]
              #rowIdx = sum(idxLen.upper)+idx.added

              # new math approach
              idx.added = idxPair[2]-idxPair[1]
              rowIdx = len/2*(len-1)-(1+len-idxPair[1])/2*(len-idxPair[1]) + idx.added
   }
   return(rowIdx)
}

util.it.triMatch2 <-
function(dipIdx, len, re.homo=FALSE){
       ## rely on the structure of tri idx pair
       ## rely on the order to the pair

       upper = .5*len*(len-1)
       if ((dipIdx)>upper){
         if(re.homo){
           return(c(rep(dipIdx-upper, 2), 0))
         }else{
           return(rep(dipIdx-upper, 2))
         }
       }

       cutoff = (len-1):1
       cutoff = cumsum(cutoff)

       block = qing.cut(dipIdx, cutPt=cutoff, cutPt.ordered = TRUE, right.include=TRUE)
       if(re.homo){
         return(c(block, dipIdx - c(0, cutoff)[block] + block, 1))
       }else{
         return(c(block, dipIdx - c(0, cutoff)[block] + block))
       }
#       if (block==1){
#         return( c(block, dipIdx+1))
#       }else{
#         return( c(block, dipIdx - cutoff[block-1] + block))
#       }

}

util.it.upTriCombIdx <-
function(idx, diag=TRUE, re.ordered = FALSE){
     ct = length(idx)
     reRow = (ct^2 - ct)/2

     if(diag){
       reRow = reRow + ct
       re = matrix(NA, ncol=2, nrow=reRow)
       if(re.ordered){
          it = 0
          for ( i in 1:ct ){
            for ( j in i:ct){
                it = it + 1
                re[it,] = range(c(idx[i], idx[j]))
            }
          }
       }else{
          it = 0
          for ( i in 1:ct ){
            for ( j in i:ct){
                it = it + 1
                re[it,] = c(idx[i], idx[j])
            }
          }
       }
       colnames(re) = c("id1", "id2")
       return(re)
     }

     re = matrix(NA, ncol=2, nrow=reRow)
     if(re.ordered){
        it = 0
        for ( i in 1:ct ){
          j = i + 1 
          while ( j <= ct){
              ## print(i)
              ## print(j)
              it = it + 1
              re[it,] = range(c(idx[i], idx[j]))
              j = j + 1
          }
        }
     }else{
        it = 0
        for ( i in 1:ct ){
          j = i + 1
          while( j <= ct){
              it = it + 1
              re[it,] = c(idx[i], idx[j])
              j = j + 1
          }
        }
     }
     colnames(re) = c("id1", "id2")
     return(re)

}

util.list.2matrix <-
function(list, byRow = TRUE, add.dimnames=FALSE) 
{
    colCt = length(list)
    rowCt = length(list[[1]])
    re = matrix(unlist(list), ncol = colCt, nrow = rowCt, byrow = FALSE)
    if(add.dimnames){
      header = dimnames( list[[1]] )[[1]]
      if(is.null(header)) header = paste("v", 1:colCt, sep="")
      rownames(re)=header
    }
    if (byRow) 
        re = t(re)
    return(re)
}

util.list.ex <-
function(nameList, varList, na.allow = TRUE, na.replace = NULL){
  ifD = FALSE

  varNames = names(varList)
  varLen = length(varList)
  exLen = length(nameList)

  if(varLen!=exLen) warning(paste("Two lists are of different length. nameList of length(", exLen, "); varList of length(", varLen,").", sep="") )
  ## create a map with key = var name, item = var position 
  varMap = data.frame(item = seq.int(from=1, to=varLen, by=1), key = varNames)

  namePos = rep(0, length=exLen)
  
  tryCatchEnv = new.env(parent=qingEnv)
  assign("namePos", namePos, envir=tryCatchEnv)
  
  ## based on this map, find the right sequence of positions to extract values
  for ( i in 1:exLen){
    tryCatch({namePos[i] = varMap$item[varMap$key == nameList[i]]},
             error = function(e){
               if(!na.allow){
                 stop(paste("var with name=[", nameList[i], "] not found in the varList.", sep=""))
               }else{
				 a = get("namePos", envir=tryCatchEnv)
				 a[i]=0
				 assign("namePos", a, envir=tryCatchEnv)
                 #namePos[i] <<- 0
               }
             })

  }

  ## assign default value
  varEx = rep(NA, length=exLen)
  if(!is.null(na.replace)){
    varEx = rep(na.replace, length=exLen)
  }
  
  for ( i in 1:exLen ){
    if(namePos[i]!=0){
      varEx[i]= unlist(varList[namePos[i]])
      if(ifD) print(paste("i=(", i, ") namePos=(", namePos[i], ")."))
    }
  }
    
  return (varEx)
}

util.list.rmByKeyVal <-
function(dataList, keys, keyRemoved){
    seqOrder = is.element(keys, keyRemoved)
    return(dataList[!seqOrder])
}

util.listMatrix.2matrix <-
function(list, rbind=TRUE){

  rowCt = length(list)

  re = NULL
  if(rbind){
      for ( i in 1:rowCt ){
         re = rbind(re, list[[i]])
      }
  }else{
      for ( i in 1:rowCt ){
        re = cbind(re, list[[i]])
      }
  }

  return(re)
}

util.matrix.2list <-
function(ma, byRow=TRUE){

  if(byRow){
    colNum = dim(ma)[2]
    rowNum = dim(ma)[1]
    ma = t(ma)
  }else{
    colNum = dim(ma)[1]
    rowNum = dim(ma)[2]
  }

  maItem = unlist(as.list(ma))
  re = NULL
  for(row in 1:rowNum){
    cur = maItem[seq.int(from=(row-1)*colNum+1, to=row*colNum, by=1)]
    re = c(re, list(cur))
  }
  return(re)
}

util.matrix.cat <-
function(data, cols, sep=""){
  len = length(cols)
  for(i in 1:len){
    if(i==1) {
      re = data[,cols[i]]
    }else{
      re = paste(re, data[,cols[i]], sep=sep)
    }
  }
  return(re)
}

util.matrix.catm <-
function(inMatrix, colVec, discIn=NULL, discVec=NULL, delimVec, digitVec, missingVec=NULL){

     anyMatrix = inMatrix[,colVec]
     mRow = dim(anyMatrix)[1]
     mCol = dim(anyMatrix)[2]
     
     re = matrix(ncol =2, nrow = mRow)

     if(!is.null(discIn)){
       re[,1]= inMatrix[,discIn]
     }else{
       re[,1]=discVec
     }
     for(row in 1:mRow) {
       curRow = NULL
       for(col in 1:mCol){
         num = anyMatrix[row, col]
         if(is.na(num)){
           if(is.null(missingVec)){
             num = "NA"
           }else{
             num = missingVec[col]
           } ## if(is.null(missingVec)){           
         }else{
           ##num = round(num, digitVec[col])
           if(digitVec[col]==0){
             num=format(num, nsmall=0)
           }else{
             num = format(num, digits=digitVec[col], nsmall=digitVec[col])
           }
         } ##if(is.na(num)){
         
         if(col == 1){
           curRow = paste(delimVec[(0+col)], num, sep="")
         }else{
           if(col == mCol){
             curRow = paste(curRow, delimVec[col], num, delimVec[1+col], sep="")
           }else{
             curRow = paste(curRow, delimVec[col], num, sep="")
           } ##if(col = mCol){
         } ## if(col = 1){
         
       } ##for(col in 1:mCol){
       re[row,2] = curRow
     } ##for(row in 1:mRow) {
     return(re)
}

util.matrix.catSave <-
function(data, cols, sep ="-"){
  colNum = dim(data)[2]
  keys = util.matrix.cat(data, cols, sep)
  data[,colNum+1]=keys
  return(data)
}

util.matrix.clone <-
function(ma, n, rowAppend=TRUE){
  rnm = dim(ma)[1]
  cnm = dim(ma)[2]
  if(rowAppend){
    t1 = replicate(n, t(ma))
    re = t(matrix(t1, nrow=cnm, byrow=FALSE))
    re
  }else{
    t1 = replicate(n, ma)
    re = matrix(t1, nrow=rnm, byrow=FALSE)
    re
  }
  return(re)
}

util.matrix.col.shuffle <-
function(ma){

  if(!is.null(dim(ma)))  return(ma)

  colNum = dim(ma)[2]
 
  filterSeq = matrix(1:colNum, ncol=2)
  filterSeq = as.vector(t(filterSeq))

  re = ma[,filterSeq]

  return(re)
}

util.matrix.col.shuffle2 <-
function(ma1, ma2){

  re = cbind(ma1, ma2)

  colNum = 1
  if(!is.null(dim(ma1))){
    colNum = dim(ma1)[2]
  }
  filterSeq = matrix(1:(2*colNum), ncol=2)
  filterSeq = as.vector(t(filterSeq))

  re = re[,filterSeq]

  return(re)
}

util.matrix.colComp <-
function(ma, values, operator="and"){
  ifD = FALSE
  vCt = length(values)

  ## if "and" operator is chosen, the original list will decrease to speed up the comparison
  matchId = 1:(dim(ma)[1])
  if(ifD) print(matchId)
  for (i in 1:vCt ){
    if(ifD) print(i)
    compIdx = which(ma[matchId,i]==values[i])
    matchId = matchId[compIdx]
    if(ifD) print(paste("compIdx=", compIdx, collapse=", ", sep=""))
    if(ifD) print(paste("matchId=", matchId, collapse=", ", sep=""))
    
    if(length(matchId)==0) return (NULL)
  }
  return(matchId)

  ## if "or" operator is chose, the original list will remain the same to capture all that meet it.
  ## not implemented
}

util.matrix.colIdx4Match <-
function(ma, val){

  rowCt = nrow(ma)

  #print("util.matrix.colIdx4Match")
  #print(dim(ma))

  matchIdx = qing.mulMatch(val, ma)
  ##matchIdx = which (ma == val)
  
  if(length(matchIdx)>0 & matchIdx[1]!=0){
    
    re = qing.cut(matchIdx,
                  cutPt = seq.int(from=rowCt, to=rowCt*ncol(ma), by=rowCt), 
                  cutPt.ordered = TRUE, right.include=TRUE)

    ##  remove because it cause memeory problems
    ##  as.numeric(cut(matchIdx, breaks = c(0, seq(rowCt, rowCt*ncol(ma), by=rowCt)), include.lowest=TRUE, right=TRUE))
    
  }else{
    
    re =  NULL
  }

  return(re)
}

util.matrix.csvText <-
function(numericMa, fileName=NULL, colNames=NULL, rowNames=NULL, digit = NULL,  keepZero=FALSE){
     mRow = dim(numericMa)[1]
     mCol = dim(numericMa)[2]
     i = 0

     lineComment = ""
     
     comm = vector( )

            if(!is.null(colNames)){
              if(!is.null(rowNames))     lineComment = " \t"
              for( col in 1:mCol ){
                lineComment = paste(lineComment, colNames[col], sep="\t")  
              }
              i = i+1
              comm[i] = lineComment
            }     

     lineComment = ""
     
     for( row in 1:mRow ){
       i = i+1
          for( col in 1:mCol ){
          
                if(col == 1){
     
                  if(!is.null(rowNames)){
                        lineComment = paste(lineComment, rowNames[row], sep="\t")  
                      }
                  }
                  tmp = numericMa[row,col]
                  if(is.numeric(tmp) & (!is.null(digit))){
                    if(!keepZero){
                      if(digit==0){
                          num = format(tmp, nsmall=0)
                      }else{
                          num = format(tmp, digits=digit, nsmall=digit)
                      }
                      lineComment = paste(lineComment, num, sep="\t")
                    }else{
                      lineComment = paste(lineComment, round(tmp, digit), sep="\t")
                    }
                  }else{
                    lineComment = paste(lineComment, tmp, sep="\t") 
                  }
     
            }  ## for( col in 1:mCol ){
          comm[i] = lineComment
          lineComment = ""
     }  ## for( row in 1:mRow ){

     if(!is.null(fileName)) write(comm, file = fileName,  append = TRUE)
     return(comm)
}

util.matrix.delCol <-
function(ma, colPos){
	rnm = dim(ma)[1]
	cnm = dim(ma)[2]
	cutPt = (colPos-1)*rnm
	if(colPos==cnm){
	    re = matrix(ma[1:cutPt], nrow = rnm, ncol=(cnm-1))
	    return(re)
	}
        if(colPos == 1){
            re = matrix(ma[(rnm+1):(rnm*cnm)], nrow = rnm, ncol=(cnm-1))
	    return(re)
        }else{	
	    re = c(ma[1:cutPt], ma[(cutPt+rnm+1):(rnm*cnm)])
	    re = matrix(re, nrow = rnm, ncol=(cnm-1))									
	    return(re)
	}																			
}

util.matrix.exByKeyInRow <-
function(ma, keyCol, keys){

    maNew = ma[ is.element( ma[,keyCol], keys), ]

    return(maNew)
}

util.matrix.insertCol <-
function(ma, insertVec, afterCol){
	rnm = dim(ma)[1]
	cnm = dim(ma)[2]
	startnm = afterCol*rnm+1
        re = NULL
	if(afterCol==cnm){
	    re = matrix(c(ma, insertVec), ncol=cnm+1)
	}else{	
	    re = c(ma[1:startnm-1], insertVec, ma[startnm:(rnm*cnm)])
	    re = matrix(re, nrow = rnm, ncol=cnm+1)									
	}
        return(re)
      }

util.matrix.merge <-
function(ma, ma2){
  maVec = as.vector(ma)
  maVec2 = as.vector(ma2)
  re = c(maVec, maVec2)
  renames = NULL
  if(is.matrix(ma)) {
    maLen = dim(ma)[2]
    renames = colnames(ma)
  }else{
    maLen = 1
    if(is.null(names(ma))){
      renames = ""
    }else{
      renames = names(ma)
    }
  }
  

  if(is.matrix(ma2)) {
    maLen2 = dim(ma2)[2]
    renames = c(renames, colnames(ma2))

  }else{
    maLen2 = 1
    if(is.null(names(ma2))){
      renames = c(renames, "")
    }else{
      renames= c(renames, names(ma2))
    }
  }    
  if(F)print(renames)
  colNum = maLen + maLen2
  if(F)print(colNum)
  re = matrix(re, ncol=colNum)
  colnames(re) = renames
  
  return(re)
}

util.matrix.rmSparseRow <-
function(ma, colChecked, minNumInCol=1){

  ifD = FALSE
  checked = ma[,colChecked]
  if(ifD) print(checked)

  len = length(checked)
  rmList = NULL
  reList = NULL
  for( i in 1:len ){
    if(ifD) print(checked[i])
    if( checked[i] <= minNumInCol ){
      rmList = c(rmList, i)
    }else{
      reList = c(reList, i)

    }
  }
  if(ifD) print(rmList)
  if(is.null(rmList)) return(NULL)
  len = length(rmList)
  re = ma
    for( j in 1:len ){
      re = t(util.matrix.delCol (t(re), rmList[j]))
      if(ifD) print(re)
      rmList = rmList - 1
      if(ifD) print(rmList)
    }

  return(list(re, reList))
}

util.sql.groupby <-
function(data, groupCols, sep="-", varCol, type=c("max", "min", "sum", "mean"), na.rm=FALSE){

  m = match(c("max", "min", "sum", "mean"), type, 0)
  matchFun = which(m==1)
  
  
  if(length(groupCols)>1){
    key = util.matrix.cat(data, cols=groupCols, sep=sep)
  }else{
    key = data[,groupCols]
  }

  keyVal = unique(key)
  grpMax=NULL
  grpMin=NULL
  grpSum=NULL
  grpMean= NULL
  varVec=NULL
  for( i in keyVal){
    filter = key==i
    varVec = unlist(data[filter, varCol])
    for(j in matchFun){
      if(j==1) grpMax = c(grpMax, max(varVec, na.rm=na.rm))
      if(j==2) grpMin = c(grpMin, min(varVec, na.rm=na.rm))
      if(j==3) grpSum = c(grpSum, sum(varVec, na.rm=na.rm))
      if(j==4) grpMean = c(grpMean, mean(varVec, na.rm=na.rm))
    }

  }
  
  re = NULL
  re = cbind(re, key=keyVal) 
  for (j in matchFun){
    if(j==1) re = cbind(re, max=grpMax)
    if(j==2) re = cbind(re, min=grpMin)
    if(j==3) re = cbind(re, sum=grpSum)
    if(j==4) re = cbind(re, mean=grpMean)
  }

  return(re)
}

util.str.2CharArray <-
function(str, len=nchar(str)){

  startSeq = 1:len
  re = lapply(startSeq, FUN = function(inStr, startSeq){
    char = substr(inStr, startSeq, startSeq)
    char
  }, inStr = str)
  re = unlist(re)
  return(re)
}

util.str.replace <-
function(str, replaced, new, replace.all=FALSE){

  pos = util.char.1stStrIdx(str, find=replaced)
  if(pos==0){
    # if not found, return the original
    return(str)
  }

  if(  (pos+nchar(replaced))<=nchar(str)){
    str.left = substr(str, pos+nchar(replaced), nchar(str))
  
    if(replace.all){
      # iteratively replace the left-over part
      str.left = util.str.replace(str.left, replaced=replaced, new=new, replace.all=replace.all)
    }
    str.re = paste(substr(str, 1, pos-1), new, str.left, sep="")
  }else{
    str.re = paste(substr(str, 1, pos-1), new, sep="")
  }
  return(str.re)

}

util.str.rmSpacePadder <-
function(str, rm=c("b", "e")){

    charArray = util.str.2CharArray(str, len=nchar(str))
    blankPos = which(charArray!=" ")
    if(length(blankPos)==0) return("")
    if(is.element("b", rm)){
       if(is.element("e", rm)){
         re = substr(str, blankPos[1], blankPos[length(blankPos)])
       }else{
         re = substr(str, blankPos[1], nchar(str))
       }
    }else{
       if(is.element("e", rm)){
         re = substr(str, 1, blankPos[length(blankPos)])
       }else{
         re = str
       }
    }
    return(re)
}

util.str.seqCutter <-
function(str, delims){
   i = 1
   j = 1
   stopS = FALSE
   re = NULL
   len = length(delims)
   curStr = str
   lastIdx = 0
   while (i<=len & !stopS){
     idx.f = util.char.1stStrIdx(curStr, delims[i])
     if (idx.f==0){
       stopS = TRUE
       return(NA)
     }else{
       # find the current one
       #print(curStr)
       #print(idx.f)
       seg = substr(curStr, 1, idx.f-1)
       if(i==len) lastIdx = nchar(curStr)
       curStr = substr(curStr, idx.f+nchar(delims[i]), nchar(curStr))
       #print(curStr)
       re = c(re, seg)
     }
     i = i +1 
   }
   # if(idx.f+nchar(delims[i])<lastIdx) re = c(re, curStr) 
   if(idx.f + 2 < lastIdx) re = c(re, curStr)
   return(re)
 }

util.str.splitAndOrNot <-
function(str, split="and", check.space=FALSE){
  split.ex = util.str.splitEx(str, split)
  ## for example "and"
  if(length(split.ex)==1){
    ## check whether no "and" is found
    if (nchar(split.ex[1])==nchar(str)) {
      re = NULL
      return(re)
    }
  }
  ## remove space
  split.ex.rm= unlist(lapply(split.ex, util.str.rmSpacePadder, rm=c("b", "e")))
  items=util.array.rmEmptyStr(split.ex.rm)

  ## if check.space, return NULL if any of the item contains space in it
  if(check.space){
    found=FALSE
    i=1
    while(i<=length(items) & !found){
      tmp.idx = util.char.1stIdx(items[i], " ")
      if(tmp.idx>0) found=TRUE
      i = i + 1
    }
    if(found) return(NULL)
  }
  ## create the list for return
  re = list(op=split, items=util.array.rmEmptyStr(split.ex.rm))
  return(re)
 
}

util.str.splitEx <-
function(str, split, case.sen=TRUE){

  if(!case.sen){
    str = tolower(str)
    split = tolower(split)
  }
  split.re = strsplit(str, split=split)
  split.re=split.re[[1]]

  len.s = nchar(split.re)
  re = split.re[len.s!=0]
  return(re)
}

util.str.tokenPicker <-
function(str, delim, indexPicked){
  tokens = unlist(strsplit(str, delim))
  return(tokens[indexPicked])
}

util.tbl.exFactorInfo <-
function(tbl){
  rnum = dim(tbl)[1]
  cnum = dim(tbl)[2]

  headList = dimnames(tbl)[[2]]
  numList  = dimnames(tbl)[[1]]
  
  num = NULL
  head = NULL
  ct = NULL
  for ( i in 1:rnum){
    for (j in 1:cnum){
      cur = tbl[i,j]

      num = c(num, numList[i])
      head = c(head, headList[j])
      ct = c(ct, cur)
        
    }
  }
  re = data.frame( row.level=num, col.level=head, ct=as.numeric(ct))
  return(re)
}

util.vec.2maIdx <-
function(rowCt, idx){
  ## assuming putted in  column by column
  rowIdx = idx%%rowCt
  colIdx = (idx-rowIdx)/rowCt+1
  if(rowIdx==0) 
    rowIdx=rowCt
    
  return(c(rowIdx, colIdx))
}

util.vec.exByKey <-
function(vecKey, vec, keysOrder){
    filter =  match(keysOrder, vecKey)

    ex = vec[filter]

    return(ex)

}

util.vec.matchVec <-
function(vec, benchMark, vecLen, benchLen){

  dup = vecLen/benchLen

  bench = rep(benchMark, dup)

  matchset = vec == bench

  re = NULL
  for( i in 0:(dup-1)){
    pos = i*benchLen + 1
    re = c(re, as.logical(min(matchset[pos:(pos+benchLen-1)])))
  }

  return(re)
}

util.vec.matchVecIdx <-
function(vec, benchMark, vecLen, benchLen){

  aa = util.vec.matchVec(vec, benchMark, vecLen, benchLen)
 
  return(which(aa))
}

util.vec.orderPair <-
function(pair){
  front = min(pair)
  back = max(pair)
  return(c(front, back))
}

util.vec.replace <-
function(vec, orignal, replaceBy){
     vec.idx = match(vec, orignal)
     re = replaceBy[vec.idx]
     return(re)
}

wrComTbl <-
function(numericMa,       fileName=NULL, colNames=NULL, rowNames=NULL, digit = NULL){
     mRow = dim(numericMa)[1]
     mCol = dim(numericMa)[2]
     i = 0

     lineComment = ""
     
     comm = vector( )

            if(!is.null(colNames)){
              if(!is.null(rowNames))     lineComment = " \t"
              for( col in 1:mCol ){
                lineComment = paste(lineComment, colNames[col], sep="\t")  
              }
              i = i+1
              comm[i] = lineComment
            }     

     lineComment = ""
     
     for( row in 1:mRow ){
       i = i+1
          for( col in 1:mCol ){
          
                if(col == 1){
     
                  if(!is.null(rowNames)){
                        lineComment = paste(lineComment, rowNames[row], sep="\t")  
                      }
                  }
                  tmp = numericMa[row,col]
                  if(is.numeric(tmp) & (!is.null(digit))){
                    lineComment = paste(lineComment, round(tmp, digit), sep="\t") 
                  }else{
                    lineComment = paste(lineComment, tmp, sep="\t") 
                  }
     
            }  ## for( col in 1:mCol ){
          comm[i] = lineComment
          lineComment = ""
     }  ## for( row in 1:mRow ){

     if(!is.null(fileName)) write(comm, file = fileName,  append = FALSE)
     return(comm)
}

