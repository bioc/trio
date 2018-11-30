trio.prepare <-
function(trio.dat, freq=NULL, blocks=NULL, logic=TRUE, ...){
   # requireNamespace("haplo.stats")
   ifD = FALSE

   alleleCode = 1:2
   
   if (!is.null(trio.dat$errors)) {     
     stop(paste("The input value for trio.dat does not contain valid trio data."))
   }
   if (is.null(trio.dat$trio)) {     
     stop(paste("The input value for trio.dat does not contain valid trio data."))
   }
   trio.datIn = trio.dat$trio
   
   trio.re = NULL
   if(!is.null(freq)){
     #print("Use frequencies provided by users.")

     freq.fromhaponly = freqbuild.haponly(hap=freq, alleleCode=alleleCode)
     
     if(logic){
       imputed = trio.impu(triodd=trio.datIn, freq=freq.fromhaponly, impu.missingOnly=FALSE)
     }else{
       imputed = trio.impu(triodd=trio.datIn, freq=freq.fromhaponly, impu.missingOnly=TRUE)
       imputed = cbind(trio.datIn[,1:2], imputed)
     }
   }else{
     if( !is.null(blocks)){
       #print("Use haplo.em() to estimate frequencies based on block sizes provided by users." )
       if( (ncol(trio.datIn)-2) != sum(blocks)){
         stop(paste("There are ", ncol(trio.datIn)-2, " SNPs in the trio data, but a total of ",
         sum(blocks), " SNPs based on the argument, blocks.", sep=""))
       }

       if(sum(blocks>=8)>0)
         stop("At least one of the block size exceeds the maximum value of 7 loci.") 
       
     }else{
#       print(paste("Use haplo.em() to estimate frequencies. ",
#                   "Because block sizes are not provided by users, will assume all SNPs are single SNPs.",
#                   sep=""))
       blocks = rep(1, times=ncol(trio.datIn)-2)
     }
     trio.re = trioFile.proc(data=trio.datIn,
                             key.prefix="ch",
                             bk.sizes=blocks, dig1Code=c(NA, 0, 1, 2),
                             dig2Code=c(0,1,2), action = c("freq estimate"), ...)


	 
     bbbb=simuHapMap.build(hapInfoFrame=trio.re$freq$hapMap.info,
                 genoInfoFrame=trio.re$freq$genoMap.info)
     
     if(is.null(trio.re$freq$hapMap.info)){
		#snp only
	 	trio.hapFreqOnly = bbbb$newDf
 	}else{
	 	trio.hapFreqOnly = bbbb$newDf[,c(1,4,5)]
 	}
     colnames(trio.hapFreqOnly)=c("key", "hap","freq")		 

 
     if(logic){     
       imputed = trio.impu(triodd=trio.datIn, freq=trio.re$freq,  impu.missingOnly=FALSE)
     }else{
       imputed = trio.impu(triodd=trio.datIn, freq=trio.re$freq,  impu.missingOnly=TRUE)
       imputed = cbind(trio.datIn[,1:2], imputed)
     }
   }

   if(logic){
      fSeq = seq(1, nrow(imputed), by=6)
      mSeq = seq(2, nrow(imputed), by=6)
      cSeq = seq(3, nrow(imputed), by=6)
      pSeq = c(fSeq, mSeq)
      
      child = imputed[-pSeq,]
      to2digit = exchangeDigit(child, cols = NULL,
                     dig1Code = c(NA, 0, 1, 2),
                     dig2Code = c(0, 2, 1), action = c("1to2"))
      bina = 2-to2digit
   
      if(!is.null(colnames(trio.datIn))){
        ## check unique
        tmpColName = paste( rep( colnames(trio.datIn)[-c(1,2)], each=2), c(".D", ".R"), sep="")
        if( length(unique(tmpColName)) != length(tmpColName)) {
          #print("Column names based on original name is not unique. Will use sequential indexes of SNPs as column names in the returned binary data.")
          colnames(bina) = paste( rep(1:ncol(child), each=2), c(".D", ".R"), sep="")
        }else{
          colnames(bina) = tmpColName
        }
      }else{
        colnames(bina) = paste( rep(1:ncol(child), each=2), c(".D", ".R"), sep="")
      }
      if(ifD){
        print("check code: trio.datIn, imputed")
        print(cbind(trio.datIn[1:6, 1:7], NA, imputed[c(1:3, 7:9), 1:5]))
        print("imputed")
        print(imputed[1:12, 1:5])
        for( i in 3: min(10,ncol(trio.datIn)-2)){
          print(paste("snp", i-2, ": imputed, NA, bina"))
          print(cbind(child[1:4, i-2], NA, bina[1:4, c((i-2-1)*2+1, (i-2)*2)]))
        }
      }

      y = rep(0, nrow(bina))
      y[ seq(1, nrow(bina), by=4)]=3
      bin = cbind(y=y, bina)
      re.fin = list(bin=bin)
      
    }else{ ## if(logic){

      if(!is.null(colnames(trio.datIn))){
        colnames(imputed) = colnames(trio.datIn)
      }else{
        colnames(imputed) = c(colnames(trio.datIn)[1:2], paste("snp", 1:(ncol(trio.datIn)-2), sep=""))
      }
      
      re.fin = list(trio=imputed)
    }


#   if(is.element("miss", op.re)){
     trio.missIdx = trioFile.proc(data=trio.datIn, dig1Code=c(NA, 0, 1, 2),
                 action = c("missingReport"))$missIdx
     # before it is row, snp
     # famid, pid, SNP number, and then row (r) and column (c) index (5 elements).
     if(!is.null(trio.missIdx)){
       if(ifD){
         print("here")
         print(dim(trio.missIdx))
         print(dim(trio.datIn))
         #print(trio.missIdx[1:10, ])
         #print(trio.datIn[1:10, 1:5])
         print(trio.missIdx)
  
         #return( trio = list(trio.datIn), miss=list(trio.missIdx))
         print(dim(   trio.datIn[trio.missIdx[1:300,1], 1:2, drop=FALSE]  ))
         print(dim(trio.missIdx[,2:1, drop=FALSE]))
         print(length(trio.missIdx[,2]+2))
       }
       trio.missNewInfo = data.frame(
         trio.datIn[trio.missIdx[,1],1:2, drop=FALSE],
         trio.missIdx[,2:1, drop=FALSE],
         trio.missIdx[,2]+2)
       trio.missNewInfo = trio.missNewInfo[order(trio.missNewInfo[,1], trio.missNewInfo[,2], trio.missNewInfo[,3]),]
       colnames(trio.missNewInfo) = c("famid", "pid", "snp", "r", "c")
       rownames(trio.missNewInfo) = NULL
       re.fin= c(re.fin, miss=list(trio.missNewInfo))
     }else{

       re.fin= c(re.fin, miss=list(NULL))

     }

#   }

#   if(is.element("freq", op.re)){
     if(!is.null(freq)) {
       #print("No estimated frequencies are returned, because users have provided it.")
       re.fin= c(re.fin, freq=list(freq))
     }else{
       re.fin= c(re.fin, freq=list(trio.hapFreqOnly))
     }
   class(re.fin) <- "trioPrepare"
   re.fin
 }

