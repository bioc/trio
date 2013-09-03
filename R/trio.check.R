trio.check <-
function(dat, is.linkage = TRUE, replace=FALSE){
  re=NULL
  snpIdxRange=NULL
  minor = 2
  if (!is.element(minor, c(1,2))) stop("The code value for minor allele must to be either 1 or 2.")
  if(is.linkage){
    if(minor==2){
      re = linkageFile.proc(data=dat, snpIdxRange=snpIdxRange,
                          action = c("formTrio1digit", "Mendelian check"),
                         # dig2Code=0:2, ## only need to estimate the frequencies
                          dig1Code=c(NA,0,1,2))
        
    }else{
      ## in SZ trio
      re = linkageFile.proc(data=dat, snpIdxRange=snpIdxRange,
                          action = c("formTrio1digit", "Mendelian check"),
                         # dig2Code=0:2, ## only need to estimate the frequencies
                          dig1Code=c(NA,2,1,0))
    }

    trio1digit = re$trio

    colnames(trio1digit) = c( colnames(dat)[1:2], paste("snp", 1:((ncol(dat)-6)/2), sep=""))
    
  }else{
    re = trioFile.proc(data=dat,  dig1Code=c(NA,0,1,2), action = c("Mendelian check"))
    trio1digit = dat
  }

  # in case of Mendelian error, replace those with NA by request
  if(replace){
    if(!is.null(re$MedErr)){
      # err
      repdd = setTrioMissingSNP(trioDf=trio1digit, cord=re$MedErr,
                 snp1digit=TRUE, missingDigit = NA)
      reo = c(trio=list(repdd), errors=list(NULL))
    }else{
      # no err
      reo = c(trio=list(trio1digit), errors=list(NULL))
    }
    
  }else{
    if(!is.null(re$MedErr)){
      # err 

      print("Found Mendelian error(s).")

      # reformat the error
      #trio, FamID, SNP, r and c
      newMedErr = data.frame(re$MedErr[,3, drop=FALSE],
                        trio1digit[re$MedErr[,1],1, drop=FALSE],
                        re$MedErr[,c(4,1,2), drop=FALSE])
      newMedErr = newMedErr[order(newMedErr[,1], newMedErr[,2], newMedErr[,3]),]
      colnames(newMedErr) = c("trio", "famid", "snp", "r", "c")
      rownames(newMedErr) = NULL
      reo = c(trio=list(NULL),
        errors=list(newMedErr), trio.err = list(trio1digit))
      
    }else{
      # no err
      reo = c(trio=list(trio1digit), errors=list(NULL))
    }
  }
  return (reo)
}

