trio.sim <-
function(freq,  interaction="1R and 2D", prev=1e-3, OR=1, pi.usr=0, n=100, rep=1,  step.save=NULL, step.load=NULL,  verbose=FALSE){
  sigType="D/R"

  if (prev >= 1 | prev <= 0) stop("Argument, prev, is out of range of (0, 1).")
  if (pi.usr >= 1 | pi.usr < 0) stop("Argument, pi.usr, is out of range of (0, 1).")

  ## NEW ADDITION to replace old code: Added By Qing 28Jan2013
  ## para = c( log(prev/(1-prev)), log(OR))
  alpha = log(prev/(1-prev))
  beta = log(OR)
  beta.new = log(exp(alpha+beta)+pi.usr) - log (exp(alpha)+pi.usr)
  alpha.new = log(exp(alpha)+pi.usr) - log(1-pi.usr)  
  para = c(alpha.new, beta.new)
  ## END of replacement
  
  dig1Code=c(3,2,1,0)
  ddF = NULL
  spSuphap.namePrefix=NULL
  alleleCode=1:2

  # print(paste("Try to simulate trio data: # datasets=", rep, "; # trio =", n, sep=""))
  infoS.s = NULL
  
  bkMap = bkMap.constr(data=freq, keyCol=1, hapLenCol=NULL, expCol=2, probCol=3, alleleCode=alleleCode)

  ## whether it contain hap and geno only (estimated by imputation)


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
                    spStrata.saveFN = step.load,
                    spStrata.name = step.save,
                    reControl =FALSE,  dig1Code=dig1Code 
                    )
  
  if(verbose) print(paste("Time used to generate ", rep,  " dataset(s).", sep=""))
  if(verbose) print(proc.time() - ptm)
  gc.e = gc()
 
  if(verbose) print("Info on memory usage:")
  if(verbose) print(gc.e)

  #print(paste("Finished generating ", rep,  " dataset(s).", sep=""))

  return(trioData)
}



