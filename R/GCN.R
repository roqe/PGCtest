#' @param HA Dataset
#' @param mc Number of cores for parallel computing, default=5.
#' @param prec The unit to construct empirical normal product pdf using in composite test, default=1e-3.
#' @import parallel
#' @import SKAT
#' @import GBJ
#' @import ACAT
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' stat=GCN(HA)

GCN=function(HA,mc=5,prec=1e-3){
  cpl=sapply(HA$M,is.null)
  Y=HA$Y
  X=HA$X
  GS=parallel::mclapply(HA$M[which(!cpl)],function(M){
    PV=apply(Y,2,function(y){
      nad=which(!is.na(y))
      PG_VCT=VCT(M=M[nad,],Y=y[nad],X=X[nad,])
      RG=pre_fit_one(M=M[nad,],Y=y[nad],X=as.data.frame(X[nad,]))
      PG_TSQ=TSQ(RG$b,RG$corrB)
      PG_GBJ=GBJ::GBJ(RG$b,RG$corrB)
      PG_GHC=GBJ::GHC(RG$b,RG$corrB)
      PG_mnP=GBJ::minP(RG$b,RG$corrB)
      PG_ACAT=ACAT::ACAT(2*pnorm(abs(RG$b),lower.tail = F))
      return(list(VCT=safe_z(PG_VCT$pv)*RG$bs,
                  TSQ=safe_z(PG_TSQ$pv)*RG$bs,
                  GBJ=safe_z(PG_GBJ$GBJ_pvalue)*RG$bs,
                  GHC=safe_z(PG_GHC$GHC_pvalue)*RG$bs,
                  mnP=safe_z(PG_mnP$minP_pvalue)*RG$bs,
                  ACAT=safe_z(PG_ACAT)*RG$bs))
    })
    return(data.table::rbindlist(PV,idcol="phenotype"))
  },mc.cores = mc,mc.preschedule = F,mc.cleanup = T)
  return(data.table::rbindlist(GS,idcol="gene"))
}

