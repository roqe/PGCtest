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

GCN=function(HA,mc=5,prec=1e-3,GC=NULL,bSNP=2){
  Y=HA$Y
  X=HA$X
  nM=names(HA$M)[!sapply(HA$M,is.null)]
  GS=parallel::mclapply(nM,function(nm){
    M=HA$M[[nm]]
    PV=lapply(names(Y),function(ny){
      if(!is.null(GC)){
        cad=which(!GC[[ny]]$gene%in%nm)[1:(nrow(GC[[1]])-1)]
        gs=cbind(GC[[ny]]$gene[cad],GC[[ny]]$leadSNP[cad])
        XX=do.call(cbind,apply(gs,1,function(gs){
          slist=strsplit(gs[2],";")[[1]]
          return(HA$M[[gs[1]]][slist])
        }))
        X=cbind(X,XX)
      }
      y=Y[[ny]]
      nad=which(!is.na(y))
      RG=pre_fit_one(M=M[nad,],Y=y[nad],X=X[nad,])
      PG_TSQ=TSQ(RG$b,RG$corrB)
      PG_GBJ=GBJ::GBJ(RG$b,RG$corrB)
      PG_GHC=GBJ::GHC(RG$b,RG$corrB)
      PG_mnP=GBJ::minP(RG$b,RG$corrB)
      PG_ACAT=ACAT::ACAT(2*pnorm(abs(RG$b),lower.tail = F))
      return(data.table(TSQ=safe_z(PG_TSQ$pv)*RG$bs,
                  GBJ=safe_z(PG_GBJ$GBJ_pvalue)*RG$bs,
                  GHC=safe_z(PG_GHC$GHC_pvalue)*RG$bs,
                  mnP=safe_z(PG_mnP$minP_pvalue)*RG$bs,
                  ACAT=safe_z(PG_ACAT)*RG$bs,
                  leadSNP=paste(unique(c(names(RG$b)[which.max(abs(RG$b))],names(RG$b)[which(abs(RG$b)>bSNP)])),collapse = ";")))
    })
    names(PV)=names(Y)
    return(data.table::rbindlist(PV,idcol="trait"))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  names(GS)=nM
  return(data.table::rbindlist(GS,idcol="gene"))
}

