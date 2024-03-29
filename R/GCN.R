#' For gene-trait associations.
#'
#' @param HA Dataset, output of sim_pleiotropy_data.
#' @param mc Number of cores for parallel computing, default=5.
#' @param GC List for adjusting co-regulating variants, the output from select_GC.
#' @param bSNP Threshold for selecting significant variants to be adjusted, default=2.
#' @param apply_TSQ Use Hotelling's T-squared statistic for gene-trait association, default=False.
#' @param apply_GBJ Use Generalized Berk-Jones for gene-trait association, default=False.
#' @param apply_GHC Use Generalized Higher Criticism for gene-trait association, default=False.
#' @param apply_mnP Use minimum p value for gene-trait association, default=False.
#' @param apply_ACAT Use ACAT for gene-trait association, default=True.
#' @import parallel
#' @import SKAT
#' @import GBJ
#' @import ACAT
#' @import data.table
#' @export
#' @examples
#' H00=sim_pleiotropy_data(hypo = "H0")
#' preR=GCN(H00,mc=40,apply_TSQ = T,apply_GBJ = T,apply_GHC = T,apply_mnP = T)
#' FR00=PGC(preR)
#' HAA=sim_pleiotropy_data(hypo = "HA",mm = 1,vv = 1,gpr = 1,rho = 0.3)
#  preR=GCN(HAA,mc=40,apply_TSQ = T,apply_GBJ = T,apply_GHC = T,apply_mnP = T)
#. GC=select_GC(preR,10)
#. aftR=GCN(HAA,mc=40,apply_TSQ = T,apply_GBJ = T,apply_GHC = T,apply_mnP = T,GC=GC)
#  FRAA=PGC(aftR)

GCN=function(HA,mc=5,GC=NULL,bSNP=2,apply_TSQ=F,apply_GBJ=F,apply_GHC=F,apply_mnP=F,apply_ACAT=T,single=F,fill=F){
  if(sum(apply_TSQ,apply_GBJ,apply_GHC,apply_mnP,apply_ACAT)==0){
    print("Please choose at least one statistic for gene-trait association.")
    return(NULL)
  }
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
          slist=paste0(gs[1],"_",strsplit(gs[2],";")[[1]])
          return(HA$M[[gs[1]]][slist])
        }))
        X=data.frame(cbind(X,XX))
      }
      y=Y[[ny]]
      nad=which(!is.na(y))
      RG=pre_fit_one(M=M[nad,],Y=y[nad],X=X[nad,])
      if(single){
        return(data.table(SNP=names(RG$b),beta_hat=RG$bh,beta_sd=RG$be))
      }else{
        if(apply_TSQ){ PG_TSQ=TSQ(RG$b,RG$corrB) }
        if(apply_GBJ){ PG_GBJ=GBJ::GBJ(RG$b,RG$corrB) }
        if(apply_GHC){ PG_GHC=GBJ::GHC(RG$b,RG$corrB) }
        if(apply_mnP){ PG_mnP=GBJ::minP(RG$b,RG$corrB) }
        if(apply_ACAT){ PG_ACAT=ACAT::ACAT(2*pnorm(abs(RG$b),lower.tail = F)) }
      }
      names(RG$b)=gsub(paste0(nm,"_"),"",names(RG$b))
      return(cbind(data.table(TSQ=ifelse(apply_TSQ,safe_z(PG_TSQ$pv)*RG$bs,NA),
                  GBJ=ifelse(apply_GBJ,safe_z(PG_GBJ$GBJ_pvalue)*RG$bs,NA),
                  GHC=ifelse(apply_GHC,safe_z(PG_GHC$GHC_pvalue)*RG$bs,NA),
                  mnP=ifelse(apply_mnP,safe_z(PG_mnP$minP_pvalue)*RG$bs,NA),
                  ACAT=ifelse(apply_ACAT,safe_z(PG_ACAT)*RG$bs,NA),
                  leadSNP=paste(unique(c(names(RG$b)[which.max(abs(RG$b))],names(RG$b)[which(abs(RG$b)>bSNP)])),collapse = ";")),t(RG$b)))
    })
    names(PV)=names(Y)
    return(data.table::rbindlist(PV,idcol="trait"))
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  names(GS)=nM
  return(data.table::rbindlist(GS,idcol="gene",fill = fill,use.names = F))
}

