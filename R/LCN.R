#' Local approach: BJ, HC, minP
#'
#' @param PS The pre-fitting results (optional when using VCT, use for determining the signs of the statistics if available)
#' @param method We provide "BJ","HC", and "minP", default="minP".
#' @param mc Number of cores for parallel computing, default=5.
#' @importFrom parallel mclapply
#' @import GBJ
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' PS=pre_stat(HA$S,HA$M,HA$Y,HA$X)
#' CSp=cn_stat(PS,trans = T)
#' BJp=LCN(CSp,"BJ")
#' HCp=LCN(CSp,"HC")
#' MPp=LCN(CSp,"minP")
#' CSm=cn_stat(PS,trans = F)
#' BJm=LCN(CSm,"BJ")
#' HCm=LCN(CSm,"HC")
#' MPm=LCN(CSm,"minP")

LCN=function(CS,method="minP",mc=5){
  LS=parallel::mclapply(CS,function(cs){
    Zcn=cs$Mcn; Zsb=cs$Msb; Zjs=cs$Mjs; SS=round(cs$MR,digits = 2)
    if(method=="BJ"){
      if(is.na(Zcn)){
        return(list(BJcn=NA,Pcn=NA,BJsb=NA,Psb=NA,BJjs=NA,Pjs=NA)) }
      cn=GBJ::BJ(Zcn,SS); pc=cn$BJ_pvalue; zc=cn$BJ
      sb=GBJ::BJ(Zsb,SS); ps=sb$BJ_pvalue; zs=sb$BJ
      js=GBJ::BJ(Zjs,SS); pj=js$BJ_pvalue; zj=js$BJ
      return(list(BJcn=zc,Pcn=pc,BJsb=zs,Psb=ps,BJjs=zj,Pjs=pj))
    }else if(method=="HC"){
      if(is.na(Zcn)){
        return(list(HCcn=NA,Pcn=NA,HCsb=NA,Psb=NA,HCjs=NA,Pjs=NA)) }
      cn=GBJ::HC(Zcn,SS); pc=cn$HC_pvalue; zc=cn$HC
      sb=GBJ::HC(Zsb,SS); ps=sb$HC_pvalue; zs=sb$HC
      js=GBJ::HC(Zjs,SS); pj=js$HC_pvalue; zj=js$HC
      return(list(HCcn=zc,Pcn=pc,HCsb=zs,Psb=ps,HCjs=zj,Pjs=pj))
    }else{
      if(is.na(Zcn)){
        return(list(MPcn=NA,Pcn=NA,MPsb=NA,Psb=NA,MPjs=NA,Pjs=NA)) }
      cn=GBJ::minP(Zcn,SS); pc=cn$minP_pvalue; zc=cn$minP
      sb=GBJ::minP(Zsb,SS); ps=sb$minP_pvalue; zs=sb$minP
      js=GBJ::minP(Zjs,SS); pj=js$minP_pvalue; zj=js$minP
      return(list(MPcn=zc,Pcn=pc,MPsb=zs,Psb=ps,MPjs=zj,Pjs=pj))
    }
  },mc.cores = mc,mc.preschedule = T,mc.cleanup = T)
  names(LS)=names(CS)
  return(data.table::rbindlist(LS,idcol = "ensg"))
}
