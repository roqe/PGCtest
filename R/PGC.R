#' @param stat Statistics evaluating the mechansisms
#' @param mc Number of cores for parallel computing, default=5.
#' @param prec The unit to construct empirical normal product pdf using in composite test, default=1e-3.
#' @import data.table
#' @export
#' @examples
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)
#' stat=GCN(HA)
#' FR=PGC(stat)

PGC=function(stt,mc=5,dcrr=T){
  RR=list()
  for(app in c("TSQ","GBJ","GHC","mnP","ACAT")){
    ZZZ=data.frame(data.table::dcast(stt,gene~trait,value.var = app))
    Q=apply((combn(ncol(ZZZ)-1,2)+1),2,function(id){
      a=unlist(ZZZ[,id[1]])
      b=unlist(ZZZ[,id[2]])
      if(dcrr){
        rh=cor(a,b)
        ss=svd(matrix(c(1,rh,rh,1),nrow=2))
        zz=ss$u%*%diag(1/sqrt(ss$d))%*%ss$v%*%rbind(a,b)
        a=zz[1,]
        b=zz[2,]
      }
      Mcn=CompTest(a,b,mc)
      Msb=Sobel(a,b)
      Mjs=JointSig(a,b)
      return(data.table::data.table(ensg=ZZZ$gene,app=app,Za=a,Zb=b,Pcn=Mcn$pp,Psb=Msb$pp,Pjs=Mjs$pp))
    })
    names(Q)=apply((combn(ncol(ZZZ)-1,2)+1),2,function(id){
      return(paste(colnames(ZZZ)[id],collapse = "_"))
    })
    RR=c(RR,Q)
  }
  return(RR)
}
