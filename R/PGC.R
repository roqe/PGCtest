#' For pleiotropic assessment.
#'
#' @param stt Estimates of gene-trait associations, output of GCN.
#' @param dcrr Apply decorrelation or not, default=False
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

PGC=function(stt,dcrr=F){
  RR=list()
  for(app in c("TSQ","GBJ","GHC","mnP","ACAT")){
    if(all(is.na(stt[[app]]))) next
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
      Mcn=CompTestER(a,b)
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
