#' @export
CompTest<-function(a, b, mc=5, prec=1e-3){
  pp=zz=rep(NA,length(a))
  nx=(is.na(a)|is.na(b))
  a=a[!nx]
  b=b[!nx]
  ab=abs(a*b)
  upb=ceiling(max(ab))*3
  B=min(upb/prec,10000)
  pdf=produce_pdf(upb=upb,B=B)
  sq=upb*1:B/B
  sumpdf=sum(pdf)
  pp0=unlist(parallel::mclapply(ab,myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  pp1=unlist(parallel::mclapply(ab/sd(a),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  pp2=unlist(parallel::mclapply(ab/sd(b),myp,sq,pdf,mc.cores = mc,mc.cleanup = T))/sumpdf
  # pp0=sapply(ab,myp,sq=sq,pdf=pdf)/sumpdf
  # pp1=sapply(ab/sd(a),myp,sq=sq,pdf=pdf)/sumpdf
  # pp2=sapply(ab/sd(b),myp,sq=sq,pdf=pdf)/sumpdf
  pp.comp<-pp1+pp2-pp0
  min_pp=min(pp.comp[pp.comp>0])
  pp.comp[pp.comp<=0]=min_pp
  pp.comp[pp.comp>=1]=0.999999
  pp[!nx]=pp.comp
  zz[!nx]=safe_z(pp.comp)*ifelse(a*b>0,1,-1)
  return(list(pp=pp,zz=zz))
}
