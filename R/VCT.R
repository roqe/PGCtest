#' @export
VCT<-function(M, Y, X=NULL, prec=1e-3){
  p<-dim(M)[2]
  n<-dim(M)[1]
  if(is.null(X)){
    X<-matrix(rep(1, n))
    q<-1
  }else{ q<-dim(X)[2]+1 }

  if(all(Y %in% 0:1)){
    fit.y0<-glm(Y~X+M, family="binomial")
  }else{ fit.y0<-lm(Y~X+M) }

  b.hat<-fit.y0$coef[(q+1):(q+p)]
  bhp=(!is.na(b.hat)&b.hat!=0&abs(b.hat)<1e+50)
  if(sum(bhp)<2){
    print(paste("no convergence, family=",ifelse(ty,"binomial","continuous")))
    b.hat=apply(M,2,function(M){
      if(ty){
        fit.y0<-glm(Y~X+M, family="binomial")
      }else{ fit.y0<-lm(Y~X+M) }
      return(fit.y0$coefficients[length(fit.y0$coefficients)])
    })
  }else{ b.hat=b.hat[which(bhp)] }

  if (all(Y %in% 0:1)) obj<-SKAT_Null_Model(Y~X, out_type="D")
  if (!all(Y %in% 0:1)) obj<-SKAT_Null_Model(Y~X, out_type="C")
  pp<-SKAT(M, obj)$p.value
  zz<-safe_z(pp)*ifelse(sum(b.hat)>0, 1, -1)
  return(list(VCT=zz,pv=pp))
}
