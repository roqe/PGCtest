#' @export
pre_fit_one=function(M,Y,X){
  if(anyNA(X)){
    print("warning: missing covariate data")
    X<-matrix(rep(1, length(Y)))
  }
  if(all(Y %in% 0:1)){
    fit.y<-glm(Y~.+M, family="binomial", data = X)
  }else{
    fit.y<-lm(Y~.+M, data = X)
  }

  fit.y.sum<-summary(fit.y)
  lengX=0
  for(x in 1:ncol(X)){
    lengX=lengX+1
    if(class(X[,x])=="factor"){
      lengX=lengX+(length(levels(X[,x]))-2)
    }
  }

  indB=(lengX+2):nrow(fit.y.sum$coefficients)
  b.hat<-fit.y.sum$coefficients[indB, 1]
  b.sd<-fit.y.sum$coefficients[indB, 2]
  b=b.hat/b.sd

  indB=which(!is.na(fit.y$coefficients))[indB]
  covB=as.matrix(vcov(fit.y)[indB,indB])

  cb.svd<-svd(covB)
  b.sum<-sum(t(b.hat)%*%cb.svd$u%*%diag(1/sqrt(cb.svd$d),ncol=length(b.hat)))
  corrB=cov2cor(covB)

  return(list(b=b,corrB=corrB,bs=ifelse(b.sum>0,1,-1)))
}
