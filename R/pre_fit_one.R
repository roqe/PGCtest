#' @export
pre_fit_one=function(M,Y,X){
  if(all(Y %in% 0:1)){
    fit.y<-glm(Y~., family="binomial", data = cbind(M,X))
  }else{
    fit.y<-lm(Y~., data = cbind(M,X))
  }

  fit.y.sum<-summary(fit.y)
  indB=2:(ncol(M)+1)
  indB=indB[which(indB%in%which(!is.na(fit.y$coefficients)))]
  b.hat<-fit.y.sum$coefficients[2:(length(indB)+1), 1]
  b.sd<-fit.y.sum$coefficients[2:(length(indB)+1), 2]
  b=b.hat/b.sd
  covB=as.matrix(vcov(fit.y)[indB,indB])
  cb.svd<-svd(covB)
  b.sum<-sum(t(b.hat)%*%cb.svd$u%*%diag(1/sqrt(cb.svd$d),ncol=length(b.hat)))
  corrB=cov2cor(covB)

  return(list(b=b,corrB=corrB,bs=ifelse(b.sum>0,1,-1),bh=b.hat,be=b.sd))
}
