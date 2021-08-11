TSQ=function(Z,cm){
  TS=t(Z)%*%ginv(cm,tol=1e-8)%*%Z
  r=rankMatrix(cm,tol=1e-8)[1]
  pv=pchisq(TS,r,lower.tail = F)
  return(list(TSQ=TS,pv=pv))
}
