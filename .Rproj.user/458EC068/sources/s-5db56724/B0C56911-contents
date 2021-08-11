para=function(hypo,sm,mm,vv,m){
  if(hypo=="H00"){
    beta=c(rep(0,m))
  }else if(hypo=="HA"){
    beta=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
  }else{
    beta=c(rep(0,(m-sm)),rnorm(sm,mm,vv))# ##equals to the length of m
  }
  return(beta)
}
