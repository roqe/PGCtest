para=function(hypo,sm,mm,vv,m){
  if(hypo=="H00"){
    beta=c(rep(0,m))
  }else{
    beta=c(rnorm(sm,mm,vv),rep(0,(m-sm)))
    #beta=c(rep(mm,sm),rep(0,(m-sm)))
  }
  return(beta)
}
