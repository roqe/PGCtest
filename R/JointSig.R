#' @export
JointSig=function(a, b){
  zz=pmin(abs(a),abs(b))*ifelse(a*b>0,1,-1)
  pp=2*pnorm(-abs(zz))
  return(list(pp=pp,zz=zz))
}
