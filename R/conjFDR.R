#' @export
conjFDR = function(a, b, p_cut = 1E-3){

  pa=pnorm(abs(a),lower.tail = F)
  pb=pnorm(abs(b),lower.tail = F)
  ndx=(pa<=p_cut)|(pb<=p_cut)
  cfdr1=cfdr2=ccfdr=rep.int(1,length(a))
  if(!any(ndx)){ return(ccfdr) }
  p1=pa[ndx]
  p2=pb[ndx]
  denoms = matrix( sapply(1:length(p1), calc.denoms, p1=p1, p2=p2), nrow = 2)
  cfdr1[ndx]=p1/denoms[1,]
  cfdr2[ndx]=p2/denoms[2,]
  ccfdr=pmax(cfdr1,cfdr2)
  cc=data.table(za=a,zb=b,cfdr1=cfdr1,cfdr2=cfdr2,ccfdr=ccfdr)

  return(cc)
}

calc.denoms = function(i, p1, p2){
  ## The edge point
  x = p1[i]
  y = p2[i]

  ## Vectors
  dd1 = p2 <= y
  dd2 = p1 <= x
  ee = dd1 & dd2

  ## Combine
  ee_n = length(which(ee))
  denom1 =  ee_n / length(which(dd2))
  denom2 =  ee_n / length(which(dd2))

  return(c(denom1, denom2))
}

