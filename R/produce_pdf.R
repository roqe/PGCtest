produce_pdf=function(upb=10,B=10000){
  return( besselK(x=upb*(1:B)/B, nu=0) )
}
