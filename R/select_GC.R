#' @export
select_GC=function(preR,top=5){
  GB=preR[order(trait,abs(mnP),decreasing = T)]
  GB$rank=rep(1:length(unique(preR$gene)),length(unique(preR$trait)))
  GG=GB[rank%in%1:(top+1)] #top 11 trait genes
  GC=split(GG,GG$trait)
  return(GC)
}
