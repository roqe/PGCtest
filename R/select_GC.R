#' @export
select_GC=function(preR,top=5){
  GBB=rbindlist(apply(preR,1,function(gtp){
    return(data.table(gene=gtp[1],trait=gtp[2],maxZ=max(abs(as.numeric(gtp[3:7])))))
  }))
  GB=GBB[order(trait,maxZ,decreasing = T)]
  GB$rank=rep(1:length(unique(preR$gene)),length(unique(preR$trait)))
  GG=GB[rank%in%1:(top+1)] #top 11 trait genes
  GC=split(GG$gene,GG$trait)
  return(GC)
}
