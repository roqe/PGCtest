#' Simulate data for analysis.
#'
#' @param hypo Hypothesis type: "H00" refers to the complete null hypothese (default);
#'             "H02" refers to the cases with disjoint signals, if rho=0 then it is a null, otherwise alternative;
#'             "HA" refer to the alternative hypothesis with connected signals.
#' @param sample_size Number of samples, default=1000.
#' @param num_variants Number of variants within one gene, default=30.
#' @param num_genes Number of genes, default=10000.
#' @param num_traits Number of outcome, default=5.
#' @param num_covariates Number of covariates, default=2.
#' @param rho The correlation coefficients between the variants, default=0.3.
#' @param har The percentage of genes with signals within one experiment, default=0.1. (affect when hypo is not "H00")
#' @param gpr The percentage of signal overlaps between genes, default=0.5. (affect when hypo is "HA")
#' @param mm The mean of the signal, default=0.1. (affect when hypo is not "H00")
#' @param vv The variance of the signal, default=0.1. (affect when hypo is not "H00")
#' @param sm Number of mediators with signal within one gene, default=5. (affect when hypo is not "H00")
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' H00=sim_mediation_data()
#' H00_ss=sim_mediation_data(sample_size=100)
#' H00_rh=sim_mediation_data(rho=0)
#' H02=sim_mediation_data(hypo="H02",har=0.1,mm=0,vv=0.05,sm=2)
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)

sim_pleiotropy_data=function(hypo="H00",sample_size=1000,num_variants=20,num_genes=10000,
                             num_traits=5,num_covariates=2,rho=0,har=0.01,gpr=1,mm=0,vv=0,sm=2){
  num_pleio_gene=num_genes*har
  X=data.frame(sapply(1:num_covariates,function(x){ return(rnorm(sample_size,mean = 1)) }))
  mi=rep(num_variants,num_genes)
  SS=matrix(0,nrow=num_variants,ncol=num_variants)
  for(i in 1:(num_variants/4)){
    for(j in 1:(num_variants/4)){
      SS[i,j]=rho
    }
  }
  diag(SS)=1
  M=lapply(mi,function(m){
    M=MASS::mvrnorm(n=sample_size,mu=rnorm(m,1),Sigma=SS)
    return(M)
  })
  Y=data.frame(sapply(1:num_traits,function(g){
    if(hypo!="HA"){ gpr=0 }
    if(g<num_traits/2){ # first cluster
      jump=(1-gpr)*(g-1)*num_pleio_gene
      MM=M[(jump+1):(jump+num_pleio_gene)]
      #print((jump+1):(jump+num_pleio_gene))
    }else{ # second cluster
      jump=(1-gpr)*(num_traits-g)*num_pleio_gene
      MM=M[(length(M)-num_pleio_gene+1-jump):(length(M)-jump)]
      #print((length(M)-num_pleio_gene+1-jump):(length(M)-jump))
    }
    MMM=sapply(MM,function(mt){
      beta=para(hypo,sm,mm,vv,ncol(mt))
      return(mt%*%beta)
    })
    return(apply(X,1,sum)+apply(MMM,1,sum)+rnorm(sample_size))
  }))
  names(X)=paste0("c",1:ncol(X))
  names(Y)=paste0("p",1:ncol(Y))
  names(M)=paste0("g",1:length(M))
  M=sapply(names(M),function(nm){
    SV=data.frame(M[[nm]])
    names(SV)=paste0(nm,"_s",1:ncol(SV))
    return(list(SV))
  })
  return(list(M=M,Y=Y,X=X))
}
