#' Simulate data for analysis.
#'
#' @param hypo Hypothesis type: "H00" refers to the complete null hypothese (default);
#'             "H02" refers to the cases with disjoint signals, if rho=0 then it is a null, otherwise alternative;
#'             "HA" refer to the alternative hypothesis with connected signals.
#' @param sample_size Number of samples, default=300.
#' @param num_variants Number of variants within one gene, default=20.
#' @param num_genes Number of genes, default=10000.
#' @param num_pleio_genes Number of pleiotropic genes, default=30.
#' @param num_traits Number of outcome, default=5.
#' @param num_covariates Number of covariates, default=2.
#' @param rho The correlation coefficients between the variants, default=0.
#' @param gpr The percentage of signal overlaps between genes, default=0.
#' @param mm The mean of the signal, default=0.
#' @param vv The variance of the signal, default=0.
#' @param sm Number of mediators with signal within one gene, default=2.
#' @export
#' @importFrom MASS mvrnorm
#' @examples
#' H00=sim_mediation_data()
#' H00_ss=sim_mediation_data(sample_size=100)
#' H00_rh=sim_mediation_data(rho=0)
#' H02=sim_mediation_data(hypo="H02",mm=0,vv=0.05,sm=2)
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)

sim_pleiotropy_data=function(hypo="H00",sample_size=300,num_variants=20,num_genes=10000,num_pleio_gene=20,
                             num_traits=5,num_covariates=2,rho=0,gpr=0,mm=0,vv=0,sm=2,bin=T){
  X=data.frame(sapply(1:num_covariates,function(x){ return(rnorm(sample_size)) }))
  names(X)=paste0("c",1:ncol(X))
  SS=matrix(0,nrow=num_variants,ncol=num_variants)
  for(i in 1:(num_variants/4)){
    for(j in 1:(num_variants/4)){
      SS[i,j]=rho
    }
  }
  diag(SS)=1
  mi=rep(num_variants,num_genes)
  M=lapply(1:length(mi),function(m){
    M=data.frame(MASS::mvrnorm(n=sample_size,mu=rnorm(mi[m],1),Sigma=SS))
    if(bin){ M=data.frame(apply(M,2,function(m){ return(as.numeric(m>1)) })) }
    names(M)=paste0("g",m,"_s",1:ncol(M))
    return(M)
  })
  names(M)=paste0("g",1:length(mi))
  print("simulation paramter --- ")
  print(paste("hypothesis:",hypo,"/ correlation coef.:",rho,"/ gene overlapping %:",gpr,"/ variants with signal:",sm))
  print(paste("smaple size:",sample_size,"/ total gene:",num_genes,"/ pleiotropic genes:",num_pleio_gene))
  print(paste0("signal: Normal(",mm,",",vv,")"))
  Y=data.frame(sapply(1:num_traits,function(p){
    if(hypo!="HA"){ gpr=0 }
    if(p<num_traits/2){ # first cluster
      jump=(1-gpr)*(p-1)*num_pleio_gene
      MM=M[(jump+1):(jump+num_pleio_gene)]
    }else{ # second cluster
      jump=(1-gpr)*(num_traits-p)*num_pleio_gene
      MM=M[(length(M)-num_pleio_gene+1-jump):(length(M)-jump)]
    }
    print(paste0("pleiotropic genes in trait ",p,": ",names(MM)[1],"-",names(MM)[length(MM)]))
    MMM=sapply(MM,function(mt){
      beta=para(hypo,sm,mm,vv,ncol(mt))
      #print(beta)
      return(as.matrix(mt)%*%beta)
    })
    return(apply(X,1,sum)+apply(MMM,1,sum)+rnorm(sample_size))
  }))
  names(Y)=paste0("p",1:ncol(Y))

  return(list(M=M,Y=Y,X=X))
}
