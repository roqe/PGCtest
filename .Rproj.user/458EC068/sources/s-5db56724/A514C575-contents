#' Simulate data for analysis.
#'
#' @param hypo Hypothesis type: "H00" refers to the complete null hypothese (default);
#'             "H02" refers to the cases with disjoint signals, if rho=0 then it is a null, otherwise alternative;
#'             "HA" refer to the alternative hypothesis with connected signals.
#' @param sample_size Number of samples, default=300.
#' @param num_mediators Number of mediators, default=30.
#' @param num_mechanisms Number of mechanisms, default=2000.
#' @param num_covariates Number of covariates, default=2.
#' @param rho The correlation coefficients between the mediators, default=0.3.
#' @param har The percentage of mechanisms with signals within one experiment, default=0.2. (affect when hypo is not "H00")
#' @param mm The mean of the signal, default=0.1. (affect when hypo is not "H00")
#' @param vv The variance of the signal, default=0.1. (affect when hypo is not "H00")
#' @param sm Number of mediators with signal within one mechanism, default=5. (affect when hypo is not "H00")
#' @param mc Number of cores for parallel computing, default=5.
#' @export
#' @importFrom parallel mclapply
#' @importFrom MASS mvrnorm
#' @examples
#' H00=sim_mediation_data()
#' H00_ss=sim_mediation_data(sample_size=100)
#' H00_rh=sim_mediation_data(rho=0)
#' H02=sim_mediation_data(hypo="H02",har=1,mm=0,vv=0.05,sm=2)
#' HA=sim_mediation_data(hypo="HA",mm=0.1,vv=0.1,sm=10)

sim_mediation_data=function(hypo="H00",sample_size=300,num_mediators=30,num_mechanisms=2000,
                            num_outcomes=20,num_covariates=2,rho=0.3,har=0.2,mm=0.1,vv=0.1,sm=5,mc=5){
  X=sapply(1:num_covariates,function(x){ return(rnorm(sample_size)) })
  mi=rep(num_mediators,num_mechanisms)
  D=parallel::mclapply(mi,function(m){
    if(runif(1)>har){ sm=0 }
    SS=matrix(rho,nrow=m,ncol=m)
    diag(SS)=1
    M=MASS::mvrnorm(n=sample_size,mu=rnorm(m),Sigma=SS)
    Y=sapply(1:num_outcomes,function(g){
      if(hypo=="H02"){
        if(runif(1)>0.5){ beta=para(hypo,sm,mm,vv,m) }else{ beta=para(hypo="HA",sm,mm,vv,m) }
      }else{
        beta=para(hypo,sm,mm,vv,m)
      }
      return(apply(X,1,sum)+M%*%beta+rnorm(sample_size))
    })
    return(list(M=M,Y=Y))
  }, mc.cores = mc, mc.set.seed = T)
  M=lapply(D,function(d){ return(d$M) })
  Y=lapply(D,function(d){ return(d$Y) })
  return(list(M=M,Y=Y,X=X))
}
