
## simulation study of LD

# assumptions:
# K = 2
# Ancestral populations are known as A and B
# Sampling sites are at a single point
# In the absence of selection, all loci are independently sorting
# p is the frequency of population B alleles
# q is the frequency of population A alleles
# n is the number of loci
# z is the number of population B alleles
# Q is the admixture proportion of population B alleles
# Hardy-Weinburg equilibrium (ie AB hybrids occur with frequence 2pq)
# All loci are autosomal

# Vague terms:
# Arbitrary number of individuals sampled from an arbitrary number of sites
# values can be true, observed, and estimated
# values can also be columns, matrices, vectors, sums, means or variances
# values are site specific
# z values when aggregated are summed across loci
# all other aggregations are averaged unless otherwise specified
# matrices and columns are only for observed values
# matrices contain observed values for each locus and individual
# columns contain observed values aggregated across loci
# vectors either contain observed values aggregated across individuals 
#   or true / estimated values for each locus
# eg the true vector of population B allele frequencies is pTvec

LD_mean_1 <- function(Q_var,Q_est,p_var,n_loci) 
  ( 2*Q_var*n_loci+p_var - Q_est * (1 - Q_est) )/(n_loci-1)

LD_mean_2 <- function(p_mat) 
  LD_mean_1(var(apply(p_mat,1,mean)),mean(p_mat),var(apply(p_mat,2,mean)),ncol(p_mat))

LD_mean_3 <- function(z_mat)
  LD_mean_2(z_mat/2)

obs_stats <- function(p_mat){
  res<- data.frame(Q_var=var(apply(p_mat,1,mean)),
                   Q_est=mean(p_mat),
                   p_var=var(apply(p_mat,2,mean)),
                   n_loci=ncol(p_mat))
  p_est=apply(p_mat,2,mean)
  res$pq_mean= mean(p_est*(1-p_est))
  res$LD_mean= with(res,LD_mean_1(Q_var,Q_est,p_var,n_loci) )
  res
}
            
## Very simple situation (A)
# assume p_true is constant across loci and loci are independent
A_res_true <- function(p_true, n_loci){
  z_true = 2*p_true*n_loci
  q_true = 1-p_true
  z_var = 2*n_loci*p_true*q_true
  Q_var = p_true*q_true / (2*n_loci*p_true)
  p_var = p_true*q_true / (2*n_loci*p_true)
  
}

A_z_mat <-  function(p_true, n_loci, ind_count)
  matrix(rbinom(n_loci*ind_count,2,p_true),ncol=n_loci)


## Some fun to be had...
library(foreach)

summary(foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
library(lattice)
xyplot(LD_mean~pq_mean, data=(foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
)
xyplot(p_var~Q_var, data=foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
xyplot(p_var~Q_var, data=foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
xyplot(p_var~Q_var, data=foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
xyplot(p_var~Q_var, data=foreach(a=1:300,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
xyplot(pq_mean~Q_est, data=foreach(a=1:300,.combine = rbind) %do% obs_stats(A_z_mat(0.4,10,10)/2) )
xyplot(LD_mean~pq_mean, data=foreach(a=1:30,.combine = rbind) %do% obs_stats(A_z_mat(0.5,10,10)/2) )
xyplot(LD_mean~pq_mean, data=foreach(a=1:300,.combine = rbind) %do% obs_stats(A_z_mat(0.5,10,10)/2) )

