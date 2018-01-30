
# Example LD site regression, 
# Assumes clinal region with K=2 with ancestral populations P_a and P_b, alternatively 1 and 2
# Names of loci 
loci = paste0("Locus_",c("A","B","C"))

# Names of sites (ie localities)
sites = paste0("Site_",c("X","Y","Z"))

# number of diagnostic loci (3 in this example)
n_loci = length(loci)

# Observed allelle (from parental population 2) frequency estimate at each site for each locus
p_obs = matrix(0,length(sites),length(loci))
dimnames(p_obs) = list(site=sites,locus=loci)
p_obs["Site_X",] = c(0.1,0.18,0.03)
p_obs["Site_Y",] = c(0.5,0.6,0.4)
p_obs["Site_Z",] = c(0.95,0.9,0.85)

# Observed allelle (from parental population 1) frequency estimate at each site for each locus
q_obs = 1-p_obs

# Number of observed allelles from parental population 2
z_obs = 2 * p_obs

# Observed Q value
Q_obs=apply(z_obs, 1, sum) / (2 * length(loci))

# Expected Q value (from cline model)
Q_est=c(0.11,0.49,0.9)
names(Q_est)<-sites

# Expected Q variance (from cline model)
Q_var=c(0.0155,0.0417,0.0151)
names(Q_var)<-sites

# site variance of allelle frequency 
p_var=apply(p_obs, 1, var)

# expected allelle frequence variance
pq_mean = apply(p_obs*q_obs,1, mean)

# reference equation 
# var(z) = 1/2n * (mean.z * (1 - mean.z) - var(p)) + 1/2 * (1 - 1/n) * mean.D_linkage.disequilibrium
# Q_var  = 1/(2*n_loci) * (Q_est * (1 - Q_est) - p_var) + 1/2 * (1 - 1/n_loci) * mean.LD
# 2*Q_var*n_loci = Q_est * (1 - Q_est) - p_var + (n_loci - 1) * mean.LD
# PQ_est = Q_est * (1 - Q_est)
site_stats=data.frame(Q_est=Q_est, Q_var=Q_var,p_var=p_var,pq_mean=pq_mean,row.names = sites)
site_stats$PQ_est = Q_est * (1 - Q_est)
site_stats$Q_var_scaled = 2*Q_var*n_loci
site_stats$mean.LD= ( 2*Q_var*n_loci+p_var - Q_est * (1 - Q_est) )/(n_loci-1)





