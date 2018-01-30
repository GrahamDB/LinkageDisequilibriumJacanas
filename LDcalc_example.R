
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
p_obs["Site_X"] = c(0.1,0.15,0.05)
p_obs["Site_Y"] = c(0.5,0.6,0.4)
p_obs["Site_Z"] = c(0.95,0.9,0.85)

# Observed allelle (from parental population 1) frequency estimate at each site for each locus
q_obs = 1-p_obs

# Number of observed allelles from parental population 2
z_obs = 2 * p_obs

# Observed Q value
Q_obs=apply(z_obs, 1, sum) / (2 * length(loci))

# Expected Q value (from cline model)
Q_est=c(0.11,0.49,0.9)

# Expected Q variance (from cline model)
Q_var=c(0.001,0.05,0.001)

# site variance of allelle frequency 
p_var=apply(p_obs, 1, var)

# expected allelle frequence variance
pq_mean = apply(p_obs*q_obs,1, mean)




