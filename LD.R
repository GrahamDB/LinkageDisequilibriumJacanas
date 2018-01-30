## Calculating Linkage Disequilibrium in order to estimate Selection

## Selection equation
# w_width = sigma_dispersal.distance * sqrt(8/s_selection)
# sigma_dispersal.distance / w_width = sqrt(s_selection/8)

## Dispersal equation
# sigma_dispersal.distance = sqrt(r_recombination.rate * mean.D_linkage.disequilibrium * w_width^2)
# sigma_dispersal.distance / w_width = sqrt(r_recombination.rate * mean.D_linkage.disequilibrium)

## Recombination rate for unlinked loci 
# r_recombination.rate = 0.5

## Selection from LD for center of cline
# sqrt(s_selection/8) = sqrt(r_recombination.rate * mean.D_linkage.disequilibrium)
# sqrt(s_selection/8)^2 = sqrt(r_recombination.rate * mean.D_linkage.disequilibrium)^2
# s_selection/8 = r_recombination.rate * mean.D_linkage.disequilibrium
# s_selection/8 = 0.5 * mean.D_linkage.disequilibrium
# s_selection = 4 * mean.D_linkage.disequilibrium

## Calculating LD from diagnostic loci for each site
# Q (STRUCTURE K = 2 admixture proportion) can be considered a hybrid index 
# based on n independent loci (discrete Mendelian markers)
# HZAR calculates Mean.Q at Distance.X along cline as well as Variance.Q

## Barton and Gale (1993) Equations - 2B pg 23
# z_hybrid.index =: Q
# var(z_hybrid.index) =: var.Q
# mean.z_hybrid.index =: mean.Q
# var(p) is the variance of allele freq across n loci
# var(z) = 1/2n * (mean.z * (1 - mean.z) - var(p)) + 1/2 * (1 - 1/n) * mean.D_linkage.disequilibrium

## Calculating LD shortcut - use allele frequencies for each locus at each site
# at each site we have n loci and mean allele frequency [p_i] for each loci 
# therefore, at each site we can calculate the variance of p_i [var(p)]
# and the mean of the allele frequency variance [mean.pq] =: mean(p_i*q_i)

# Need allele frequencies of diagnostic loci - calculate from observed samples at each site