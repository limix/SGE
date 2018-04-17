arg = as.numeric(commandArgs(trailingOnly = TRUE))

library(MASS)
load('data_dir/simulations/polygenicNdgeQTL/image_simulations.RData')

nb_sims=100

overall_cov = liste_effects [arg,'var_Ad'] * scaled_K + liste_effects [arg,'corr_Ads']*sqrt(liste_effects [arg,'var_Ad'])*sqrt(liste_effects [arg,'var_As']) * (scaled_crossK%*%t(Z) + Z %*% t(scaled_crossK)) + liste_effects [arg,'var_As'] * Z%*%scaled_K_cm%*%t(Z) + liste_effects [arg,'var_Ed']*scaled_env + liste_effects [arg,'corr_Eds']*sqrt(liste_effects [arg,'var_Ed'])*sqrt(liste_effects [arg,'var_Es']) * (scaled_env_cross%*%t(Z) + Z %*% t(scaled_env_cross)) + liste_effects [arg,'var_Es'] * Z%*%scaled_env_cm%*%t(Z) + liste_effects [arg,'var_C'] * scaled_cage_mat
n=dim(overall_cov)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%overall_cov%*%p))
#denom/num = sample var of overall_cov
total_var = denom/num
my_scaling_factor_overall_cov=1/total_var
scaled_overall_cov = overall_cov * my_scaling_factor_overall_cov
#so sample variance of scaled_overall_cov is 1

simulations=t(mvrnorm(n=nb_sims,mu=rep(0,dim(scaled_overall_cov)[1]),scaled_overall_cov))
rownames(simulations) = rownames(scaled_overall_cov)
random_only_vars = apply(simulations,FUN = var, MAR = 2)
mean_random_only_var_sims = mean(random_only_vars)
#mean_random_only_var_sims: 0.99

id_dgeQTLs = colnames(dosages)[sample(1:dim(dosages)[2],size = nb_sims, replace = F)]
std_dosages = apply(dosages[,id_dgeQTLs],FUN = scale, MAR = 2, center = T, scale = T)
rownames(std_dosages) = rownames(dosages)
motch = match(rownames(scaled_overall_cov),rownames(std_dosages))
if (any(is.na(motch))) stop('pb motch')
std_dosages = std_dosages[motch,]
#divide liste_effects [arg,'var_dgeQTL'] by 100 as now simulations have sample variance 1
dgeQTL_var = liste_effects [arg,'var_dgeQTL'] / 100
dgeQTL = std_dosages * sqrt(dgeQTL_var)
dgeQTL_vars = apply(dgeQTL,FUN = var, MAR = 2)
mean_dgeQTL_vars = mean(dgeQTL_vars)
#mean_dgeQTL_vars: 0.5

simulations = simulations + dgeQTL
colnames(simulations) = paste('simulation_le',arg,'_',id_dgeQTLs,sep='')
vars = apply(simulations,FUN = var, MAR = 2)
mean_var_sims = mean(vars)

mafs = apply(dosages[,id_dgeQTLs],FUN = mean, MAR = 2)
w = which(mafs>0.5)
mafs[w] = 1 - mafs[w]
names(mafs) = id_dgeQTLs

total_var_fixedNrandom = total_var + liste_effects [arg,'var_dgeQTL']

target_effects = liste_effects[arg,]
target_effects[c('var_Ad','var_As','var_Ed','var_Es','var_C')] = target_effects[c('var_Ad','var_As','var_Ed','var_Es','var_C')] / total_var_fixedNrandom * 100

targetted_target_effects = liste_effects[arg,]

#targetted_target_effects is targets initially entered in liste_effects
#target_effects scales them by total_var to have proportion of phenotypic variance they actually explain in simulations

save(simulations,file = paste('data_dir/simulations/polygenicNdgeQTL/simulations/simulations',arg,'.RData',sep=''))
save(total_var, total_var_fixedNrandom, mean_random_only_var_sims, mean_var_sims, targetted_target_effects,target_effects,mafs, file = paste('data_dir/simulations/polygenicNdgeQTL/targets/targets',arg,'.RData',sep=''))
