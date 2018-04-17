arg = as.numeric(commandArgs(trailingOnly = TRUE))

library(MASS)
load('data_dir/simulations/image_simulations.RData')

nb_sims=500

overall_cov = liste_effects [arg,'var_Ad'] * scaled_K + liste_effects [arg,'corr_Ads']*sqrt(liste_effects [arg,'var_Ad'])*sqrt(liste_effects [arg,'var_As']) * (scaled_crossK%*%t(Z) + Z %*% t(scaled_crossK)) + liste_effects [arg,'var_As'] * Z%*%scaled_K_cm%*%t(Z) + liste_effects [arg,'var_Ed']*scaled_env + liste_effects [arg,'corr_Eds']*sqrt(liste_effects [arg,'var_Ed'])*sqrt(liste_effects [arg,'var_Es']) * (scaled_env_cross%*%t(Z) + Z %*% t(scaled_env_cross)) + liste_effects [arg,'var_Es'] * Z%*%scaled_env_cm%*%t(Z) + liste_effects [arg,'var_C'] * scaled_cage_mat
n=dim(overall_cov)[1]
vecti=matrix(1,ncol=1,nrow=n)
p=diag(n)-(vecti%*%t(vecti))/n
num=(n-1)
denom=sum(diag(p%*%overall_cov%*%p))
my_scaling_factor_overall_cov=num/denom
scaled_overall_cov = overall_cov * my_scaling_factor_overall_cov
total_var = 1/my_scaling_factor_overall_cov

simulations=t(mvrnorm(n=nb_sims,mu=rep(0,dim(scaled_overall_cov)[1]),scaled_overall_cov))
colnames(simulations) = paste('simulation',paste(liste_effects [arg,],collapse = '_'), 1:dim(simulations)[2],sep='_')
rownames(simulations) = rownames(scaled_overall_cov)

target_effects = liste_effects[arg,]
target_effects[c('var_Ad','var_As','var_Ed','var_Es','var_C')] = target_effects[c('var_Ad','var_As','var_Ed','var_Es','var_C')] / total_var * 100

targetted_target_effects = liste_effects[arg,]

#targetted_target_effects is targets initially entered in liste_effects
#target_effects scales them by total_var to have proportion of phenotypic variance they actually explain in simulations

save(simulations,file = paste('data_dir/simulations/simulations/simulations',arg,'.RData',sep=''))
save(total_var, targetted_target_effects,target_effects,file = paste('data_dir/simulations/targets/targets',arg,'.RData',sep=''))
