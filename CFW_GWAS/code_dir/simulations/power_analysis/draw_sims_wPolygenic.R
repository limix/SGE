library(MASS)

# code below is messy because I added new simulation sets separately. Important is that hard coding below needs to be consistent with sims_design.RData and efficient_map_SGE_single_snp.py
row_sims_design = as.numeric(commandArgs(trailingOnly = TRUE))
if (61 <= row_sims_design & row_sims_design <= 72) {
	gen_model = 'colMax'
} else if (73 <= row_sims_design & row_sims_design <= 84) {
	gen_model = 'add'
}

load('data_dir/simulations/power_analysis/sims_design.RData')

if (gen_model == 'colMax'){
	new_gen_model = '_colMax'
	new_row_sims_design = row_sims_design-dim(sims_design)[1]+12
} else if (gen_model == 'add') {
	new_gen_model = NULL
	new_row_sims_design = row_sims_design-dim(sims_design)[1]
} else {
	new_gen_model = NULL
	new_row_sims_design = row_sims_design
}
gs = sims_design[new_row_sims_design,'gs']

nb_sims=300

#load("/homes/abaud/CFW/output/reproduce/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/all_VCs.RData")
#all_VCs = all_VCs[all_VCs[,'var_Ad']>0.1 & all_VCs[,'var_As']>0.1,]
#medians = apply(all_VCs[,-c(1,12)],FUN=median,MAR=2)
#medians
var_Ad = 0.17
var_As = 0.17
corr_Ads = 0.65
var_Ed = 0.19
var_Es = 0.15
corr_Eds = -0.8
var_C = 0.25

load(paste('data_dir/simulations/power_analysis/image_simulations_gs',gs,'.RData',sep=''))
overall_cov = var_Ad * scaled_K + corr_Ads*sqrt(var_Ad)*sqrt(var_As) * (scaled_crossK%*%t(Z) + Z %*% t(scaled_crossK)) + var_As * Z%*%scaled_K_cm%*%t(Z) + var_Ed*scaled_env + corr_Eds*sqrt(var_Ed)*sqrt(var_Es) * (scaled_env_cross%*%t(Z) + Z %*% t(scaled_env_cross)) + var_Es * Z%*%scaled_env_cm%*%t(Z) + var_C * scaled_cage_mat
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


load(paste('data_dir/simulations/power_analysis/dosages/dir_soc_dosages_gs',gs,new_gen_model,'.RData', sep=''))
#loads dosages and social_dosage_sumUnstd
#social_dosage_sumUnstd range between 0 and 2*(gs-1) (as it was calculated using colSums function in build_social_genotypes.R)
all(rownames(dosages) == rownames(simulations))
#TRUE
all(rownames(dosages) == rownames(social_dosage_sumUnstd))
#TRUE

#if want to simulate from additive model: use social_dosage_sumUnstd as loaded
#if want to simulate from max model, use colMax-derived social_dosage_sumUnstd
#if want to simulate from proportional model:
if (gen_model == 'prop') {
	social_dosage_sumUnstd = social_dosage_sumUnstd/(gs - 1)
}

load(paste('data_dir/simulations/power_analysis/SNPsCanUse/SNPsCanUse',new_row_sims_design,new_gen_model,'.RData', sep=''))
#loads SNPsCanUse which is a subset of colnames(dosages) 
if (any(!SNPsCanUse %in%  colnames(dosages))) stop('pb overlap SNPs')

SNPs2use = sample(SNPsCanUse, size = nb_sims, replace = F)
sgeQTLs = social_dosage_sumUnstd[,paste(SNPs2use,'_gs',gs,sep='')] * sims_design[new_row_sims_design,'sgeQTL_a']
dgeQTLs = dosages[,SNPs2use] * sims_design[new_row_sims_design,'dgeQTL_a']

simulations = simulations + sgeQTLs + dgeQTLs
colnames(simulations) = paste('simulation_sd',row_sims_design,'_',SNPs2use,sep='')

save(simulations,file = paste('data_dir/simulations/power_analysis/simulations/simulations_sd',row_sims_design,'.RData',sep=''))
