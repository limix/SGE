load('output_dir/effect_sizes_additive_model/betas_conditioning/allelic_effect/betas.RData')
head(indiv_SGE_QTLs_betas)
#            V1                    V2  sig_beta1
#Bioch.ALAT_BWcorr chr14_98852062_social 0.03445159
add_indiv_SGE_QTLs_betas = indiv_SGE_QTLs_betas
load('output_dir/effect_sizes_additive_model/betas_conditioning/allelic_effect/betas.RData')
head(indiv_SGE_QTLs_betas)
#            V1                    V2  sig_beta1
#Bioch.ALAT_BWcorr chr14_98852062_social 0.03445159
prop_indiv_SGE_QTLs_betas = indiv_SGE_QTLs_betas
all(add_indiv_SGE_QTLs_betas[,1] == prop_indiv_SGE_QTLs_betas[,1])
all(add_indiv_SGE_QTLs_betas[,2] == prop_indiv_SGE_QTLs_betas[,2])
#both TRUE
indiv_SGE_QTLs_betas = cbind(add_indiv_SGE_QTLs_betas, prop_indiv_SGE_QTLs_betas[,3])
colnames(indiv_SGE_QTLs_betas) = c('Measure','sgeQTL','add_allelic_effect','prop_allelic_effect')

#to get allele dosages
load('output_dir/effect_sizes/DGE_sig_QTLs_genos.RData')
#dosages btw 0 and 2 but means between 0 and 1 as these are number of minor alleles
#load('output_dir/effect_sizes/SGE_sig_QTLs_social_genos.RData')
#social_dosage_sumUnstd btw 0 and 4 but means between 0 and 2 as number of minor alleles for additive model
load('output_dir/effect_sizes_additive_model/SGE_sig_QTLs_social_genos.RData')
all(rownames(dosages) == rownames(social_dosage_sumUnstd))
#[1] TRUE

load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')
#to get sample size (confirm with all_VCs) and realised minor allele dosage
library(rhdf5)
null_covars_dir = 'output_dir/null_covars/pruned_dosages_include_DGE_IGE_IEE_cageEffect/'
indiv_DGE_QTLs_betas$sample_size = NA
#MAD is mean allele dosage
indiv_DGE_QTLs_betas$MAD = NA
for (i in 1:dim(indiv_DGE_QTLs_betas)[1]) {
	pheno = indiv_DGE_QTLs_betas[i,1]
	mice = h5read(paste(null_covars_dir,pheno,'.h5',sep=''), '/sampleID')
	indiv_DGE_QTLs_betas[i,'sample_size'] = length(mice)
	if (all_VCs[pheno,'sample_size'] != indiv_DGE_QTLs_betas[i,'sample_size']) stop('pb sample size')
	motch = match(mice, rownames(social_dosage_sumUnstd))
	if (any(is.na(motch))) stop('pb')
	MAD = mean(dosages[motch,indiv_DGE_QTLs_betas[i,2]])
	indiv_DGE_QTLs_betas[i,'MAD'] = MAD
}
indiv_SGE_QTLs_betas$sample_size = NA
indiv_SGE_QTLs_betas$MAD = NA
for (i in 1:dim(indiv_SGE_QTLs_betas)[1]) {
	pheno = indiv_SGE_QTLs_betas[i,1]
	mice = h5read(paste(null_covars_dir,pheno,'.h5',sep=''), '/sampleID')
	indiv_SGE_QTLs_betas[i,'sample_size'] = length(mice)
	if (all_VCs[pheno,'sample_size'] != indiv_SGE_QTLs_betas[i,'sample_size']) stop('pb sample size')
	motch = match(mice, rownames(social_dosage_sumUnstd))
	if (any(is.na(motch))) stop('pb')
	MAD = mean(social_dosage_sumUnstd[motch,indiv_SGE_QTLs_betas[i,2]])
	indiv_SGE_QTLs_betas[i,'MAD'] = MAD
}

# get logPs, 
effect = 'SGE'
load(paste('output_dir/pvalues_pointwise_conditional/',effect,'/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
SGE_QTLs = all_QTLs
effect = 'DGE'
load(paste('output_dir/pvalues_pointwise_conditional/',effect,'/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
DGE_QTLs = all_QTLs
head(DGE_QTLs)
#                   measure      marker chr      pos     cumpos
#   FACS.CD45posCD3posCD8pos merged_peak  17 34370716 2238018246
#  FACS.CD45posCD3negDX5pos merged_peak  17 34284719 2237932249
#          FACS.CD3posCD4pos merged_peak  17 34887946 2238535476
motch = match(paste(indiv_DGE_QTLs_betas[,1],sub('_direct','',indiv_DGE_QTLs_betas[,2]),sep='_'),paste(DGE_QTLs[,'measure'],'_chr',DGE_QTLs[,'chr'],'_',DGE_QTLs[,'pos'],sep=''))
any(is.na(motch))
#FALSE
DGE_QTLs = DGE_QTLs[motch,]
indiv_DGE_QTLs_betas$logP = DGE_QTLs$logP
motch = match(paste(indiv_SGE_QTLs_betas[,1],sub('_social','',indiv_SGE_QTLs_betas[,2]),sep='_'),paste(SGE_QTLs[,'measure'],'_chr',SGE_QTLs[,'chr'],'_',SGE_QTLs[,'pos'],sep=''))
any(is.na(motch))
#FALSE
SGE_QTLs = SGE_QTLs[motch,]
indiv_SGE_QTLs_betas$logP = SGE_QTLs$logP

af_indiv_SGE_QTLs_betas = indiv_SGE_QTLs_betas
af_indiv_DGE_QTLs_betas = indiv_DGE_QTLs_betas
load('output_dir/effect_sizes_additive_model/betas_conditioning/var_expl/betas.RData')
ve_indiv_SGE_QTLs_betas = indiv_SGE_QTLs_betas
ve_indiv_DGE_QTLs_betas = indiv_DGE_QTLs_betas

all(af_indiv_SGE_QTLs_betas[,1] == ve_indiv_SGE_QTLs_betas[,1])
all(af_indiv_DGE_QTLs_betas[,1] == ve_indiv_DGE_QTLs_betas[,1])
#both TRUE

af_indiv_SGE_QTLs_betas$var_expl = ve_indiv_SGE_QTLs_betas[,'sig_beta1']
indiv_SGE_QTLs_betas = af_indiv_SGE_QTLs_betas
af_indiv_DGE_QTLs_betas$var_expl = ve_indiv_DGE_QTLs_betas[,'sig_beta1']
af_indiv_DGE_QTLs_betas$allelic_effect = af_indiv_DGE_QTLs_betas$sig_beta1
af_indiv_DGE_QTLs_betas = af_indiv_DGE_QTLs_betas[,-which(colnames(af_indiv_DGE_QTLs_betas) == 'sig_beta1')]
indiv_DGE_QTLs_betas = af_indiv_DGE_QTLs_betas

save(indiv_SGE_QTLs_betas,indiv_DGE_QTLs_betas,file = 'output_dir/effect_sizes_additive_model/betas_conditioning/all_betas.RData')
