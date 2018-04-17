load('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData')
SGE_sig_QTLs = all_QTLs
SGE_sig_QTLs$id = paste('chr',SGE_sig_QTLs[,'chr'],'_',SGE_sig_QTLs[,'pos'],sep='')
load('output_dir/pvalues_pointwise_conditional/DGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData')
DGE_sig_QTLs = all_QTLs
DGE_sig_QTLs$id = paste('chr',DGE_sig_QTLs[,'chr'],'_',DGE_sig_QTLs[,'pos'],sep='')

all_measures = unique(c(SGE_sig_QTLs[,'measure'],DGE_sig_QTLs[,'measure']))
length(all_measures)
#[1] 72

load('data_dir/phenotypes/all_pheno_data.RData')
phenotypes_bc = as.matrix(phenotypes_bc)
motch = match(all_measures,colnames(phenotypes_bc))
any(is.na(motch))
phenotypes_bc = phenotypes_bc[,motch]

motch = match(all_measures,rownames(my_model_menu))
any(is.na(motch))
covariates_strings = my_model_menu[motch,'my_covs']

library(rhdf5)
fid='data_dir/CFWmice_effect_sizes_sig.h5'
h5createFile(fid)

######## phenotypes first
h5createGroup(fid,"phenotypes")
h5createGroup(fid,"phenotypes/row_header")
max(nchar(rownames(phenotypes_bc)))
#15
h5createDataset(file=fid,dataset="phenotypes/row_header/sample_ID",dims=dim(phenotypes_bc)[1],storage.mode='character',size=20)
h5write(obj=rownames(phenotypes_bc),file=fid,name="phenotypes/row_header/sample_ID")
max(nchar(cages),na.rm=T)
#7
h5createDataset(file=fid,dataset="phenotypes/row_header/cage",dims=length(cages),storage.mode='character',size=10)
h5write(obj=cages,file=fid,name="phenotypes/row_header/cage")
#h5createDataset(file=fid,dataset="phenotype/matrix",dims=dim(phenotypes))


h5write(obj=phenotypes_bc,file=fid,name="phenotypes/matrix")

h5createGroup(fid,"phenotypes/col_header")
max(nchar(colnames(phenotypes_bc)))
#29
h5createDataset(file=fid,dataset="phenotypes/col_header/phenotype_ID",dims=dim(phenotypes_bc)[2],storage.mode='character',size=35)
h5write(obj=colnames(phenotypes_bc),file=fid,name="phenotypes/col_header/phenotype_ID")

SGE_covariates_strings = covariates_strings
for (pheno in colnames(phenotypes_bc)) {
	w_sig_SGE = which(SGE_sig_QTLs[,'measure']== pheno)
	if (length(w_sig_SGE) == 0) next
	SGE_covariates_strings[pheno]= paste(c(SGE_covariates_strings[pheno],paste(SGE_sig_QTLs[w_sig_SGE,'id'],'direct',sep='_'),paste(SGE_sig_QTLs[w_sig_SGE,'id'],'social',sep='_')),collapse = ',')
# not conditioning
#	SGE_covariates_strings[pheno]= paste(c(SGE_covariates_strings[pheno],paste(SGE_sig_QTLs[w_sig_SGE,'id'],'social',sep='_')),collapse = ',')
}

DGE_covariates_strings = covariates_strings
for (pheno in colnames(phenotypes_bc)) {
	w_sig_DGE = which(DGE_sig_QTLs[,'measure']== pheno)
	if (length(w_sig_DGE) == 0) next
	DGE_covariates_strings[pheno]= paste(c(DGE_covariates_strings[pheno],paste(DGE_sig_QTLs[w_sig_DGE,'id'],'direct',sep='_'),paste(DGE_sig_QTLs[w_sig_DGE,'id'],'social',sep='_')),collapse = ',')
# not conditioning
#	DGE_covariates_strings[pheno]= paste(c(DGE_covariates_strings[pheno],paste(DGE_sig_QTLs[w_sig_DGE,'id'],'direct',sep='_')),collapse = ',')
}


h5createGroup(fid,"phenotypes/col_header/covariatesUsed")
max(nchar(SGE_covariates_strings))
#122
h5createDataset(file=fid,dataset="phenotypes/col_header/covariatesUsed/SGE",dims=length(SGE_covariates_strings),storage.mode='character',size=130)
h5write(obj=SGE_covariates_strings,file=fid,name="phenotypes/col_header/covariatesUsed/SGE")

max(nchar(DGE_covariates_strings))
#279
h5createDataset(file=fid,dataset="phenotypes/col_header/covariatesUsed/DGE",dims=length(DGE_covariates_strings),storage.mode='character',size=285)
h5write(obj=DGE_covariates_strings,file=fid,name="phenotypes/col_header/covariatesUsed/DGE")


##### covariates now with QTLs as covariates
h5createGroup(fid,"covariates")
h5createGroup(fid,"covariates/row_header")
max(nchar(rownames(cov_data)))
#15
h5createDataset(file=fid,dataset="covariates/row_header/sample_ID",dims=dim(cov_data)[1],storage.mode='character',size=20)
h5write(obj=rownames(cov_data),file=fid,name="covariates/row_header/sample_ID")

load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')
#ref allele, dosages arbitrariy between 0 and 1.
motch = match(unique(c(SGE_sig_QTLs[,'id'],DGE_sig_QTLs[,'id'])),colnames(dosages))
any(is.na(motch))
#FALSE
dosages = dosages[,motch]
motch = match(rownames(cov_data),rownames(dosages))
any(is.na(motch))
#FALSE
dosages = dosages[motch,]
colnames(dosages) = paste(colnames(dosages),'direct',sep = '_')
means = apply(dosages, FUN = mean, MAR = 2)
w = which(means >0.5)
if (length(w)!=0) {
	dosages[,w] = 1 - dosages[,w]
}
#now working with minor alleles rather than refs
dosages = dosages*2 #for interpretability: number of minor alleles
dir.create('output_dir/effect_sizes')
save(dosages,file = 'output_dir/effect_sizes/DGE_sig_QTLs_genos.RData')

load('data_dir/dosages/pruned_dosages/unstd_social_dosages_sumUnstd.RData')
# refs alleles, social dosages btw 0 and 2
motch = match(unique(c(SGE_sig_QTLs[,'id'],DGE_sig_QTLs[,'id'])),colnames(social_dosage_sumUnstd))
any(is.na(motch))
#FALSE
social_dosage_sumUnstd = social_dosage_sumUnstd[,motch]
motch = match(rownames(cov_data),rownames(social_dosage_sumUnstd))
any(is.na(motch))
#FALSE
social_dosage_sumUnstd = social_dosage_sumUnstd[motch,]
colnames(social_dosage_sumUnstd) = paste(colnames(social_dosage_sumUnstd),'social',sep = '_')
means = apply(social_dosage_sumUnstd, FUN = mean, MAR = 2)
w = which(means >1)
if (length(w)!=0) {
	social_dosage_sumUnstd[,w] = 2 - social_dosage_sumUnstd[,w]
}
#now working with minor alleles. *2 to have number of minor alleles (additive model for SGE)
social_dosage_sumUnstd = social_dosage_sumUnstd*2
save(social_dosage_sumUnstd,file = 'output_dir/effect_sizes/SGE_sig_QTLs_social_genos.RData')

cov_data = cbind(cov_data,dosages,social_dosage_sumUnstd)

h5write(obj=cov_data,file=fid,name="covariates/matrix")
h5createGroup(fid,"covariates/col_header")
max(nchar(colnames(cov_data)))
#24
h5createDataset(file=fid,dataset="covariates/col_header/covariate_ID",dims=dim(cov_data)[2],storage.mode='character',size=30)
h5write(obj=colnames(cov_data),file=fid,name="covariates/col_header/covariate_ID")

#### GRM now
load('data_dir/dosages/pruned_dosages/GRM.RData')

h5createGroup(fid,"GRM")
h5createGroup(fid,"GRM/pruned_dosages")

h5createGroup(fid,"GRM/pruned_dosages/row_header")
max(nchar(rownames(GRM)))
#15
h5createDataset(file=fid,dataset="GRM/pruned_dosages/row_header/sample_ID",dims=dim(GRM)[1],storage.mode='character',size=20)
h5write(obj=rownames(GRM),file=fid,name="GRM/pruned_dosages/row_header/sample_ID")

h5write(obj=GRM,file=fid,name="GRM/pruned_dosages/matrix")

####subsets
load('data_dir/dosages/close_pairs_GRM.RData')
motch = match(all_cagemates_close_pairs, rownames(GRM))
any(is.na(motch))
#FALSE
include = rownames(GRM)[-motch]
length(include)
#[1] 1812

h5createGroup(file=fid,'subsets')
max(nchar(include))
#15
h5createDataset(file=fid,dataset="subsets/include",dims=length(include),storage.mode='character',size=20)
h5write(obj=include,file=fid,name="subsets/include")






