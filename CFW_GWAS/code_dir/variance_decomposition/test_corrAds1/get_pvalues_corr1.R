load('output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect/significance_SGE_VC.RData')

load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/constrained1/all_VCs_constrained1.RData')

unique_phenos = unique(all_VCs[,'trait'])
corr1 = NULL
for (pheno in unique_phenos) {
	ff = all_VCs[all_VCs[,'trait']==pheno & abs(all_VCs[,'corr_Ads'])!=1,,drop=F]
	lr = all_VCs[(all_VCs[,'trait']==pheno)& abs(all_VCs[,'corr_Ads'])==1,,drop=F]
	pv_chi2dof2 = 1 - pchisq(2*(-ff$LML+lr$LML),df=1)
	corr1 = rbind(corr1,c(unlist(ff[,'var_Ad']),unlist(ff[,'var_As']),unlist(ff[,'corr_Ads']),pv_chi2dof2))
}
corr1 = as.data.frame(corr1)
rownames(corr1) = unique_phenos
colnames(corr1) = c('var_Ad','var_As','corr DGE SGE','corr P value')
corr1[,'corr DGE SGE'] = round(corr1[,'corr DGE SGE'],digits = 2)

library(gap)
pdf('plots_dir/QQplot_pvalues_corrsAds1.pdf',bg='white')
qqunif(corr1[,'corr P value'],ci=T,pch = 16,col='black',las=1)
dev.off()


corr_threshold = 0.05
reliable_corrs = which(corr1[,'var_As']>corr_threshold & corr1[,'var_Ad']>corr_threshold)

#Q value    
library(qvalue)
corr1[,'corr Q value'] = NA
corr1[reliable_corrs,'corr Q value'] = qvalue(corr1[reliable_corrs,'corr P value'], lambda = 0)$qvalue
corr1[is.na(corr1[,'corr Q value']),'corr DGE SGE'] = NA
corr1[is.na(corr1[,'corr Q value']),'corr P value'] = NA
save(corr1,corr_threshold,file = 'output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/corr1_pvalues.RData')



