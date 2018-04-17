
load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')
all_VCs_full = all_VCs
load('output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect/all_VCs.RData')
all_VCs_null = all_VCs
    
inter=intersect(rownames(all_VCs_null),rownames(all_VCs_full))
all_VCs_full=all_VCs_full[match(inter,rownames(all_VCs_full)),]
all_VCs_null=all_VCs_null[match(inter,rownames(all_VCs_null)),]
    
pv_chi2dof2 = 1 - pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=2)
names(pv_chi2dof2) = inter

library(gap)
pdf('plots_dir/QQplot_pvalues_SGE_VC.pdf',bg='white')
qqunif(pv_chi2dof2,ci=T,col='black',las=1)
dev.off()

#Q value    
library(qvalue)
qv_chi2dof2 = qvalue(pv_chi2dof2)$qvalue

save(pv_chi2dof2,qv_chi2dof2, file = 'output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect/significance_SGE_VC.RData')




load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')
all_VCs_full = all_VCs
load('output_dir/VD/pruned_dosages_include_IGE_IEE_cageEffect/all_VCs.RData')
all_VCs_null = all_VCs
    
inter=intersect(rownames(all_VCs_null),rownames(all_VCs_full))
all_VCs_full=all_VCs_full[match(inter,rownames(all_VCs_full)),]
all_VCs_null=all_VCs_null[match(inter,rownames(all_VCs_null)),]
    
pv_chi2dof2 = 1 - pchisq(2*(-all_VCs_full$LML+all_VCs_null$LML),df=2)
names(pv_chi2dof2) = inter

mean(all_VCs_full[pv_chi2dof2<0.05,'var_Ad'])
#15%



