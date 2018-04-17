load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/corr1_pvalues.RData')
mean(corr1[,'corr DGE SGE'], na.rm = T)
load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')
load('output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect/significance_SGE_VC.RData')
load('data_dir/phenotypes/all_pheno_data.RData')

motch = match(rownames(all_VCs),names(pv_chi2dof2))
any(is.na(motch))
#FALSE 
pv_chi2dof2 = pv_chi2dof2[motch]

motch = match(rownames(all_VCs),rownames(corr1))
any(is.na(motch))
#FALSE 
corr1 = corr1[motch,]

motch = match(rownames(all_VCs),names(qv_chi2dof2))
any(is.na(motch))
#FALSE 
qv_chi2dof2 = qv_chi2dof2[motch]

motch = match(rownames(all_VCs),my_model_menu[,'Name'])
any(is.na(motch))
#FALSE 
my_model_menu = my_model_menu[motch,]

model_menu = read.delim('data_dir/phenotypes/TableS1_phenotypes.txt',as.is=T,check.name = F)
model_menu = model_menu[1:200,]
motch = match(sub('_BWcorr','',rownames(corr1)), model_menu[,'Name'])
any(is.na(motch))
#TRUE - osteoporosis manually changed later
model_menu = model_menu[motch,]

save(pv_chi2dof2, all_VCs, qv_chi2dof2, my_model_menu, model_menu, corr1, file = 'output_dir/VD/pruned_dosages_include_DGE_IEE_cageEffect/Table1.RData')

toWrite=cbind(model_menu[,'Measure'],names(pv_chi2dof2), my_model_menu[,'my_covs'],all_VCs[,'sample_size'],paste(round(all_VCs[,'var_As']*100,digits = 0),'+/-',round(all_VCs[,'STE_As']*100,digits = 0),sep=''),paste(round(all_VCs[,'var_Ad']*100,digits=0),'+/-',round(all_VCs[,'STE_Ad']*100,digits =0),sep=''),signif(pv_chi2dof2,digits = 2),signif(qv_chi2dof2, digits = 2),paste(round(all_VCs[,'corr_Ads'],digits = 2),'+/-',round(all_VCs[,'STE_corr_Ads'],digits = 2)), signif(corr1[,'corr P value'],digits = 2), signif(corr1[,'corr Q value'],digits = 2))
colnames(toWrite)=c('Phenotype','Name','Covariates','Sample size','Aggregate contribution SGE (%)','Aggregate contribution DGE (%)','P value aggregate SGE','Q value aggregate SGE','Correlation  DGE SGE', 'P value rho != 1', 'Q value rho != 1')

ordered_toWrite = toWrite[order(toWrite[,'Name']),]
write.table(ordered_toWrite,file = 'plots_dir/Supplementary_Table1.xls',sep='\t',quote = F,col.names = T,row.names = F)

cols2exclude = c('Covariates','Sample size','Q value aggregate SGE', 'Q value rho != 1')
keep = which(!is.na(corr1[,'corr DGE SGE']) | pv_chi2dof2<0.05)
toWrite = toWrite[keep,-which(colnames(toWrite) %in% cols2exclude)]
pv_chi2dof2 = pv_chi2dof2[keep]
toWrite[is.na(toWrite[,'P value rho != 1']),'Correlation  DGE SGE'] = NA
ordered_toWrite = toWrite[order(pv_chi2dof2),]
write.table(ordered_toWrite,file = 'plots_dir/Table1.xls',sep='\t',quote = F,col.names = T,row.names = F)


