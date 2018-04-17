effect = 'DGE'
pvalue_dir = paste('output_dir/pvalues_pointwise_conditional/',effect,'/pruned_dosages_include_DGE_IGE_IEE_cageEffect/',sep='')
load(paste(pvalue_dir,'QTLs_alpha0.01.RData',sep=''))
uniq_phenos_real = unique(all_QTLs[,'measure'])
length(uniq_phenos_real)
#170
pheno2run = read.delim('data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',as.is = T, header = F)[,1]
length(pheno2run)
#170
uniq_phenos_real = intersect(uniq_phenos_real,pheno2run)
length(uniq_phenos_real)
#170

perms_dir = paste('output_dir/pvalues_perms_pointwise_conditional/',effect,'/pruned_dosages_include_DGE_IGE_IEE_cageEffect/',sep='')
liste = list.files(perms_dir)
liste = liste[grep('_QTLs_alpha',liste)]
liste = sub("_QTLs_alpha.*.RData",'',liste,perl = T)
length(liste)
#170
liste = intersect(liste,pheno2run)
length(liste)
#170

phenos = intersect(uniq_phenos_real, liste)
length(phenos)
#170

chosen_FDR = 0.1

all_QTLs$FP = NA
all_QTLs$P = NA
all_QTLs$FDR = NA
all_nb_perms = c()
for (pheno in phenos) {
	load(paste(perms_dir,pheno,'_QTLs_alpha0.01.RData',sep=''))	

	nb_perms = length(unique(all_QTLs_perms[,'perm']))
	all_nb_perms = c(all_nb_perms,nb_perms)

	w = which(all_QTLs[,'measure'] == pheno)

	#all_QTLs does not need to be ordered
	for (i in w) {
		alpha = all_QTLs[i,'pvalue']
		P = sum(all_QTLs[w,'pvalue']<=as.numeric(alpha))
		FP = sum(all_QTLs_perms[,'pvalue']<=as.numeric(alpha)) / nb_perms
		all_QTLs[i,'FP'] = FP
		all_QTLs[i,'P'] = P
		all_QTLs[i,'FDR'] = FP / P
	}
}
table(all_nb_perms)
#SGE and DGE: all phenos with 100 perms

all_QTLs = all_QTLs[!is.na(all_QTLs[,'FDR']),]
all_QTLs = all_QTLs[order(all_QTLs[,'FDR']),]
all_QTLs = all_QTLs[all_QTLs[,'FDR']<=chosen_FDR,]
save(all_QTLs,file = paste(pvalue_dir,'QTLs_FDR',chosen_FDR,'.RData',sep=''))
all_QTLs[order(all_QTLs[,'chr'],all_QTLs[,'pos']),]

