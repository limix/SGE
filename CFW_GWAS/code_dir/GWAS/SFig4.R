load(paste('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData',sep=''))
SGE_QTLs = all_QTLs

library(rhdf5)
pvalues_dir = 'output_dir/pvalues_pointwise_conditional/DGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect'
files = list.files(pvalues_dir)
files = files[grep('h5',files)]
length(files)
pheno2run = read.delim('data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt',as.is = T, header = F)[,1]
inter = intersect(paste(pheno2run,'.h5',sep=''),files)
length(inter)
#170
my_get = function(nome) {
	ret = my_h5[[nome]]
	names(ret) = paste(nome,1:length(ret),sep='_')
	return(ret)
}

pvalues_vars_in_sgeQTLs_same_pheno = c()
pvalues_vars_not_in_sgeQTLs_same_pheno = c()
for (file in inter) {
	pheno_name = sub('.h5','',file,fixed=T)
	print(pheno_name)

	my_h5 = try(h5read(paste(pvalues_dir,'/',pheno_name,'.h5',sep=''),'/'))
	if (inherits(my_h5,'try-error')) next
	if (length(my_h5)!=4 || all(is.na(my_h5$pvalues))) next
	my_h5 = as.data.frame(my_h5)

	g = grep(pheno_name,SGE_QTLs[,'measure'])
	if (length(g) != 0) {
		for (k in g) {
			is_in = rep(FALSE, dim(my_h5)[1])
			is_in[which(my_h5$chr == SGE_QTLs[k,'chr'] & my_h5$pos >= SGE_QTLs[k,'ci_starts'] & my_h5$pos <= SGE_QTLs[k,'ci_stops'])] = NA
			is_in[which(my_h5$chr == SGE_QTLs[k,'chr'] & abs(my_h5$pos -SGE_QTLs[k,'pos'])<100000)] = TRUE
			pvalues_vars_in_sgeQTLs_same_pheno = c(pvalues_vars_in_sgeQTLs_same_pheno,my_h5[which(is_in),'pvalues'])
			pvalues_vars_not_in_sgeQTLs_same_pheno = c(pvalues_vars_not_in_sgeQTLs_same_pheno,sample(my_h5[which(!is_in),'pvalues'],size = length(which(is_in))))
		}
	}
}
length(pvalues_vars_not_in_sgeQTLs_same_pheno)
length(pvalues_vars_in_sgeQTLs_same_pheno)
#both 1895

library(gap)
jpeg('plots_dir/Suppl_Figure4.jpeg')
qqunif(pvalues_vars_in_sgeQTLs_same_pheno, ci = T, col = 'black', main = 'dgeGWAS P values, variants at SGE loci for the same phenotype', las = 1, cex.lab = 1.5)
#qqunif(pvalues_vars_not_in_sgeQTLs_same_pheno, ci = T, col = 'black',  main = 'dgeGWAS P values, variants not at sgeQTLs for the same phenotype', las = 1, cex.lab = 1.5)
dev.off()
