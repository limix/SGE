library(parallel)
library(rhdf5)
library(gap)

my_f = function(file) {
	h5 = try(h5read(file, '/'))
	if (inherits(h5, 'try-error')) return(NULL)
	id_file = strsplit(file,'_')[[1]]
	combin = sub('le','',id_file[length(id_file)-2])
	id_file = paste(sub('chr','',id_file[length(id_file)-1]),sub('.h5','',id_file[length(id_file)]),sep='_')
	h5 = unlist(h5)
	id_h5 = paste(h5[1],h5[2],sep='_')
	return(c(h5,id_file,id_h5,combin))
}

for (cond in c('not_conditioning','conditioning')) {
	dir = paste('output_dir/pvalues_pointwise_conditional_simulations/polygenicNdgeQTL/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/',cond,'/',sep='')
	sge_files = list.files(dir)
	sge_files = sge_files[grep('h5',sge_files)]
	res_sge = mclapply(paste(dir,sge_files,sep=''), my_f, mc.cores = 14)
	lengths = lapply(res_sge,length)
	w = which(lengths == 0)
	if (length(w)!=0) {
		res_sge = res_sge[-w]
		sge_files = sge_files[-w]
	}
	res_sge = do.call('rbind',res_sge)
	all(res_sge[,4] == res_sge[,5])
	#TRUE so did test snp that was used to simulate dgeQTL
	res_sge = as.data.frame(res_sge,stringsAsFactors = F)
	for (i in c(1:3,6)) {
		res_sge[,i] = as.numeric(res_sge[,i])
	}
	print(dim(res_sge))
	assign(paste('res_sge',cond,sep='_'),res_sge)
}

pdf('plots_dir/polygenicNdgeQTL/fig_need4conditioning.pdf', width = 12)
par(mfrow = c(1,2), mar = c(5,5,1,1))
#6th row corresponds to var_Ad and var_As 20 and var_dgeQTL 20 (ie 16 var expl - high in real data = where problem is)
qqunif(res_sge_not_conditioning[res_sge_not_conditioning[,6] == 6,'pvalues'], col = 'black', pch = 16, ci = T, main = '', las = 1, cex.lab = 1.5)
qqunif(res_sge_conditioning[res_sge_conditioning[,6] == 6,'pvalues'], col = 'black', pch = 16, ci = T, main = '', las = 1, cex.lab = 1.5)
dev.off()

