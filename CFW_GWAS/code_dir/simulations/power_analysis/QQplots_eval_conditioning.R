library(parallel)
library(rhdf5)
library(gap)

my_f = function(file) {
	h5 = try(h5read(file, '/'))
	if (inherits(h5, 'try-error')) return(NULL)
	id_file = strsplit(file,'_')[[1]]
	combin = sub('sd','',id_file[length(id_file)-2])
	id_file = paste(id_file[length(id_file)-1],sub('.h5','',id_file[length(id_file)]),sep='_')
	h5 = unlist(h5)
	id_h5 = sub('_gs.*','',h5[1])
	return(c(h5,id_file,id_h5,combin))
}

logP_th = 5

root_dir = 'output_dir/pvalues_pointwise_conditional_simulations/power_analysis/SGE/single_snp/pruned_dosages_include_DGE_IGE_IEE_cageEffect'
for (sge_single_snp_dir in c('_conditioning','_not_conditioning')) {
	print(sge_single_snp_dir)
	sge_files = list.files(paste(root_dir,sge_single_snp_dir, sep=''))
	if (length(sge_files) == 0) next
	sge_files = sge_files[grep('h5',sge_files)]
	res_sge = mclapply(paste(root_dir,sge_single_snp_dir,'/',sge_files,sep=''), my_f, mc.cores = 14)
	lengths = lapply(res_sge,length)
	w = which(lengths == 0)
	if (length(w)!=0) {
		res_sge = res_sge[-w]
		sge_files = sge_files[-w]
	}
	res_sge = do.call('rbind',res_sge)
	all(res_sge[,3] == res_sge[,4])
	#TRUE so did test snp that was used to simulate dgeQTL
	res_sge = as.data.frame(res_sge,stringsAsFactors = F)
	colnames(res_sge) = c('social_geno_SNP_id','pvalues','id1','id2','row_sims_design')
	for (i in c(2,5)) {
		res_sge[,i] = as.numeric(res_sge[,i])
	}
	res_sge$logP = -log10(res_sge[,'pvalues'])
	dim(res_sge)
	#5076 5
	print(sort(table(res_sge[,'row_sims_design'])))

	load('data_dir/simulations/power_analysis/sims_design.RData')
	sims_design = sims_design[1:60,]
	sims_design = rbind(sims_design,sims_design[sims_design[,'gs']==3,],sims_design[sims_design[,'gs']==3,])
	sims_design$detecteds = NA
	sims_design$nb_detecteds = NA
	sims_design$nb_sims = NA
	sims_design$sim_type = 'prop'
	for (row_sims_design in 1:dim(sims_design)[1]) {
		w = which(res_sge[,'row_sims_design'] == row_sims_design)
		sims_design[row_sims_design,'nb_sims'] = length(w)
		if (length(w) == 0) next
		sims_design[row_sims_design,'detecteds'] = sum(res_sge[w,'logP']>logP_th) / length(w)
		sims_design[row_sims_design,'nb_detecteds'] = sum(res_sge[w,'logP']>logP_th)
		if (61 <= row_sims_design & row_sims_design <=72) sims_design[row_sims_design,'sim_type'] = 'colMax'
		if (73 <= row_sims_design & row_sims_design <=84) sims_design[row_sims_design,'sim_type'] = 'add'
	}

	sims_design$gs = factor(sims_design$gs,c("2","3","4","5","6"))
	sims_design$MAF = factor(sims_design$MAF,c("low","mid","high"))
	sims_design$corsP = factor(sims_design$corsP,c("low","high"))

	assign(paste('res_sge',sge_single_snp_dir,sep=''),res_sge)
	assign(paste('sims_design',sge_single_snp_dir,sep=''),sims_design)
}
w1 = which(sims_design_not_conditioning[,'sim_type'] == 'add' & sims_design_not_conditioning[,'gs'] == 3 & sims_design_not_conditioning[,'sgeQTL_a'] != 0 & sims_design_not_conditioning[,'dgeQTL_a'] == 0 & sims_design_conditioning[,'corsP'] == 'low')
w2 = which(sims_design_not_conditioning[,'sim_type'] == 'add' & sims_design_not_conditioning[,'gs'] == 3 & sims_design_not_conditioning[,'sgeQTL_a'] != 0 & sims_design_not_conditioning[,'dgeQTL_a'] == 0 & sims_design_not_conditioning[,'corsP'] == 'high')
w3 = which(sims_design_not_conditioning[,'sim_type'] == 'add' & sims_design_not_conditioning[,'gs'] == 3 & sims_design_not_conditioning[,'sgeQTL_a'] != 0 & sims_design_not_conditioning[,'dgeQTL_a'] != 0 & sims_design_conditioning[,'corsP'] == 'low')
w4 = which(sims_design_not_conditioning[,'sim_type'] == 'add' & sims_design_not_conditioning[,'gs'] == 3 & sims_design_not_conditioning[,'sgeQTL_a'] != 0 & sims_design_not_conditioning[,'dgeQTL_a'] != 0 & sims_design_not_conditioning[,'corsP'] == 'high')

summary = c(mean(sims_design_not_conditioning[w1,'detecteds']),mean(sims_design_conditioning[w1,'detecteds']))
summary = cbind(summary, c(mean(sims_design_not_conditioning[w2,'detecteds']),mean(sims_design_conditioning[w2,'detecteds'])))
summary = cbind(summary, c(mean(sims_design_not_conditioning[w3,'detecteds']),mean(sims_design_conditioning[w3,'detecteds'])))
summary = cbind(summary, c(mean(sims_design_not_conditioning[w4,'detecteds']),mean(sims_design_conditioning[w4,'detecteds'])))
rownames(summary) = c('not_conditioning','conditioning')
colnames(summary) = c('noDGE_lowCorsP','noDGE_highCorsP','DGE_lowCorsP','DGE_highCorsP')

#Suppl figure 3a-d (a-b for additive model, c-d for proportional model)
pdf(paste('plots_dir/power_analysis/QQplot_sims_power_analysis_eval_conditioning_logP_th',logP_th,'_add_model.pdf',sep=''), width = 12, height = 8)
#pdf(paste('plots_dir/power_analysis/QQplot_sims_power_analysis_eval_conditioning_logP_th',logP_th,'_model.pdf',sep=''), width = 12, height = 8)
par(mfrow = c(1,2), mar = c(5,5,1,1))
barplot(summary[,1:2], beside = T, las = 2, xaxt="n", ylab = 'Power', ylim = c(0,0.85), cex.lab = 1.5)
axis(1,at = c(2,5), labels = c('No local DGE,\nlow genotypic corr','No local DGE,\nhigh genotypic corr'), cex.axis = 1.5, tick = FALSE, padj = 1)
legend(x = 'top',legend = c('Not conditioning', 'Conditioning'), fill = c('black','grey'), bty = 'n', cex = 1.5)

barplot(summary[,3:4], beside = T, las = 2, xaxt="n", ylab = 'Power', ylim = c(0,0.85), cex.lab = 1.5)
axis(1,at = c(2,5), labels = c('Local DGE,\nlow genotypic corr','Local DGE,\nhigh genotypic corr'), cex.axis = 1.5, tick = FALSE, padj = 1)
legend(x = 'top',legend = c('Not conditioning', 'Conditioning'), fill = c('black','grey'), bty = 'n', cex = 1.5)
dev.off()


