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
sge_single_snp_dir = '_conditioning'
root_dir = paste('output_dir/pvalues_pointwise_conditional_simulations/power_analysis/SGE/single_snp/pruned_dosages_include_DGE_IGE_IEE_cageEffect',sge_single_snp_dir,sep='')
sge_files = list.files(root_dir)
if (length(sge_files) == 0) next
sge_files = sge_files[grep('h5',sge_files)]
res_sge = mclapply(paste(root_dir,'/',sge_files,sep=''), my_f, mc.cores = 14)
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

sims_design = sims_design[sims_design[,'gs'] <= 3 & sims_design[,'dgeQTL_a'] == 0,]
for (i in seq(1,dim(sims_design)[1],by = 2)) {
	sims_design[i,'nb_detecteds'] = sims_design[i,'nb_detecteds'] + sims_design[i+1,'nb_detecteds']
	sims_design[i,'nb_sims'] = sims_design[i,'nb_sims'] + sims_design[i+1,'nb_sims']
}
sims_design = sims_design[seq(1,dim(sims_design)[1],by = 2),]
sims_design$detecteds = sims_design$nb_detecteds / sims_design$nb_sims
#1 2 3 correspond to gs 2 (equivalently DGE); 10 11 12 correspond to SGE gs3 additive model, 4 5 6 to SGE gs3 proportional model, 7 8 9 to SGE gs3 max model
sims_design = sims_design[c(1,2,3,10,11,12,4,5,6,7,8,9),]
#freq
#within generative model (add, prop, max)
#within gs
# figure 3a:
pdf(paste('plots_dir/power_analysis/paper_figure_logP_th',logP_th,'.pdf',sep=''), width = 8, height = 6)
par(mar = c(6,5,1,1))
barplot(sims_design[sims_design[,'sim_type'] != 'colMax', 'detecteds'], las = 2, xaxt="n", ylab = 'Power', ylim = c(0,1), col = 'black', cex.lab = 1.5, space = c(0,0.1,0.1,rep(c(0.8,0.1,0.1),2)))
axis(1,at = c(1.7,5.7,9.7), labels = c('DGE or\nSGE 1 cage mate','SGE additive\n2 cage mates','SGE proportional\n2 cage mates'), cex.axis = 1.5, tick = FALSE, padj = 1)
dev.off()

