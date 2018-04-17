load('data_dir/simulations/power_analysis/fake_cage_cm.RData')
load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')

motch = match(rownames(fake_cage_cm),rownames(dosages))
any(is.na(motch))
#[1] FALSE
dosages = dosages[motch,]
dosages = dosages*2
#now varies between 0 and 2, ie ref allele dosage

mean_dosages = apply(dosages, FUN = mean, MAR = 2)
names(mean_dosages) = colnames(dosages)
mean_dosages = mean_dosages/2
w = which(mean_dosages>0.5)
mean_dosages[w] = 1 - mean_dosages[w]
MAFs = mean_dosages
low_MAFs = names(MAFs)[MAFs<0.05]
length(low_MAFs)
mid_MAFs = names(MAFs)[MAFs>0.225 & MAFs<0.275]
length(mid_MAFs)
high_MAFs = names(MAFs)[MAFs>0.45]
length(high_MAFs)
#resp number of SNPs: 67,246, 31,335, and 19,623
save(MAFs, low_MAFs, mid_MAFs, high_MAFs,file= 'data_dir/simulations/power_analysis/MAFs.RData')

std_dosages = scale(dosages,center = T, scale = T)

###

#run script for gs in 2,3,4,6
gs = 3
library(parallel)
social_dosage_sumUnstd = matrix (NA,ncol = dim(dosages)[2],nrow = dim(dosages)[1]) 
rownames(social_dosage_sumUnstd) = rownames(dosages)
colnames(social_dosage_sumUnstd) = paste(colnames(dosages),'_gs',gs,sep='')
for (mouse in rownames(dosages)) {
	first_cage = fake_cage_cm[mouse,paste('cage_gs',gs,sep='')]
	cagemates = rownames(fake_cage_cm)[which(fake_cage_cm[,paste('cage_gs',gs,sep='')]==first_cage)]
	cagemates = cagemates[cagemates!=mouse]

	if (length(cagemates) != (gs-1)) stop('pb')

	# create 2 social_dosage_sumUnstd: one for additive/proportional models, and one for maximum model (not in paper)
	# add/prop:
	social_dosage_sumUnstd[mouse,] = colSums(dosages[cagemates,,drop=F])
	#max
	#social_dosage_sumUnstd[mouse,] = apply(dosages[cagemates,,drop=F], FUN = max, MAR = 2)
}
std_social_dosage_sumUnstd = scale(social_dosage_sumUnstd,center = T, scale = T)

#calculate genotypic correlations now
corsP = rep(NA,dim(dosages)[2])
corsE = rep(NA,dim(dosages)[2])
for (i in 1:dim(dosages)[2]) {
	test = cor.test(social_dosage_sumUnstd[,i],dosages[,i], method = 'spearman')
	corsP[i] = -log10(test$p.value)
	corsE[i] = test$estimate
}
names(corsP) = colnames(dosages)
names(corsE) = colnames(dosages)

high_corsP = names(corsP)[corsP > 2]
low_corsP = names(corsP)[corsP < 0.05]

range(corsE[corsP > 2])
range(corsE[corsP < 0.05])

save(corsP, corsE, low_corsP, high_corsP, file= paste('data_dir/simulations/power_analysis/genotypic_correlations_gs',gs,'.RData', sep=''))
	
#now will save reduced version of direct and social dosages with 1,800 mice and only SNPs that are in
SNPs_to_save = c(intersect(low_MAFs,low_corsP),intersect(low_MAFs,high_corsP),intersect(mid_MAFs,low_corsP),intersect(mid_MAFs,high_corsP),intersect(high_MAFs,low_corsP),intersect(high_MAFs,high_corsP))
length(SNPs_to_save)

dosages = dosages[,SNPs_to_save]
social_dosage_sumUnstd = social_dosage_sumUnstd[,paste(SNPs_to_save,'_gs',gs,sep='')]
assign(paste('std_dosages_gs',gs,sep=''), std_dosages[,SNPs_to_save])
assign(paste('std_social_dosage_sumUnstd_gs',gs,sep=''), std_social_dosage_sumUnstd[,paste(SNPs_to_save,'_gs',gs,sep='')])
#_colMax manually added if that's what run
save(list = c(paste('std_dosages_gs',gs,sep=''),paste('std_social_dosage_sumUnstd_gs',gs,sep='')), file= paste('data_dir/simulations/power_analysis/dosages/dir_soc_std_dosages_gs',gs,'.RData', sep=''))
save(dosages,social_dosage_sumUnstd, file= paste('data_dir/simulations/power_analysis/dosages/dir_soc_dosages_gs',gs,'.RData', sep=''))

