#Supplementary Figure 2a-c

### load direct and social genotypes and exclude close pairs as done in the analysis
load('data_dir/dosages/pruned_dosages/unstd_social_dosages_sumUnstd.RData')
load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')
load('data_dir/dosages/close_pairs_GRM.RData')
motch = match(all_cagemates_close_pairs,rownames(social_dosage_sumUnstd))
any(is.na(motch))
#FALSE
social_dosage_sumUnstd = social_dosage_sumUnstd[-motch,]
motch = match(all_cagemates_close_pairs,rownames(dosages))
any(is.na(motch))
#FALSE
dosages = dosages[-motch,]
all(colnames(social_dosage_sumUnstd) == colnames(dosages))
# TRUE

### calculate correlation estimates and P values
corsP = rep(NA,dim(dosages)[2])
corsE = rep(NA,dim(dosages)[2])
for (i in 1:dim(dosages)[2]) {
	test = cor.test(social_dosage_sumUnstd[,i],dosages[,i], method = 'spearman')
	corsP[i] = test$p.value
	corsE[i] = test$estimate
}
names(corsP) = colnames(dosages)
names(corsE) = colnames(dosages)

corsP_perm = rep(NA,dim(dosages)[2])
corsE_perm = rep(NA,dim(dosages)[2])
for (i in 1:dim(dosages)[2]) {
	test = cor.test(dosages[,i],sample(social_dosage_sumUnstd[,i], replace = F), method = 'spearman')
	corsP_perm[i] = test$p.value
	corsE_perm[i] = test$estimate
}

save(dosages, social_dosage_sumUnstd, corsE, corsE_perm, corsP, corsP_perm, file = 'data_dir/dosages/pruned_dosages/corsE_SFigs.RData')

library(gap)
r = range(c(corsE^2,corsE_perm^2))
r 0.04

jpeg('plots_dir/genotypes_correlation/SupplFig2a.jpeg', height = 550)
hist(corsE^2, xlab = parse(text = 'R^2'), breaks = seq(0,0.04, by = 0.002), ylim = c(0,500), freq = F, col = rgb(1,0,0,0.5), cex.lab = 1.5, las = 1, cex.axis = 1.2, main = '')
hist(corsE_perm^2, col = rgb(0,0,0,0.3), breaks = seq(0,0.04, by = 0.002), add = T, freq = F)
legend(x = 'topright', legend = c('Unpermuted social genotypes','Permuted social genotypes'), fill = c(rgb(1,0,0,0.5), rgb(0,0,0,0.3)), bty = 'n', cex = 1.5)
dev.off()

jpeg('plots_dir/genotypes_correlation/SupplFig2b.jpeg')
qq2 = qqunif(corsP_perm, col = 'white')
qq1 = qqunif(corsP, col = 'white', ci = T, las = 1, cex.axis = 1.5, cex.lab = 1.5)
points(qq1$x, qq1$y, col = rgb(1,0,0,0.5))
points(qq2$x, qq2$y, col = rgb(0,0,0,0.3))
legend(x = 'topleft', legend = c('Unpermuted social genotypes','Permuted social genotypes'), fill = c(rgb(1,0,0,0.5), rgb(0,0,0,0.3)), bty = 'n', cex = 1.5)
dev.off()



#
MAFs_direct = apply(dosages,FUN = mean, MAR = 2)
range(MAFs_direct)
#0 to 1
w = which(MAFs_direct>0.5)
MAFs_direct[w] = 1 - MAFs_direct[w]

MAFs_social = apply(social_dosage_sumUnstd,FUN = mean, MAR = 2)
range(MAFs_social)
#0 to 2
MAFs_social = MAFs_social/2
w = which(MAFs_social>0.5)
MAFs_social[w] = 1 - MAFs_social[w]

jpeg('plots_dir/genotypes_correlation/SupplFig2c.jpeg')
par(mar = c(5,5,1,1))
smoothScatter(MAFs_direct, -log10(corsP), xlab = 'MAF', ylab = '-logP', las = 1, cex = 3, cex.axis = 1, cex.lab = 1.5)
#smoothScatter(MAFs_direct, corsE^2, xlab = 'MAF', ylab = parse(text = 'R^2'), las = 1, cex = 3, cex.axis = 1, cex.lab = 1.5)
dev.off()


# Supplementary Figure 3e

all_pos = do.call('rbind',sapply(names(corsP),strsplit,'_'))
all_pos = as.data.frame(all_pos,stringsAsFactors=F)
all_pos[,1] = as.numeric(sub('chr','',all_pos[,1]))
all_pos[,2] = as.numeric(all_pos[,2])
colnames(all_pos) = c('chr','pos')

#load('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_alpha0.01.RData')
load('output_dir/pvalues_pointwise_conditional/SGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData')
SGE_QTLs = all_QTLs
load('output_dir/pvalues_pointwise_conditional/DGE/pruned_dosages_include_DGE_IGE_IEE_cageEffect/QTLs_FDR0.1.RData')
DGE_QTLs = all_QTLs

all_pos$in_FDR_SGE_QTL = FALSE
for (i in 1:dim(SGE_QTLs)[1]) {
	w = which(all_pos[,'chr'] == SGE_QTLs[i,'chr'] & all_pos[,'pos'] >= SGE_QTLs[i,'ci_starts'] & all_pos[,'pos'] <= SGE_QTLs[i,'ci_stops'])
	all_pos[w,'in_FDR_SGE_QTL'] = NA
	w = which(all_pos[,'chr'] == SGE_QTLs[i,'chr'] & all_pos[,'pos'] >= (SGE_QTLs[i,'pos']-100000) & all_pos[,'pos'] <= (SGE_QTLs[i,'pos']+100000))
	all_pos[w,'in_FDR_SGE_QTL'] = TRUE
}
all_pos$in_FDR_DGE_QTL = FALSE
for (i in 1:dim(DGE_QTLs)[1]) {
	w = which(all_pos[,'chr'] == DGE_QTLs[i,'chr'] & all_pos[,'pos'] >= DGE_QTLs[i,'ci_starts'] & all_pos[,'pos'] <= DGE_QTLs[i,'ci_stops'])
	all_pos[w,'in_FDR_DGE_QTL'] = NA
	w = which(all_pos[,'chr'] == DGE_QTLs[i,'chr'] & all_pos[,'pos'] >= (DGE_QTLs[i,'pos']-100000) & all_pos[,'pos'] <= (DGE_QTLs[i,'pos']+100000))
	all_pos[w,'in_FDR_DGE_QTL'] = TRUE
}

library(gap)
jpeg('plots_dir/genotypes_correlation/SupplFig3e.jpeg')
#par(mfrow = c(1,2))
qq2 = qqunif(corsP[which(all_pos$in_FDR_SGE_QTL | all_pos$in_FDR_DGE_QTL)])
qq1 = qqunif(corsP[which(!all_pos$in_FDR_SGE_QTL & !all_pos$in_FDR_DGE_QTL)], col = 'white', ci =T,las = 1, cex.lab = 1.5, cex.axis = 1.5)
#qqunif(corsP[which(!all_pos$in_FDR_SGE_QTL & !all_pos$in_FDR_DGE_QTL)], col = rgb(0,0,0,0.3), ci =T, main = 'Variants outside QTLs', ylim = c(0,6), xlim = c(0,4))
#qqunif(corsP[which(all_pos$in_FDR_SGE_QTL | all_pos$in_FDR_DGE_QTL)], col = rgb(1,0,1,0.5), ci =T, main = 'Variants at sgeQTL or dgeQTL', ylim = c(0,6), xlim = c(0,4))
#plot(qq1$x, qq1$y, xlab = '-log10(Expected', ylab = '-log10(Observed', col = rgb(0,0,0,0.3), las = 1, cex.lab = 1.5, cex.axis = 1.5)
points(qq1$x, qq1$y, col = 'black')
points(qq2$x, qq2$y, col = rgb(1,0,1,0.5))
legend(x = 'topleft',legend = c('Variants outside SGE and DGE loci','Variants at SGE or DGE locus'), fill = c(rgb(0,0,0,0.3), rgb(1,0,1,0.5)), cex = 1.5, bty = 'n', border = 'white')
#abline(0,1)
dev.off()






