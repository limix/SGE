load('output_dir/VD/pruned_dosages_include_DGE_IGE_IEE_cageEffect/w_STE/all_VCs.RData')
load('output_dir/effect_sizes_additive_model/betas_conditioning/var_expl/betas.RData')

r = range(c(indiv_DGE_QTLs_betas[,'sig_beta1'],indiv_SGE_QTLs_betas[,'sig_beta1']))
pdf('plots_dir/effect_sizes/Fig3b.pdf', width = 7)
par(mar = c(5,5,1,1))
hist(indiv_DGE_QTLs_betas[,'sig_beta1'], main = '', col = 'black', xlim = c(0,40), breaks = seq(0,40, by = 2), las = 1, xlab = 'Proportion of phenotypic variance explained (%)', ylab = 'Count', cex.lab = 1.5, cex.axis = 1.5)
hist(indiv_SGE_QTLs_betas[,'sig_beta1'], col = rgb(1,0,0,0.8), xlim = c(0,40), breaks = seq(0,40, by = 2), cex.lab = 1.5, cex.axis = 1.5, add = T)
legend(x = 'top', legend = c('Genome-wide significant SGE','Genome-wide significant DGE'), fill = c('red','black'), bty = 'n', border = 'white', cex = 1.5)
dev.off()


pdf('plots_dir/effect_sizes/Fig3cd.pdf', width = 8)
par(mar = c(5,8,1,1))
plot(all_res_DGE[,'herit_DGE']*100, all_res_DGE[,'var_expl_DGE_QTLs']*100, col = 'black', xlim = c(0,50), ylim = c(0,50), pch = 16, las = 1, xlab = 'Variance explained in aggregate', ylab = 'Variance explained by\n all genome-wide significant associations\nfor a given phenotype', cex.lab = 1.5, cex.axis = 1.5)
abline(0,1,col = 'darkgrey')
points(all_res_SGE[,'herit_SGE']*100, all_res_SGE[,'var_expl_SGE_QTLs']*100, pch = 16, col = 'red')
legend(x = 'topleft', legend = c('SGE','DGE'), fill = c('red','black'), bty = 'n', border = 'white', cex = 1.5)
abline(v = 5, col = 'darkgrey', lty = 5)

props_DGE = all_res_DGE[,'var_expl_DGE_QTLs']*100 / all_res_DGE[,'herit_DGE']
props_DGE = props_DGE[all_res_DGE[,'herit_DGE']>0.05]
props_SGE = all_res_SGE[,'var_expl_SGE_QTLs']*100 / all_res_SGE[,'herit_SGE']
props_SGE = props_SGE[all_res_SGE[,'herit_SGE']>0.05]

hist(props_DGE, main = '', col = 'black', ylim = c(0,20), xlim = c(0,100), breaks = seq(0,100, by = 10),  freq = T, axes = F, xlab = '', ylab = '')
hist(props_SGE, col = rgb(1,0,0,0.8), breaks = seq(0,100, by = 10), add = T, freq = T)
axis(1,las = 1, cex.lab = 1.5, cex.axis = 1.5)
axis(2,las = 1, cex.lab = 1.5, cex.axis = 1)
mtext(side = 1, text = 'Proportion of genetic variance explained (%)', padj = 3, cex = 1.5)
mtext(side = 2, text = 'Count', padj = -4, cex = 1.5)
legend(x = 'topright', legend = c('SGE','DGE'), fill = c('red','black'), bty = 'n', border = 'white', cex = 1.5)
dev.off()
