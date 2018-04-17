load('output_dir/effect_sizes_additive_model/betas_conditioning/all_betas.RData')
#betas are additive betas^2

#colour for var_expl
N=4
hot=heat.colors(N)
hot=colorRampPalette(c('red','yellow'))(N)
hot=hot[seq(length(hot),1,by=(-1))]
plotted = indiv_DGE_QTLs_betas[,'var_expl']
r=range(plotted)
r = c(0, 3)
cutted=cut(plotted,breaks=seq(r[1]-0.000000001,r[2]+0.000000001,length.out=N),labels=F)
cols_DGE=hot[cutted]
cols_DGE[plotted>3] = 'red4'

hot=heat.colors(N)
hot=colorRampPalette(c('red','yellow'))(N)
hot=hot[seq(length(hot),1,by=(-1))]
plotted = indiv_SGE_QTLs_betas[,'var_expl']
r=range(plotted)
r = c(0, 3)
cutted=cut(plotted,breaks=seq(r[1]-0.000000001,r[2]+0.000000001,length.out=N),labels=F)
cols_SGE_add=hot[cutted]


#cex for sample size
N=6
r = c(800, 1700)
plotted = indiv_DGE_QTLs_betas[,'sample_size']
cutted=cut(plotted,breaks=seq(r[1],r[2],length.out=N),labels=F)
cexes_DGE=seq(0.6,1.4,length = N)[cutted]
plotted = indiv_SGE_QTLs_betas[,'sample_size']
cutted=cut(plotted,breaks=seq(r[1],r[2],length.out=N),labels=F)
cexes_SGE=seq(0.6,1.4,length = N)[cutted]

load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages.RData')
MAF = apply(dosages,FUN = mean, MAR = 2)
range(MAF)
#0 to 1
w = which(MAF>0.5)
MAF[w] = 1 - MAF[w]
#0 to 0.5
dens_MAF_all = density(MAF, from = 0, to = 0.5)

#MAD is always number of minor alleles (additive SGE model)!
#max(indiv_SGE_QTLs_betas[,'MAD'])
#[1] 1.50781 when it could vary btw 0 and 2
#max(indiv_DGE_QTLs_betas[,'MAD'])
#[1] 0.9971054 when it could vary btw 0 and 1

#for Suppl Figure 6 figure, want to plot genotype variance for DGE and SGE under two models
#the genotype variance is 2p(1-p) for DGE, 2Np(1-p) for SGE under the additive model, and 2Np(1-p)/N2 for SGE under the proportional model
motch_SGE = match(sub('_social','',indiv_SGE_QTLs_betas[,2]),names(MAF))
MAF_SGE_QTLs = MAF[motch_SGE]
SGE_QTLs_vars_add = 2 * 2 * MAF_SGE_QTLs * (1-MAF_SGE_QTLs)
SGE_QTLs_vars_prop = MAF_SGE_QTLs * (1-MAF_SGE_QTLs)

dens_SGE_QTLs_vars_add = density(SGE_QTLs_vars_add, from = 0, to = 1)
dens_SGE_QTLs_vars_prop = density(SGE_QTLs_vars_prop, from = 0, to = 0.25)

#allelic effect (betas squared) for additive SGE model
r = range(indiv_SGE_QTLs_betas[,'add_allelic_effect'])
dens_SGE_QTLs_all_eff = density(indiv_SGE_QTLs_betas[,'add_allelic_effect'], from = r[1], to = r[2])

r_SGE_QTLs_vars = range(c(dens_SGE_QTLs_vars_add$y,dens_SGE_QTLs_vars_prop$y))
pdf('plots_dir/effect_sizes/paper_marginals.pdf', width = 13)
par(mfrow = c(1,2))
plot(dens_SGE_QTLs_vars_add$x, dens_SGE_QTLs_vars_add$y, type = 'l', xlab = 'Social genotype variance', ylab = 'Density', las = 1, col = 'red', xlim = c(0,1), ylim = r_SGE_QTLs_vars)
lines(dens_SGE_QTLs_vars_prop$x, dens_SGE_QTLs_vars_prop$y, col = 'orange')
legend(x = 'topright', legend = c('Additive model','Proportional model'), fill = c('red','orange'), bty = 'n', border = 'white') 

plot(dens_SGE_QTLs_all_eff$x, dens_SGE_QTLs_all_eff$y, type = 'l', xlab = 'Social allelic effect', ylab = 'Density', las = 1, col = 'red')

dev.off()
