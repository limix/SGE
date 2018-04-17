cd code_dir/GWAS/effect_sizes


#LIMIX only calculates betas for covariates so include these QTLs as "covariates" to null model with no (other) variant
#creates new HDF5 file specifically for this analysis
add_sig_QTLs_as_covariates_H5.R
# 1 - 72 for sig

qsub_get_effect_sizes.sh #calls exp_varianceDecomp.py

explore_betas_sig.R #gets effect sizes of SGE and DGE QTLs; run with type = 'allelic_effect' and type = 'var_expl'

gather_af_ve.R # brings together allelic effects and variances explained, as well as minor allele dosages and sample sizes

plot_af.R # plots Suppl Figure 6
plot_ve.R # plots Figure 3b-d
