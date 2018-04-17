cd code_dir/variance_decomposition

## important:
## in all exp_varianceDecomp.py scripts, always set random effects as desired and also the calc_ste argument to DirIndirVD 
##

#initial variance decomposition
qsub_varianceDecomp.sh to run exp_varianceDecomp_init.py (with DGE = "DGE", IGE = "IGE", IEE = "IEE", cageEffect = "cageEffect",  calc_ste=True)
mine_results_DGE_IGE.R # with path = "/homes/abaud/CFW/output/reproduce/pruned_dosages_include_DGE_IGE_IEE_cageEffect" and constrained = NULL
# shows that no collider effect due to having body weight as covariate for many phenotypes
code_dir/create_input/phenotypes/phenos_to_run.R

#now without IGE to calculate P values for aggregate contribution of SGE from LRT
qsub_varianceDecomp.sh to run exp_varianceDecomp.py (with DGE = "DGE", IGE = None, IEE = "IEE", cageEffect = "cageEffect", calc_ste=False)
mine_results_DGE_IGE.R # with path = "/homes/abaud/CFW/output/reproduce/pruned_dosages_include_DGE_IEE_cageEffect" and constrained = NULL
calc_VC_pvalues.R

test_corrAds1/qsub_varianceDecomp_constrained.sh
mine_results_DGE_IGE.R # with path = "/homes/abaud/CFW/output/reproduce/pruned_dosages_include_DGE_IGE_IEE_cageEffect" and constrained = 1
test_corrAds1/get_pvalues_corr1.R

write_VD.R # Table1.xls and Supplementary_Table1.xls
