#Simulations 1: "polygenic_only"; for Supplementary Figure 1
cd code_dir/simulations/polygenic_only
simulate.R #builds covariance matrices and list of effect sizes (based on estimated averages)
qsub_draw_sims.sh #calls draw_sims.R; draws 500 simulations from mvnorm; saves one RData file per combination of effect sizes (67 combinations total) with simulations, targetted_target_effects and target_effects; simulations is #mice x 500
create_h5.R #writes simulations in HDF5 files; samples columns so that different combinations get done at the same time
qsub_varianceDecomp.sh #calls exp_varianceDecomp.py
mine_results_simulations.R #builds all_VCs.RData
compare_results_sims_targets.R # boxplots for Supplementary Figure 1


#Simulations 2: "polygenicNdgeQTL"; for Supplementary Figure 2d and 2e
cd code_dir/simulations/polygenicNdgeQTL
simulate.R #builds covariance matrices and list of effect sizes (high and low polyegenic effects; 4 QTL effect sizes so 8 combinations total)
qsub_draw_sims.sh #calls draw_sims.R; draws 100 simulations from mvnorm; saves one RData file per combination of effect sizes (8 combinations total) with simulations, targetted_target_effects and target_effects; simulations is #mice x 100
create_h5.R #writes simulations in HDF5 files; samples columns so that different combinations get done at the same time
qsub_varianceDecomp.sh #calls exp_varianceDecomp.py
qsub_map_SGE.sh #to launch map_SGE_single_snp.py
need4conditioning.R #QQ plots for Supplementary Figure 2d and 2e (correspond to a single set of parameters)


#Simulations 3: "power_analysis"; for Figure 3a and Supplementary Figure 3a-d 
cd code_dir/simulations/power_analysis
build_cages.R #randomly assigns mice to cages of 2 to 6 mice
build_social_genotypes.R # builds social genotypes based on the random cage assignments; lists variants with low (MAF<0.05) mid (0.225<MAF<0.275) or high (MAF>0.45) MAF and low (-logP < 0.05) or high (-logP > 2) correlation between direct and social genotypes
simulate.R #builds list of parameters (group size 2 to 6, high or low P value of direct/social genotype correlation, effect sizes of local SGE 0 or 0.2, effect sizes of local DGE 0 or 1); lists appropriate variants for each set of parameters; builds covariance matrices - those are the same for all paramater sets)
qsub_draw_sims.sh #calls draw_sims_wPolygenic.R; draws 300 simulations from mvnorm; uses median values for random effects; draws from proportional, additive or maximum models for SGE (latter not discussed in paper)
create_h5.R #writes simulations in HDF5 files; samples columns so that different combinations get done at the same time
qsub_varianceDecomp.sh #calls exp_varianceDecomp.py
qsub_map_SGE.sh #to launch efficient_map_SGE_single_snp.py
QQplots_eval_conditioning.R # Suppl figure 3a-d
QQplots_conditioning_figure.R # barplots for Figure 3a 



