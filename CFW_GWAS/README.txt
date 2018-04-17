The code provided in this folder is relevant to the following manuscript:
Genome-wide association study of social genetic effects on 170 phenotypes in laboratory mice.
Amelie Baud, Francesco Paolo Casale, Jerome Nicod, Oliver Stegle
https://doi.org/10.1101/302349



This code should be seen as a guided path to reproducing the analyses of the paper and learn how to use LIMIX for SGE analyses. Do not expect well-written code that 'just runs'... I would be happy to help anyone keen to use it though: abaud@ebi.ac.uk.



When running the analyses, my folders were structured as follows:
- code_dir with all the code
- data_dir for primary data (which I retrieved from http://wp.cs.ucl.ac.uk/outbredmice/) and secondary data
- output_dir for tertiary data 
- plots_dir for results (plots, tables etc.)

In this github repository, only code_dir folder is present. You can download the primary data from their original website and derive subsequent data and results from the code.
Note that non-SGE-specific LIMIX functions need to be downloaded from http://github.com/limix/limix.



To perform analyses, follow the following order
1) create_input/Master_script_create_input.R
2) variance_decomposition/MasterScript_VD.R
3) GWAS/MasterScript_GWAS.R
4) simulations/Master_script_simulations.R
These MasterScripts are not R scripts, they are text files better viewed in a code editor hence the .R extension


Code for each display item in manuscript:
Figure 1: none
Figure 2: code_dir/GWAS/porcupine_plot.R
Figure 3a: code_dir/simulations/power_analysis/QQplots_conditioning_figure.R
Figure 3b-d: code_dir/GWAS/effect_sizes/plot_ve.R
Table 1: code_dir/variance_decomposition/write_VD.R
Table 2: code_dir/GWAS/write_sgeQTLs.R
Supplementary Table 1: code_dir/variance_decomposition/write_VD.R
Supplementary Table 2: code_dir/GWAS/write_dgeQTLs.R
Supplementary Table 3: none
Supplementary Figure 1: code_dir/simulations/polygenic_only/compare_results_sims_targets.R
Supplementary Figure 2a-c: code_dir/create_input/dosages/Supplementary_Figures_2ac_3e.R
Supplementary Figure 2d-e: code_dir/simulations/polygenicNdgeQTL/need4conditioning.R
Supplementary Figure 3a-d: code_dir/simulations/power_analysis/QQplots_eval_conditioning.R
Supplementary Figure 3e: code_dir/create_input/dosages/Supplementary_Figures_2ac_3e.R
Supplementary Figure 4: code_dir/GWAS/SFig4.R
Supplementary Figure 5: code_dir/GWAS/fine_map/paper_fine_maps.R
Supplementary Figure 6: code_dir/GWAS/effect_sizes/plot_af.R
