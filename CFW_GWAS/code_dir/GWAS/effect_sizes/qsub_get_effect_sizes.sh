#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -J "size[1-72]"

python code_dir/GWAS/effect_sizes/exp_varianceDecomp.py $LSB_JOBINDEX pruned_dosages include DGE # with task = 'effect_sizes_sig', std_genos = True (to calculate variance explained) or False (to get betas)
python code_dir/GWAS/effect_sizes/exp_varianceDecomp.py $LSB_JOBINDEX pruned_dosages include SGE # with task = 'effect_sizes_sig', std_genos = True or False

#submit with bsub < code_dir/GWAS/effect_sizes/get_effect_sizes.sh
