#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -J "VD[1-170]"

python code_dir/variance_decomposition/whole_sample/exp_varianceDecomp.py $LSB_JOBINDEX pruned_dosages include
# 1-271 for exp_variance_decomp_init.py
# 1- 170 for exp_variance_decomp.py

#submit with bsub < code_dir/variance_decomposition/whole_sample/qsub_varianceDecomp.sh
