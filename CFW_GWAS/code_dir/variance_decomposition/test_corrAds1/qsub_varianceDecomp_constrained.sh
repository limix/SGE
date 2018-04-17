#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -J "cVD[1-170]"
#BSUB -M 20000
#BSUB -R "rusage[mem=20000]"
#BSUB -q highpri

python code_dir/variance_decomposition/whole_sample/test_corrAds1/exp_varianceDecomp_constrained.py $LSB_JOBINDEX pruned_dosages include

#submit with bsub < code_dir/variance_decomposition/whole_sample/test_corrAds1/qsub_varianceDecomp_constrained.sh
