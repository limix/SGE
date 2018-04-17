#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 10000
#BSUB -R "rusage[mem=10000]"
#BSUB -J "simsVD[1-300]"

python code_dir/simulations/power_analysis/exp_varianceDecomp.py $LSB_JOBINDEX pruned_dosages include

#submit with bsub < code_dir/simulations/power_analysis/qsub_varianceDecomp.sh
