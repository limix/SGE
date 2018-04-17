#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 30000
#BSUB -R "rusage[mem=30000]"
#BSUB -J "map_SGE_sims[1-800]"

python code_dir/simulations/polygenicNdgeQTL/map_SGE_single_snp.py $LSB_JOBINDEX pruned_dosages include False

#submit with bsub < code_dir/simulations/polygenicNdgeQTL/qsub_map_SGE.sh

