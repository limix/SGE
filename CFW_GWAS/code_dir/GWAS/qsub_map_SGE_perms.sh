#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 25000
#BSUB -R "rusage[mem=25000]"
#BSUB -J "map_SGEperms[1-170]"

python code_dir/GWAS/map_SGE_perms.py $LSB_JOBINDEX pruned_dosages include

#submit with bsub < code_dir/GWAS/qsub_map_SGE_perms.sh
