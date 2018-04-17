#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 40000
#BSUB -R "rusage[mem=40000]"
#BSUB -J "hp_map_DGEperms[1-170]"
#BSUB -q highpri 

python code_dir/GWAS/map_DGE_perms.py $LSB_JOBINDEX pruned_dosages include

#submit with bsub < code_dir/GWAS/qsub_map_DGE_perms.sh
