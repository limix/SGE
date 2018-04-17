#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -J "QTL_perms[1-170]"
#BSUB -M 20000
#BSUB -R "rusage[mem=20000]"

phenotype=`head -$LSB_JOBINDEX data_dir/phenotypes/colnames_phenotypes_bc_to_run.txt | tail -1`
Rscript --vanilla code_dir/GWAS/call_QTLs_permutations.R SGE $phenotype

#submit with bsub < code_dir/GWAS/qsub_call_QTLs_permutations.sh
