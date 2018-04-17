# GWAS of real and permuted data + call genome-wide significant associations (QTLs)
cd code_dir/GWAS/
qsub_map_SGE.sh calls map_SGE.py
qsub_map_DGE.sh calls map_DGE.py
call_QTLs.R # makes list of would be QTLs at P 0.01 for DGE and SGE

qsub_map_SGE_perms.sh calls map_SGE_perms.py
qsub_map_DGE_perms.sh calls map_DGE_perms.py
qsub_call_QTLs_permutations.sh # for DGE then SGE; calls call_QTLs_permutations.R #identifies loci associated at nominal P < 0.01 (merge neighbouring loci when necessary)
list_QTLs_per_pheno_FDRbis.R #to call SGE and DGE QTLs at per-phenotype FDR 10%

SFig4.R #plots Suppl_Figure4.jpeg

# effect sizes of genome-wide significant associations
code_dir/GWAS/effect_sizes/MasterScript_effect_sizes.sh # master script for effect size analysis

write_sgeQTLs.R #Table 2
write_dgeQTLs.R #Suppl Table 2

# fine-map (ie zoom in on genome-wide significant locus and plot association for all variants rather than LD-pruned ones only)
code_dir/GWAS/fine_map/MasterScript_fine_map.sh # master script for zoom in plots (Suppl Figure 5)

#porcupine plot (Figure 2)
porcupine_plot.R

