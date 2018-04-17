#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 15000
#BSUB -R "rusage[mem=15000]"
#BSUB -J "draw_sims[61-72]"
Rscript --no-restore --no-save code_dir/simulations/power_analysis/draw_sims_wPolygenic.R $LSB_JOBINDEX

#submit with bsub < code_dir/simulations/power_analysis/qsub_draw_sims.sh
#61-72 for gs = 3 colMax
#73-84 for gs = 3 add
