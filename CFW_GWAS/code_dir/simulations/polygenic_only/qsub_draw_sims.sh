#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 10000
#BSUB -R "rusage[mem=10000]"
#BSUB -J "draw_sims[1-67]"

Rscript --no-restore --no-save code_dir/simulations/draw_sims.R $LSB_JOBINDEX

#submit with bsub < code_dir/simulations/qsub_draw_sims.sh
