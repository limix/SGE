#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -M 15000
#BSUB -R "rusage[mem=15000]"
#BSUB -J "draw_sims[1-8]"

Rscript --no-restore --no-save code_dir/simulations/polygenicNdgeQTL/draw_sims.R $LSB_JOBINDEX

#submit with bsub < code_dir/simulations/polygenicNdgeQTL/qsub_draw_sims.sh
