#BSUB -M 40000
#BSUB -R "rusage[mem=40000]"
#BSUB -o output_dir/output/%J_%I.out
#BSUB -e output_dir/error/%J_%I.err
#BSUB -J "dos[1-19]"

Rscript code_dir/create_input/dosages/get_web_dosages_ready.R $LSB_JOBINDEX

#submit with bsub < code_dir/create_input/dosages/qsub_get_web_dosages_ready.sh
