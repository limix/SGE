chr=$2
beg=$4
end=$5
phenotype=$1
peak=$3

locus="chr"$chr"_"$beg"_"$end

Rscript --vanilla code_dir/GWAS/fine_map/get_all_local_variants.R $chr $beg $end
Rscript --vanilla code_dir/GWAS/fine_map/add_local_dosages_to_h5.R $locus
python code_dir/GWAS/fine_map/fine_map_SGE.py $phenotype $locus
python code_dir/GWAS/fine_map/fine_map_DGE.py $phenotype $locus
Rscript --vanilla code_dir/GWAS/fine_map/paper_fine_maps.R $phenotype $locus
