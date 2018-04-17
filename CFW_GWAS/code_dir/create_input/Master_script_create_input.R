# input
# reformat Supplementary_Table3.xlsx: remove top lines and Sleep column, save as tab delimited file with extension xls
imputed_cage_info = read.delim('data_dir/Supplementary_Table3.xls',as.is = T, header = T)
# rest is from http://wp.cs.ucl.ac.uk/outbredmice/

# 1) Cage information and filtering out of some mice
#identify mice that are in cages where not all mice have been genotyped, mice with NA for all cages, and mice that are singly housed from the beginning
#those will be excluded from all subsequent analyses
code_dir/create_input/subsets/globally_excluded_mice.R


# 2) Phenotypes, covariates
cd code_dir/create_input/phenotypes
# when  grouping was not stable across the whole phenotyping pipeline (ie when cage changed): turn phenotype of mice and all their cage mates (before and after event) to NA for phenotypes aftern a cage change (ie cage turns to NA or changes)
# note that Sleep Bioch and Haem were assayed after all mice singly housed. Will nevertheless investigate "persistent effect of cage mates" ie effect of the cage mates they had before being singly housed
create_phenotypes_w_NA_for_chged_grouping.R #globally_excluded mice excluded at this point
#decide on important covariates
check_covariates.R
#boxcox phenotypes and recode covariates as numerical
prep_all_pheno_data.R
#create HDF5 file
create_HDF5.R


# 3) Genotypes and GRM
cd code_dir/create_input/dosages
qsub_get_web_dosages_ready.sh #starts get_web_dosages_ready.R
merge_web_dosages.R #merge all chromosomes and exclude globally excluded mice 
build_GRM_from_studyDosages.R
add_GRM_to_h5.R
add_direct_dosages_to_h5.R
add_social_dosages_to_h5.R # social genotype is sum of direct genotypes of cage mates
Supplementary_Figures_2ac_3e.R # correlation between direct and social genotypes

# 4) subset
# remove closely related mice (based on GRM) and their cage mates. done by creating a /subset in HDF5
code_dir/create_input/dealw_close_pairs.R





