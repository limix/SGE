imputed_cage_info = read.delim('data_dir/Supplementary_Table3.xls',as.is = T, header = T)

# first will want to exclude mice whose cage is always NA (which is equivalent to first cage is NA)
# for QC purposes distinguish between different kinds of NA patterns in cage info (internal NAs, all NAs, remaining NAs - see examples below)
all_nas_mice = c()
remaining_nas_mice = c()
internal_nas_mice = c()
for (i in 1:dim(imputed_cage_info)[1]) {
	w = which(is.na(imputed_cage_info[i,]))
	if (length(w)==dim(imputed_cage_info)[2]) all_nas_mice = c(all_nas_mice,rownames(imputed_cage_info)[i]) else {
		if (length(w)!=0) {
			if ((length(w)!=length(seq(w[1],dim(imputed_cage_info)[2]))) | any(w!=seq(w[1],dim(imputed_cage_info)[2]))) internal_nas_mice = c(internal_nas_mice,rownames(imputed_cage_info)[i])  else {
				remaining_nas_mice = c(remaining_nas_mice, rownames(imputed_cage_info)[i])
			}
		}
	}
}
imputed_cage_info[internal_nas_mice,]
#               cage_WH cage_EPM cage_OFT cage_PAS cage_Neo cage_SPPI cage_FC
#Q_CFW-SW/33.0h 3845727  3845727  3845727  3845727  3845727   3845727 3845727
#               cage_PST cage_Hypoxia cage_Cardio
#Q_CFW-SW/33.0h  3845727         <NA>     3860027
# indeed cannot impute
length(internal_nas_mice)
#[1] 1 - - those will not be globally excluded but will turn their and their cage mates' genotypes to NA from time NA appears on

imputed_cage_info[all_nas_mice[1],]
#                cage_WH cage_EPM cage_OFT cage_PAS cage_Neo cage_SPPI cage_FC
#Q_CFW-SW/115.0c    <NA>     <NA>     <NA>     <NA>     <NA>      <NA>    <NA>
#                cage_PST cage_Hypoxia cage_Cardio
#Q_CFW-SW/115.0c     <NA>         <NA>        <NA>
length(all_nas_mice)
#64 mice without any cage info - those will be globally excluded

imputed_cage_info[remaining_nas_mice[1],]
#                cage_WH cage_EPM cage_OFT cage_PAS cage_Neo cage_SPPI cage_FC
#Q_CFW-SW/100.0b 3895273     <NA>     <NA>     <NA>     <NA>      <NA>    <NA>
#                cage_PST cage_Hypoxia cage_Cardio
#Q_CFW-SW/100.0b     <NA>         <NA>        <NA>
length(remaining_nas_mice)
#[1] 73 where cage info becomes missing at some point - those will not be globally excluded but will turn their and their cage mates' genotypes to NA from time NA appears on


# second there may also be mice which are with N!=2 cage mates from beginning on (cage changes will be dealth with in create_phenotypes_w_NA_for_chged_grouping.R)
not3_mice = c()
for (i in 1:dim(imputed_cage_info)[1]) {
	if (is.na(imputed_cage_info[i,1])) next
	mouse = rownames(imputed_cage_info)[i]
	cagemates = rownames(imputed_cage_info)[which(imputed_cage_info[,1]==imputed_cage_info[i,1])]
	cagemates = cagemates[cagemates!=mouse]
	#if this mouse is alone in (non NA) cage, won't have social partners.
	if (length(cagemates) != 2) not3_mice = c(not3_mice,mouse)
}
length(not3_mice)
#[1] 121


# third exclude all mice of cages where some cage mates have not been genotyped
unique_cages = unique(imputed_cage_info[!rownames(imputed_cage_info) %in% c(all_nas_mice,not3_mice),1])
all_cagemates = rownames(imputed_cage_info)[imputed_cage_info[,1] %in% unique_cages]

load('data_dir/dosages/pruned_dosages/unstd_my_final_dosages_chr19.RData')
all_genotyped_mice = rownames(dosages)
non_genotyped_cagemates = all_cagemates[which(!all_cagemates %in% all_genotyped_mice)]
pb_cages = unique(imputed_cage_info[non_genotyped_cagemates,1])
w_non_geno_cm = rownames(imputed_cage_info)[imputed_cage_info[,1] %in% pb_cages]
length(w_non_geno_cm)
#[1] 63

# globally_excluded will be excluded as both focal mice and cage mates
globally_excluded = c(all_nas_mice,not3_mice,w_non_geno_cm)
length(globally_excluded)
#[1] 248
save(globally_excluded,file = 'data_dir/subsets/globally_excluded.RData')



