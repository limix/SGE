load('data_dir/phenotypes/phenotypes_wNA.RData')

covariates = read.delim('data_dir/phenotypes/CFW_covariates.txt',as.is=T)
#str(covariates,list.len=ncol(covariates))
rownames(covariates) = covariates[,1]

#exclude globally excluded from covariates table (already excluded from phenotype table)
load('data_dir/subsets/globally_excluded.RData')
motch = match(globally_excluded,rownames(covariates))
covariates = covariates[-na.exclude(motch),]

load('data_dir/phenotypes/imputed_cage_info.RData')
motch = match(rownames(phenotypes),rownames(imputed_cage_info))
any(is.na(motch))
#FALSE
imputed_cage_info = imputed_cage_info[motch,]

all(rownames(phenotypes)==rownames(covariates))
#TRUE

# add some covariates to covariates table
## first
Weight.Diss = phenotypes[,'Weight.Diss']
Weight.Hypo = phenotypes[,'Weight.Hypo']
Weight.Startle = phenotypes[,'Weight.Startle']
Tibia.Length = phenotypes[,'Tibia.Length']
BMC.osteoporosis = phenotypes[,'BMC.osteoporosis']
covariates = cbind(covariates,Weight.Diss,Weight.Hypo,Weight.Startle,Tibia.Length,BMC.osteoporosis)


### strategy will be to fit cage effect and no other random effect. cage is (mostly) nested within batch, Neuro.Staining.batch etc. as shown below
### ideal model would fit cage nested within batch etc but not possible. as a consequence 'cage effects' will soak up batch effects etc. model mismatch still possible, but it should be anticonservative with regard to SGE (as non-SGE signal remains unaccounted for by the model) as I believe there is covariate that would create fake SGE that would not be accounted for by cage and/or DGE.
uniq_cages = unique(imputed_cage_info[,1])
for (cage in uniq_cages) {
	cagemates = rownames(imputed_cage_info)[which(imputed_cage_info[,1]==cage)]
	l = length(na.exclude(unique(covariates[cagemates,'Neuro.Staining.batch'])))
	if (l >1) print(covariates[cagemates,'Neuro.Staining.batch'])
}
#boxplot(phenotypes[,"Neuro.Ki67"] ~ covariates[,'Neuro.Staining.batch'])
#a few batches for which Ki67 staining is different
#boxplot(phenotypes[,"Neuro.DCX"] ~ covariates[,'Neuro.Staining.batch'])
#looks ok
for (cage in uniq_cages) {
	cagemates = rownames(imputed_cage_info)[which(imputed_cage_info[,1]==cage)]
	l = length(na.exclude(unique(covariates[cagemates,'Batch'])))
	if (l >1) stop('pb')
}
#cage mates always from same batch as expected

#conclusion: Sex and Batch for all. Diss BW whenever Jerome used it for one of the measures of the assay. 
#Tibia.Length. Experimenter. Neuro.Staining.batch Neuro.Nb_slices_Ki67 Neuro.Nb_slices_DCX Hypoxia.Chamber Hypoxia.Time.Hour.Minute
#not day of year nor month nor year. 

#download Suppl Table 1 from Jerome's NG paper and save as Tab Delimited file
### see data_dir/phenotypes/model_menu.txt (downloaded from NG paper) for guidance on which covariates to use
model_menu = read.delim('data_dir/phenotypes/TableS1_phenotypes.txt',as.is=T,check.name = F)
model_menu[201:dim(model_menu)[1],]
model_menu = model_menu[1:200,]
rownames(model_menu) = model_menu[,'Name']

colnames(phenotypes)[!(colnames(phenotypes) %in% rownames(model_menu))]
# "BMC.osteoporosis"
cates = unlist(lapply(strsplit(rownames(model_menu),'.',fixed=T),'[',1))
uniq_cates = unique(cates)
# [1] "Weight"       "Adrenals"     "Muscles"      "Tibia"        "Diss"        
# [6] "BMC"          "WH"           "Neuro"        "Bioch"        "Haem"        
#[11] "Serotonin"    "Micronucleus" "FACS"         "Hypoxia"      "Cardio"      
#[16] "PST"          "EPM"          "FC"           "OFT"          "PAS"         
#[21] "Sleep"        "SPPI"         "Neo"         

#do not include Sex as this will be included for all
chosen_covs = list()
chosen_covs[['Weight']] = NA
chosen_covs[['Adrenals']] = NA
chosen_covs[['Adrenals_BWcorr']] = 'Weight.Hypo'
chosen_covs[['Muscles']] = NA
chosen_covs[['Muscles_BWcorr']] = 'Weight.Hypo'
chosen_covs[['Tibia']] = NA
chosen_covs[['Tibia_BWcorr']] = 'Weight.Hypo'
chosen_covs[['Diss']] = NA
chosen_covs[['Diss_BWcorr']] = 'Weight.Hypo'
chosen_covs[['BMC']] = 'BMC.osteoporosis'
chosen_covs[['WH']] = NA
chosen_covs[['Bioch']] = NA
chosen_covs[['Bioch_BWcorr']] = 'Weight.Diss'
chosen_covs[['Haem']] = NA
chosen_covs[['Serotonin']] = NA
chosen_covs[['Serotonin_BWcorr']] = 'Weight.Hypo'
chosen_covs[['Micronucleus']] = NA
chosen_covs[['FACS']] = NA
chosen_covs[['Hypoxia']] = c('Hypoxia.Chamber','Hypoxia.Time.Hour.Minute')
chosen_covs[['Hypoxia_BWcorr']] = c('Weight.Hypo','Hypoxia.Chamber','Hypoxia.Time.Hour.Minute')
chosen_covs[['Cardio']] = NA
chosen_covs[['Cardio_BWcorr']] = 'Weight.Hypo'
chosen_covs[['PST']] = NA
chosen_covs[['EPM']] = NA
chosen_covs[['EPM_BWcorr']] = 'Weight.Startle'
chosen_covs[['FC']] = 'FC.Box.Cue'
chosen_covs[['OFT']] = NA
chosen_covs[['PAS']] = NA
chosen_covs[['Sleep']] = NA
chosen_covs[['Sleep_BWcorr']] = 'Weight.Hypo'
chosen_covs[['SPPI']] = NA
chosen_covs[['SPPI_BWcorr']] = 'Weight.Startle'
chosen_covs[['Neo']] = NA
chosen_covs[['Neo_BWcorr']] = 'Weight.Startle'

my_model_menu = NULL
for (cate in names(chosen_covs)) {
	w = which(unlist(sapply(cates,grepl,cate)))
	for (k in w) {
		pheno = model_menu[k,'Name']
		if (grepl('_BWcorr',cate)) pheno = paste(pheno,'_BWcorr',sep='')
		my_model_menu = rbind(my_model_menu,c(pheno,paste(c('Sex',na.exclude(chosen_covs[[cate]])),collapse=',')))
	}
}

my_model_menu = rbind(my_model_menu,c('Neuro.Ki67','Sex,Neuro.Nb_slices_Ki67'))
my_model_menu = rbind(my_model_menu,c('Neuro.Ki67_BWcorr','Sex,Weight.Hypo,Neuro.Nb_slices_Ki67'))
my_model_menu = rbind(my_model_menu,c('Neuro.DCX','Sex,Neuro.Nb_slices_DCX'))
my_model_menu = rbind(my_model_menu,c('Neuro.DCX_BWcorr','Sex,Weight.Hypo,Neuro.Nb_slices_DCX'))

my_model_menu = rbind(my_model_menu,c('BMC.osteoporosis','Sex'))

rownames(my_model_menu) = my_model_menu[,1]
colnames(my_model_menu) = c('Name','my_covs')

to_add = my_model_menu[grep('_BWcorr',my_model_menu[,'Name']),'Name']
for (add in to_add) {
	phenotypes = cbind(phenotypes,phenotypes[sub('_BWcorr','',add)])
	colnames(phenotypes)[dim(phenotypes)[2]] = add
}


check = unname(my_model_menu[,'Name'][!my_model_menu[,'Name'] %in% colnames(phenotypes)])
check[!grepl('_BWcorr$',check)]
my_model_menu = my_model_menu[my_model_menu[,'Name'] %in% colnames(phenotypes),]

all(colnames(phenotypes)[!colnames(phenotypes) %in% my_model_menu[,'Name']])
# TRUE

#couple of final checks
all(rownames(covariates)==rownames(phenotypes))
# TRUE
covariates_strings = my_model_menu[,'my_covs']
names(covariates_strings) = my_model_menu[,'Name']
unique_covs = unique(unlist(sapply(covariates_strings,FUN = strsplit,',')))
all(unique_covs %in% colnames(covariates))
#TRUE


save(phenotypes,covariates,my_model_menu,file = 'data_dir/phenotypes/my_model_menu.RData')
write.table(my_model_menu, file = 'plots_dir/my_model_menu.xls',sep='\t',quote = F, row.names = T, col.names =T)




