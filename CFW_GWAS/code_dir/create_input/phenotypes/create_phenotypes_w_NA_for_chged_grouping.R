phenotypes = read.delim('data_dir/phenotypes/CFW_measures.txt',as.is=T)
#str(phenotypes,list.len=ncol(phenotypes))
rownames(phenotypes) = phenotypes[,1]
phenotypes = phenotypes[,-1]

# from Jerome's Evernotes: should remove mice where there was an issue with Sleep.Light
covariates = read.delim('data_dir/phenotypes/CFW_covariates.txt',as.is=T)
all(rownames(phenotypes) == covariates[,1])
#TRUE
pb_sleep_mice = covariates[which(covariates[,'Sleep.Light']==1),1]
sleep_cols = colnames(phenotypes)[grep('Sleep',colnames(phenotypes))]
phenotypes[pb_sleep_mice,sleep_cols] = NA

#exclude globally excluded
load('data_dir/subsets/globally_excluded.RData')
motch = match(globally_excluded,rownames(phenotypes))
phenotypes = phenotypes[-na.exclude(motch),]
dim(phenotypes)
#[1] 1869  200

imputed_cage_info = read.delim('data_dir/Supplementary_Table3.xls',as.is = T, header = T)

motch=match(rownames(phenotypes),rownames(imputed_cage_info))
any(is.na(motch))
#FALSE
imputed_cage_info = imputed_cage_info[motch,]
all(rownames(imputed_cage_info)== rownames(phenotypes))
#TRUE

dim(imputed_cage_info)
#[1] 1869   10
colnames(imputed_cage_info)
# [1] "cage_WH"      "cage_EPM"     "cage_OFT"     "cage_PAS"     "cage_Neo"    
# [6] "cage_SPPI"    "cage_FC"      "cage_PST"     "cage_Hypoxia" "cage_Cardio" 


####now create phenotype table with appropriate NAs
#need to match phenotypes to a cage info so that when a cage info is problematic, the relevant phenotypes are changed to NA
match_liste = list()
match_liste[["cage_WH"]] = c("WH.Ears_Area")
match_liste[["cage_EPM"]] = c("EPM.ClosedArms.Distance","EPM.ClosedArms.Entries","EPM.ClosedArms.Time","EPM.OpenArms.Distance",
	"EPM.OpenArms.Entries","EPM.OpenArms.Time","EPM.Total.Distance")
match_liste[["cage_OFT"]] = c("OFT.Boli","OFT.Centre.Time","OFT.Centre.Distance","OFT.Centre.Entries","OFT.Arena.Distance",
	"OFT.Periphery.Time","OFT.Periphery.Distance","OFT.Open.Time","OFT.Open.Distance")
match_liste[["cage_PAS"]] = c("PAS.First5","PAS.Last10","PAS.Total_Activity")
match_liste[["cage_Neo"]] = c("Neo.Latency")
match_liste[["cage_SPPI"]] = c("Weight.Startle","SPPI.ln_pa","SPPI.ln_pb","SPPI.ln_pc","SPPI.Habituation",
	"SPPI.ppReactivity","SPPI.pReactivity","SPPI.pc_average_pA","SPPI.pc_average_pB","SPPI.pc_average_pC","SPPI.pc_average_ABC",
	"SPPI.slope_pA","SPPI.slope_pB","SPPI.slope_pC","SPPI.slpPPI_average")
match_liste[["cage_FC"]] = c("FC.Training.Baseline","FC.Cue.Baseline","FC.Training.UnconditionedFreeze",
	"FC.Training.UnconditionedFreeze.Corrected","FC.Context.Freeze","FC.Context.Freeze.Corrected","FC.Cue.MeanFreeze",
	"FC.Cue.MeanFreeze.Corrected")
match_liste[["cage_PST"]] = c("PST.Immobility.First2min","PST.Immobility.Last4min")
match_liste[["cage_Hypoxia"]] = c("Hypoxia.MV_Baseline","Hypoxia.MV_AHR","Hypoxia.MV_HVD","Hypoxia.MV_Undershoot",
	"Hypoxia.MV_Off_response","Hypoxia.MV_SHR","Hypoxia.MV_NR","Hypoxia.f_Baseline","Hypoxia.f_AHR","Hypoxia.f_HVD",
	"Hypoxia.f_Undershoot","Hypoxia.f_Off_response","Hypoxia.f_SHR","Hypoxia.f_NR","Hypoxia.TV_Baseline","Hypoxia.TV_AHR",
"Hypoxia.TV_HVD","Hypoxia.TV_Undershoot","Hypoxia.TV_Off_response","Hypoxia.TV_SHR","Hypoxia.TV_NR","Weight.Hypo")
match_liste[["cage_Cardio"]] = c("Weight.ECG","Cardio.ECG.Heart_Rate","Cardio.ECG.PR_main","Cardio.ECG.PR_peak",
	"Cardio.ECG.P_Duration","Cardio.ECG.QRS_main","Cardio.ECG.QRS_peak","Cardio.ECG.QT_main","Cardio.ECG.QTcorr_main",
	"Cardio.ECG.QTcorr_peak","Cardio.ECG.JT_Interval","Cardio.ECG.Tpeak_Tend")
#not necessary
match_liste[["cage_Sleep"]] = c("Adrenals.Adrenals_g","Weight.Average","Bioch.Albumin","Bioch.ALP",                                
"Bioch.ALAT","Bioch.Amylase","Bioch.ASAT","Bioch.Calcium","Bioch.Chloride","Bioch.CreatinineEnzymatic",                
"Bioch.FreeFattyAcid","Bioch.Glucose","Bioch.Glycerol","Bioch.HDL","Bioch.Iron","Bioch.LDH","Bioch.LDL",
"Bioch.Phosphorous","Bioch.Potassium","Bioch.Sodium","Bioch.Tot.Billirubin","Bioch.Tot.Cholesterol","Bioch.Tot.Protein",
"Bioch.Triglycerides","Bioch.Urea","FACS.CD3pos","FACS.CD45posCD3posCD4pos","FACS.CD45posCD3posCD8pos","FACS.CD45posCD3negCD19pos",
"FACS.CD45posCD3negDX5pos","FACS.CD3posCD4pos","FACS.CD3posCD8pos","FACS.CD3posCD4CD8Ratio","FACS.CD3posCD4posCD44pos",
"FACS.CD3posCD8posCD44pos","FACS.CD3posCD44posCD4CD8Ratio","FACS.CD3posCD44negCD4CD8Ratio","Haem.WBCP","Haem.WBCB","Haem.RBC",
"Haem.measHGB","Haem.HCT","Haem.MCV","Haem.MCH","Haem.MCHC","Haem.CHCM","Haem.RDW","Haem.HDW","Haem.MPV","Haem.PDW","Haem.PCT",
"Haem.NEUT_percent","Haem.LYM_percent","Haem.MONO_percent","Haem.EOS_percent","Haem.LUC_percent","Haem.BASO_percent",
"Haem.abs_neuts","Haem.abs_lymphs","Haem.abs_mono","Haem.abs_eos","Haem.abs_lucs","Haem.abs_basos","Haem.PLT","Haem.Large_PLT",
"Muscles.TA.g","Muscles.EDL.g","Muscles.Gast.g","Muscles.Plant.g","Muscles.Sol.g","Tibia.Length","BMC.Mode","BMC.Median",
"BMC.Mean","BMC.Area","BMC.StdDev","BMC.Min","BMC.Max","BMC.Perim","BMC.Width","BMC.Height","BMC.Skew","BMC.Kurt","BMC.Mode.N",
"BMC.Median.N","BMC.Mean.N","BMC.StdDev.N","BMC.Min.N","BMC.Max.N","BMC.Skew.N","BMC.Kurt.N","Weight.Diss","Weight.BMI.body","Weight.BMI.tibia",
"Diss.Body.Length","Diss.Tail.Length","Diss.Brain.Weight","Micronucleus.Mn.NCE","Neuro.Ki67","Neuro.DCX","Serotonin.Serotonin_nM",
"Sleep.s24h","Sleep.s12h_L","Sleep.s12h_D","Sleep.sDif_LD","Sleep.VAR","Sleep.Percent_wake_over_17min","Sleep.Max_wake_episode.h",
"Sleep.short_sleep","Sleep.long_sleep","Sleep.long_sleep.percent_total","Sleep.longest_sleep.min","Sleep.VAR_1h","Sleep.VAR_12hL",
"Sleep.VAR_12hD","Sleep.VAR_24h","Sleep.sleep_L_onset","Sleep.sleep_D_onset","Sleep.T_max","Sleep.Ampl")                               

length(match_liste)
#[1] 11 (but remember cage Sleep in not in imputed_Cage_info)

# explore cage switches
#when mouse has been seen in more than one cage, print cage across time. this does not print mice with missing cage info
switching_mice = c()
for (i in 1:dim(imputed_cage_info)[1]) {
	if (length(na.exclude(unique(unlist(imputed_cage_info[i,]))))>1) {
		switching_mice = c(switching_mice,i)
	}
}
length(switching_mice)
# 15
head(imputed_cage_info[switching_mice,])
#true cage switches

#phenotypes need to be set to NA in following cases:
#cages missing from one point on
#one single cage missing for one mouse above
#cage switches

#for check purposes - see below
save_phenotypes = phenotypes

check_cage_in_row = function(row,cage) {
	return(cage %in% row)
}
#reminder
all(rownames(imputed_cage_info)== rownames(phenotypes))
#TRUE

count = 0
for (i in 1:dim(phenotypes)[1]) {
	#reminder imputed_cage_info does not have cage_Sleep column so cage change or cage going missing at that point won't be a problem
	eq=(imputed_cage_info[i,]==imputed_cage_info[i,1])
	pbs = which(is.na(eq) | !eq)
	if (length(pbs)==0) next else count = count + 1

	#turn all posterior phenotypes (including sleep etc) to NA - for focal mouse and its cage mates before and after the problem
	phenotypes[i,unlist(match_liste[pbs[1]:length(match_liste)])] = NA
	pbic_cages_before = unique(unlist(imputed_cage_info[i,1:(pbs[1]-1)]))
	for (pbic_cage_before in pbic_cages_before) {
		befores = which(apply(imputed_cage_info, FUN = check_cage_in_row,MAR = 1, cage = pbic_cage_before))
#		phenotypes[befores,unlist(match_liste[pbs[1]:length(match_liste)])] = (-99999)
		phenotypes[befores,unlist(match_liste[pbs[1]:length(match_liste)])] = NA
	}
	pbic_cages_after = unique(unlist(imputed_cage_info[i,pbs[1]:dim(imputed_cage_info)[2]]))
	for (pbic_cage_after in na.exclude(pbic_cages_after)) {
		afters = which(apply(imputed_cage_info, FUN = check_cage_in_row,MAR = 1, cage = pbic_cage_after))
#		phenotypes[afters,unlist(match_liste[pbs[1]:length(match_liste)])] = (-99999)
		phenotypes[afters,unlist(match_liste[pbs[1]:length(match_liste)])] = NA
	}

	# to double check print out before and after changing to NA
	#if (!all(is.na(save_phenotypes[i,unlist(match_liste[pbs[1]:length(match_liste)])]))) {
	#	print(paste('changing2',rownames(phenotypes)[i]))
	#	print(cbind(unlist(save_phenotypes[i,]),unlist(phenotypes[i,])))
	#	print(imputed_cage_info[i,])
	#}
	#looks correct
}
#count 71

#osteoporosis as per Jerome Nicod's evernotes
phenotypes$BMC.osteoporosis = as.numeric(phenotypes[,'BMC.Mode'] > 210)
### remove BMC.*.N traits (ie traits limited to rats without osteoporosis) as now fitting osteoporosis as a covariate
colnames(phenotypes)[grep('^BMC.*N$',colnames(phenotypes),perl = T)]
phenotypes = phenotypes[,-grep('^BMC.*N$',colnames(phenotypes),perl = T)]

g = grep('Corrected$',colnames(phenotypes))
remove = sub('.Corrected','',colnames(phenotypes)[g])
phenotypes = phenotypes[,-match(remove,colnames(phenotypes))]

save(phenotypes,file='data_dir/phenotypes/phenotypes_wNA.RData')


