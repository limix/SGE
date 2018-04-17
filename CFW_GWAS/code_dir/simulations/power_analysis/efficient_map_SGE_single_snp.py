# high mem usage

import sys
import os
import h5py
import scipy as sp
sys.path.insert(0,'code_dir/LIMIX')
# runs GWAS without estimating the covariance structure: uses the total covariance matrix provided (built from DGE SGE DEE SEE and cage effects)
import noNoise_qtl as qtl
from social_data_CFW import SocialData
from dirIndirVD_wSS import DirIndirVD
import h5py
import pdb
import re

#first fetch trait name
task = 'VD_CFW_sims_power_analysis'

kinship_type = sys.argv[2]
subset = sys.argv[3]
data = SocialData(task,kinship_type,subset)
conditioning = sys.argv[4]

DGE = "DGE"
IGE = "IGE"
IEE = "IEE"
cageEffect = "cageEffect"

#import genotypes. direct and social geno already matched on rows and cols
input_file_name = 'data_dir/simulations/power_analysis/simulations/simulations.hdf5'
input_file = h5py.File(input_file_name,'r')

geno  = input_file['direct_pruned_dosages']
saved_geno_matrix = geno['matrix'][:].T
geno_sample_ID = geno['row_header']['sample_ID'][:]
geno_SNP_id = geno['col_header']['id'][:]
#chr1_7529341...
social_geno  = input_file['social_pruned_dosages']
saved_social_geno_matrix = social_geno['matrix'][:].T
social_geno_sample_ID = social_geno['row_header']['sample_ID'][:]
#checked that same as geno_sample_ID
social_geno_SNP_id = social_geno['col_header']['id'][:]
#chr1_7529341_gs2... !! not in same order as geno_SNP_id (not same length)
input_file.close()

if conditioning == 'True':
	pvalues_file_dir = "".join(['output_dir/pvalues_pointwise_conditional_simulations/power_analysis/SGE/single_snp/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect,'conditioning']))])
else:
	pvalues_file_dir = "".join(['output_dir/pvalues_pointwise_conditional_simulations/power_analysis/SGE/single_snp/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect,'not_conditioning']))])

if not os.path.exists(pvalues_file_dir):
	os.makedirs(pvalues_file_dir)

col_init=int(sys.argv[1])
boy = 24 #by
mox = 7200
col_start=(col_init-1)*boy+1-1
col_end=col_init*boy
if col_end > mox:
    col_end =mox
all_cols=range(col_start,col_end)

for col in all_cols:

	doto = data.get_data(col)
	trait = doto['trait']
	print trait
	#e.g. simulation_sd22_chr13_17471073
	#########

	#### CAREFUL HARD CODING HERE - ASSUMING sims_design has 60 rows with first 12 for gs2, then gs3 etc. rows 61 to 85 are for max and additive SGE models (max model not in paper)
	#### needs to match simulate.sh and draw_sims_wPolygenic.R
	########
	row_sd = int(re.sub('sd','',trait.split('_')[1]))
	if row_sd in range(1,13):
		gs = 2
	elif row_sd in range(13,25):
		gs = 3
	elif row_sd in range(25,37):
		gs = 4
	elif row_sd in range(37,49):
		gs = 5
	elif row_sd in range(49,61):
		gs = 6
	elif row_sd in range(61,85):
		gs = 3

	snp_id = "_".join([trait.split('_')[2],trait.split('_')[3]])
	ouca = sp.where(geno_SNP_id == snp_id)[0]
	if len(ouca) == 0:
		continue
	idx_snp_geno = ouca[0]
	social_snp_id = "".join([snp_id,'_gs',str(gs)])
	idx_snp_social_geno = sp.where(social_geno_SNP_id == social_snp_id)[0][0]

	pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'.h5'])
	if os.path.exists(pvalues_file_name):
		continue

	covar_outfile_dir = "".join(["output_dir/null_covars_simulations/power_analysis/",kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
	covar_outfile_name = "".join([covar_outfile_dir,"/",trait,".h5"])
	
	if not os.path.exists(covar_outfile_name):
		continue
	covar_outfile = h5py.File(covar_outfile_name,'r')
	saved = {}
	saved['sampleID'] = covar_outfile['sampleID'][:]
	saved['pheno'] = covar_outfile['pheno'][:]
	saved['covs'] = covar_outfile['covs'][:]
	saved['covar_mat'] = covar_outfile['covar_mat'][:]
	covar_outfile.close()

	#match genotypes with the rest
	Imatch = sp.nonzero(saved['sampleID'][:,sp.newaxis]==geno_sample_ID)
	print "Imatch"+str(len(Imatch[0]))
	saved['sampleID'] = saved['sampleID'][Imatch[0]]
	saved['pheno'] = saved['pheno'][Imatch[0]]
	saved['covs'] = saved['covs'][Imatch[0],:]
	saved['covar_mat'] = saved['covar_mat'][Imatch[0],:][:,Imatch[0]]
	geno_matrix = saved_geno_matrix[Imatch[1],]
	social_geno_matrix = saved_social_geno_matrix[Imatch[1],:]

	pvalues_file = h5py.File(pvalues_file_name,'w')

	if conditioning == 'True':
		lmm = qtl.test_lmm(snps=social_geno_matrix[:,idx_snp_social_geno],pheno=saved['pheno'],covs = sp.concatenate([saved['covs'],geno_matrix[:,idx_snp_geno][:,sp.newaxis]], axis=1),K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)
	else:
		lmm = qtl.test_lmm(snps=social_geno_matrix[:,idx_snp_social_geno],pheno=saved['pheno'],covs = saved['covs'],K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)
	pvalues = lmm.getPv()[0,0]

	pvalues_file.create_dataset(name = 'id',data = social_snp_id)
	pvalues_file.create_dataset(name = 'pvalues',data = pvalues)
	pvalues_file.close()


