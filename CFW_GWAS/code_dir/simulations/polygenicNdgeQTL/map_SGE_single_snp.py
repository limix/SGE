#mem usage 20GB

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
task = 'VD_CFW_sims_polygenicNdgeQTL'

col=int(sys.argv[1])
kinship_type = sys.argv[2]
subset = sys.argv[3]
conditioning = sys.argv[4]

data = SocialData(task,kinship_type,subset)
doto = data.get_data(col)
trait = doto['trait']
print trait
#e.g. simulation_le4_chr3_111250435

snp_chr = int(re.sub('chr','',trait.split('_')[2]))
snp_pos = int(trait.split('_')[3])

DGE = "DGE"
IGE = "IGE"
IEE = "IEE"
cageEffect = "cageEffect"

#try opening pvalues file early so that lmm doesnt run if file opening is going to fail...
if conditioning == 'True':
	pvalues_file_dir = "".join(['output_dir/pvalues_pointwise_conditional_simulations/polygenicNdgeQTL/SGE/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect])),'/conditioning'])
elif conditioning == 'False':
	pvalues_file_dir = "".join(['output_dir/pvalues_pointwise_conditional_simulations/polygenicNdgeQTL/SGE/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect])),'/not_conditioning'])

if not os.path.exists(pvalues_file_dir):
	os.makedirs(pvalues_file_dir)
pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'.h5'])
if os.path.exists(pvalues_file_name):
	sys.exit(0)

covar_outfile_dir = "".join(["output_dir/null_covars_simulations/polygenicNdgeQTL/",kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
covar_outfile_name = "".join([covar_outfile_dir,"/",trait,".h5"])
covar_outfile = h5py.File(covar_outfile_name,'r')
saved = {}
saved['sampleID'] = covar_outfile['sampleID'][:]
saved['pheno'] = covar_outfile['pheno'][:]
saved['covs'] = covar_outfile['covs'][:]
saved['covar_mat'] = covar_outfile['covar_mat'][:]
covar_outfile.close()

#import genotypes. direct and social geno already matched on rows and cols
input_file_name = 'data_dir/CFWmice.h5'
input_file = h5py.File(input_file_name,'r')
geno  = input_file['direct_pruned_dosages']
geno_matrix = geno['matrix'][:].T
geno_sample_ID = geno['row_header']['sample_ID'][:]
social_geno  = input_file['social_pruned_dosages']
social_geno_matrix = social_geno['matrix'][:].T
position = {
	"chr" : geno['col_header']['chr'][:],
	"pos"   : geno['col_header']['pos'][:],
	"cumpos" : geno['col_header']['cumpos'][:],
}
idx_snp = sp.where(sp.logical_and(position['chr'] == snp_chr,position['pos'] == snp_pos))[0][0]
input_file.close()

#match genotypes with the rest
Imatch = sp.nonzero(saved['sampleID'][:,sp.newaxis]==geno_sample_ID)
print "Imatch"+str(len(Imatch[0]))
saved['sampleID'] = saved['sampleID'][Imatch[0]]
saved['pheno'] = saved['pheno'][Imatch[0]]
saved['covs'] = saved['covs'][Imatch[0],:]
saved['covar_mat'] = saved['covar_mat'][Imatch[0],:][:,Imatch[0]]
geno_matrix = geno_matrix[Imatch[1],]
geno_sample_ID=geno_sample_ID[Imatch[1]]
social_geno_matrix = social_geno_matrix[Imatch[1],:]

pvalues_file = h5py.File(pvalues_file_name,'w')

if conditioning == 'True':
	lmm = qtl.test_lmm(snps=social_geno_matrix[:,idx_snp],pheno=saved['pheno'],covs = sp.concatenate([saved['covs'],geno_matrix[:,idx_snp][:,sp.newaxis]], axis=1),K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)
elif conditioning == 'False':
	lmm = qtl.test_lmm(snps=social_geno_matrix[:,idx_snp],pheno=saved['pheno'],covs = saved['covs'],K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)

pvalues = lmm.getPv()[0,0]

pvalues_file.create_dataset(name = 'chr',data = position['chr'][idx_snp])
pvalues_file.create_dataset(name = 'pos',data = position['pos'][idx_snp])
pvalues_file.create_dataset(name = 'pvalues',data = pvalues)
pvalues_file.close()

