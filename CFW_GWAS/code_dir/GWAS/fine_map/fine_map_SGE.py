#high mem usage

import sys
import os
import h5py
import scipy as sp
sys.path.insert(0,'code_dir/LIMIX')
# will map under EMMAX model
import noNoise_qtl as qtl
from social_data_CFW import SocialData
from dirIndirVD_wSS import DirIndirVD
import h5py
import pdb
import csv

#first fetch trait name
task = 'VD_CFW'
kinship_type = 'pruned_dosages'
subset = 'include'

pheno_name = sys.argv[1]
print "Pheno should be "+pheno_name
locus = sys.argv[2]

phenos_to_run_file = open('data_dir/phenotypes/colnames_phenotypes_bc.txt')
phenos_to_run = csv.reader(phenos_to_run_file)
for i, row in enumerate(phenos_to_run):
	if row[0] == pheno_name:
		col = i
		break

data = SocialData(task,kinship_type,subset)
doto = data.get_data(col)
trait = doto['trait']
print "Pheno is "+trait

#then fetch pickled object
DGE = 'DGE'
IGE = 'IGE'
IEE = 'IEE'
cageEffect = 'cageEffect'
covar_outfile_dir = "".join(["output_dir/null_covars/",kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
covar_outfile_name = "".join([covar_outfile_dir,"/",trait,".h5"])
covar_outfile = h5py.File(covar_outfile_name,'r')
saved = {}
saved['sampleID'] = covar_outfile['sampleID'][:]
saved['pheno'] = covar_outfile['pheno'][:]
saved['covs'] = covar_outfile['covs'][:]
saved['covar_mat'] = covar_outfile['covar_mat'][:]
covar_outfile.close()

#import genotypes
input_file_name = 'data_dir/CFWmice_locals.h5'
input_file = h5py.File(input_file_name,'r')
geno_sample_ID = input_file['local_social_dosages']['row_header']['sample_ID'][:]
social_geno  = input_file['local_social_dosages'][locus]
social_geno_matrix = social_geno['matrix'][:].T
position = {
	"id" : social_geno['col_header']['id'][:],
	"chr" : social_geno['col_header']['chr'][:],
	"pos"   : social_geno['col_header']['pos'][:],
	"hwe_fail"   : social_geno['col_header']['hwe_fail'][:],
	"info_fail"   : social_geno['col_header']['info_fail'][:],
	"pruned_in"   : social_geno['col_header']['pruned_in'][:],
}
geno  = input_file['local_direct_dosages'][locus]
geno_matrix = geno['matrix'][:].T

input_file.close()

#match genotypes with the rest
Imatch = sp.nonzero(saved['sampleID'][:,sp.newaxis]==geno_sample_ID)
print "Imatch"+str(len(Imatch[0]))
saved['sampleID'] = saved['sampleID'][Imatch[0]]
saved['pheno'] = saved['pheno'][Imatch[0]]
saved['covs'] = saved['covs'][Imatch[0],:]
saved['covar_mat'] = saved['covar_mat'][Imatch[0],:][:,Imatch[0]]
social_geno_matrix = social_geno_matrix[Imatch[1],:]
geno_matrix = geno_matrix[Imatch[1],:]
geno_sample_ID=geno_sample_ID[Imatch[1]]


#try opening pvalues file early so that lmm doesnt run if file opening is going to fail...
pvalues_file_dir = "".join(['output_dir/fine_map_files/local_pvalues/SGE/',kinship_type,'_',subset,'_DGE_IGE_IEE_cageEffect'])
if not os.path.exists(pvalues_file_dir):
	os.makedirs(pvalues_file_dir)
pvalues_file_name = "".join([pvalues_file_dir,'/',trait,'_',locus,'.h5'])
pvalues_file = h5py.File(pvalues_file_name,'w')

pvalues = [-1] * social_geno_matrix.shape[1]
lmm = qtl.test_lmm(snps=social_geno_matrix[:,0:1],pheno=saved['pheno'],covs = sp.concatenate([saved['covs'],geno_matrix[:,0:1]], axis=1),K=saved['covar_mat'],NumIntervalsDelta0=2,searchDelta=False)
for i in xrange(social_geno_matrix.shape[1]):
	covs = sp.concatenate([saved['covs'],geno_matrix[:,i:i+1]], axis=1)
	snp = social_geno_matrix[:,i:i+1]
	lmm.setCovs(covs)
	lmm.setSNPs(snp)
	lmm.process()
	pvalues[i] = lmm.getPv()[0,0]

pvalues_file.create_dataset(name = 'chr',data = position['chr'])
pvalues_file.create_dataset(name = 'pos',data = position['pos'])
pvalues_file.create_dataset(name = 'pvalues',data = pvalues)
pvalues_file.create_dataset(name = 'id',data = position['id'])
pvalues_file.create_dataset(name = 'hwe_fail',data = position['hwe_fail'])
pvalues_file.create_dataset(name = 'info_fail',data = position['info_fail'])
pvalues_file.create_dataset(name = 'pruned_in',data = position['pruned_in'])
pvalues_file.close()

