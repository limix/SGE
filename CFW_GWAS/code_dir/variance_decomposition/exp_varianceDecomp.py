import os
import sys
sys.path.insert(0,'code_dir/LIMIX')
from social_data_CFW import SocialData
from dirIndirVD_wSS import DirIndirVD
import pdb
import csv
import h5py

if __name__=='__main__':
    
    phenos_to_run_file = open('data_dir/phenotypes/phenos_to_run.txt')
    phenos_to_run = csv.reader(phenos_to_run_file)
    for i, row in enumerate(phenos_to_run):
        if i==(int(sys.argv[1])-1):
            col = int(row[0])-1
            break

    task = 'VD_CFW'
   
    kinship_type = sys.argv[2]
    
    subset = sys.argv[3]

    DGE = None
    IGE = "IGE"
    IEE = "IEE"
    cageEffect = "cageEffect"
    
    data = SocialData(task,kinship_type,subset)
    doto = data.get_data(col)
    
    trait=doto['trait']
    print 'trait is ' + trait
    
    covar_outfile_dir = "".join(["output_dir/null_covars/",kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect]))])
    if not os.path.exists(covar_outfile_dir):
        os.makedirs(covar_outfile_dir)
    covar_outfile_name = "".join([covar_outfile_dir,"/",trait,".h5"])
    if os.path.exists(covar_outfile_name):
        sys.exit(0)

    VD_outfile_dir = "".join(['output_dir/VD/',kinship_type,'_',"_".join(filter(None, [subset, DGE, IGE, IEE, cageEffect])),'/'])
    if not os.path.exists(VD_outfile_dir):
        os.makedirs(VD_outfile_dir)
    VD_outfile_name="".join([VD_outfile_dir,'/',trait,'.txt'])
   
    vc = DirIndirVD(doto['pheno'],doto['pheno_ID'],doto['covs'],doto['covs_ID'], doto['covariates_names'],doto['kinship_full'],doto['kinship_full_ID'],doto['cage_full'],doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'])

    rv=vc.getOutput()

    toWrite=(trait,rv['sample_size'],rv['sample_size_cm'],rv['var_Ad'],rv['var_As'],rv['STE_Ad'],rv['STE_As'],rv['total_var'],rv['corr_Ads'],rv['STE_corr_Ads'],rv['corr_params'],rv['conv'],rv['LML'],rv['var_Ed'],rv['var_Es'],rv['corr_Eds'],rv['var_C'])
    VD_outfile=open(VD_outfile_name,'w')
    VD_outfile.write("\t".join(str(e) for e in toWrite)+'\n')
    VD_outfile.close()

    toSave = vc.getToSave()
    toSave_file = h5py.File(covar_outfile_name,'w')
    toSave_file.create_dataset(name = 'sampleID',data = toSave['sampleID'])
    toSave_file.create_dataset(name = 'pheno',data = toSave['pheno'])
    toSave_file.create_dataset(name = 'covs',data = toSave['covs'])
    toSave_file.create_dataset(name = 'covar_mat',data = toSave['covar_mat'])
    toSave_file.close()


