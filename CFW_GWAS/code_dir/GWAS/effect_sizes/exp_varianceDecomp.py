import os
import sys
sys.path.insert(0,'code_dir/LIMIX')
from social_data_CFW import SocialData
from dirIndirVD_wSS import DirIndirVD
import pdb
import csv
import numpy as sp

if __name__=='__main__':
    
    task = 'VD_CFW_effect_sizes_sig'
   
    col=int(sys.argv[1])-1

    kinship_type = sys.argv[2]
    
    subset = sys.argv[3]

    effect = sys.argv[4]

    DGE = "DGE"
    IGE = "IGE"
    IEE = "IEE"
    cageEffect = "cageEffect"
    
    ### that's where QTLs are added as covariates
    data = SocialData(task,kinship_type,subset,effect)
    doto = data.get_data(col)
    
    trait=doto['trait']
    print 'trait is ' + trait

    std_genos = False    
    vc = DirIndirVD(doto['pheno'],doto['pheno_ID'],doto['covs'],doto['covs_ID'],doto['covariates_names'],doto['kinship_full'],doto['kinship_full_ID'],doto['cage_full'],doto['cage_full_ID'], DGE = (DGE is not None), IGE = (IGE is not None), IEE = (IEE is not None), cageEffect = (cageEffect is not None), calc_ste=False, subset_IDs = doto['subset_IDs'], std_genos = std_genos)

    rv=vc.getOutput()

    toWrite=sp.concatenate(((trait,rv['sample_size'],rv['sample_size_cm']),rv['covariates_names'],rv['covs_betas'][:,],(rv['var_Ad'],rv['var_As'],rv['total_var'],rv['conv'],rv['LML'])))

    if 'effect_sizes_sig' in task: 
        if std_genos:
            betas_file_dir = "".join(['/homes/abaud/CFW/output/reproduce/effect_sizes/betas/var_expl/',effect,'/'])
        elif not std_genos:
            betas_file_dir = "".join(['/homes/abaud/CFW/output/reproduce/effect_sizes/betas/allelic_effect/',effect,'/'])
    
    if not os.path.exists(betas_file_dir):
        os.makedirs(betas_file_dir)
    betas_file_name = "".join([betas_file_dir,trait,'.txt'])
    betas_file=open(betas_file_name,'w')
    betas_file.write("\t".join(str(e) for e in toWrite)+'\n')
    betas_file.close()

