import scipy as sp
import h5py
#re for regular expressions
import re
import pdb

class SocialData():
    
    def __init__(self, task = None, kinship_type = "", subset = None, effect = None, chr = None):
        assert task is not None, 'Specify task!'
        self.task=task
        self.kinship_type=kinship_type
        self.subset = subset
        self.effect = effect
        self.chr = chr
        self.load()
    
    
    def load(self):
        
        if 'VD_CFW_sims' in self.task:
            if self.task == 'VD_CFW_sims_polygenicNdgeQTL':
                in_file = '/homes/abaud/CFW/data/reproduce/simulations/polygenicNdgeQTL/simulations/simulations.hdf5'
            elif self.task == 'VD_CFW_sims_polygenic_only':
                in_file = '/homes/abaud/CFW/data/reproduce/simulations/polygenic_only/simulations/simulations.hdf5'
            elif self.task == 'VD_CFW_sims_power_analysis':
                in_file = '/homes/abaud/CFW/data/reproduce/simulations/power_analysis/simulations/simulations.hdf5'
            f = h5py.File(in_file,'r')
    
            self.measures = f['simulations']['cols_measures']['measures'][:]
            self.all_pheno = f['simulations']['array'][:].T
            self.pheno_ID = f['simulations']['rows_subjects']['outbred'][:]
            self.all_covs2use = None
            f.close()

            in_file = '/homes/abaud/CFW/data/reproduce/CFWmice.h5'
            f = h5py.File(in_file,'r')
            self.kinship_full = f['GRM'][self.kinship_type]['matrix'][:]
            self.kinship_full_ID = f['GRM'][self.kinship_type]['row_header']['sample_ID'][:]
            self.cage_full = f['phenotypes']['row_header']['cage'][:]
            self.cage_full_ID = f['phenotypes']['row_header']['sample_ID'][:]
            self.subset_IDs = f['subsets'][self.subset][:]
            f.close()

        elif 'VD_CFW' in self.task:
            if 'effect_sizes_sig' in self.task: 
                in_file = '/homes/abaud/CFW/data/reproduce/CFWmice_effect_sizes_sig.h5'
            else:
                in_file = '/homes/abaud/CFW/data/reproduce/CFWmice.h5'

            f = h5py.File(in_file,'r')
    
            self.measures = f['phenotypes']['col_header']['phenotype_ID'][:]
            self.all_pheno = f['phenotypes']['matrix'][:].T
            self.pheno_ID = f['phenotypes']['row_header']['sample_ID'][:]
            self.all_covs = f['covariates']['matrix'][:].T
            self.covs_ID = f['covariates']['row_header']['sample_ID'][:]
            self.covariates = f['covariates']['col_header']['covariate_ID'][:]
            self.cage_full = f['phenotypes']['row_header']['cage'][:]
            self.cage_full_ID = f['phenotypes']['row_header']['sample_ID'][:]

            if 'effect_sizes' in self.task: 
                if self.effect == 'SGE':
                    self.all_covs2use = f['phenotypes']['col_header']['covariatesUsed']['SGE']
                elif self.effect == 'DGE':
                    self.all_covs2use = f['phenotypes']['col_header']['covariatesUsed']['DGE']
            else:
                self.all_covs2use = f['phenotypes']['col_header']['covariatesUsed']

            if self.chr is not None:
                self.kinship_full = f['GRM'][self.kinship_type][''.join(['chr',str(chr)])],['matrix'][:]
            else:    
                self.kinship_full = f['GRM'][self.kinship_type]['matrix'][:]
            self.kinship_full_ID = f['GRM'][self.kinship_type]['row_header']['sample_ID'][:]
 
            if self.subset is None:
                self.subset_IDs = None
            else:
                self.subset_IDs = f['subsets'][self.subset][:]

        else:
            print "Nothing done: task unknown!"


    def get_data(self,col,col_MT = None):
        
        self.trait = self.measures[col]

        self.pheno = self.all_pheno[:,col]
        self.trait_MT = None
        # not used for paper analyses
        self.track_trait = None

        if self.all_covs2use is None:
            self.covs = None
            self.covs_ID = None
            covariates_names = None
        else:
            covs2use = self.all_covs2use[col].split(',')
            Ic = sp.zeros(self.covariates.shape[0],dtype=bool)
            for cov in covs2use:
                Ic = sp.logical_or(Ic,self.covariates==cov)
            covariates_names = self.covariates[Ic]
            print covariates_names
            self.covs = self.all_covs[:,Ic]
    
        return {'trait' : self.trait,
                'trait_MT' : self.trait_MT,
                'pheno' : self.pheno,
                'pheno_ID' : self.pheno_ID,
                'track_trait' : self.track_trait,
                'covs' : self.covs,
                'covs_ID' : self.covs_ID,
                'covariates_names' : covariates_names,
                'kinship_type' : self.kinship_type,
                'kinship_full' : self.kinship_full,
                'kinship_full_ID' : self.kinship_full_ID,
                'cage_full' : self.cage_full,
                'cage_full_ID' : self.cage_full_ID,
                'subset_IDs' : self.subset_IDs}






