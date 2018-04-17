# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.
"""
qtl.py contains wrappers around C++ Limix objects to streamline common tasks in GWAS.
"""

import numpy as np
import limix
import scipy as sp
import limix.deprecated as dlimix
import time

class lmm:
	def __init__(self, snps, pheno, K=None, covs=None, test='lrt', NumIntervalsDelta0=2, NumIntervalsDeltaAlt=100, Ldeltamax0=(-4.), searchDelta=False, 
verbose=None):
		"""
		Univariate fixed effects linear mixed model test for all SNPs

		If phenotypes have missing values, then the subset of individuals used for each phenotype column
		will be subsetted

		Args:
			snps:   [N x S] np.array of S SNPs for N individuals
			pheno:  [N x P] np.array of P phenotype sfor N individuals
			K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
							If not provided, then linear regression analysis is performed
			covs:   [N x D] np.array of D covariates for N individuals
			test:   'lrt' for likelihood ratio test (default) or 'f' for F-test
			NumIntervalsDelta0:     number of steps for delta optimization on the null model (100)
			NumIntervalsDeltaAlt:   number of steps for delta optimization on the alt. model (100), requires searchDelta=True to have an effect.
			searchDelta:     Carry out delta optimization on the alternative model? if yes We use NumIntervalsDeltaAlt steps
			verbose: print verbose output? (False)
		"""
		print 'my Ldeltamax0 is:' 
		print Ldeltamax0 

		#create column of 1 for fixed if nothing provide
		if len(pheno.shape)==1:
			pheno = pheno[:,sp.newaxis]

		self.verbose = dlimix.getVerbose(verbose)
		self.snps = snps
		self.pheno = pheno
		self.K = K
		self.covs = covs
		#two lines below say don't estimate covariance, instead use provided covariance matrix
		self.searchDelta = searchDelta
		self.Ldeltamax0 = Ldeltamax0
		self.test = test
		self.NumIntervalsDelta0 = NumIntervalsDelta0
		self.NumIntervalsDeltaAlt = NumIntervalsDeltaAlt
		self.verbose = verbose
		self.N       = self.pheno.shape[0]
		self.P       = self.pheno.shape[1]
		self.Iok     = ~(np.isnan(self.pheno).any(axis=1))
		if self.K is None:
			self.searchDelta=False
			self.K = np.eye(self.snps.shape[0])
		if self.covs is None:
			raise Exception('covs should have at least column of 1s')
			self.covs = np.ones((self.snps.shape[0],1))

		self._lmm = None
		#run
		self.verbose = verbose
		self.process()

	def process(self):
		t0 = time.time()
		if self._lmm is None:
			self._lmm = limix.deprecated.CLMM()
			self._lmm.setK(self.K)
			self._lmm.setSNPs(self.snps)
			self._lmm.setPheno(self.pheno)
			self._lmm.setCovs(self.covs)
			self._lmm.setLdeltamax0(self.Ldeltamax0)

			if self.test=='lrt':
				self._lmm.setTestStatistics(self._lmm.TEST_LRT)
			elif self.test=='f':
				self._lmm.setTestStatistics(self._lmm.TEST_F)
			else:
				print((self.test))
				raise NotImplementedError("only f and lrt are implemented")
			#set number of delta grid optimizations?
			self._lmm.setNumIntervals0(self.NumIntervalsDelta0)
			if self.searchDelta:
				self._lmm.setNumIntervalsAlt(self.NumIntervalsDeltaAlt)
			else:
				self._lmm.setNumIntervalsAlt(0)

		if not np.isnan(self.pheno).any():
			#process
			self._lmm.process()
			self.pvalues = self._lmm.getPv()
			self.beta_snp = self._lmm.getBetaSNP()
			self.beta_ste = self._lmm.getBetaSNPste()
			self.ldelta_0 = self._lmm.getLdelta0()
			self.ldelta_alt = self._lmm.getLdeltaAlt()
			self.NLL_0 = self._lmm.getNLL0()
			self.NLL_alt = self._lmm.getNLLAlt()
		else:
			if self._lmm is not None:
				raise Exception('cannot reuse a CLMM object if missing variables are present')
			else:
				raise Exception('should never get here')
				self._lmm = limix.deprecated.CLMM()
			#test all phenotypes separately
			raise Exception('should really never get here')
			self.pvalues = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.beta_snp = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.beta_ste = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.ldelta_0 = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.ldelta_alt = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.NLL_0 = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.NLL_alt = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			self.test_statistics = np.zeros((self.phenos.shape[1],self.snps.shape[1]))
			for ip in np.arange(self.phenos.shape[1]):
				pheno_ = self.phenos[:,ip]
				i_nonz = ~(pheno_.isnan())

				self._lmm.setK(self.K[i_nonz,i_nonz])
				self._lmm.setSNPs(self.snps[i_nonz])
				self._lmm.setPheno(pheno_[i_nonz,np.newaxis])
				self._lmm.setCovs(self.covs[i_nonz])
				self._lmm.process()
				self.pvalues[ip:ip+1] = self._lmm.getPv()
				self.beta_snp[ip:ip+1] = self._lmm.getBetaSNP()
				self.beta_ste[ip:ip+1] = self._lmm.getBetaSNPste()
				self.ldelta_0[ip:ip+1] = self._lmm.getLdelta0()
				self.ldelta_alt[ip:ip+1] = self._lmm.getLdeltaAlt()
				self.NLL_0[ip:ip+1] = self._lmm.getNLL0()
				self.NLL_alt[ip:ip+1] = self._lmm.getNLLAlt()
				self.test_statistics[ip:ip+1] = self._lmm.getTestStatistics()
				pass
		if self._lmm.getTestStatistics() == self._lmm.TEST_LRT and self.test != "lrt":
			raise NotImplementedError("only f and lrt are implemented")
		elif self._lmm.getTestStatistics() == self._lmm.TEST_F and self.test != "f":
			raise NotImplementedError("only f and lrt are implemented")

		if self._lmm.getTestStatistics() == self._lmm.TEST_F:
			self.test_statistics = (self.beta_snp*self.beta_snp)/(self.beta_ste*self.beta_ste)
		if self._lmm.getTestStatistics() == self._lmm.TEST_LRT:
			self.test_statistics = 2.0 * (self.NLL_0 - self.NLL_alt)
		t1=time.time()

		if self.verbose:
			print(("finished GWAS testing in %.2f seconds" %(t1-t0)))

	def setCovs(self,covs):
		self._lmm.setCovs(covs)

	def setSNPs(self,snps):
		self._lmm.setSNPs(snps)

	def getBetaSNP(self):
		return self.beta_snp

	def getPv(self):
		"""
		Returns:
			[P x S] np.array of P-values
		"""
		return self.pvalues

def test_lm(snps,pheno, covs=None, test='lrt',verbose=None):
	"""
	Univariate fixed effects linear model test for all SNPs
	(wrapper around test_lmm, using identity kinship)

	If phenotypes have missing values, then the subset of individuals used for each phenotype column
	will be subsetted

	Args:
		snps:   [N x S] np.array of S SNPs for N individuals
		pheno:  [N x 1] np.array of 1 phenotype for N individuals
		covs:   [N x D] np.array of D covariates for N individuals
		test:   'lrt' for likelihood ratio test (default) or 'f' for F-test
		verbose: print verbose output? (False)

	Returns:
		limix LMM object
	"""
	lm = test_lmm(snps=snps,pheno=pheno,K=None,covs=covs, test=test,verbose=verbose, NumIntervalsDelta0=100,NumIntervalsDeltaAlt=100,searchDelta=False)
	return lm

def test_lmm(snps,pheno,K=None,covs=None, test='lrt',NumIntervalsDelta0=100,NumIntervalsDeltaAlt=100,searchDelta=False,verbose=None):
	"""
	Univariate fixed effects linear mixed model test for all SNPs

	If phenotypes have missing values, then the subset of individuals used for each phenotype column
	will be subsetted

	Args:
		snps:   [N x S] np.array of S SNPs for N individuals
		pheno:  [N x 1] np.array of 1 phenotype for N individuals
		K:      [N x N] np.array of LMM-covariance/kinship koefficients (optional)
						If not provided, then linear regression analysis is performed
		covs:   [N x D] np.array of D covariates for N individuals
		test:   'lrt' for likelihood ratio test (default) or 'f' for F-test
		NumIntervalsDelta0:     number of steps for delta optimization on the null model (100)
		NumIntervalsDeltaAlt:   number of steps for delta optimization on the alt. model (100), requires searchDelta=True to have an effect.
		searchDelta:     Carry out delta optimization on the alternative model? if yes We use NumIntervalsDeltaAlt steps
		verbose: print verbose output? (False)

	Returns:
		LMM object
	"""
	lmm_ = lmm(snps=snps, pheno=pheno, K=K, covs=covs, test=test, NumIntervalsDelta0=NumIntervalsDelta0, NumIntervalsDeltaAlt=NumIntervalsDeltaAlt, searchDelta=searchDelta, verbose=verbose)
	return lmm_


