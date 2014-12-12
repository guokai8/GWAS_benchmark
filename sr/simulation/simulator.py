import numpy as np
import scipy as sp
import pdb, sys, pickle
#import matplotlib.pylab as plt
import datetime
import scipy.stats
from optparse import OptionParser
import pandas
import cPickle
import os
import scipy.stats as st
import external.fastlmmc as fc
import fastlmm.util.util as util
import logging

#class snps(object):

#class bed_snps(snps):

#    def __init__(self, bedfilename):

#     def generate_snps(self, sample_size, population_index=None, snp_index = None):
#        """
        
#        """
    
#     def sample_frequencies(self):
#         raise Exception("not implemented for this class")  
      

class sim_snps(object):
    """
    Generates genotypes with a certain MAF and optionally with population structure.
    In case of no population structure, they are sampled from a binomial,
    otherwise from a Beta-Binomial (Balding and Nichols, 1995).
    """
    def __init__(self, num_snps, num_differentiated = np.array([0]), MAF_ancestral = np.array([0.05,0.5]), Fst = np.array([[0.1,0.1]]), diploid = True, quiet = True, p_ancestral=None):
        self.differentiated = None
        self.num_snps=num_snps
        self.quiet = quiet
        self.MAF_ancestral = MAF_ancestral
        self.Fst = Fst                                  #dimensions: number populations times number SNP groups
        self.num_differentiated = num_differentiated    #dimensions: number of SNP groups
        self.p_ancestral=p_ancestral

        #set differentiated SNPs
        i_0 = 0
        self.differentiated = np.zeros((len(self.num_differentiated),self.num_snps),dtype = 'bool')
        
        for i_group, num_diff in enumerate(self.num_differentiated):
            self.differentiated[i_group,np.arange(i_0,i_0+num_diff)]=True
            i_0+=num_diff
        
        if diploid:
            self.chr_copies = 2
        else:
            self.chr_copies = 1

        #set the allele frequencies
        self.sample_frequencies()



    def sample_frequencies(self):
        '''
        sample the allele frequencies for all ancestral SNPs and all populations
        ancestral frequencies are sampled uniformly according to MAFs
        population frequencies from a Beta distribution (Balding and Nichols, 1995).
        '''
        #sample ancestral allele frequencies, len(p)=num_snps
        if self.p_ancestral is None:
            self.p_ancestral = np.random.uniform(self.MAF_ancestral[0], self.MAF_ancestral[1], self.num_snps)
        
        #alphas: shape number(populations) times number(snps) and copy over ancestral frequencies
        self.alphas = np.zeros((len(self.Fst[0]),self.p_ancestral.shape[0]))
        self.alphas[:,:]=self.p_ancestral
        
        #re-sample differentiated population frequencies
        for i_snpgroup in xrange(len(self.num_differentiated)):
            for i_population, F in enumerate(self.Fst[i_snpgroup]):
                i_snps = self.differentiated[i_snpgroup]
                p_anc = self.p_ancestral[i_snps]
                if F == 0.0: 
                    self.alphas[i_population,i_snps] = p_anc #special treatment if no population structure
                else:
                    self.alphas[i_population,i_snps] = np.random.beta(p_anc*(1.0-F)/F,(1.0-p_anc)*(1.0-F)/F, p_anc.shape[0])
        pass


    def generate_snps(self, sample_size, population_index=None, snp_index = None):
        """
        Generates genotypes with a certain MAF and optionally with population structure.
        In case of no population structure, they are sampled from a binomial,
        otherwise from a Beta-Binomial (Balding and Nichols, 1995).
        """
        if snp_index is None:
            snp_index = np.arange(self.num_snps)
        if population_index is None:
            if not self.quiet:
                print "Simulating SNPs from ancestral frequencies..."
            #generate from ancestral frequencies
            pgen = self.p_ancestral[snp_index]
        else:        
            if not self.quiet:
                print ("Simulating SNPs from population %i" % population_index)

            #generate from population frequencies    
            pgen = self.alphas[population_index,snp_index]
        snps = np.zeros((sample_size,pgen.shape[0]),dtype='int8')
        for i in xrange(self.chr_copies):
            #sample each allele
            rand = np.random.random((sample_size,pgen.shape[0]))
            snps[rand<pgen]+=1
        #snps_ = np.random.binomial(self.chr_copies,pgen,(sample_size,pgen.shape[0]))
        return snps

    def gen_snps_familystruct(self, snps_parents=None):
        pass

    def mate(self, snps_parents, i_parent, num_children_per_couple):
        '''
        given a set of snps for the parents, mate the individuals
        '''
        assert self.chr_copies==2, "assert self.chr_copies==2"
        num_snps = snps_parents[0].shape[1]
        num_trios = i_parent.shape[0]
        num_children = num_trios*num_children_per_couple #10 children per couple
        snps = np.zeros((num_children,num_snps),dtype='int8')
        
        #import pdb; pdb.set_trace()
        for i in xrange(len(snps_parents)):
            #sample each allele
            for j in xrange(num_children_per_couple):
                rand = np.random.random((num_trios,num_snps))
                snps[j*num_trios:(j+1)*num_trios][rand<0.5*snps_parents[i]]+=1
        return snps

    def generate_trios(self, snps_parents=None, num_trios=0, population_percentages=None, snp_index=None, num_children_per_couple=1):
        '''
        generate a single set of trios
        '''    
        assert self.chr_copies==2, "assert self.chr_copies==2"
        num_trios = int(num_trios)
        if snps_parents is None:
            #sample parents
            sample_size = 2*num_trios
            snps_parents = self.generate_snps_populations(sample_size=sample_size,population_percentages=population_percentages,snp_index=snp_index)
        else:
            sample_size = snps_parents.shape[0]
        num_snps = snps_parents.shape[1]
        assert sample_size>=2*num_trios, "sample_size>=2*num_trios"
        potential_parents = np.random.permutation(sample_size)
        i_done=0
        i_parent = np.empty([num_trios,self.chr_copies],dtype = 'int64')
        snps_parents_sampled = []
        for i in xrange(self.chr_copies):
            #sample each allele
            i_parent[:,i] = potential_parents[i_done:i_done+num_trios]
            snps_parents_sampled.append(snps_parents[i_parent[:,i]])
            i_done+=num_trios
        snps = self.mate(snps_parents_sampled, i_parent, num_children_per_couple)
        return snps, i_parent

    def generate_snps(self, sample_size, snps=None, population_index=None, snp_index = None):
        """
        Generates genotypes with a certain MAF and optionally with population structure.
        In case of no population structure, they are sampled from a binomial,
        otherwise from a Beta-Binomial (Balding and Nichols, 1995).
        """
        if snp_index is None:
            snp_index = np.arange(self.num_snps)
        if population_index is None:
            if not self.quiet:
                print "Simulating SNPs from ancestral frequencies..."
            #generate from ancestral frequencies
            pgen = self.p_ancestral[snp_index]
        else:        
            if not self.quiet:
                print ("Simulating SNPs from population %i" % population_index)

            #generate from population frequencies    
            pgen = self.alphas[population_index,snp_index]
        if snps is None:
            snps = np.zeros((sample_size,pgen.shape[0]),dtype='int8')
        else:
            snps[:]=0
        for i in xrange(self.chr_copies):
            #sample each allele
            rand = np.random.random((sample_size,pgen.shape[0]))
            snps[rand<pgen]+=1
        #snps_ = np.random.binomial(self.chr_copies,pgen,(sample_size,pgen.shape[0]))
        return snps

    def generate_snps_populations(self, sample_size=1000, population_percentages=None, snp_index=None):
        n_pop=self.alphas.shape[0]
        if snp_index is None:
            X=np.zeros((sample_size,self.num_snps),dtype='int8')
        else:
            X=np.zeros((sample_size,self.p_ancestral[snp_index].shape[0]),dtype='int8')
        if population_percentages is None:
            if not self.quiet:
                print "Simulating SNPs from all populations in equal proportions..."
            #create evenly distributed populations
            population_percentages = np.ones(n_pop)/n_pop
        if np.absolute(population_percentages.sum()-1.0)>0.0001:
            raise Exception("np.absolute(n_pop.sum()-1.0)>0.0001")
        n_done = 0
        for i_pop in xrange(n_pop):
            n_population = int(population_percentages[i_pop]*sample_size)
            if i_pop ==n_pop-1:
                n_population = X.shape[0] - n_done                
            X[n_done:n_done+n_population,:] = self.generate_snps(sample_size=n_population,population_index=i_pop,snp_index=snp_index)
            n_done += n_population
        return X


class sim_pheno(object):
    def __init__(self, num_causal=1, i_causal=None, W=None, W_covariates=None, noise_var=0.1, genetic_var=0.1, num_phenotypes=1, perc_causal_differentiated=np.array([0.5]), weight_distribution='normal', quiet=True, covariates_var=0.3, sim_snps=None, noise_distribution='normal'):
        #parameters needed for the liability threshold model
        self.W = W                          #snp effects
        self.W_covariates = W_covariates    #covariates effects
        self.i_causal = i_causal
        self.num_causal = num_causal
        self.quiet = quiet
        self.sim_snps = sim_snps            #snp simulation class
        self.num_causal = num_causal
        self.perc_causal_differentiated = perc_causal_differentiated
        self.num_phenotypes=num_phenotypes
        self.weight_distribution=weight_distribution
        self.noise_distribution=noise_distribution
        self.noise_var=noise_var,
        
        self.genetic_var=genetic_var
        
        #generate indices of causal SNPs
        if self.i_causal is None:
            self.generate_i_causal()

        #generate SNP effects
        if self.W is None:
            self.W = self.generateW(num_causal=self.num_causal, num_phenotypes=self.num_phenotypes, weight_distribution=self.weight_distribution,genetic_var=self.genetic_var)

        #generate covariates effects
        if self.W_covariates is None:
            pass
        pass


    def generate_i_causal(self):
        '''
        generate the indices of the causal SNPs
        Randomly sample causal SNPs
        '''
        diff =  self.sim_snps.differentiated
        n_diff = self.sim_snps.num_differentiated
        n_snps = self.sim_snps.num_snps
        n_groups = diff.shape[0]
        n_causal= self.num_causal
        i_causal = np.zeros(n_snps,dtype='bool')
        
        n_causal_tot=0
        #first non-differentiated causal SNPs
        i_notdiff=~diff.any(0)
        n_diff_g=int((1.0-self.perc_causal_differentiated.sum())*n_causal)
        if (n_snps-n_diff.sum())<n_diff_g:
                raise Exception("n_snps-n_diff.sum())<n_diff_g")
        int_index_g=np.where(i_notdiff)[0]
        perm = np.random.permutation(int_index_g.shape[0])
        i_causal[int_index_g[perm[:n_diff_g]]]=True
        n_causal_tot+=n_diff_g
        #then for each differentiatedness group mark differentiated causal SNPs
        for i_group in xrange(n_groups):
            n_diff_g=int(self.perc_causal_differentiated[i_group]*n_causal)
            if n_diff[i_group]<n_diff_g:
                raise Exception("n_diff[i_group]<n_diff_g")
            #pdb.set_trace()
            int_index_g=np.where(diff[i_group])[0]
            perm = np.random.permutation(int_index_g.shape[0])
            i_causal[int_index_g[perm[:n_diff_g]]]=True
            n_causal_tot+=n_diff_g
            pass
        #correct the number of causal SNPs
        n_diff=n_causal_tot-n_causal
        
        if n_diff>0:#too many causals
            int_c = np.where(i_causal)[0]
            perm = np.random.permutation(int_c.shape[0])
            i_causal[int_c[perm[:n_diff]]]=False
        elif n_diff<0:#too few causals
            int_c = np.where(~i_causal)[0]
            perm = np.random.permutation(int_c.shape[0])
            i_causal[int_c[perm[:-n_diff]]]=True
        self.i_causal=i_causal
        

    def __str__(self):
        message = "Simulation generated on %s \n\n" % self.generated_at
        message += "Dimensions: \n\t %d SNPs (%d causal), %d phenotype(s) \n\n" % (self.sim_snps.num_snps, self.num_causal, self.num_phenotypes)
        #message += "SNP info: \n\t MAF = %.3f, %d chromosome copies \n\n" % (self.sim_phen.MAF_, self.sim_phen.chr_copies)
        #if self.pop_struct:
        #    message += "Population structure: \n\t %d differentiated causal SNPs (%d total differentiated SNPs), F_st = %.3f\n\n" % (int(self.num_causal*self.perc_causal_differentiated), self.num_differentiated, self.Fst)
        message += "Variance components: \n\t noise=%.4f, genetic=%.4f" % (self.noise_var, self.genetic_var.mean())
        if self.pheno_transform != None:
            message += "\n\nPhenotype transformation: \n\t %s, parameter = %.2f" % (self.pheno_transform, self.transform_param)
        return message

    def generateW(self, num_causal=1, num_phenotypes=1, weight_distribution='normal', genetic_var=0.1, dof=1, mean = 0.0):
        '''
        generate weights
        '''
        if weight_distribution == 'Student':
            W = sp.stats.t.rvs(dof, mean, genetic_var, size=(num_causal,num_phenotypes))
        elif weight_distribution == 'normal':
            W = np.random.randn(num_causal, num_phenotypes) * np.sqrt(genetic_var) + mean
        else:
            raise Exception("not implemented")
        return W

    def generate_phenotype(self, X_causal=None, X=None, covariates = None, add_noise=True, add_covariates=False, add_assoc=True, add_interactions=False, freq_causal=None,freq=None,make_binary=False):
        """
        Generates the phenotype by summing together all the individual variance components.
        Note: the phenotype is not guaranteed to be centered  and scaled, because that
        functionality does not belong here. It's much better if the individual methods
        implement it.

        """
        if X_causal is None:
            X_causal=X[:,self.i_causal]
        N=X_causal.shape[0]
        
        Y = np.zeros((N, self.num_phenotypes))
        Z={}
        if add_assoc:
            if freq_causal is None and freq is not None:
                freq_causal=freq[self.i_causal]
            Z['assoc'] = self.add_associations(X_causal = X_causal, freq_causal=freq_causal)
        if add_noise:
            Z['noise'] = self.add_noise(X=X_causal)
        if add_interactions:
            raise NotImplementedError("interactions")
        if add_covariates:
            raise NotImplementedError("covariates")    

        for k in Z.keys():
            Y += Z[k]
            if not self.quiet:
                print "Variance due to %s = %.3f" % (k, Z[k].var())

        if make_binary:
            print "The mean (which should be about zero) is {0}".format(Y.mean())
            lessThanMeanVector = Y < 0 #!! make sure roughly 50/50
            print "The # of 1 is {0} of {1}".format(lessThanMeanVector.sum(),len(Y))
            Y[:] = 1
            Y[lessThanMeanVector] = 0

        return Y,Z

    def add_associations(self, X_causal, freq_causal=None):
        """
        Simulates associations.
        """

        if freq_causal is not None:
            if not self.quiet:
                print "standardizing SNPs by ancestral freqs..."
            X_causal=standardize_freq(snps=X_causal,freq=freq_causal)
        if not self.quiet:
            print "Adding associations..."
        XW = np.dot(X_causal, self.W)
        return XW

    def add_noise(self,X=None,dof=1,loc=0.0):
        """
        Adds Gaussian/Student noise
        """
        if not self.quiet:
            print "Adding noise..."
        noise_std = np.sqrt(self.noise_var)
        if self.noise_distribution=='normal':
            Z = noise_std*sp.randn(X.shape[0], self.num_phenotypes)
        elif self.noise_distribution=='Student':
            Z= sp.stats.t.rvs(dof, loc, self.noise_var, size=(X.shape[0],num_phenotypes))
        else:
            raise Exception("not implemented")
        return Z


    def get_phenotype_transform(self, Y, pheno_transform = 'liability', transform_param = None, prevalence = 0.5, population_percentages=None, num_samples=3000):
        N=Y.shape[0]
        if pheno_transform != None:
            if pheno_transform == "exp_ay":
                Y_transformed = np.zeros_like(Y)
                Y_transformed[:N/5.0] = np.exp( Y[:N/5.0] * transform_param)
            elif pheno_transform == "exp_root":
                Y_transformed = np.exp(Y)**(1.0/transform_param)
            elif pheno_transform == "rounding":
                if type(self.transform_param)==np.ndarray:
                    Y_transformed = np.zeros_like(Y)
                    for i in xrange(Y.shape[1]):
                        Y_transformed[:,i]=Y[:,i].round(transform_param)
                else:
                    Y_transformed = Y.round(transform_param)
            elif pheno_transform =='liability':
                if transform_param is None:
                    transform_param = self.estimate_liab_threshold(num_samples=num_samples, prevalence=prevalence, sim_snps=self.sim_snps, population_percentages=population_percentages)
                Y_transformed = np.array(Y>transform_param, dtype = np.float)
            return Y_transformed,transform_param
        else:
            return Y,transform_param
    
    def estimate_liab_threshold(self,  sim_snps, prevalence=0.5, num_samples=3000, population_percentages=None ):
        #TODO pass add_covariates etc.
        print "liab estimation only works with default params"
        X_causal = sim_snps.generate_snps_populations(sample_size=num_samples, population_percentages=population_percentages, snp_index=self.i_causal)
        Y,Z = self.generate_phenotype(X_causal=X_causal,add_noise = self.add_noise,add_assoc=self.add_associations,freq=sim_snps.p_ancestral)
        Y.sort(0)
        threshold =  Y[int(num_samples*(1.0-prevalence))]        
        return threshold

    def write_plink(self, X, Y, basefilename='test',i_SNP=None,make_bed = False,merge_obs=False):
        import pandas
        import plink_write as pw
        sample_names =  np.array(['s_%d' % i for i in range(Y.shape[0])], dtype=str)
        family_names = np.array(['f_%d' % i for i in range(Y.shape[0])], dtype=str)[:, None]
        paternal_ids = np.zeros_like(family_names,dtype=str)
        maternal_ids = np.zeros_like(family_names,dtype=str)
        sex = np.zeros_like(family_names,dtype=str)
        fam_pid_mid_sex = np.concatenate((family_names, paternal_ids, maternal_ids, sex), axis=1)
        F = pandas.DataFrame(data=fam_pid_mid_sex, index=sample_names, columns=['Family', 'Paternal ID', 'Maternal ID', 'Sex'])
        if merge_obs:
            snp_names=np.empty(X.shape[1],dtype='|S15')
            for i in range(X.shape[1]):
                if i_SNP['causal'][i]:
                    snp_names[i] = 'snp_%d_c' %i
                else:
                    snp_names[i] = 'snp_%d' %i
            G = pandas.DataFrame(data=X, index=sample_names, columns=snp_names)
            inter_fn=basefilename
            pw.write_plink(G,F,intermediate_file_name=inter_fn,final_file_name=basefilename,make_bed=make_bed)
        else:
            if i_SNP is not None:
                for k in i_SNP.keys():
                    i_SNP[k]
                    snp_names = ['snp_%s_%d' % (k,i) for i in range(X[:,i_SNP[k]].shape[1])]
                    G = pandas.DataFrame(data=X[:,i_SNP[k]], index=sample_names, columns=snp_names)
                    basefilename_k=basefilename+'_'+k
                    #inter_fn=basefilename_k+'_intermediate'
                    inter_fn=basefilename_k
                    pw.write_plink(G,F,intermediate_file_name=inter_fn,final_file_name=basefilename_k,make_bed=make_bed)
            else:
                snp_names = ['snp_%d' % i for i in range(X.shape[1])]
                G = pandas.DataFrame(data=X, index=sample_names, columns=snp_names)
                #inter_fn=basefilename+'_intermediate'
                inter_fn=basefilename
                pw.write_plink(G,F,intermediate_file_name=inter_fn,final_file_name=basefilename,make_bed=make_bed)
        pheno_array = np.empty((Y.shape[0],Y.shape[1]+2),dtype='|S20')
        pheno_array[:,0]=family_names[:,0]
        pheno_array[:,1]=sample_names
        pheno_array[:,2:]=Y
        filename_phen = basefilename+'.phen.txt'
        np.savetxt(filename_phen, pheno_array, delimiter='\t', fmt='%s')
            
        
        

def scale_K(K, X, verbose=False, return_K = False):
    """scale covariance K such that it explains unit variance"""
    c = sp.sum((sp.eye(len(K)) - (1.0 / len(K)) * sp.ones(K.shape)) * sp.array(K))
    scalar = (len(K) - 1) / c
    if verbose:
        print 'Kinship scaled by: %0.4f' % scalar
        print scalar
    K = scalar * K
    X = X*np.sqrt(scalar)
    if return_K:
        return K
    return X

def standardize_freq(snps,freq):
    X=snps
    X=(X-2.0*freq)/np.sqrt(2.0*freq*(1.0-freq))
    return X

class gwas_dataset(object):
    '''
    a dataset containing SNPs, phenotypes, etc.
    '''
    def __init__(self, snps=None, pheno=None, covariates=None, transformation=None, pheno_transformed=None):
        self.snps = snps                                #can be a reader or a Matrix
        self.covariates = covariates                    
        self.pheno = pheno
        self.pheno_transformed = pheno_transformed
        self.transformation = transformation




def generate_data(options,args):
    num_snps = options.numSnps+options.csnps_hidden #csnps not rel

    #set random seed
    np.random.seed(options.randomseed)

    num_children = 10 #!!!cmk make this param  sib_per_family
    num_trios = int(options.numIndividuals*options.fracSibs/(2 * num_children)) #!!trio is a misnomer because mom+dad+10 kids
    num_samples = options.numIndividuals-options.numIndividuals*options.fracSibs
    num_causal_obs = options.csnps
    num_causal_hidden = options.csnps_hidden
    num_causal = num_causal_obs+num_causal_hidden
    assert num_causal<=num_snps,"num_causal<=num_snps"
    num_phenos = options.num_phen
    num_differentiated = np.array([int(options.numSnps*options.diff),num_causal_hidden])
    Fst = np.array([[options.fst,options.fst],[np.NaN,np.NaN]])
    
    quiet=False
    assert options.pop_perc<=1.0 and options.pop_perc>=0.0,"    assert options.pop_prec<=1.0 and options.pop_perc>=0.0"
    pop_perc = np.array([options.pop_perc, 1.0-options.pop_perc])
    #prevalence = 0.3
    #transform = 'liability'
    maf=options.minFreq
    transform_param = None
    
    perc_causal_differentiated= np.array([0,0]) #!!!cmk remove options.diffCause*(1.0*num_causal_obs)/num_causal,(1.0*num_causal_hidden)/num_causal])
   
    assert options.h2<=1.0 and options.h2 >=0.0,"assert h2<=1.0 and hidden_var >=0.0"
    if num_causal_hidden==0:
        options.var_hidden=0.0
    assert options.var_hidden<1.0 and options.var_hidden >=0.0,"assert hiddenvar<1.0 and hidden_var >=0.0"
    noise_var = (1.0-options.var_hidden)*(1.0-options.h2)

    

    num_trios_pop= pop_perc*num_trios
    #new snp simulator
    simsnps = sim_snps(num_snps, num_differentiated = num_differentiated, MAF_ancestral = np.array([maf,0.5]), Fst = Fst, quiet = quiet, p_ancestral=None)
    
    
    simphen=sim_pheno(num_causal=num_causal, i_causal=None, W=None, W_covariates=None, noise_var=noise_var, genetic_var=1.0, num_phenotypes=num_phenos, perc_causal_differentiated=perc_causal_differentiated, weight_distribution='normal', quiet=quiet, covariates_var=0.3, sim_snps=simsnps)

    i_causal_diff = simsnps.differentiated[:,simphen.i_causal]
    

    #SNP indicators for output writing
    i_SNPs={}
    #generate vector of causal observed SNPs
    i_SNPs['causal_obs'] = simphen.i_causal.copy()
    i_SNPs['causal_obs'][simsnps.differentiated[1]]=False
    #generate indicator vector of causal hidden SNPs (hidden confounders)
    if 1: #!!!cmk remove if
        i_SNPs['causal_hidden'] = simsnps.differentiated[1].copy()
    if 0:
        i_snps['observed'] = ~i_SNPs['causal_hidden']
    #generate vector of non-causal SNPs
    i_SNPs['noncausal']=~simphen.i_causal
    #i_SNPs['differentiated']=simsnps.differentiated
    
    #gen_var = np.array([genetic_var/num_causal_obs, hidden_var/num_causal_hidden])
    gen_var = np.zeros(num_causal)
    gen_var[~i_causal_diff.any(0)]=0#!!!np.sqrt(options.h2/num_causal_obs)
    gen_var[i_causal_diff[0]]=0#!!!np.sqrt(options.h2/num_causal_obs)
    if num_causal_hidden>0:
        gen_var[i_causal_diff[1]]=0#!!!np.sqrt(options.var_hidden*(1.0-options.h2)/num_causal_hidden)
    simphen.W=(simphen.W.T*gen_var).T


    snps_pop=[]
    i_parent_pop=[]
    nonchild_index_list = []
    nonchild_start = 0
    for i_pop in xrange(len(pop_perc)):
        #import pdb; pdb.set_trace()
        snps=simsnps.generate_snps(int(num_samples*pop_perc[i_pop]),population_index = i_pop, snp_index = None)
        nonchild_index_list = nonchild_index_list + range(nonchild_start,nonchild_start+len(snps))
        #import pdb; pdb.set_trace()
        snps_kids,i_parent = simsnps.generate_trios(snps_parents=snps, num_trios=num_trios_pop[i_pop], population_percentages=None, snp_index=None, num_children_per_couple=num_children)
        nonchild_start += len(snps) + len(snps_kids)
        snps_pop.append(np.concatenate([snps,snps_kids],0))
        i_parent_pop.append(i_parent)
    #import pdb; pdb.set_trace()
    snps_all = np.concatenate(snps_pop,0)
    
    snps_kids,i_parent = simsnps.generate_trios(snps_parents=snps_all, num_trios=num_trios, population_percentages=None, snp_index=None, num_children_per_couple=num_children)
    snps_all = np.concatenate([snps_all,snps_kids],0)
    Y,Z=simphen.generate_phenotype(X=snps_all,freq=simsnps.p_ancestral,make_binary=options.make_binary)
    return snps_all,Y,simsnps,simphen,i_SNPs,nonchild_index_list
    #Y_bin, transform_param = simphen.get_phenotype_transform(Y=Y, pheno_transform = 'liability', transform_param = transform_param, prevalence = prevalence, population_percentages=np.array(pop_perc))

def run_select(options,eigvecs_fn,strategy):
    res = []
    import fastlmm.feature_selection.feature_selection_example as fs
    #clusterize this?
    strategy_=strategy
    k_values=None

    if strategy=='lmm_ins':
        k_values=[346456734]
        strategy_='insample_cv'
    elif strategy=='lmm_oos':
        k_values=[346456734]
        strategy_='lmm_full_cv'
    elif strategy=='lin':
        k_values=[0]
        strategy_='insample_cv'

    for i_pc,pc_fn in enumerate(eigvecs_fn): 
        cov_fn = pc_fn
        output_basefilename = "%s_%s"%(strategy,pc_fn)
        output_prefix = os.path.join(options.outdir,output_basefilename)
        debug_prefix = os.path.join(options.outdir, strategy+'_%s'%options.outdir)

        if not os.path.exists(debug_prefix):
            os.makedirs(debug_prefix)

        picklefile = output_prefix+".pickle"
        
        if (options.recompute) or (not os.path.exists(output_prefix+".pickle")): 
            myres = fs.runselect(bed_fn=options.bfile,
                                 pheno_fn=options.bfile+"."+options.suff_pheno+".txt",
                                 strategy=strategy_, select_by_ll=True,
                                 output_prefix=output_prefix, num_pcs=0, 
                                 random_state=3, num_snps_in_memory=5000, cov_fn=cov_fn,
                                 k_values=k_values, delta_values=None,num_folds=options.num_folds,
                                 penalty=options.penalty)
            myres['output_prefix']=output_prefix
            res.append(myres)
            cPickle.dump(res[i_pc],open(picklefile,'wb'))
        else:
            #print "WARNING: stuff disabled"
            res.append(cPickle.load(open(picklefile,'r')))
    return res

def plot_sim(snps,simsnps,Y):
        snps_std = standardize_freq(snps,simsnps.p_ancestral)
        K = snps_std.dot(snps_std.T)/snps_std.shape[1]
    
        import pylab as pl
        pl.ion()
        pl.figure();pl.imshow(K)
        pl.figure();pl.imshow(Y.dot(Y.T))
        #pl.figure();pl.imshow(Y_bin.dot(Y_bin.T))
def compute_pcfiles(options,recomputePC=False):
    numpc = np.arange(3)
    eigvecs_fn = []
    exists = True
    for i_pc,pc in enumerate(numpc):
        from fastlmm.util.computePC import getEigvecs_fn
        eigvecs_fn.append(getEigvecs_fn(options.bfile,pc))
        exists=exists and os.path.exists(eigvecs_fn[i_pc])# and os.path.exists("%s.pca.eigvals.%i.txt"%(options.bfile,pc))
    if (options.recompute) or (not exists) or recomputePC:
        from fastlmm.util.computePC import computePC
        computePC(file=options.bfile, filepath = None, numpc = numpc)
    else:
        print "PCs already computed"
    return eigvecs_fn

#!!move this to more central place
def compute_auto_pcs(snpreader, useChrom, cutoff, pc_file_name=None, k_values = np.arange(0,5+1,1)):
    snpMatrix = snpreader.read()
    snps = snpMatrix['snps']
    import fastlmm.util.standardizer as stdizer
    snps = stdizer.Unit().standardize(snps)
    
    #use vertex cut to find just parents
    rrm = snps.dot(snps.T) / snps.shape[1]
    import simulation.VertexCut as vc
    removeList = vc.VertexCut().work(rrm,cutoff) #These are the indexes of the IIDs to remove
    logging.info("removing {0} of {1} iids".format(len(removeList), rrm.shape[0]))
    keepList = [x for x in range(rrm.shape[0]) if x not in set(removeList)]
    nofam = snps[keepList,:]
    
    #learn # of pcs and generate on nonchild view of data
    import simulation.mainpca as mainpca
    geno = {'snps':nofam,'pos':snpMatrix['pos']}
    bestNumPCs = mainpca.FindBestPCCount(geno, predictSNPs=True,useChrom=useChrom, k_values=k_values) # run the search
    logging.info('best num PCs: {0}'.format(bestNumPCs))
    
    logging.info("computing svd...")
    Xnofam = nofam #in the original order
    import time
    t0 = time.time()
    import scipy.linalg as la
    Utr,Str,Vtr = la.svd(Xnofam, full_matrices=False)
    t1 = time.time()
    logging.info("done after %.4f seconds" % (t1 - t0))
    from simulation.pca import PCA, fast_logdet
    pca = PCA(n_components = bestNumPCs)
    pca._fit(Xnofam,Utr,Str,Vtr)
    
    #apply those pcs to all the data (i.e. transform)
    Xfam = snps
    logging.info('Projecting individuals to PCs space...')
    X_fit = pca.transform(Xfam)
    
    #write the file out
    if pc_file_name is not None:
        logging.info('writing results to file...')
        with open(pc_file_name, 'w') as f:
            iid = snpMatrix['iid']
            for i_ind in xrange(X_fit.shape[0]):
                f.write("{0} {1} ".format(iid[i_ind][0],iid[i_ind][1]))
                f.write(' '.join([str(pc) for pc in X_fit[i_ind, :]]))
                f.write('\n')

    return bestNumPCs, X_fit

def compute_pcfiles_auto(options,recomputePC=False): #!!
    from fastlmm.util.computePC import getEigvecs_fn
    fn = "%s_pc%s.vecs"%(options.bfile,"auto")
    if (not options.recompute) and (os.path.exists(fn)) and not recomputePC:
        return [fn]

    #input is options.bfile
    import fastlmm.pyplink.snpreader as sr
    snpreader = sr.Bed(options.bfile)
    bestNumPCs, X_fit = compute_auto_pcs(snpreader, False, options.vertex_cutoff, fn)


    # return its name in a singlton list
    return [fn]

    #old code...
    #eigvecs_fn = []
    ##exists = True
    #from fastlmm.util.computePC import getEigvecs_fn
    #eigvecs_fn.append(getEigvecs_fn(options.bfile,pc))
    #exists=exists and os.path.exists(eigvecs_fn[i_pc])# and os.path.exists("%s.pca.eigvals.%i.txt"%(options.bfile,pc))
    #if (options.recompute) or (not exists) or recomputePC:
    #    from fastlmm.util.computePC import computePC
    #    computePC(file=options.bfile, filepath = None, numpc = numpc)
    #else:
    #    print "PCs already computed"
    #return eigvecs_fn


def parseArgs():
    parser = OptionParser()
    parser.add_option('--bfile', metavar='bfile', help='output bfile name (default = out)', default="")
    parser.add_option('--csnps', metavar='csnps', type=int, default=100, help='number of causal SNPs')
    parser.add_option('--csnps_hidden', metavar='csnps_hidden', type=int, default=100, help='number of hidden causal SNPs')
    parser.add_option('--fst', metavar='fst',  type=float, default=0.025, help='Fst distance between the two populations at observed SNPs')
    parser.add_option('--h2', metavar='h2', type=float, default=0.5, help='trait heritability')
    parser.add_option('--var_hidden', metavar='var_hidden', type=float, default=0.5, help='fraction of noise variance that is due to hidden SNPs')
    parser.add_option('--numSnps', metavar='numSnps', type=int, default=50000, help='number of observed SNPs')
    parser.add_option('--numIndividuals', metavar='numIndividuals', type=int, default=2000, help='number of individuals')
    parser.add_option('--minFreq', metavar='minFreq', type=float, default=0.1, help='minimum minor allele frequency')
    parser.add_option('--fracSibs', metavar='fracSibs', type=float, default=0.3, help='fraction of siblings in the data')
    parser.add_option('--diffCause', metavar='diffCause', type=float, default=1.0, help='fraction of observed causal SNPs that are differented: 1.0 if all causal SNPs are differentiated, 0.0 if none')
    parser.add_option('--diff', metavar='diff', type=float, default=1.0, help='fraction of observed (causal or not) SNPs that are differented: 1.0 if all SNPs are differentiated, 0.0 if none')
    parser.add_option('--pop_perc', metavar='pop_perc', type=float, default=0.5, help='fraction of individuals from population 1 (rest 2)')
    parser.add_option('--num_phen', metavar='num_phen', type=int, default=1, help='number of phenotypes being generated')
    parser.add_option('--recompute', metavar='recompute', type=int, default=0, help='If intermediate files are present do not recompute previous steps' )
    parser.add_option('--suff_pheno',metavar='suff_pheno',type=str,default='phen',help='suffix of the pheno file')
    parser.add_option('--suff_pickle',metavar='suff_pickle',type=str,default='pickle',help='suffix of the pickle file')
    parser.add_option('--randomseed',metavar='randomseed',type=int,default=1,help='random seed')
    parser.add_option('--outdir',metavar='outdir',default="out",help="output directory")
    parser.add_option('--short_fn', metavar='short_fn', type=int, default=0, help='generate a random unique filename')
    parser.add_option('--num_folds', metavar='num_folds', type=int, default=10, help='number of train test splits for feature selection')
    parser.add_option('--penalty', metavar='penalty', type=float, default=0.0, help='regression penalty on fixed effects (e.g. PCs)')
    parser.add_option('--make_binary', metavar='make_binary', type=int, default=0, help='Should the phenotype be boolean?')
    parser.add_option('--vertex_cutoff', metavar='vertex_cutoff', type=float, default=0.1, help='Cutoff removing related iids')
    (options, args) = parser.parse_args()
    return (options,args)

def estimate_lambda(pv):
    """estimate lambda form a set of PV"""
    LOD2 = sp.median(st.chi2.isf(pv,1))
    L = (LOD2/0.456)
    return (L)

def sim_main(options, args):
    util.create_directory_if_necessary(options.outdir)
    if (len(options.bfile)==0):
        if options.short_fn==0:
            options.bfile = os.path.join(options.outdir,options.bfile)
            #make unique bfile
            options.bfile = "N%i_snps%i_causal%i_h2%.2f_sibs%.2f_pop%.2f_Fst%.4f_FstH%.4f_varh%.2f_seed%i"%(options.numIndividuals,options.numSnps,options.csnps,options.h2,options.fracSibs,options.pop_perc,options.fst,options.fst_hidden,options.var_hidden,options.randomseed)
        else:
            options.bfile = "N%iS%ic%ih%.2fs%.2fp%.2fF%.4fFH%.4fv%.2f_%i"%(options.numIndividuals,options.numSnps,options.csnps,options.h2,options.fracSibs,options.pop_perc,options.fst,options.fst_hidden,options.var_hidden,options.randomseed)
    recomputePC=False
    
    if (not options.recompute) and os.path.exists(options.bfile+".bed") and os.path.exists(options.bfile+".bim") and os.path.exists(options.bfile+".fam"):
        logging.info( "no sim needed")
        # make sure bed files are not corrupted
        try:
            import fastlmm.pyplink.snpreader as sr
            snp_reader = sr.Bed(options.bfile)
            tmp_data = snp_reader.read()
            
        except Exception, detail:
            print "problem reading bed file", options.bfile
            extensions = [".bed", ".fam", ".bim", ".ped", ".map"]
            # clean things up
            for ext in extensions:
                tmp_fn = options.bfile + ext
                if os.path.exists(tmp_fn):
                    os.remove(tmp_fn)
            raise
    else:
        logging.info( "sim needed. options.recompute={0}. os.path.exists({1})={2}. .path.exists({3})={4} os.path.exists({5})={6}".format(
            options.recompute, options.bfile+".bed", os.path.exists(options.bfile+".bed"),
            options.bfile+".bim", os.path.exists(options.bfile+".bim"),  options.bfile+".fam", os.path.exists(options.bfile+".fam")))
        snps,Y,simsnps,simphen,i_SNPs,nonchild_index_list = generate_data(options,args)
            #write out all SNPs in different files
        #pdb.set_trace()
        #plot simulations
        if 0: plot_sim(snps,simsnps,Y)
        
        if 0: simphen.write_plink(snps,Y,options.bfile,i_SNPs,make_bed=True)
    
        #write out all observed SNPs
        if 1:
            i_SNPs_obs={}
            i_SNPs_obs['causal']=i_SNPs['causal_obs'][~i_SNPs['causal_hidden']]
            snps_obs=snps[:,~i_SNPs['causal_hidden']]
            simphen.write_plink(snps_obs,Y,options.bfile,i_SNPs_obs,make_bed=True,merge_obs=True)
            
            # make sure bed files are not corrupted
            try:
                import fastlmm.pyplink.snpreader as sr
                snp_reader = sr.Bed(options.bfile)
                tmp_data = snp_reader.read()
                
            except Exception, detail:
                print "problem reading bed file", options.bfile
                extensions = [".bed", ".fam", ".bim", ".ped", ".map"]
                # clean things up
                for ext in extensions:
                    tmp_fn = options.bfile + ext
                    if os.path.exists(tmp_fn):
                        os.remove(tmp_fn)
                raise
                
        snps = None #save memory by thrwoing away snps
        #sim['causal']=i_SNPs_obs['causal']


    #compute PCs
    #eigvecs_fn = compute_pcfiles(options, recomputePC) #!!
    eigvecs_fn = compute_pcfiles_auto(options,recomputePC) #!!

    #run selection
    res_select_oosample=[]
    res_select_insample=[]
    res_select_lmm_oosample=[]
    res_select_lmm_insample=[]
    res_select_linreg=[]
    if 1: #!!
        strategy = "lmm_full_cv"
        res_select_oosample += run_select(options,eigvecs_fn,strategy)
    if 1: #!!
        strategy = "insample_cv"
        res_select_insample += run_select(options,eigvecs_fn,strategy)
    if 1: #!!
        strategy = "lmm_oos"
        res_select_lmm_oosample += run_select(options,eigvecs_fn,strategy)    
    if 1: #!!
        strategy = "lmm_ins"
        res_select_lmm_insample += run_select(options,eigvecs_fn,strategy)    
    if 1: #!!
        strategy = "lin"
        res_select_linreg += run_select(options,eigvecs_fn,strategy)    
    fastlmm_runs = []
    outfiles = []
    #run fastlmmc insample
    
    runindiv = 1
    index_out = []
    pc_file_out = []
    objective_sel=[]
    k_sel=[]
    for i_pc, pc_file in enumerate(eigvecs_fn):

        if "_pcauto." in pc_file:
            i_pcString = "auto"
        else:
            continue #ckw: this is for the revised version of the experiments that doesn't use PCp
            i_pcString = str(i_pc)

        fastlmmcstr = fc.executable+" -mpheno 1 -bfile "+options.bfile+" -pheno "+options.bfile+".phen.txt"+" -covar "+pc_file
        
        if 1:#adding true causals
            
            # we just added output files
            outfile = pc_file+".lmm_true_causals"
            if options.recompute or not os.path.exists(outfile):

                from fastlmm.pyplink.snpreader.Bed import Bed
                snp_reader = Bed(options.bfile)
                snp_reader.run_once()
                extract_fn = options.bfile + ".extract"
                extract_file = open(extract_fn, "w")
                for snp_rs in [c for c in snp_reader.rs if "_c" in c]:
                    extract_file.write(snp_rs + "\n")
                extract_file.close()
                print "extract_fn", extract_fn
                fstr = fc.executable+" -mpheno 1 -bfile "+options.bfile+" -pheno "+options.bfile+".phen.txt" + " -excludebygeneticdistance 0 -bfilesim "+options.bfile+" -out "+outfile + " -extractsim "+extract_fn
                
                print "(lmm_fs_cond_full) file does not exist: %s" % (outfile)
                #assert os.path.exists(outfile), "(lmm_true_causals) file does not exist: %s" % (outfile)
                os.system(fstr)

            index_out.append('lmm_true_causals_%s'%i_pcString)
            pc_file_out.append(pc_file)
            outfiles.append(outfile)
            
            objective_sel.append(0.0)
            k_sel.append(0)
            #objective_sel.append(res_select_linreg[i_pc]['best_obj'])
            #k_sel.append(res_select_linreg[i_pc]['best_k'])
            
        if 1:#the new run
            
            # we just added output files
            outfile = pc_file+".lmm_fs_cond_full"
            if options.recompute or not os.path.exists(outfile):
                print "(lmm_fs_cond_full) file does not exist: %s" % (outfile)
                from synthetic_2k import compute_fs, compute_gwas

                fs_fn = pc_file+".meta.pickle"
                phen_fn = options.bfile+".phen.txt"

                input_tuple_fs = options.bfile, phen_fn, pc_file, fs_fn                
                compute_fs(input_tuple_fs)

                tmp_dir = "D:/%s/%s/" % (os.environ["userdomain"], os.environ["username"])
                tmp_fn = tmp_dir + pc_file + "_big_tmp_file"
                input_tuple_gwas = options.bfile, fs_fn, pc_file, phen_fn, tmp_dir, tmp_fn
                compute_gwas(input_tuple_gwas)

            index_out.append('lmm_fs_cond_full_%s'%i_pcString)
            pc_file_out.append(pc_file)
            outfiles.append(outfile)
            
            objective_sel.append(0.0)
            k_sel.append(0)
                #objective_sel.append(res_select_linreg[i_pc]['best_obj'])
                #k_sel.append(res_select_linreg[i_pc]['best_k'])
            
        if 0:#the new run (subsampled)
            
            # we just added output files
            for subsample in [2,8]:
                outfile = pc_file+".lmm_fs_cond_full_%i" % (subsample)
                assert os.path.exists(outfile), "(lmm_fs_cond_full) file does not exist: %s" % (outfile)
                index_out.append('lmm_fs_cond_full_%s_%i'% (i_pcString, subsample))
                pc_file_out.append(pc_file)
                outfiles.append(outfile)
            
                objective_sel.append(0.0)
                k_sel.append(0)
                #objective_sel.append(res_select_linreg[i_pc]['best_obj'])
                #k_sel.append(res_select_linreg[i_pc]['best_k'])
            

        if 1:#linreg !!
            fstr = fastlmmcstr
            outfile = pc_file+".lin"
            index_out.append('linreg_%s'%i_pcString)
            pc_file_out.append(pc_file)
            outfiles.append(outfile)
            fstr = fstr + " -linreg -out "+outfile
            #objective_sel.append(0.0)
            #k_sel.append(0)
            objective_sel.append(res_select_linreg[i_pc]['best_obj'])
            k_sel.append(res_select_linreg[i_pc]['best_k'])
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)
        if 1:#full insample  !!
            fstr = fastlmmcstr
            outfile = pc_file+".lmm_ins"
            fstr = fstr + " -excludebygeneticdistance 0 -bfilesim "+options.bfile+" -out "+outfile
            outfiles.append(outfile)
            index_out.append('lmm_ins_%s'%i_pcString)
            pc_file_out.append(pc_file)
            objective_sel.append(res_select_lmm_insample[i_pc]['best_obj'])
            k_sel.append(res_select_lmm_insample[i_pc]['best_k'])
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)    #run fastlmmc oos
        if 1:#full oos  !!
            fstr = fastlmmcstr
            outfile = pc_file+".lmm_oos"
            fstr = fstr + " -excludebygeneticdistance 0 -bfilesim "+options.bfile+" -out "+outfile
            fstr = fstr + " -logdelta %.4f"%np.log(res_select_lmm_oosample[i_pc]['best_delta'])
            outfiles.append(outfile)
            index_out.append('lmm_oos_%s'%i_pcString)
            pc_file_out.append(pc_file)
            objective_sel.append(res_select_lmm_oosample[i_pc]['best_obj'])
            k_sel.append(res_select_lmm_oosample[i_pc]['best_k'])
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)    #run fastlmmc oos
        if 1:#insample !!
            fstr = fastlmmcstr
            outfile = pc_file+".ins"
            outfiles.append(outfile)
            index_out.append('insample_%s'%i_pcString)
            pc_file_out.append(pc_file)
            objective_sel.append(res_select_insample[i_pc]['best_obj'])
            k_sel.append(res_select_insample[i_pc]['best_k'])
            fstr = fstr + " -out "+outfile
            if res_select_insample[i_pc]['best_k']<=0:
                fstr =fstr+' -linreg'
            else:
                if res_select_insample[i_pc]['best_k']<=10000:
                    extractfile = res_select_insample[i_pc]['output_prefix']+'_best_snps.txt'
                    np.savetxt(extractfile,res_select_insample[i_pc]['best_snps'],fmt='%s')
                    fstr = fstr+' -extractsim '+extractfile
                fstr = fstr + " -bfilesim "+options.bfile
                fstr = fstr + " -excludebygeneticdistance 0"
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)
        if 1:#oosample #!!
            fstr = fastlmmcstr
            index_out.append('oosample_%s'%i_pcString)
            pc_file_out.append(pc_file)
            objective_sel.append(res_select_oosample[i_pc]['best_obj'])
            k_sel.append(res_select_oosample[i_pc]['best_k'])
            outfile = pc_file+".oos"
            outfiles.append(outfile)
            fstr = fstr + " -out "+outfile
            if res_select_oosample[i_pc]['best_k']<=0:
                fstr =fstr+' -linreg'
            else:
                if res_select_oosample[i_pc]['best_k']<=10000:
                    extractfile = res_select_oosample[i_pc]['output_prefix']+'_best_snps.txt'
                    np.savetxt(extractfile,res_select_oosample[i_pc]['best_snps'],fmt='%s')
                    fstr = fstr+' -extractsim '+extractfile
                fstr = fstr + " -logdelta %.4f"%np.log(res_select_oosample[i_pc]['best_delta'])
                fstr = fstr + " -bfilesim "+options.bfile
                fstr = fstr + " -excludebygeneticdistance 0"
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)

        if 0:#true causals
            from fastlmm.pyplink.snpreader.Bed import Bed
            snp_reader = Bed(options.bfile)
            snp_reader.run_once()
            extract_fn = options.bfile + ".extract"
            extract_file = open(extract_fn, "w")
            for snp_rs in [c for c in snp_reader.rs if "_c" in c]:
                extract_file.write(snp_rs + "\n")
            extract_file.close()
            print "extract_fn", extract_fn
            fstr = fastlmmcstr
            outfile = pc_file+".lmm_true_causals"
            fstr = fstr + " -excludebygeneticdistance 0 -bfilesim "+options.bfile+" -out "+outfile + " -extractsim "+extract_fn
            outfiles.append(outfile)
            index_out.append('lmm_true_causals_%s'%i_pcString)
            pc_file_out.append(pc_file)
            objective_sel.append(0.0)
            k_sel.append(0)
            if options.recompute or not os.path.exists(outfile):
                if runindiv:
                    os.system(fstr)
                else:
                    fastlmm_runs.append(fstr)    #run fastlmmc oos

    #run fastlmmc
    if not runindiv:
        for run in os.system(fastlmm_runs):
            os.system(run)

    pv_thresholds = [1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8]
    sim_params = ['h2', 'h2 hidden', '% trios', '#causal', 'Fst']
    cols = sum([sim_params, ['Power @ %.1e' % i for i in pv_thresholds], ['Type I @ %.1e' % i for i in pv_thresholds],
                ['lambda', 'k', 'obj'],['Power corr @ %.1e' % i for i in pv_thresholds], 
                ['Type I corr @ %.1e' % i for i in pv_thresholds], ['numPCs']], [])
    results_df = pandas.DataFrame(index=np.unique(index_out), columns=cols)
    
    #TODO get objectives for 0 and all SNPs

    #load output
    for j, outputfile in enumerate(outfiles):
        logging.info("enumerating {0} of {1}, outfile {2}".format(j,len(outfiles),outputfile))
        #try:
        output = np.loadtxt(outputfile,dtype = str,usecols=(0,1,2,3,4),comments=None) #The usecols is needed to avoid thanye #INDs in other columns that cause np.loadtext to fail
        #except exception:
        #    logging.critical("bad file: " + outputfile)
        #    os.rename(outputfile,outputfile+".bad.txt")
        #    continue
        #import pdb; pdb.set_trace()
        results = pandas.DataFrame(data=np.array(output[1:,1:],dtype=float), index=output[1:,0], columns=output[0,1:])
        results = results.sort(columns='Pvalue', ascending=True)
        i_causal = np.zeros(len(results.index),dtype=bool)
        for i,snp in enumerate(results.index):
            i_causal[i]=snp[-2:]=='_c'
        
        lambda_gc  = estimate_lambda(results['Pvalue'].flatten())
        results_df.T[index_out[j]]['lambda'] = lambda_gc
        #import pdb; pdb.set_trace()
        print "INDEX:", index_out[j]
        n_causal = i_causal.sum()
        #compute power and Type-1 error
        from sklearn.metrics import roc_curve
        fpr, tpr, thresholds = roc_curve(i_causal, results['Pvalue'].flatten())

        if 1:
            power = np.zeros_like(pv_thresholds)
            t1err = np.zeros_like(pv_thresholds)
            power_corr = np.zeros_like(pv_thresholds)
            t1err_corr = np.zeros_like(pv_thresholds)
            pvcorr=results['Pvalue']
            pvcorr=st.chi2.sf(st.chi2.isf(pvcorr,1)/lambda_gc,1)
            for i_t, t in enumerate(pv_thresholds):
                #compute uncorrected power and T1
                i_lower = results['Pvalue']<t
                power[i_t] =  i_causal[i_lower].sum()/(1.0*(n_causal))
                t1err[i_t] = (~i_causal[i_lower]).sum()/(1.0*(len(i_causal)-n_causal))
                results_df.T[index_out[j]]['Power @ %.1e' % t] = power[i_t]
                results_df.T[index_out[j]]['Type I @ %.1e' % t] = t1err[i_t]
                #compute GC corrected Power and T1
                i_lower_corr = pvcorr<t
                power_corr[i_t] =  i_causal[i_lower_corr].sum()/(1.0*(n_causal))
                t1err_corr[i_t] = (~i_causal[i_lower_corr]).sum()/(1.0*(len(i_causal)-n_causal))
                results_df.T[index_out[j]]['Power corr @ %.1e' % t] = power_corr[i_t]
                results_df.T[index_out[j]]['Type I corr @ %.1e' % t] = t1err_corr[i_t]

        results_df.T[index_out[j]]['h2'] = options.h2
        results_df.T[index_out[j]]['h2 hidden'] = options.var_hidden
        results_df.T[index_out[j]]['% trios'] = options.fracSibs
        results_df.T[index_out[j]]['#causal'] = options.csnps
        results_df.T[index_out[j]]['Fst'] = options.fst
        results_df.T[index_out[j]]['k'] = k_sel[j]
        results_df.T[index_out[j]]['obj'] = objective_sel[j]

        # read one line of the file, to get the # of pcs
        pc_file = pc_file_out[j]
        with open(pc_file, 'r') as f:
            first_line = f.readline()
            numPCs = len(first_line.split(" ")) - 2 

        results_df.T[index_out[j]]['numPCs'] = numPCs

        if 0:
            import pylab as pl
            pl.ion()
            pl.figure()
            pl.plot(pv_thresholds,power)
            pl.title("power")
            pl.figure()
            pl.plot(pv_thresholds,t1err)
            pl.title("t1 error")
    return results_df

if __name__ == '__main__':
    (options,args) = parseArgs()
    #testing options
    #perform simulations
    results_df = sim_main(options)


