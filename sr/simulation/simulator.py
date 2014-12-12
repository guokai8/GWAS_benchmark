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
    
    
    #!!!cmk may be important only because it uses the random number generator
    state = np.random.get_state()
    simphen=sim_pheno(num_causal=num_causal, i_causal=None, W=None, W_covariates=None, noise_var=noise_var, genetic_var=1.0, num_phenotypes=num_phenos, perc_causal_differentiated=perc_causal_differentiated, weight_distribution='normal', quiet=quiet, covariates_var=0.3, sim_snps=simsnps)
    np.random.set_state(state)


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
    return snps_all




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

