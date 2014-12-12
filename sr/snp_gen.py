import sys
import os
import numpy as np
import scipy as sp
import logging
import pysnptools.util as pstutil
from pysnptools.snpreader import SnpData

def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, maf_high=.5, seed=0,sibs_per_family=10,freq_pop_1=.5): #!!!cmk move this to the front
    """
    #!!!cmk fill in with docs including example

    fst is Degree of Population Structure, e.g. [0.005, 0.01, 0.05, 0.1] !!!cmk was self.FSTs
        fst=0 is a special case #!!!is it?
    dfr is Degree of Family Relatedness, the fraction of individuals belonging to a family [0.0, 0.5, 0.6, 0.7, 0.8, 0.9], e.g. #!!!cmk was self.fracSibs
    iid_count !!!cmk was numIndividuals
    sid_count == numSnps
    MAF is Minor allele frequency
    freq_pop_1 !!cmk was caseFrac ????

    """
    val = _generate_data(num_snps=sid_count,randomseed=seed,fracSibs=dfr,numIndividuals=iid_count,num_children=sibs_per_family,pop_perc=freq_pop_1,maf_low=maf_low,maf_high=maf_high,fst=fst)
    iid = np.array([["i_{0}".format(iid_index),"f_{0}".format(iid_index)] for iid_index in xrange(val.shape[0])])
    sid=np.array(["snp_{0}".format(sid_index) for sid_index in xrange(val.shape[1])])
    pos = np.array(list([sid_index,0,0] for sid_index in xrange(len(sid)))) # every snp has position 0,0 on its own chrom

    snpdata = SnpData(iid, sid, pos, val, 
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4})".format(fst, dfr, iid_count, sid_count, maf_low) #!!!cmk is this up-to-date?
                      )

    return snpdata

class _sim_snps(object):
    """
    Generates genotypes with a certain MAF and optionally with population structure.
    In case of no population structure, they are sampled from a binomial,
    otherwise from a Beta-Binomial (Balding and Nichols, 1995).
    """
    def __init__(self, num_snps, MAF_ancestral, fst):
        self.num_snps=num_snps
        self.MAF_ancestral = MAF_ancestral
        self.fst = fst                                  #dimensions: number populations times number SNP groups

        #set the allele frequencies
        self.sample_frequencies()



    def sample_frequencies(self):
        '''
        sample the allele frequencies for all ancestral SNPs and all populations
        ancestral frequencies are sampled uniformly according to MAFs
        population frequencies from a Beta distribution (Balding and Nichols, 1995).
        '''
        #sample ancestral allele frequencies, len(p)=num_snps
        p_ancestral = np.random.uniform(self.MAF_ancestral[0], self.MAF_ancestral[1], self.num_snps)
        
        #alphas: shape number(populations) times number(snps) and copy over ancestral frequencies
        self.alphas = np.zeros([2,p_ancestral.shape[0]]) #"2" is the number of populations
        self.alphas[:,:]=p_ancestral
        
        for i_population in xrange(2): #"2" is the number of populations
            if self.fst == 0.0: 
                self.alphas[i_population,:] = p_ancestral #special treatment if no population structure
            else:
                self.alphas[i_population,:] = np.random.beta(p_ancestral*(1.0-self.fst)/self.fst,(1.0-p_ancestral)*(1.0-self.fst)/self.fst, p_ancestral.shape[0])
        pass



    def generate_trios(self, snps_parents, num_trios=0, population_percentages=None, snp_index=None, num_children_per_couple=1):
        '''
        generate a single set of trios
        '''    
        num_trios = int(num_trios)
        sample_size = snps_parents.shape[0]
        num_snps = snps_parents.shape[1]
        assert sample_size>=2*num_trios, "sample_size>=2*num_trios"
        potential_parents = np.random.permutation(sample_size)
        i_done=0
        i_parent = np.empty([num_trios,2],dtype = 'int64') #"2" for diploid
        snps_parents_sampled = []
        for i in xrange(i_parent.shape[1]):
            #sample each allele
            i_parent[:,i] = potential_parents[i_done:i_done+num_trios]
            snps_parents_sampled.append(snps_parents[i_parent[:,i]])
            i_done+=num_trios
        snps = _mate(snps_parents_sampled, i_parent, num_children_per_couple)
        return snps, i_parent

    def generate_snps(self, sample_size, population_index):
        """
        Generates genotypes with a certain MAF and optionally with population structure.
        In case of no population structure, they are sampled from a binomial,
        otherwise from a Beta-Binomial (Balding and Nichols, 1995).
        """
        snp_index = np.arange(self.num_snps)
        logging.info("Simulating SNPs from population %i" % population_index)

        #generate from population frequencies    
        pgen = self.alphas[population_index,snp_index]

        snps = np.zeros((sample_size,pgen.shape[0]),dtype='int8')

        for i in xrange(2): #"2" for diploid
            #sample each allele
            rand = np.random.random((sample_size,pgen.shape[0]))
            snps[rand<pgen]+=1
        return snps

def _generate_data(num_snps, randomseed,fracSibs,numIndividuals,num_children,pop_perc,maf_low,maf_high,fst):
    
    #set random seed
    np.random.seed(randomseed) #!!!cmk how many places is the seed set?

    num_trios = int(numIndividuals*fracSibs/(2 * num_children)) #!!trio is a misnomer because mom+dad+10 kids
    num_samples = numIndividuals-numIndividuals*fracSibs
    
    assert 0 <= pop_perc and pop_perc <=1.0,"assert 0 <= pop_perc and pop_perc <=1.0"
    pop_perc = np.array([pop_perc, 1.0-pop_perc])
    

    num_trios_pop= pop_perc*num_trios
    simsnps = _sim_snps(num_snps, MAF_ancestral = np.array([maf_low,maf_high]), fst = fst)

    snps_pop=[]
    i_parent_pop=[]
    nonchild_index_list = []
    nonchild_start = 0
    for i_pop in xrange(2): #"2" is the number of populations
        snps=simsnps.generate_snps(int(num_samples*pop_perc[i_pop]),population_index = i_pop)
        nonchild_index_list = nonchild_index_list + range(nonchild_start,nonchild_start+len(snps))
        snps_kids,i_parent = simsnps.generate_trios(snps_parents=snps, num_trios=num_trios_pop[i_pop], population_percentages=None, snp_index=None, num_children_per_couple=num_children)
        nonchild_start += len(snps) + len(snps_kids)
        snps_pop.append(np.concatenate([snps,snps_kids],0))
        i_parent_pop.append(i_parent)
    snps_all = np.concatenate(snps_pop,0)
    
    snps_kids,i_parent = simsnps.generate_trios(snps_parents=snps_all, num_trios=num_trios, population_percentages=None, snp_index=None, num_children_per_couple=num_children)
    snps_all = np.concatenate([snps_all,snps_kids],0)
    return snps_all

def _mate(snps_parents, i_parent, num_children_per_couple):
    '''
    given a set of snps for the parents, mate the individuals
    '''
    num_snps = snps_parents[0].shape[1]
    num_trios = i_parent.shape[0]
    num_children = num_trios*num_children_per_couple #10 children per couple
    snps = np.zeros((num_children,num_snps),dtype='int8')
        
    for i in xrange(len(snps_parents)):
        #sample each allele
        for j in xrange(num_children_per_couple):
            rand = np.random.random((num_trios,num_snps))
            snps[j*num_trios:(j+1)*num_trios][rand<0.5*snps_parents[i]]+=1
    return snps

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    snpdata = snp_gen(fst=.1,dfr=.1,iid_count=100,sid_count=1000)
    print snpdata


    import doctest
    doctest.testmod()


    print "done"

