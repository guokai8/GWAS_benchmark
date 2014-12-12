import sys
import os
import numpy as np
import scipy as sp
import logging
import pysnptools.util as pstutil
from pysnptools.snpreader import SnpData

def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, maf_high=.5, seed=0,sibs_per_family=10,freq_pop_0=.5): #!!!cmk move this to the front
    """
    #!!!cmk fill in with docs including example

    fst is Degree of Population Structure, e.g. [0.005, 0.01, 0.05, 0.1] !!!cmk was self.FSTs
        fst=0 is a special case #!!!is it?
    dfr is Degree of Family Relatedness, the fraction of individuals belonging to a family [0.0, 0.5, 0.6, 0.7, 0.8, 0.9], e.g. #!!!cmk was self.fracSibs
    iid_count !!!cmk was iid_count
    sid_count == numSnps
    MAF is Minor allele frequency
    freq_pop_0 !!cmk was caseFrac ????

    """
    assert 0 <= freq_pop_0 and freq_pop_0 <=1.0,"assert 0 <= freq_pop_0 and freq_pop_0 <=1.0"
    freq_pops = np.array([freq_pop_0, 1.0-freq_pop_0])

    np.random.seed(seed)

    iid_solo_count = iid_count-iid_count*dfr
    family_count = int(iid_count*dfr/(2 * sibs_per_family)) #!!trio is a misnomer because mom+dad+10 kids

    alphas = _sample_frequencies(sid_count, maf_low, maf_high, fst = fst, population_count=len(freq_pops))

    snps_pop=[]
    i_parent_pop=[]
    nonchild_index_list = []
    nonchild_start = 0
    for population_index, freq_pop in enumerate(freq_pops): #"2" is the number of populations
        logging.info("Simulating SNPs from a population %i" % population_index)
        alpha = alphas[population_index,:]

        snps_parents=_generate_snps(alpha, int(iid_solo_count*freq_pop), sid_count)
        nonchild_index_list = nonchild_index_list + range(nonchild_start,nonchild_start+len(snps_parents))
        snps_kids,i_parent = _generate_family(snps_parents=snps_parents, family_count=int(freq_pop*family_count), num_children_per_couple=sibs_per_family)
        nonchild_start += len(snps_parents) + len(snps_kids)
        snps_pop.append(np.concatenate([snps_parents,snps_kids],0))
        i_parent_pop.append(i_parent)
    val = np.concatenate(snps_pop,0)
    
    snps_kids,i_parent = _generate_family(snps_parents=val, family_count=family_count, num_children_per_couple=sibs_per_family)
    val = np.concatenate([val,snps_kids],0)

    iid = np.array([["i_{0}".format(iid_index),"f_{0}".format(iid_index)] for iid_index in xrange(val.shape[0])])
    sid=np.array(["snp_{0}".format(sid_index) for sid_index in xrange(val.shape[1])])
    pos = np.array(list([sid_index,0,0] for sid_index in xrange(len(sid)))) # every snp has position 0,0 on its own chrom

    snpdata = SnpData(iid, sid, pos, val, 
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4})".format(fst, dfr, iid_count, sid_count, maf_low) #!!!cmk is this up-to-date?
                      )

    return snpdata


def _generate_snps(alpha, sample_size, sid_count):
    """
    Generates genotypes with a certain MAF and optionally with population structure.
    In case of no population structure, they are sampled from a binomial,
    otherwise from a Beta-Binomial (Balding and Nichols, 1995).
    """
    #generate from population frequencies    
    snps = np.zeros((sample_size,sid_count),dtype='int8')
    for i in xrange(2): #"2" for diploid
        #sample each allele
        rand = np.random.random((sample_size,sid_count))
        snps[rand<alpha]+=1
    return snps

def _generate_family(snps_parents, family_count, num_children_per_couple):
    '''
    generate a single set of family members
    '''    
    sample_size = snps_parents.shape[0]
    sid_count = snps_parents.shape[1]
    assert sample_size>=2*family_count, "sample_size>=2*family_count"
    potential_parents = np.random.permutation(sample_size)
    i_done=0
    i_parent = np.empty([family_count,2],dtype = 'int64') #"2" for diploid
    snps_parents_sampled = []
    for i in xrange(i_parent.shape[1]):
        #sample each allele
        i_parent[:,i] = potential_parents[i_done:i_done+family_count]
        snps_parents_sampled.append(snps_parents[i_parent[:,i]])
        i_done+=family_count
    snps = _mate(snps_parents_sampled, i_parent, num_children_per_couple)
    return snps, i_parent

def _sample_frequencies(sid_count, maf_low, maf_high, fst, population_count):
    '''
    sample the allele frequencies for all ancestral SNPs and all populations
    ancestral frequencies are sampled uniformly according to MAFs
    population frequencies from a Beta distribution (Balding and Nichols, 1995).
    '''
    #sample ancestral allele frequencies
    ancestral = np.random.uniform(maf_low, maf_high, sid_count)
        
    #alphas: shape number(populations) times number(snps) and copy over ancestral frequencies
    alphas = np.empty([population_count,sid_count])
    for population_index in xrange(population_count):
        if fst == 0.0: 
            alphas[population_index,:] = ancestral #special treatment if no population structure
        else:
            alphas[population_index,:] = np.random.beta(ancestral*(1.0-fst)/fst,(1.0-ancestral)*(1.0-fst)/fst, sid_count)
    return alphas

def _mate(snps_parents, i_parent, num_children_per_couple):
    '''
    given a set of snps for the parents, mate the individuals
    '''
    sid_count = snps_parents[0].shape[1]
    family_count = i_parent.shape[0]
    num_children = family_count*num_children_per_couple
    snps = np.zeros((num_children,sid_count),dtype='int8')
        
    for i in xrange(len(snps_parents)):
        #sample each allele
        for j in xrange(num_children_per_couple):
            rand = np.random.random((family_count,sid_count))
            snps[j*family_count:(j+1)*family_count][rand<0.5*snps_parents[i]]+=1
    return snps

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    snpdata = snp_gen(fst=.1,dfr=.1,iid_count=100,sid_count=1000)
    print snpdata


    import doctest
    doctest.testmod()


    print "done"

