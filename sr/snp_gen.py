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

    np.random.seed(seed)

    iid_solo_count = iid_count-iid_count*dfr
    family_count = int(iid_count*dfr/(2 * sibs_per_family))

    ancestral = np.random.uniform(maf_low, maf_high, sid_count)     #sample ancestral allele frequencies

    snp_list=[]
    for population_index, freq_pop in enumerate([freq_pop_0, 1.0-freq_pop_0]):
        logging.info("Simulating SNPs from a population %i" % population_index)
        snps_parents=_generate_snps(ancestral, fst, int(iid_solo_count*freq_pop), sid_count)
        snp_list.append(snps_parents)

        snp_list.append(_generate_kids(parent_snps=snps_parents, family_count=int(freq_pop*family_count), sibs_per_family=sibs_per_family))

    snp_list.append(_generate_kids(parent_snps=np.concatenate(snp_list), family_count=family_count, sibs_per_family=sibs_per_family))
    val = np.concatenate(snp_list)

    iid = np.array([["i_{0}".format(iid_index),"f_{0}".format(iid_index)] for iid_index in xrange(val.shape[0])])
    sid = np.array(["snp_{0}".format(sid_index) for sid_index in xrange(val.shape[1])])
    pos = np.array(list([sid_index,0,0] for sid_index in xrange(len(sid)))) # every snp has position 0,0 on its own chrom

    snpdata = SnpData(iid, sid, pos, val, 
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4})".format(fst, dfr, iid_count, sid_count, maf_low) #!!!cmk is this up-to-date?
                      )

    return snpdata


def _generate_snps(ancestral, fst, sample_size, sid_count):
    """
    Generates genotypes with a certain MAF and optionally with population structure.
    In case of no population structure, they are sampled from a binomial,
    otherwise from a Beta-Binomial (Balding and Nichols, 1995).
    """
    if fst == 0.0: 
        alpha = ancestral #special treatment if no population structure
    else:
        alpha = np.random.beta(ancestral*(1.0-fst)/fst,(1.0-ancestral)*(1.0-fst)/fst, sid_count)


    #generate from population frequencies    
    snps = np.zeros((sample_size,sid_count),dtype='int8') #.zeros not .empty because will be adding 1's to it
    for i in xrange(2): #"2" for diploid
        #sample each allele
        rand = np.random.random((sample_size,sid_count))
        snps[rand<alpha]+=1
    return snps


def _generate_kids(parent_snps, family_count, sibs_per_family): #!!!cmk should it be sibs, kids, or children
    '''
    generate a single set of family members
    '''    
    parent_count, sid_count = parent_snps.shape
    assert parent_count>=2*family_count, "sample_size>=2*family_count"


    parent_permutation = np.random.permutation(parent_count)
    snps = np.zeros((family_count*sibs_per_family,sid_count),dtype='int8')
    for copy_index in xrange(2):#"2" for diploid
        #sample each allele
        sample = parent_snps[parent_permutation[copy_index*family_count:(copy_index+1)*family_count],:]
        for kid_index in xrange(sibs_per_family):
            rand = np.random.random((family_count,sid_count))
            snps[kid_index*family_count:(kid_index+1)*family_count][rand<0.5*sample]+=1
    return snps


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    snpdata = snp_gen(fst=.1,dfr=.1,iid_count=100,sid_count=1000)
    print snpdata


    import doctest
    doctest.testmod()


    print "done"

