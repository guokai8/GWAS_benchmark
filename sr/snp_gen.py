import sys
import os
import numpy as np
import scipy as sp
import logging
import pysnptools.util as pstutil
from pysnptools.snpreader import SnpData

def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, maf_high=.5, seed=0,sibs_per_family=10,freq_pop_0=.5):
    """Generates a random :class:`.SnpData`

    :param fst: Degree of Population Structure, e.g. 0 (a special case), 0.005, 0.01, 0.05, 0.1
    :type fst: float

    :param dft: Degree of Family Relatedness, the fraction of individuals belonging to a family, ie. fracSibs, e.g. 0.0, 0.5, 0.6, 0.7, 0.8, 0.9
    :type dft: float

    :param iid_count: The number of individuals to generate. Because of rounding the actual number may be less.
    :type iid_count: int

    :param sid_count: The number of snps to generate.
    :type sid_count: int

    :param maf_low: (default .05) lower bound of uniformly-generated Minor allele frequency
    :type maf_low: float

    :param maf_high: (default .5) upper bound of uniformly-generated Minor allele frequency
    :type maf_high: float

    :param seed: (default 0) Random seed
    :type seed: int

    :param sibs_per_family: (default 10) number of siblings in each family
    :type sibs_per_family: int

    :param freq_pop_0: (default .5) Fraction of individuals in population 0 (the rest will be in population 1)
    :type freq_pop_0: float

    :rtype: :class:`.SnpData`

    :Example:

    >>> snpdata = snp_gen(fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
    >>> print snpdata.iid_count, snpdata.sid_count #because of rounding got 190 individuals
    190, 20
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
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4}, maf_high={5}, seed={6}, sibs_per_family={7}, freq_pop_0={8})"
                      .format(fst, dfr, iid_count, sid_count, maf_low, maf_high, seed, sibs_per_family, freq_pop_0)
                      )

    if snpdata.iid_count != iid_count:
        logging.warn("Because of rounding the actual number of iids is {0} rather than the requested {1}".format(snpdata.iid_count, iid_count))

    return snpdata


def _generate_snps(ancestral, fst, iid_count, sid_count):
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
    snps = np.zeros((iid_count,sid_count),dtype='int8') #.zeros not .empty because will be adding 1's to it
    for i in xrange(2): #"2" for diploid
        #sample each allele
        rand = np.random.random((iid_count,sid_count))
        snps[rand<alpha]+=1
    return snps


def _generate_kids(parent_snps, family_count, sibs_per_family): #!!! should it be sibs, kids, or children
    '''
    generate a single set of family members
    '''    
    parent_count, sid_count = parent_snps.shape
    assert parent_count>=2*family_count, "parent_count>=2*family_count"


    parent_permutation = np.random.permutation(parent_count)
    snps = np.zeros((family_count*sibs_per_family,sid_count),dtype=np.float64)
    for copy_index in xrange(2):#"2" for diploid
        sample = parent_snps[parent_permutation[copy_index*family_count:(copy_index+1)*family_count],:]         #sample each allele
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

