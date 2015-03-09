import sys
import os
import numpy as np
import scipy as sp
import logging
import pysnptools.util as pstutil
from pysnptools.snpreader import SnpData, Bed, Dat



def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, maf_high=.5, seed=0,sibs_per_family=10,freq_pop_0=.5,chr_count=None):
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

    :param chr_count: (default one chromosome per SNP) Number of chromosomes to which SNPs should be assigned. The SNPs will
    be assigned as evenly as possible. Chromosome names are integers starting with 1. SNP positions within a chromosome are sequential
    integers starting with 1.
    :type chr_count: int

    :rtype: :class:`.SnpData`
    :Example:

    >>> snpdata = snp_gen(fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
    >>> print int(snpdata.iid_count), int(snpdata.sid_count) #because of rounding got 190 individuals
    190 20

    """
    assert 0 <= freq_pop_0 and freq_pop_0 <=1.0,"assert 0 <= freq_pop_0 and freq_pop_0 <=1.0"

    if seed is not None:
        np.random.seed(int(seed % sys.maxint))

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

    if chr_count is None:
        chr_count = len(sid)

    assert len(sid) == 0 or chr_count > 0, "chr_count must be at least 1 (unless sid_count is 0)"
    sid_per_chrom = int(sp.ceil(float(len(sid))/max(1,chr_count)))
    pos = np.array(list([1+sid_index//sid_per_chrom, 1+sid_index%sid_per_chrom, 1+sid_index%sid_per_chrom] for sid_index in xrange(len(sid))))
    if len(sid) == 0: #make it work when no sids are wanted
        pos = pos.reshape(len(sid),3)

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


from pysnptools.snpreader import SnpReader

def encode_snp(entry):
    if entry == 0:
        return "A A"
    elif entry == 1:
        return "A C"
    elif entry == 2:
        return "C C"

def write_tped(snpdata, basefilename):
    #\\bobd02\Public\PLink\x64\Plink.exe --noweb --tfile test --make-bed --out test2
    
    
    SnpReader._write_fam(snpdata, basefilename, remove_suffix="tped")
    #SnpReader._write_map_or_bim(snpdata, basefilename, remove_suffix="dat", add_suffix="map")

    snpsarray = snpdata.val
    with open(basefilename,"w") as dat_filepointer:
        for sid_index, sid in enumerate(snpdata.sid):
            if sid_index % 1000 == 0:
                logging.info("Writing snp # {0} to file '{1}'".format(sid_index, basefilename))
            dat_filepointer.write("1 {0} {1} {1} ".format(sid ,sid_index)) #use "j" and "n" as the major and minor allele
            row = snpsarray[:,sid_index]
            dat_filepointer.write(" ".join((encode_snp(i) for i in row)) + "\n")
    logging.info("Done writing " + basefilename)
    
    

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    #snpdata = snp_gen(fst=0.005, dfr=0.5, iid_count=100, sid_count=1000,chr_count=11)
    #Bed.write(snpdata, "25k")
    #write_tped(snpdata, "test.tped")

    import doctest
    doctest.testmod()


    print "done"

