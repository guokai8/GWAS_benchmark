import logging
import pysnptools.util as pstutil
from pysnptools.snpreader import SnpData
import numpy as np
import time
#import md5
import random


def write_plink(X): #!!!cmk rename
    import pandas
    #import plink_write as pw
    sample_names =  np.array(['s_%d' % i for i in range(X.shape[0])], dtype=str)
    family_names = np.array(['f_%d' % i for i in range(X.shape[0])], dtype=str)[:, None]
    paternal_ids = np.zeros_like(family_names,dtype=str)
    maternal_ids = np.zeros_like(family_names,dtype=str)
    sex = np.zeros_like(family_names,dtype=str)
    fam_pid_mid_sex = np.concatenate((family_names, paternal_ids, maternal_ids, sex), axis=1)
    #F = pandas.DataFrame(data=fam_pid_mid_sex, index=sample_names, columns=['Family', 'Paternal ID', 'Maternal ID', 'Sex'])
    snp_names=np.empty(X.shape[1],dtype='|S15')
    for i in range(X.shape[1]):
        snp_names[i] = 'snp_%d' %i
    #G = pandas.DataFrame(data=X, index=sample_names, columns=snp_names)
    return sample_names, family_names, snp_names

def snp_gen(fst, dfr, iid_count, sid_count, maf_low=.05, seed=0): #!!!cmk move this to the front
    """
    #!!!cmk fill in with docs including example

    fst is Degree of Population Structure, e.g. [0.005, 0.01, 0.05, 0.1] !!!cmk was self.FSTs
        fst=0 is a special case #!!!is it?
    dfr is Degree of Family Relatedness, the fraction of individuals belonging to a family [0.0, 0.5, 0.6, 0.7, 0.8, 0.9], e.g. #!!!cmk was self.fracSibs
    iid_count !!!cmk was numIndividuals
    sid_count == numSnps
    MAF is Minor allele frequency
    freq_pop_1 !!cmk was caseFrac ????

    #!!!cmk do we still want maf_high?
    #!!!cmk, freq_pop_1=.5, 


    """
    #!!!cmk want??? t0 = time.time()
    #random.seed(seed) #!!!cmk still needed?
    #np.random.seed(seed)

    #self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
    #self.fracSibs=np.array([0.0,0.05,0.1,0.2])
    #self.h2s = np.arange(0.1, 0.7, 0.1)
    #self.var_hidden = np.arange(0.0, 1.0, 0.3)
    #self.num_causal = np.array([10, 50, 100, 500, 1000])
    #self.FSTs, self.fracSibs, self.h2s, self.var_hidden, self.num_causal
    #fst, fracSibs, h2, hh2, causal = params


    # freq_pop_1=.5 #!!!cmk these three are being ignored

    import sr.simulation.simulator as sim

    #fst, fracSibs, h2, hh2, causal = params
    options, args = sim.parseArgs()
    # here we just call parseargs because I'm not sure whether
    # the default params are set in the init or not. It's just for safety
    # we override the important stuff anyway
    options.csnps_hidden = 0
    options.h2 = .5  #!!!cmk needed?  #not rel
    options.fracSibs = dfr
    options.csnps = 0 #!!!cmk what is this?  #not rel
    options.minFreq = maf_low
    options.fst = fst

    options.var_hidden = 0 #!!!cmk needed? #not rel
    options.short_fn=None#!!!cmk remove
    options.numIndividuals=iid_count
    options.numSnps = sid_count
    options.num_folds=None #!!!cmk needed? 10? #not rel
    options.randomseed=seed
    options.penalty=None #!!!!cmk needed? 0 #not rel

    snps = sim.generate_data(options,args)


    sample_names, family_names, sid = write_plink(snps)
    iid = np.array(list(zip(sample_names,family_names[:,0])))
    pos = np.array(list([i,0,0] for i in xrange(len(sid)))) # every snp has position 0,0 on its own chrom

    snpdata = SnpData(iid, sid, pos, snps, 
                      parent_string="snp_gen(fst={0}, dfr={1}, iid_count={2}, sid_count={3}, maf_low={4})".format(fst, dfr, iid_count, sid_count, maf_low) #!!!cmk is this up-to-date?
                      )

    return snpdata


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    snpdata = snp_gen(fst=.1,dfr=.1,iid_count=100,sid_count=1000)
    print snpdata


    import doctest
    doctest.testmod()


    print "done"

