import numpy as np
import scipy as sp
import logging
import doctest

import unittest
import os.path
import os
import time

from GWAS_benchmark import snp_gen
from pysnptools.snpreader import Bed
import pysnptools.util as pstutil
   
class TestSnpGen(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    @staticmethod
    def is_same(snpdata0, snpdata1): #!!! should this be an equality _eq_ operator on snpdata?
        if not (np.array_equal(snpdata0.iid,snpdata1.iid) and 
                  np.array_equal(snpdata0.sid, snpdata1.sid) and 
                  np.array_equal(snpdata0.pos, snpdata1.pos)):
            return False

        try:
            np.testing.assert_equal(snpdata0.val, snpdata1.val)
        except:
            return False
        return True

    def gen_and_compare(self, output_file, **kwargs):
        gen_snpdata = snp_gen(**kwargs)
        #pstutil.create_directory_if_necessary(self.currentFolder + "/tempdir/" + output_file,isfile=True)
        #Bed.write(gen_snpdata, self.currentFolder + "/tempdir/" + output_file)  #comment out
        ref_snpdata = Bed(self.currentFolder + "/expected/" + output_file).read()
        assert TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Failure on "+output_file
        return gen_snpdata
        #!!! Ped doesn't seem to round trip well
        #!!! Hdf5 doesn't seem to round trip well


    def test_gen1(self):
        self.gen_and_compare("gen1", fst=0,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2(self):
        self.gen_and_compare("gen2", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2b(self):
        """
        Test that different seed produces different result
        """
        gen_snpdata = self.gen_and_compare("gen2b", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
        ref_snpdata = Bed(self.currentFolder + "/expected/gen2").read()
        assert not TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Expect different seeds to produce different results"

    def test_gen3(self):
        self.gen_and_compare("gen3", fst=.1,dfr=0,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen4(self):
        self.gen_and_compare("gen4", fst=.1,dfr=.01,iid_count=200,sid_count=20,maf_low=.1,seed=5)

    def test_gen5(self):
        gen_snpdata = self.gen_and_compare("gen5", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,maf_high=.4, seed=5)
        ref_snpdata = Bed(self.currentFolder + "/expected/gen2").read()
        assert not TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Expect different seeds to produce different results"

    def test_gen6(self):
        gen_snpdata = self.gen_and_compare("gen6", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,sibs_per_family=5)
        ref_snpdata = Bed(self.currentFolder + "/expected/gen2").read()
        assert not TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Expect different seeds to produce different results"

    def test_gen7(self):
        gen_snpdata = self.gen_and_compare("gen7", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,freq_pop_0=.75)
        ref_snpdata = Bed(self.currentFolder + "/expected/gen2").read()
        assert not TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Expect different seeds to produce different results"

    def test_gen8(self):
        self.gen_and_compare("gen8a", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=3)
        self.gen_and_compare("gen8b", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=4)
        self.gen_and_compare("gen8c", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5,chr_count=6)

    def test_gensmall(self):
        #Just checking that doesn't generate errors
        for iid_count in [10, 5, 3, 2, 1, 0]:
            for sid_count in [0, 10, 5, 3, 2, 1]:
                for chr_count in [30, 10, 5, 3, 2, 1, 0]:
                    if chr_count == 0 and sid_count > 0:
                        continue # not break
                    logging.debug("{0}, {1}, {2}".format(iid_count, sid_count, chr_count))
                    snpdata = snp_gen(fst=.1,dfr=.5,iid_count=iid_count,sid_count=sid_count,maf_low=.05,seed=6,chr_count=chr_count)
                    assert snpdata.iid_count <= iid_count
                    assert snpdata.sid_count == sid_count
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,0]) <= chr_count
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,1]) <= int(max(1,np.ceil(float(sid_count) / chr_count)))
                    assert len(snpdata.pos) == 0 or max(snpdata.pos[:,2]) <= int(max(1,np.ceil(float(sid_count) / chr_count)))

    def test_doc_test(self):
        import GWAS_benchmark
        import sys
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/..")
        result = doctest.testmod(sys.modules['GWAS_benchmark.snp_gen'])
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__

def getTestSuite():
    """
    set up test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpGen))
    return test_suite

if __name__ == '__main__':

    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)
