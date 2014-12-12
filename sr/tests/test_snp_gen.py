import numpy as np
import scipy as sp
import logging
import doctest

import unittest
import os.path
import os
import time

from sr.snp_gen import snp_gen #!!!cmk shorten namespace so not snp_gen.snp_gen
from pysnptools.snpreader import Dat
import pysnptools.util as pstutil
   
class TestSnpGen(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    @staticmethod
    def is_same(snpdata0, snpdata1): #!!! should this be an equality _eq_ operator on snpdata?
        result = (np.array_equal(snpdata0.iid,snpdata1.iid) and 
                  np.array_equal(snpdata0.sid, snpdata1.sid) and 
                  np.array_equal(snpdata0.iid, snpdata1.iid) and 
                  np.array_equal(snpdata0.pos, snpdata1.pos) and
                  np.array_equal(snpdata0.val, snpdata1.val))
        return result

    def gen_and_compare(self, output_file, **kwargs):
        gen_snpdata = snp_gen(**kwargs)
        #pstutil.create_directory_if_necessary(elf.currentFolder + "/tempdir/" + output_file,isfile=True) #!!!cmk move this to Ped.write?
        Dat.write(gen_snpdata, self.currentFolder + "/tempdir/" + output_file) #!!!cmk comment out
        ref_snpdata = Dat(self.currentFolder + "/expected/" + output_file).read()
        assert TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Failure on "+output_file
        return gen_snpdata
        #!!!cmk Ped doesn't seem to round trip well
        #!!!cmk Hdf5 doesn't seem to round trip well


    def test_gen1(self):
        self.gen_and_compare("gen1.dat", fst=0,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2(self):
        self.gen_and_compare("gen2.dat", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen2b(self):
        """
        Test that different seed produces different result
        """
        gen_snpdata = self.gen_and_compare("gen2b.dat", fst=.1,dfr=.5,iid_count=200,sid_count=20,maf_low=.05,seed=6)
        ref_snpdata = Dat(self.currentFolder + "/expected/gen2.dat").read()
        assert not TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Expect different seeds to produce different results"

    def test_gen3(self):
        self.gen_and_compare("gen3.dat", fst=.1,dfr=0,iid_count=200,sid_count=20,maf_low=.05,seed=5)

    def test_gen4(self):
        self.gen_and_compare("gen4.dat", fst=.1,dfr=.01,iid_count=200,sid_count=20,maf_low=.1,seed=5)


def getTestSuite():
    """
    set up test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSnpGen))
    return test_suite

if __name__ == '__main__':

    suites = getTestSuite()
    if True:
        # TestFeatureSelection().test_aaa_hdf5_speed()
        r = unittest.TextTestRunner(failfast=False)
        r.run(suites)
    else: #Cluster test run
        task_count = 500
        runner = HPC(task_count, 'RR1-N13-09-H44',r'\\msr-arrays\Scratch\msr-pool\Scratch_Storage6\Redmond',
                     remote_python_parent=r"\\msr-arrays\Scratch\msr-pool\Scratch_Storage6\REDMOND\carlk\Source\carlk\july_7_14\pythonpath",
                     update_remote_python_parent=True,
                     min=150,
                     priority="AboveNormal",mkl_num_threads=1)
        runner = Local()
        #runner = LocalMultiProc(taskcount=4,mkl_num_threads=5)
        #runner = LocalInParts(1,2,mkl_num_threads=1) # For debugging the cluster runs
        #runner = Hadoop2(100, mapmemory=8*1024, reducememory=8*1024, mkl_num_threads=1, queue="default")
        distributable_test = DistributableTest(suites,"temp_test")
        print runner.run(distributable_test)


