import unittest
import os.path
import sys
import logging

from GWAS_benchmark.semisynth_simulations import run_simulation
from fastlmm.util.runner import Local
import numpy as np
from GWAS_benchmark.semisynth_simulations import generate_discrete_ascertained
from pysnptools.snpreader import Bed
import pysnptools.util as pstutil
from GWAS_benchmark.tests.test_snp_gen import TestSnpGen
import doctest
from GWAS_benchmark.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs



def compare_np_arrays(lhs, rhs):
    np.testing.assert_array_almost_equal(lhs, rhs, decimal=4, verbose=True)

def compare_lists(lhs, rhs):

    assert len(lhs) == len(rhs)
    for idx in xrange(len(lhs)):
        compare_nested(lhs[idx], rhs[idx])

def compare_dicts(lhs, rhs):
    
    #assert (lhs.keys() == rhs.keys()).all()
    assert len(lhs) == len(rhs)
    
    for key in lhs.keys():
        assert key in rhs.keys()
        compare_nested(lhs[key], rhs[key])

def compare_nested(lhs, rhs):

    assert type(lhs) == type(rhs)
    
    if type(lhs) == list or type(lhs) == tuple:
        compare_lists(lhs, rhs)
    elif type(lhs) == np.ndarray:
        compare_np_arrays(lhs, rhs)
    elif type(lhs) == dict:
        compare_dicts(lhs, rhs)
    elif type(lhs) == float:
        np.testing.assert_approx_equal(lhs, rhs)
    else:
        assert lhs == rhs
   
class TestSemiSynth(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        self.currentFolder = os.path.dirname(os.path.realpath(__file__))

    from GWAS_benchmark.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs

    def test_lmm(self):
        self.run_sim_and_compare("lmm", execute_lmm)

    def test_lr(self):
        self.run_sim_and_compare("lr", execute_linear_regression)

    def test_dual_fs(self):
        self.run_sim_and_compare("dual_fs", execute_dual_fs)

    def test_fs(self):
        self.run_sim_and_compare("fs", execute_fs)


    def run_sim_and_compare(self, name, method):
        logging.info('in test_all')
        import fastlmm.util.runner as runner

        currentFolder = os.path.dirname(os.path.realpath(__file__))
        snp_fn = os.path.realpath(currentFolder + "/../../data/mouse/alldata")
        out_prefix = currentFolder + "/tempdir/mouse_"

    
        description = "test_run_{0}".format(name)
        runner = Local()
    
        num_causals = 500
        num_repeats = 1
        num_pcs = 5
        
        expected_prefix = currentFolder + "/expected/"
        methods = [method]
        combine_output = run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, description, runner, plot_fn="out.png", seed=42)
        from fastlmm.util.pickle_io import load
        filename = "%s%s.bzip" % (expected_prefix, name)
        co = load(filename)
        compare_nested(combine_output, co)

    def compare(self, gen_snpdata, pheno, output_file):
        assert len(pheno) % 2 == 0 # even number of iids
        assert np.all(pheno[:len(pheno)//2] == 0) #first half are controls
        assert np.all(pheno[len(pheno)//2:] == 1) #second half are cases
        #pstutil.create_directory_if_necessary(self.currentFolder + "/tempdir/" + output_file,isfile=True) #comment out
        #Bed.write(gen_snpdata, self.currentFolder + "/tempdir/" + output_file)  #comment out
        ref_snpdata = Bed(self.currentFolder + "/expected/" + output_file).read()
        assert TestSnpGen.is_same(gen_snpdata, ref_snpdata), "Failure on "+output_file


    def test_generate_discrete_ascertained_1(self):
        for prevalence in [.5, .01,.1]:
            for sid_count in [2, 200,20,1,0]:
                for iid_count in [2,1000, 100,10,1,0]:
                    #print prevalence, sid_count, iid_count
                    snp_args = {"fst":.1,"dfr":.5,"sid_count":sid_count,"maf_low":.05}
                    phenotype_args = {"causals":sid_count//10,"genetic_var":0.5, "noise_var":0.5}
                    #print prevalence, iid_count,snp_args,phenotype_args
                    snps,pheno = generate_discrete_ascertained(prevalence=prevalence,iid_count=iid_count,seed=5,snp_args=snp_args,phenotype_args=phenotype_args)
                    output_file = "gda.{0}_{1}_{2}".format(prevalence,sid_count,iid_count)
                    self.compare(snps,pheno,output_file)

    def test_doc_test(self):
        import GWAS_benchmark.semisynth_simulations
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__))+"/..")
        result = doctest.testmod(GWAS_benchmark.semisynth_simulations)
        os.chdir(old_dir)
        assert result.failed == 0, "failed doc test: " + __file__



def getTestSuite():
    """
    set up test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSemiSynth))
    return test_suite

if __name__ == '__main__':
    suites = getTestSuite()
    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)

