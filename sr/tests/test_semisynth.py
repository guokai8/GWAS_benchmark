import numpy as np
import scipy as sp
import logging
import doctest

import unittest
import os.path
import os
import time

from sr.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs
from sr.semisynth_simulations import run_simulation
from fastlmm.util.runner import Local

   
class TestSemiSynth(unittest.TestCase):     


    @classmethod
    def setUpClass(self):
        pass

    def test_gen1(self):
        import fastlmm.util.runner as runner

        snp_fn = os.path.join(os.path.dirname(os.path.realpath(runner.__file__)),"../../../tests/datasets/mouse/alldata")
        out_prefix = "tempdir/mouse_"

        num_causals = 500
        num_repeats = 10
        num_pcs = 5
    
        methods = [execute_fs, execute_linear_regression]
    
        run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, "test_gen1", Local())


def getTestSuite():
    """
    set up test suite
    """
    
    test_suite = unittest.TestSuite([])
    test_suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestSemiSynth))
    return test_suite

if __name__ == '__main__':

    suites = getTestSuite()
    if True:
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


