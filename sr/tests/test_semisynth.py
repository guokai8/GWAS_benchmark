import unittest
import os.path

from sr.semisynth_simulations import run_simulation
from fastlmm.util.runner import Local

import numpy as np


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
        pass

    def test_all(self):
        import fastlmm.util.runner as runner

        currentFolder = os.path.dirname(os.path.realpath(__file__))
        snp_fn = os.path.realpath(currentFolder + "/../../data/mouse/alldata")
        out_prefix = currentFolder + "/tempdir/mouse_"

    
        description = "test_run"
        runner = Local()
    
        num_causals = 500
        num_repeats = 1
        num_pcs = 5
        
        # make this a tuple of function and kwargs
        from sr.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs
        methods_dict = {"lmm": execute_lmm, "lr": execute_linear_regression, "dual_fs": execute_dual_fs, "fs": execute_fs}
        
        for name, method in methods_dict.items():
            expected_prefix = currentFolder + "/expected/"
            methods = [method]
            combine_output = run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, description, runner, plot_fn="out.png", seed=42)
            from fastlmm.util.pickle_io import load
            co = load("%s%s.bzip" % (expected_prefix, name))
            
            compare_nested(combine_output, co)



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

