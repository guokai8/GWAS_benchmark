import logging
import unittest
import GWAS_benchmark.tests.test_snp_gen
import GWAS_benchmark.tests.test_semisynth

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    suites = unittest.TestSuite([GWAS_benchmark.tests.test_semisynth.getTestSuite(), GWAS_benchmark.tests.test_snp_gen.getTestSuite()])
    suites.debug

    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)

    logging.info("done with testing")
