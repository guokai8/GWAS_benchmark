import logging
import unittest
import sr.tests.test_snp_gen


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)

    suites = unittest.TestSuite([sr.tests.test_snp_gen.getTestSuite()])
    suites.debug

    r = unittest.TextTestRunner(failfast=False)
    r.run(suites)

    logging.info("done with testing")
