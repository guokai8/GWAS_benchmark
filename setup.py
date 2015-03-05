"""
file to set up python package, see http://docs.python.org/2/distutils/setupscript.html for details.
"""


import platform
import os
import sys
import shutil

from distutils.core import setup
from distutils.extension import Extension
from distutils.command.clean import clean as Clean

try:
	import numpy
except Exception:
	print "numpy needed for installation, please install numpy first"
	sys.exit()


def readme():
    with open('README.md') as f:
       return f.read()

#python setup.py sdist bdist_wininst upload
setup(
    name='GWAS_benchmark',
    version='0.1.0',
    description='GWAS benchmark',
    long_description=readme(),
    keywords='gwas bioinformatics benchmarks',
    url="http://research.microsoft.com/en-us/um/redmond/projects/mscompbio/",
    author='MSR',
    author_email='fastlmm@microsoft.com',
    license='Apache 2.0',
    packages=[
        "GWAS_benchmark/tests",
        "GWAS_benchmark",
	],
    requires = ['numpy', 'scipy', 'pandas', 'scikit-learn', 'matplotlib', 'pysnptools', 'fastlmm', 'fastcluster'],
  )

