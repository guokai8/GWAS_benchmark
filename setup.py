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


class CleanCommand(Clean):
    description = "Remove build directories, and compiled files (including .pyc)"

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('.'):
            for filename in filenames:
                if (   filename.endswith('.so')
                    or filename.endswith('.pyd')
                    or filename.endswith('.pyc')
                                ):
                    tmp_fn = os.path.join(dirpath, filename)
                    print "removing", tmp_fn
                    os.unlink(tmp_fn)

# set up macro
if platform.system() == "Darwin":
    macros = [("__APPLE__", "1")]
elif "win" in platform.system().lower():
    macros = [("_WIN32", "1")]
else:
    macros = [("_UNIX", "1")]

ext = []

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
	package_data={
                 },
    requires = ['numpy', 'scipy', 'pandas', 'sklearn', 'matplotlib', 'pysnptools', 'fastlmm', 'fastcluster'],
    #zip_safe=False,
    # extensions
    cmdclass = {'build_ext': build_ext, 'clean': CleanCommand},
    ext_modules = ext,
	include_dirs = [numpy.get_include()]
  )

