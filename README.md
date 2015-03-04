## GWAS_benchmark
-------------------------------------

This python code can be used to benchmark or evaluate GWAS algorithms.
  
If you use this code, please cite:

* C. Widmer, C. Lippert, O. Weissbrod, N. Fusi, C. Kadie, R. Davidson, J. Listgarten, and D. Heckerman, Further Improvements to Linear Mixed Models for Genome-Wide Association Studies, _Scientific Reports_ **4**, 6874, Nov 2014 (doi:10.1038/srep06874).

See this website for related software:  
http://research.microsoft.com/en-us/um/redmond/projects/MicrosoftGenomics/

Our documentation (including live examples) is available as ipython notebook:
http://nbviewer.ipython.org/github/MicrosoftGenomics/GWAS_benchmark/blob/master/GWAS_benchmark/simulation.ipynb

(To start ipython notebook locally, type `ipython notebook` at the command line.)

This code contains the following modules:

* semisynth_experiments: the core module for generating synthetic phenotypes based on real snps, running different methods for GWAS and evaluating them all within one pipeline

* cluster_data: module to compute and visualize a hierarchical clustering of GWAS data to get an understanding of its structure (population structure, family structure)

* split_data_helper: helper module for splitting SNPs by chromosome

For testing purposes a small data set is provided at `data/mouse` (see the `README` file within that directory for the data license).

An example run to compute type I error rate on the mouse data using 10 causal SNPs can be executed by running `python run_simulation.py`.

We recommend running this example on a cluster computer as this simulation is computationally demanding. An example result plot (of type I error) is provided in the results directory.

Further, we use the ipython-notebook to demonstrate some of the functionality of the hierarchical clustering module: 
http://nbviewer.ipython.org/github/MSRCompBio/GWAS_benchmark/blob/master/cluster_data.ipynb

### Quick install:


If you have pip installed, installation is as easy as:

```
pip install GWAS_benchmark
```


### Detailed Package Install Instructions:


fastlmm has the following dependencies:

python 2.7

Packages:

* numpy
* scipy
* matplotlib
* pandas
* scikit.learn (sklearn)
* fastcluster
* fastlmm
* pysnptools
* optional: [statsmodels -- install only required for logistic-based tests, not the standard linear LRT]


#### (1) Installation of dependent packages

We highly recommend using a python distribution such as 
Anaconda (https://store.continuum.io/cshop/anaconda/) 
or Enthought (https://www.enthought.com/products/epd/free/).
Both these distributions can be used on linux and Windows, are free 
for non-commercial use, and optionally include an MKL-compiled distribution
for optimal speed. This is the easiest way to get all the required package
dependencies.


#### (2) Installing from source

Go to the directory where you copied the source code for fastlmm.

On linux:

At the shell, type: 
```
sudo python setup.py install
```

On Windows:

At the OS command prompt, type 
```
python setup.py install
```


### For developers (and also to run regression tests)

When working on the developer version, just set your PYTHONPATH to point to the directory
above the one named GWAS_benchmark in the source code. For e.g. if GWAS_benchmark is 
in the [somedir] directory, then in the unix shell use:
```
export PYTHONPATH=$PYTHONPATH:[somedir]
```
Or in the Windows DOS terminal, one can use: 
```
set PYTHONPATH=%PYTHONPATH%;[somedir]
```
(or use the Windows GUI for env variables).

#### Running regression tests

From the directory tests at the top level, run:
```
python test.py
```
This will run a
series of regression tests, reporting "." for each one that passes, "F" for each
one that does not match up, and "E" for any which produce a run-time error. After
they have all run, you should see the string "............" indicating that they 
all passed, or if they did not, something such as "....F...E......", after which
you can see the specific errors.

Note that you must set your PYTHONPATH as described above to run the 
regression tests, and not "python setup.py install".
