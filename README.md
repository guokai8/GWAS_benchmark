This python code can be used to benchmark or evaluate GWAS algorithms.
  
If you use this code, please cite:

* C. Widmer, C. Lippert, O. Weissbrod, N. Fusi, C. Kadie, R. Davidson, J. Listgarten, and D. Heckerman, Further Improvements to Linear Mixed Models for Genome-Wide Association Studies, _Scientific Reports_ **4**, 6874, Nov 2014 (doi:10.1038/srep06874).

This code contains the following modules:

* semisynth_experiments: the core module for generating synthetic phenotypes based on real snps, running different methods for GWAS and evaluating them all within one pipeline

* cluster_data: module to compute and visualize a hierarchical clustering of GWAS data to get an understanding of its structure (population structure, family structure)

* split_data_helper: helper module for splitting SNPs by chromosome

For testing purposes a small data set is provided at `data/mouse` (see the `README` file within that directory for the data license).

An example run to compute type I error rate on the mouse data using 10 causal SNPs can be executed by running `python run_simulation.py`.

We recommend running this example on a cluster computer as this simulation is computationally demanding. An example result plot (of type I error) is provided in the results directory.

Further, we use the ipython-notebook to demonstrate some of the functionality of the hierarchical clustering module: 
http://nbviewer.ipython.org/github/MSRCompBio/GWAS_benchmark/blob/master/cluster_data.ipynb

To start ipython notebook locally, type `ipython notebook` at the command line.
