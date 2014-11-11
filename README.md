This is the python code used for the experiments in the SR-paper "Further Improvements to Linear Mixed Models for Genome-Wide Association Studies".

If this is useful to you, please cite:
C. Widmer, C. Lippert, O. Weissbrod, N. Fusi, C. Kadie, R. Davidson, J. Listgarten, and D. Heckerman. Further Improvements to Linear Mixed Models for Genome-Wide Association Studies.


It contains the following modules:
- semisynth_experiments: the core module for generating synthetic phenotypes based on real snps, running different methods for GWAS and evaluating them all within one pipeline
- cluster_data: module to compute and visualize a hierarchical clustering of GWAS data to get an understanding of its structure (population structure, family structure)
- split_data_helper: helper module for splitting SNPs by chromosome

For testing purposes a small data set is provided under (see the README file within that directory for the data license):
data/mouse

An example run to compute T1 error on the mouse data using 10 causal SNPs can be executed by running: 
python run_simulation.py

We recommend running this example on a cluster computer as this simulation is computationally demanding. An example result plot (of T1-error) is provided in the results directory.

Further, we use the ipython-notebook to demonstrate some of the functionality of the hierarchical clustering module. To start ipython notebook type:
ipython notebook
