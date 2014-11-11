"""
module to perform semi-synthetic simulations:
- take real snps
- simulate phenotypes
- perform GWAS with different methods
- measure performance
"""

from sr import LeaveTwoChrOutSimulation
from fastlmm.util.runner import Local, Hadoop2


def main():
 
    
    #snp_fn = "data/toydata.5chrom"
    snp_fn = "data/mouse/alldata"
    out_prefix = "results/mouse_"

    #queue = "shared"
    #runner = Hadoop2(1000, mapmemory=40*1024, reducememory=90*1024, mkl_num_threads=4, queue=queue)
    #print "using snps", snp_fn
    runner = Local()

    num_causals = 10
    num_repeats = 1000
    sc = LeaveTwoChrOutSimulation(snp_fn, out_prefix)
    sc.run(num_causals, num_repeats, "mouse_", runner)

if __name__ == "__main__":
    main()

