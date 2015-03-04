"""
module to perform semi-synthetic simulations:
- take real snps
- simulate phenotypes
- perform GWAS with different methods
- measure performance
"""

from GWAS_benchmark.semisynth_simulations import run_simulation
from fastlmm.util.runner import Local


def main():
    
    snp_fn = "data/mouse/alldata"
    out_prefix = "results/mouse_"

    description = "test_run"
    runner = Local()

    num_causals = 500
    num_repeats = 10
    num_pcs = 5
    
    # make this a tuple of function and kwargs
    from GWAS_benchmark.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs
 
    for name, method in {"lmm": execute_lmm, "lr": execute_linear_regression, "dual_fs": execute_dual_fs, "fs": execute_fs}.items():
        methods = [method]
        combine_output = run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, description, runner, plot_fn=None, seed=42)


if __name__ == "__main__":
    main()
    