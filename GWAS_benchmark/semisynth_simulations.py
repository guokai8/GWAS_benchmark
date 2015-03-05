"""
module to perform semi-synthetic simulations:
- take snps
- simulate phenotypes
- perform GWAS with different methods
- measure performance
"""

import logging
import os
import time
import sys

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab
import pylab

import fastlmm.association.gwas_eval as gw
from fastlmm.util.pickle_io import save, load
from fastlmm.util.runner import Local, Hadoop2, LocalMultiProc
from fastlmm.util import distributed_map
from pysnptools.standardizer import DiagKtoN
import split_data_helper
import semisynth_simulations


class LeaveTwoChrOutSimulation():

    def __init__(self, snp_fn, out_prefix):


        self.force_recompute = False

        #self.base_path = base_path
        self.snp_fn = snp_fn

        from pysnptools.snpreader import Bed
        self.snp_reader = Bed(snp_fn)
        
        self.eigen_fn = self.snp_fn + "_pcs.pickle"

        self.out_prefix = out_prefix


    def precompute_pca(self):
        """
        compute pcs
        """

        logging.info("computing PCA on train set")
        t0 = time.time()
        
        if not os.path.isfile(self.eigen_fn) or self.force_recompute:

            G = self.snp_reader.read(order='C').standardize().val
            G.flags.writeable = False
            chr1_idx, chr2_idx, rest_idx = split_data_helper.split_chr1_chr2_rest(self.snp_reader.pos)

            G_train = G.take(rest_idx, axis=1)

            from sklearn.decomposition import PCA
            pca = PCA()
            pcs = pca.fit_transform(G_train)

            logging.info("saving eigendecomp to file %s" % self.eigen_fn)
            
            eig_dec = {"pcs": pcs}
            save(self.eigen_fn, eig_dec)


            logging.info("time taken for pc computation: " + str(time.time()-t0))
        else:
            logging.info("pc file already exists: %s" % (self.eigen_fn))


    def run(self, methods, num_causal, num_repeats, num_pcs, description, runner, seed=None, plot_fn=None):
        
        
        self.precompute_pca()

        input_files = [self.snp_fn + ext for ext in [".bed", ".fam", ".bim"]] + [self.eigen_fn]
        input_args = [(methods, self.snp_fn, self.eigen_fn, num_causal, num_pcs, seed, sim_id) for sim_id in range(num_repeats)]
        output_list = distributed_map.d_map(semisynth_simulations.compute_core, input_args, runner, input_files=input_files)

        ############################################
        results_fn = "%s_results.runs_%i.causals_%i.pickle.bzip" % (description, num_repeats, num_causal)
        reduced_results_fn = results_fn.replace("runs", "reduced.runs")

        save(results_fn, output_list)

        
        methods = output_list[0][0].keys()
        arg_list = [(method, results_fn) for method in methods]

        #reduce_runner = Hadoop(len(methods), mapmemory=90*1024, reducememory=90*1024, mkl_num_threads=1, queue="shared")
        reduce_runner = Local()
        combine_output = distributed_map.d_map(semisynth_simulations.combine_results, arg_list, reduce_runner, input_files=[results_fn])
        
        save(reduced_results_fn, combine_output)
        title = "%i causal, %i repeats" % (num_causal, num_repeats)
        visualize_reduced_results(methods, combine_output, title=title, plot_fn=plot_fn)

        return combine_output



def simulate_ascertained(methods, prevalence, iid_count, num_causal, num_repeats, description, snp_args, phenotype_args, runner=Local(), seed=None, plot_fn=None):
    """
    run a synthetic simulation using ascertained data
    
    :param methods: A list of functions implementing methods to be compared.
    :type methods: list<function>
    
    :param prevalence: Prior probability of a case, e.g. .1
    :type prevalence: a float between 0.0 and 1.0 (exclusive)
       
    :param iid_count: The number of individuals to generate.
    :type iid_count: int
     
    :param num_causal: The number causal SNPs in the simulation.
    :type num_causal: int

    :param num_repeats: The number of repeats in the simulation.
    :type num_repeats: int

    :param description: Short description string of experiment (for output)
    :type description: str
    
    :param num_repeats: The number of repeats in the simulation.
    :type num_repeats: int

    :param snp_args: arguments for an internal call to :func:`GWAS_benchmark.snp_gen`. Do not include
    'iid_count' or 'seed'
    :type snp_args: dictionary

    :param phenotype_args: arguments for an internal call to :func:`.generate_phenotype`. Do not include
    'snp_count' or 'seed'
    :type phenotype_args: dictionary

    :param runner: a Runner object (e.g. Local, Hadoop, HPC)
    :type runner: Runner

    :param seed: a random seed to control random number generation
    :type seed: int

    :param plot_fn: filename under which to save the output figure
    :type plot_fn: str

    """    

    
    input_args = [(methods, num_causal, prevalence, iid_count, snp_args, phenotype_args, seed, sim_id) for sim_id in range(num_repeats)]
    output_list = distributed_map.d_map(semisynth_simulations.compute_core_ascertained, input_args, runner)


    ############################################
    results_fn = "%s_ascertained_results.runs_%i.causals_%i.pickle.bzip" % (description, num_repeats, num_causal)
    reduced_results_fn = results_fn.replace("runs", "reduced.runs")

    save(results_fn, output_list)

    
    methods = output_list[0][0].keys()
    arg_list = [(method, results_fn) for method in methods]

    combine_output = distributed_map.d_map(semisynth_simulations.combine_results, arg_list, Local(), input_files=[results_fn])
    
    save(reduced_results_fn, combine_output)
    title = "%i causal, %i repeats" % (num_causal, num_repeats)
    visualize_reduced_results(methods, combine_output, title=title, plot_fn=plot_fn)

    return combine_output



#TODO: create a different subclass for ascertained data
# call the other subclass existing phenotype

def visualize_reduced_results(methods, combine_output, title="", plot_fn=None):
        """
        set up plots: T1-error, ROC, PRC
        """

        t0 = time.time()

        fig = pylab.figure()
        fig.set_size_inches(26,7)
        for mi, method in enumerate(methods):
            o = combine_output[mi]
            pylab.subplot(131)
            gw.draw_t1err_curve(o["t1err"][0], o["t1err"][1], method, o["num_trials"])

            pylab.subplot(132)
            draw_roc_curve(o["roc"][0], o["roc"][1], o["roc"][2], method)

            pylab.subplot(133)
            gw.draw_prc_curve(o["prc"][0], o["prc"][1], o["prc"][2], method)


        pylab.suptitle(title)

        print "time taken to draw figure", time.time()-t0

        if plot_fn is None:
            print "showing figure!"
            pylab.show()
        else:
            print "saving figure!"
            pylab.savefig(plot_fn, dpi=100)


def combine_results(input_tuple):
    """
    compute performance statistics from p-values of method
    """
    
    method, results_fn = input_tuple

    logging.info("reading file: %s" % results_fn)
    output_list = load(results_fn)

    p_values_all = []
    mask_all = []    

    p_values_all = []
    p_values_chr1 = []
    p_values_chr2 = []
    mask_all = []

    t0 = time.time()
    logging.info("concatenating p-values")
    for result, idx in output_list:
        causals_chr2_idx = np.intersect1d(idx["chr2_idx"], idx["causal_idx"])

        assert len(result[method]) == len(idx["chr1_idx"]) + len(idx["chr2_idx"])

        p_vals_t1_err = result[method][idx["chr1_idx"]]
        p_vals_power = result[method][causals_chr2_idx]

        p_values_chr1.extend(p_vals_t1_err)
        p_values_chr2.extend(p_vals_power)

        p_values_all.extend(p_vals_t1_err)
        p_values_all.extend(p_vals_power)
                
        mask_t1_err = np.zeros(len(idx["chr1_idx"]), dtype=np.bool)
        mask_power = np.ones(len(causals_chr2_idx), dtype=np.bool)

        mask_all.extend(mask_t1_err)
        mask_all.extend(mask_power)
    
    logging.info("done concatenating p-values (%s)" % (str(time.time()-t0)))
    result = {}

    t0 = time.time()
    result["roc"] = gw.compute_roc_data(np.array(mask_all, dtype=np.bool), -np.array(p_values_all))
    logging.info("computed roc in (%s)" % (str(time.time()-t0)))

    t0 = time.time()
    result["prc"] = gw.compute_prc_data(np.array(mask_all, dtype=np.bool), -np.array(p_values_all))
    logging.info("computed prc in (%s)" % (str(time.time()-t0)))

    t0 = time.time()
    result["t1err"] = gw.compute_t1err_data(np.array(p_values_chr1), np.zeros(len(p_values_chr1), dtype=np.bool))
    logging.info("computed t1err in (%s)" % (str(time.time()-t0)))

    t0 = time.time()
    result["power"] = gw.compute_power_data(np.array(p_values_chr2), np.ones(len(p_values_chr2), dtype=np.bool))
    logging.info("computed power in (%s)" % (str(time.time()-t0)))

    result["method"] = method
    result["num_trials"] = len(p_values_chr1)

    return result



def generate_phenotype(snp_data, causals, genetic_var, noise_var, seed=None):
    """
    generate phenotype given genotype

    'causals' can be either an array of indexes to the causal snps or the number of causal snps desired.
    """

    if seed is not None:
        np.random.seed(int(seed % sys.maxint))
    
    try:
        num_causal = len(causals)
        causal_idx = causals
    except:
        num_causal = causals
        #the "if..else" is a work around because the linux version of np.random.choice doesn't like to select zero items from an empty list. We need to call random in either case so that the random seed ends up in the expected state
        causal_idx = np.random.choice(sp.arange(snp_data.sid_count if num_causal>0 else 1),size=num_causal,replace=False)

    num_phenotypes = 1
    mean = 0.0
    X = snp_data.val.copy()
    X.flags.writeable = False
    X_causal = X[:,causal_idx]
    X_causal = 1./np.sqrt(X_causal.shape[1]) * X_causal
    W = np.random.randn(num_causal, num_phenotypes) * np.sqrt(genetic_var) + mean
    XW = np.dot(X_causal, W)
    noise_std = np.sqrt(noise_var)
    Z = noise_std*sp.randn(X_causal.shape[0], num_phenotypes)
    y = XW + Z
    y = y[:,0]

    return y


def generate_discrete_ascertained(prevalence, iid_count, snp_args, phenotype_args, seed=0):
    """
    Generate discrete ascertained data. Internally, case will be generated at the requested
    prevalence. Before returning, however, the control will randomly sampled so 
    that in the returned data, case and control have number of examples.

    :param prevalence: Prior probability of a case, e.g. .1
    :type prevalence: a float between 0.0 and 1.0 (exclusive)

    :param iid_count: The number of examples desired in the returned data. Because of
    rounding during data generate the actual number may be lower. Of this happens,
    a warning will be shown.
    :type iid_count: int

    :param snp_args: arguments for an internal call to :func:`GWAS_benchmark.snp_gen`. Do not include
    'iid_count' or 'seed'
    :type snp_args: dictionary

    :param phenotype_args: arguments for an internal call to :func:`.generate_phenotype`. Do not include
    'snp_count' or 'seed'
    :type phenotype_args: dictionary

    :param seed: a random seed to control random number generation
    :type seed: int

    :rtype: a :class:`pysnptools.snpreader.SnpData' of genotype data and a nparray of 0,1 phenotype values.

    :Example:

    >>> snp_args = {"fst":.1,"dfr":.5,"sid_count":200,"maf_low":.05}
    >>> phenotype_args = {"causals":10,"genetic_var":0.5, "noise_var":0.5}
    >>> snps,pheno = generate_discrete_ascertained(prevalence=.1,iid_count=100,seed=5,snp_args=snp_args,phenotype_args=phenotype_args)
    >>> print int(snps.val.shape[0]),int(snps.val.shape[1]),int(len(pheno))
    98 200 98

    """
    assert 0<prevalence and prevalence <= .5, "Expect prevalence to be between 0.0 (exclusive) and .5 (inclusive)"
    assert int(iid_count) == iid_count and iid_count >= 0, "Expect iid_count to be a non-negative integer"

    # generate more examples than we ultimately want
    iid_count2 = int(float(iid_count) / 2.0 / prevalence)
    from GWAS_benchmark import snp_gen
    snp2 = snp_gen(iid_count=iid_count2, seed=seed, **snp_args)
    pheno2 = generate_phenotype(snp_data=snp2, seed=seed, **phenotype_args)

    # Sort by snps by pheno2 value
    snps2_sorted = snp2[pheno2.argsort(),:]

    # we want the top snp_count*prevalence for cases
    # and a random sample, of the same size, from the rest for control
    case_count = int(snps2_sorted.iid_count * prevalence)
    case_index = range(-1,-(case_count+1),-1) # e.g. if case_count is 3, then -1,-2,-3
    control_count = case_count

    if control_count + case_count != iid_count:
        logging.warn("iid_count is {0} instead of {1} because of rounding".format(control_count + case_count, iid_count))

    if seed is not None:
        np.random.seed(int(seed % sys.maxint))
    #print "gda", snps2_sorted.iid_count,case_count,control_count

    #the "if..else" is a work around because the linux version of np.random.choice doesn't like to select zero items from an empty list. We need to call random in either case so that the random seed ends up in the expected state
    control_index = np.random.choice(np.arange(snps2_sorted.iid_count-case_count if control_count > 0 else 1), control_count, replace=False)
    
    snp_final = snps2_sorted[np.concatenate((control_index,case_index)),:].read()
    pheno_final = np.zeros(control_count+case_count)
    pheno_final[control_count:]=1

    return snp_final, pheno_final





def compute_core(input_tuple):
    """
    Leave-two-chromosome-out evaluation scheme:
    Chr1: no causals, used for T1-error evaluation
    Chr2: has causals, not conditioned on, used for power evaluation
    Rest: has causals, conditioned on
    
      T1   Pow  [     cond     ] 
    ===== ===== ===== .... =====
            x x   x x      xx
    
    """
    
    
    
    methods, snp_fn, eigen_fn, num_causal, num_pcs, seed, sim_id = input_tuple
    
    # partially load bed file
    from pysnptools.snpreader import Bed
    snp_reader = Bed(snp_fn)

    # determine indices for generation and evaluation
    ##################################################################
    chr1_idx, chr2_idx, rest_idx = split_data_helper.split_chr1_chr2_rest(snp_reader.pos)
    
    causal_candidates_idx = np.concatenate((chr2_idx, rest_idx))
    # only compute t1-error (condition on all chr with causals on them)
    #causal_candidates_idx = rest_idx
    test_idx = np.concatenate((chr1_idx, chr2_idx))
    
    if seed is not None:
        np.random.seed(int(seed % sys.maxint))
    
    causal_idx = np.random.permutation(causal_candidates_idx)[0:num_causal]
    
    # generate phenotype
    ###################################################################
    genetic_var = 0.5
    noise_var = 0.5

    y = generate_phenotype(Bed(snp_fn).read(order='C').standardize(), causal_idx, genetic_var, noise_var)
    y.flags.writeable = False


    ############### only alter part until here --> modularize this


    # load pcs
    ###################################################################
    logging.info("loading eigendecomp from file %s" % eigen_fn)
    eig_dec = load(eigen_fn)
    G_pc = eig_dec["pcs"]
    G_pc.flags.writeable = False

    G_pc_ = G_pc[:,0:num_pcs]
    G_pc_norm = DiagKtoN(G_pc_.shape[0]).standardize(G_pc_.copy())
    G_pc_norm.flags.writeable = False
    

    # run feature selection
    #########################################################

    # generate pheno data structure
    pheno = {"iid": snp_reader.iid, "vals": y, "header": []}
    covar = {"iid": snp_reader.iid, "vals": G_pc_norm, "header": []}
    
    # subset readers
    G0 = snp_reader[:,rest_idx]
    test_snps = snp_reader[:,test_idx]
    
    result = {}
    fs_result = {}

    # additional methods can be defined and included in the benchmark
    for method_function in methods:
        result_, fs_result_ = method_function(test_snps, pheno, G0, covar)
        result.update(result_)
        fs_result.update(fs_result_)
    
    # save indices
    indices = {"causal_idx": causal_idx, "chr1_idx": chr1_idx, "chr2_idx": chr2_idx, "input_tuple": input_tuple, "fs_result": fs_result}
    #test_idx
    
    return result, indices


def compute_core_ascertained(input_tuple):
    """
    Leave-two-chromosome-out evaluation scheme:
    Chr1: no causals, used for T1-error evaluation
    Chr2: has causals, not conditioned on, used for power evaluation
    Rest: has causals, conditioned on
    
      T1   Pow  [     cond     ] 
    ===== ===== ===== .... =====
            x x   x x      xx
    
    """
    
    methods, num_causal, prevalence, iid_count, snp_args, phenotype_args, seed, sim_id = input_tuple

    
    
    # determine indices for generation and evaluation
    ##################################################################
    chr1_idx, chr2_idx, rest_idx = range(0,1000), range(1000, 2000), range(2000, 10000)
    
    causal_candidates_idx = np.concatenate((chr2_idx, rest_idx))
    # only compute t1-error (condition on all chr with causals on them)
    test_idx = np.concatenate((chr1_idx, chr2_idx))
    
    if seed is not None:
        np.random.seed(seed)
    
    causal_idx = np.random.permutation(causal_candidates_idx)[0:num_causal]
    
    
    
    # generate phenotype
    ###################################################################
    #y = generate_phenotype(Bed(snp_fn).read(order='C').standardize(), causal_idx, genetic_var, noise_var)
    #y.flags.writeable = False

    phenotype_args["causals"] = causal_idx
    
    #import pdb; pdb.set_trace()
    
    snp_reader, y = generate_discrete_ascertained(prevalence, iid_count, snp_args, phenotype_args, seed=seed)


    # run feature selection
    #########################################################

    # generate pheno data structure
    pheno = {"iid": snp_reader.iid, "vals": y, "header": []}
    covar = {"iid": snp_reader.iid, "vals": np.ones((len(y),1)), "header": []}
    
    # subset readers
    G0 = snp_reader[:,rest_idx]
    test_snps = snp_reader[:,test_idx]
    
    result = {}
    fs_result = {}

    # additional methods can be defined and included in the benchmark
    for method_function in methods:
        result_, fs_result_ = method_function(test_snps, pheno, G0, covar)
        result.update(result_)
        fs_result.update(fs_result_)
    
    # save indices
    indices = {"causal_idx": causal_idx, "chr1_idx": chr1_idx, "chr2_idx": chr2_idx, "input_tuple": input_tuple, "fs_result": fs_result}
    #test_idx
    
    return result, indices



def draw_roc_curve(fpr, tpr, roc_auc, label):
    """
    draw semi-log-scaled ROC curve
    """
    
    if len(fpr) > 1000:
        sub_idx = [int(a) for a in np.linspace(0, len(fpr)-1, num=1000, endpoint=True)]
        fpr, tpr = fpr[sub_idx], tpr[sub_idx]


    #pylab.semilogx(fpr, tpr, label='%s (area = %0.4f)' % (label, roc_auc))
    pylab.semilogx(fpr, tpr, label=label)
    #pylab.plot([0, 1], [0, 1], 'k--')
    pylab.xlim([0.0, 1.0])
    pylab.ylim([0.0, 1.0])
    #pylab.xlabel('False Positive Rate')
    pylab.xlabel('type I error', fontsize="large")
    #pylab.ylabel('True Positive Rate (Power)')
    pylab.ylabel('power', fontsize="large")
    #pylab.title('Receiver operating characteristic example')
    pylab.grid(True)
    pylab.legend(loc="lower right")



def run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, description, runner, plot_fn=None, seed=None):
    sc = LeaveTwoChrOutSimulation(snp_fn, out_prefix)
    return sc.run(methods, num_causals, num_repeats, num_pcs, "mouse_", runner, seed=seed, plot_fn=plot_fn)
    

def run_simulation_ascertained():

    snp_args = {"fst": 0.2, "dfr": 0.1, "sid_count": 10000}
    phenotype_args = {"genetic_var": 0.5, "noise_var": 0.5}
        # make this a tuple of function and kwargs
    from GWAS_benchmark.methods import execute_lmm, execute_linear_regression
    methods = [execute_lmm] #, execute_linear_regression]
    
    prevalence = 0.2
    num_causal = 20
    num_repeats = 50
    iid_count= 500
    description = "ascertained"
    
    simulate_ascertained(methods, prevalence, iid_count, num_causal, num_repeats, description, snp_args, phenotype_args) 


def main():
    logging.basicConfig(level=logging.INFO)
    
    
    #snp_fn = "data/toydata.5chrom"
    snp_fn = "data/mouse/alldata"
    out_prefix = "results/mouse_"

    description = "test_run"
    queue = "shared"
    #runner = Hadoop2(200, mapmemory=40*1024, reducememory=90*1024, mkl_num_threads=4, queue=queue)
    print "using snps", snp_fn
    #runner = LocalMultiProc(20)
    runner = Local()

    num_causals = 500
    num_repeats = 3
    num_pcs = 5
    
    # make this a tuple of function and kwargs
    from GWAS_benchmark.methods import execute_lmm, execute_linear_regression, execute_dual_fs, execute_fs
    methods = [execute_fs, execute_linear_regression]
    
    run_simulation(snp_fn, out_prefix, methods, num_causals, num_repeats, num_pcs, description, runner)
    

if __name__ == "__main__":
    run_simulation_ascertained()
    #main()
