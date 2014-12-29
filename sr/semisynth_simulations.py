"""
module to perform semi-synthetic simulations:
- take real snps
- simulate phenotypes
- perform GWAS with different methods
- measure performance
"""

import logging
import os
import numpy as np
import scipy as sp
import pandas as pd
        
import time

import fastlmm.association.gwas_eval as gw
import pylab

from fastlmm.util.pickle_io import save, load
from fastlmm.association.LocoGwas import FastGwas
from fastlmm.association import single_snp
from fastlmm.util.runner import Local, Hadoop2, LocalMultiProc
from fastlmm.util.util import argintersect_left
from fastlmm.util import distributed_map
from fastlmm.feature_selection.feature_selection_two_kernel import FeatureSelectionInSample


from pysnptools.standardizer import DiagKtoN

import split_data_helper

import semisynth_simulations



def create_feature_selection_distributable(snp_reader, phen_fn, pc_fn, num_pcs_kernel, output_prefix, snp_idx, cov_fn=None, include_all=True):

        snp_reader_subset = snp_reader[:,snp_idx]

        from fastlmm.feature_selection import FeatureSelectionStrategy
        import fastlmm.feature_selection.PerformSelectionDistributable as psd

        # set up paramters
        num_folds = 10
        random_state = 42
        num_snps_in_memory = 1000000

        ##############################
        num_steps_delta = 7
        num_steps_k = 7
        num_steps_mix = 7

        # log_2 space and all SNPs
        
        k_values = [int(k) for k in np.logspace(0, 10, base=2, num=num_steps_k, endpoint=True)]
        if include_all:
            k_values.append(len(snp_idx))
        delta_values = np.logspace(-5, 10, endpoint=True, num=num_steps_delta, base=np.exp(1))
        mix_values = np.linspace(0.0, 1.0, num=num_steps_mix, endpoint=True)

        if pc_fn is None:
            assert num_pcs_kernel == 0
            logging.info("feature selection: no PCs specified, disabling loop over mixing parameter")
            mix_values = np.array([0.0])

        strategy = "insample_cv"
        select_by_ll = True

        # go!
        feature_selector = FeatureSelectionStrategy(snp_reader_subset, phen_fn, num_folds, random_state=random_state, num_snps_in_memory=num_snps_in_memory, interpolate_delta=False, cov_fn=cov_fn)
        perform_selection_distributable = psd.PerformSelectionDistributable(feature_selector, k_values, delta_values, strategy, output_prefix, select_by_ll)

        return perform_selection_distributable


class LeaveTwoChrOutSimulation():

    def __init__(self, snp_fn, out_prefix):


        self.random_state = 42
        self.force_recompute = False
        self.mindist = 50

        #self.base_path = base_path
        self.snp_fn = snp_fn

        from pysnptools.snpreader import Bed
        self.snp_reader = Bed(snp_fn)
        
        self.cache_dir =  "data/"

        self.eigen_fn = self.cache_dir + "pcs.pickle"
        self.pc_prefix = self.cache_dir + "pcs"

        self.phen_string = "phen_sanity"
        self.phen_prefix = self.cache_dir + self.phen_string

        self.out_prefix = out_prefix

        self.simulator = None
        self.pc_selector = None
        self.feature_selector = None
        self.gwas = None

        # get from pc file
        self.S = None
        self.U = None

        self.p_values = None


    def precompute_pca(self):
        """
        compute pcs
        """
        


        logging.info("computing PCA on train set")
        t0 = time.time()
        assert os.path.exists(self.cache_dir)

        if not os.path.isfile(self.eigen_fn) or self.force_recompute:

            G = self.snp_reader.read(order='C').standardize().val
            G.flags.writeable = False
            chr1_idx, chr2_idx, rest_idx = split_data_helper.split_chr1_chr2_rest(snp_reader.pos)

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


    def run(self, methods, num_causal, num_repeats, num_pcs, description, runner, plot_fn=None):
        
        assert os.path.exists(self.cache_dir), "path does not exist %s" % (self.cache_dir)
        
        self.precompute_pca()

        input_files = [self.snp_fn + ext for ext in [".bed", ".fam", ".bim"]] + [self.eigen_fn]
        input_args = [(methods, self.snp_fn, self.eigen_fn, num_causal, num_pcs, sim_id) for sim_id in range(num_repeats)]
        output_list = distributed_map.d_map(semisynth_simulations.compute_core, input_args, runner, input_files=input_files)

        ############################################
        # power
        #indices = {"causal_idx": causal_idx, "chr1_idx": chr1_idx, "chr2_idx": chr2_idx}

        
        results_fn = "results/%s_results.runs_%i.causals_%i.pickle.bzip" % (description, num_repeats, num_causal)
        reduced_results_fn = results_fn.replace("runs", "reduced.runs")

        save(results_fn, output_list)

        
        methods = output_list[0][0].keys()
        arg_list = [(method, results_fn) for method in methods]

        #reduce_runner = Hadoop(len(methods), mapmemory=90*1024, reducememory=90*1024, mkl_num_threads=1, queue="shared")
        reduce_runner = Local()
        combine_output = distributed_map.d_map(semisynth_simulations.combine_results, arg_list, reduce_runner, input_files=[results_fn])
        
        save(reduced_results_fn, combine_output)
        title = "%i causal, %i repeats" % (num_causal, num_repeats)
        visualized_reduced_results(methods, combine_output, title=title, plot_fn=plot_fn)


def visualized_reduced_results(methods, combine_output, title="", plot_fn=None):

        t0 = time.time()

        fig = pylab.figure()
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
            pylab.show()
        else:
            fig.set_size_inches(12.5,12.5)
            pylab.savefig(plot_fn, dpi=100)


def combine_results(input_tuple):

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



def generate_phenotype(snp_data, causal_idx, genetic_var, noise_var):
    """
    generate phenotype given genotype
    """
    
    
    num_causal = len(causal_idx)
    
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
    
    methods, snp_fn, eigen_fn, num_causal, num_pcs, sim_id = input_tuple
    
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
    causal_idx = np.random.permutation(causal_candidates_idx)[0:num_causal]
    
    # generate phenotype
    ###################################################################
    genetic_var = 0.5
    noise_var = 0.5

    y = generate_phenotype(Bed(snp_fn).read(order='C').standardize(), causal_idx, genetic_var, noise_var)
    y -= y.mean()
    y /= y.std()
    y.flags.writeable = False


    # load pcs
    ###################################################################
    logging.info("loading eigendecomp from file %s" % eigen_fn)
    eig_dec = load(eigen_fn)
    G_pc = eig_dec["pcs"]
    G_pc.flags.writeable = False

    G_pc_ = G_pc[:,0:num_pcs]
    G_pc_norm = 1./np.sqrt(compute_kernel_diag_from_G(G_pc_) / float(G_pc_.shape[0])) * G_pc_
    G_pc_norm2 = DiagKtoN(G_pc_.shape[0]).standardize(G_pc_.copy())
    
    np.testing.assert_array_almost_equal(G_pc_norm, G_pc_norm2)
    G_pc_norm.flags.writeable = False
    

    # run feature selection
    #########################################################
    delta = None
    
    
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


def execute_lmm(test_snps, pheno, G0, covar):
    
    result = {}
    fs_result = {}
    
    result["full"] = single_snp(test_snps, pheno, G0=G0, covar=covar).sort(["Chr", "ChrPos"])["PValue"].as_matrix()


    #TODO: implement linear regression with same interface as single_snp
    # linear regression with causals as covariates
    #from fastlmm.inference.linear_regression import f_regression_cov
    #_, result["linreg"] = f_regression_cov(G_test.copy(), y.copy(), np.ones((len(y),1)))
    #_, result["linreg_cov_pcs"] = f_regression_cov(G_test.copy(), y.copy(), G_pc_norm.copy())
    
    return result, fs_result
    
    
def execute_fs_methods():
    
    print "TODO"    
    """
    # invoke GWAS

    G = snp_reader.read().standardize().val

    # full kernel
    # causal snps
    G_train_unnorm = G.take(rest_idx, axis=1)
    G_train_unnorm.flags.writeable = False
    G_train = DiagKtoN(G_train_unnorm.shape[0]).standardize(G_train_unnorm.copy())
    #G_train = 1./np.sqrt(G_train_unnorm.shape[1]) * G_train_unnorm
    G_train.flags.writeable = False

    G_test = G.take(test_idx, axis=1)
    G_test.flags.writeable = False
    
    
    gwas = FastGwas(G_train, G_test, y, delta=delta, cov=np.hstack((covar['vals'].copy(),np.ones((len(snp_reader.iid), 1)))), mixing=0.0)
    gwas.run_gwas()
    result["full_old"] = gwas.p_values_F
    
    # fs conditioned on full kernel
    select = FeatureSelectionInSample(max_log_k=7, order_by_lmm=True)

    fs_result["insample_cond_full"] = select.run_select(G_train_unnorm, G_train_unnorm, y, cov=G_pc_norm)
    best_k, fs_idx, best_mix, best_delta = fs_result["insample_cond_full"]
    print "best_k:", best_k, ", best_mix:", best_mix
    G_fs = G_train_unnorm.take(fs_idx, axis=1)
    G_fs *= 1./np.sqrt(best_k)
    G_fs.flags.writeable = False

    gwas = FastGwas(G_train, G_test, y, delta=delta, train_pcs=G_fs, mixing=best_mix)
    gwas.run_gwas()
    result["full_fs_low"] = gwas.p_values
    
    
    # fs unconditioned
    ########################
    out_fn = "tmp_pheno_%i.txt" % (sim_id)
    out_data = pd.DataFrame({"id1": snp_reader.iid[:,0], "id2": snp_reader.iid[:,1], "y": y})
    out_data.to_csv(out_fn, sep=" ", header=False, index=False)
    
    from pysnptools.snpreader import Bed
    fsd = create_feature_selection_distributable(Bed(snp_fn), out_fn, None, 0, "fs_out", snp_idx=rest_idx, include_all=True)
    fs_result["result_uncond_all"] = Local().run(fsd)
    best_k, best_delta, best_obj, best_snps = fs_result["result_uncond_all"]
    int_snp_idx = argintersect_left(snp_reader.sid[rest_idx], best_snps)
    fs_idx = np.array(rest_idx)[int_snp_idx]
    
    G_fs = G.take(fs_idx, axis=1)
    G_fs = 1./np.sqrt(G_fs.shape[1]) * G_fs
    G_fs.flags.writeable = False
    
    # fs all
    gwas = FastGwas(G_fs, G_test, y, delta=delta, train_pcs=None, mixing=0.0)
    gwas.run_gwas()
    result["fs_all"] = gwas.p_values

    # fs cov pcs
    gwas = FastGwas(G_fs, G_test, y, delta=delta, train_pcs=None, mixing=0.0, cov=G_pc_norm)
    gwas.run_gwas()
    result["fs_all_pcs_cov"] = gwas.p_values
    
    """ 

def draw_roc_curve(fpr, tpr, roc_auc, label):
    
    if len(fpr) > 1000:
        sub_idx = [int(a) for a in np.linspace(0, len(fpr)-1, num=1000, endpoint=True)]
        fpr, tpr = fpr[sub_idx], tpr[sub_idx]

    import pylab
    #pylab.semilogx(fpr, tpr, label='%s (area = %0.4f)' % (label, roc_auc))
    pylab.semilogx(fpr, tpr, label=label)
    #pylab.plot([0, 1], [0, 1], 'k--')
    pylab.xlim([0.0, 1.0])
    pylab.ylim([0.0, 1.0])
    #pylab.xlabel('False Positive Rate')
    pylab.xlabel('type I error',fontsize="large")
    #pylab.ylabel('True Positive Rate (Power)')
    pylab.ylabel('power',fontsize="large")
    #pylab.title('Receiver operating characteristic example')
    pylab.grid(True)
    pylab.legend(loc="lower right")

def compute_kernel_diag_from_G(G):
    # diag(K) = diag(G^TG) = \sum_{i,j) G_{i,j}^2
    return (G**2).sum()


def main():
    logging.basicConfig(level=logging.INFO)
    
    
    #snp_fn = "data/toydata.5chrom"
    snp_fn = "data/mouse/alldata"
    out_prefix = "results/mouse_"

    queue = "shared"
    #runner = Hadoop2(200, mapmemory=40*1024, reducememory=90*1024, mkl_num_threads=4, queue=queue)
    print "using snps", snp_fn
    #runner = LocalMultiProc(20)
    runner = Local()

    num_causals = 500
    num_repeats = 5
    num_pcs = 5
    
    # make this a tumple of function and kwargs
    methods = {execute_lmm}
    
    sc = LeaveTwoChrOutSimulation(snp_fn, out_prefix)
    sc.run(methods, num_causals, num_repeats, num_pcs, "mouse_", runner)

if __name__ == "__main__":
    main()
