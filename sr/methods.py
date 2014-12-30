"""
examples for method definitions
"""

import uuid
import logging
import numpy as np
import pandas as pd

from fastlmm.util.runner import Local, Hadoop2, LocalMultiProc
from fastlmm.util import distributed_map
from fastlmm.association.LocoGwas import FastGwas
from fastlmm.association import single_snp
from fastlmm.feature_selection.feature_selection_two_kernel import FeatureSelectionInSample
from fastlmm.util.util import argintersect_left
from pysnptools.standardizer import DiagKtoN


def execute_lmm(test_snps, pheno, G0, covar):
    
    result = {}
    fs_result = {}
    
    result["full"] = single_snp(test_snps, pheno, G0=G0, covar=covar).sort(["Chr", "ChrPos"])["PValue"].as_matrix()

    return result, fs_result


def execute_linear_regression(test_snps, pheno, G0, covar):
    """
    implementation of linear regression with and without covariates
    """
    
    result = {}
    fs_result = {}
    
    # linear regression with causals as covariates
    from fastlmm.inference.linear_regression import f_regression_cov
    G_test = test_snps.read().standardize().val
    _, result["linreg"] = f_regression_cov(G_test.copy(), pheno["vals"].copy(), np.ones((len(pheno["vals"]),1)))
    _, result["linreg_cov_pcs"] = f_regression_cov(G_test.copy(), pheno["vals"].copy(), covar["vals"].copy())
    
    return result, fs_result
    
    
def execute_dual_fs(test_snps, pheno, G0, covar):
    """
    implementation of dual-kernel feature selection
    """
    
    result = {}
    fs_result = {}
    
    
    # extract data
    G_test = test_snps.read().standardize().val
    G_train_unnorm = G0.read().standardize().val
    
    # fs conditioned on full kernel
    select = FeatureSelectionInSample(max_log_k=7, order_by_lmm=True)
    fs_result["insample_cond_full"] = select.run_select(G_train_unnorm, G_train_unnorm, pheno["vals"], cov=covar["vals"])
    best_k, fs_idx, best_mix, best_delta = fs_result["insample_cond_full"]
    print "best_k:", best_k, ", best_mix:", best_mix

    # set up foreground kernel
    G1 = G0[:,fs_idx]
    
    result["full_fs_low"] = single_snp(test_snps, pheno, G0=G0, covar=covar, G1=G1, mixing=best_mix).sort(["Chr", "ChrPos"])["PValue"].as_matrix()

    return result, fs_result

    
def execute_fs(test_snps, pheno, G0, covar):
    """
    run feature selection
    """
    
    result = {}
    fs_result = {}
    
    # fs unconditioned
    ########################
    tmp_uuid = str(uuid.uuid4())[0:13]
    out_fn = "tmp_pheno_%s.txt" % (tmp_uuid)
    out_data = pd.DataFrame({"id1": G0.iid[:,0], "id2": G0.iid[:,1], "y": pheno["vals"]})
    out_data.to_csv(out_fn, sep=" ", header=False, index=False)
    
    # write out covariates
    items = [
                ('id1', G0.iid[:,0]),
                ('id2', G0.iid[:,1]), 
            ]
    
    items += [("pc_%i" % i, covar["vals"][:,i]) for i in xrange(covar["vals"].shape[1])]
    cov_df = pd.DataFrame.from_items(items)
    cov_fn = "tmp_cov_%s.txt" % (tmp_uuid)
    cov_df.to_csv(cov_fn, sep=" ", header=False, index=False)
    
    #TODO: fix include_all!!
    fsd = create_feature_selection_distributable(G0, out_fn, None, 0, "fs_out", include_all=False, cov_fn=cov_fn)
    fs_result["result_uncond_all"] = Local().run(fsd)
    best_k, best_delta, best_obj, best_snps = fs_result["result_uncond_all"]
    fs_idx = argintersect_left(G0.sid, best_snps)
    
    G_fs = G0[:,fs_idx]
    
    result["fs_all"] = single_snp(test_snps, pheno, G0=G_fs).sort(["Chr", "ChrPos"])["PValue"].as_matrix()
    result["fs_all_cov"] = single_snp(test_snps, pheno, G0=G_fs, covar=covar).sort(["Chr", "ChrPos"])["PValue"].as_matrix()

    return result, fs_result


def create_feature_selection_distributable(snp_reader, phen_fn, pc_fn, num_pcs_kernel, output_prefix, cov_fn=None, include_all=True):

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
            k_values.append(snp_reader.sid_count)
        delta_values = np.logspace(-5, 10, endpoint=True, num=num_steps_delta, base=np.exp(1))

        if pc_fn is None:
            assert num_pcs_kernel == 0
            logging.info("feature selection: no PCs specified, disabling loop over mixing parameter")

        strategy = "insample_cv"
        select_by_ll = True

        # go!
        feature_selector = FeatureSelectionStrategy(snp_reader, phen_fn, num_folds, random_state=random_state, num_snps_in_memory=num_snps_in_memory, interpolate_delta=False, cov_fn=cov_fn)
        perform_selection_distributable = psd.PerformSelectionDistributable(feature_selector, k_values, delta_values, strategy, output_prefix, select_by_ll)

        return perform_selection_distributable
    
