from simulation.pca import PCA, fast_logdet
import numpy as np
import fastlmm.pyplink.plink as plink
import fastlmm.util.preprocess as up
import scipy.linalg as la
import os
import sys
import cPickle
import time
import logging


def loadBed(bed_fn):
	geno = plink.readBED_quick(bed_fn, order='F')
	X = geno['snps']
	betaNotUnitVariance = False
	betaA = 1.0
	betaB = 1.0
	blocksize = 1000
	if X.dtype == np.float32:
		# wrapped C-reader will deliver float32 encoding
		from pysnptools.pysnptools.snpreader import wrap_plink_parser
		X = wrap_plink_parser.standardize(X,betaNotUnitVariance=betaNotUnitVariance,betaA=betaA,betaB=betaB)
		X = np.array(X, dtype=np.float64, order="C")
	else:
		X = up.standardize_block(X, blocksize=blocksize)
		X = np.array(X, dtype=np.float64, order="C")		
		
	geno['snps'] = X
	return geno
	
def FindBestPCCount(geno, predictSNPs, useChrom, k_values = np.arange(0,5+1,1)):
    randomstate = 1

    X = geno['snps']
    #k_values = np.arange(0,np.min([k_max + 1,X.shape[0],X.shape[1]]), 1)
    if max(k_values) >= X.shape[0] or max(k_values) >= X.shape[1]: raise Exception("The number of PCs search should be less than the # of rows and also the # of cols in the matrix")
    
    if predictSNPs: #Predict SNPs in instead of individuals
        X = X.T
    
    if useChrom :
        if not predictSNPs : raise Exception("if useChrom, then should also predictSNPs")
        pos = geno['pos']
        chromSet = set([p[0] for p in pos])
        folds = []
        for chr in chromSet:
            chrInd = [i for i in xrange(len(pos)) if pos[i][0] == chr]
            notChrInd = [i for i in xrange(len(pos)) if pos[i][0] != chr]
            folds.append([notChrInd, chrInd])
        n_folds = len(chromSet)
        if (n_folds < 2): raise Exception('Bed file has only one chromosome')
    else:
    	from sklearn.cross_validation import KFold
    	n_folds = 10
    	folds = KFold(X.shape[0], n_folds = n_folds, shuffle=True, random_state=randomstate)
    	
    scores = np.zeros((k_values.shape[0],n_folds))
    for i_fold, [train_idx,test_idx] in enumerate(folds):
    	logging.info("computing svd...")
    	t0 = time.time()	
    	#Utr,Str,Vtr = la.svd(X[train_idx],full_matrices=False)
    	Utr,Str,Vtr = None, None, None
    	t1 = time.time()
    	logging.info("done after %.4f seconds" % (t1 - t0))
    	logging.info('test set size: {0}'.format(len(test_idx)))
    	for i_k, k in enumerate(k_values):
    		pca = PCA(n_components = k)
    		Utr,Str,Vtr = pca._fit(X[train_idx],Utr,Str,Vtr)
    	
    		scores[i_k,i_fold] = pca.score(X[test_idx])
    		logging.info("{0},{1},{2}".format(i_fold, k, scores[i_k,i_fold]))

    normalizedMean = np.zeros(k_values.shape[0])
    for i_k in xrange(k_values.shape[0]):
        kMean = 0.0
        for i_fold, [train_idx,test_idx] in enumerate(folds):
            kMean += scores[i_k,i_fold]
        kMean /= float(n_folds)
        normalizedMean[i_k] = kMean	
    logging.info('normalized Means: {0}'.format(normalizedMean))
    
    bestNumPCs = k_values[normalizedMean.argmax()]
    return bestNumPCs


if __name__ == "__main__":
    # The script accepts two input .bed files and outputs a text file in a
    # covariates format.
    # The first .bed file corresponds to a sample with family members excluded, and
    # the second refers to
    # the same sample but with the family members retained.  The output file is the
    # results of projecting the
    # individuals from the second file to PC space using a trained PPCA model.


    # Inputs:
    # bed_fn_nofam: The base name of the .bed file close family members excluded
    # bed_fn_withfam: The plink file with family members included
    # out_basefn: The base name of the output covariates file
    # k_max - the maximum number of PCs to test.  Currently it's hard coded as 5.
    bed_fn_nofam = 'MS_fixdist_dup'
    bed_fn_withfam = 'MS_fixdist_dup'
    out_basefn = "MS_PCs"

    demo = False
    saveplot = False
    k_max = 10
    randomstate = 1


    # The script first takes the reduced data set and finds the optimal #PCs for
    # predicting SNPs, using cross-validation on the chromosomes.
    geno = loadBed(bed_fn_nofam)
    bestNumPCs = FindBestPCCount(geno, predictSNPs=True,useChrom=True)
    print 'best num PCs:', bestNumPCs


    # Afterwards, it trains a PPCA model on the reduced data set, using the optimal
    # #PCs,
    # but using SNPs as features this time (with all individuals and all SNPs
    # included in the training set).
    # i.e.  train a PCA model with the required number of PCs on all the data
    print "computing svd..."
    if predictSNPs: X = X.T
    t0 = time.time()
    #Utr,Str,Vtr = la.svd(X, full_matrices=False)
    Utr,Str,Vtr = None, None, None
    t1 = time.time()
    print "done after %.4f seconds" % (t1 - t0)
    pca = PCA(n_components = bestNumPCs)
    Utr,Str,Vtr = pca._fit(X,Utr,Str,Vtr)


    # Next, it loads the .bed file corresponding to the full data set (with family
    # members included),
    # and projects all individuals there to PCs space, using the trained PPCA model
    # from the previous stage.
    # i.e.  Now load the full bed file and project all individuals to the PCs space
    geno = loadBed(bed_fn_withfam)
    X = geno['snps']
    print 'Projecting individuals to PCs space...'
    #if predictSNPs: X = X.T
    X_fit = pca.transform(X)
    #if predictSNPs: X_fit = X_fit.T

    # Finally, it writes a text file in the FastLMM covariates format, with the PCs
    # of the projected individuals.
    print 'writing results to file...'
    f = open(out_basefn + '.cov', 'w')
    for i_ind in xrange(X_fit.shape[0]):	
        f.write(geno['iid'][i_ind][0] + ' ' + geno['iid'][i_ind][1] + ' ')
        f.write(' '.join([str(pc) for pc in X_fit[i_ind, :]]))
        f.write('\n')
    f.close()

    if saveplot:
        import pylab as pl
        pl.figure()
        pl.plot(k_values,scores.mean(1))
        pl.ylabel("xval log likelihood")
        pl.xlabel("number PCs")
        pl.title("best number PCs: %i Loglik:%.2f" % (k_values[scores.mean(1).argmax()], scores.mean(1).max()))
        pl.grid()
        pl.plot(k_values[scores.mean(1).argmax()],scores.mean(1).max(),"ro")
        pl.savefig(out_basefn + ".pdf")
        pl.close("all")
