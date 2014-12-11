import fastlmm.util.runner as run
import numpy as np
import simulator as sim
import itertools
import numpy as np
import cPickle
import sys

class ParamSweepX(object):
    """
    This is a super simple class that shows how to use the util.runner library.
    """

    def __init__(self,tempdirectory, short_fn=0,output_prefix="out",numIndividuals=2000,num_folds=10, seed=5,penalty=0.0):

        # 3 or 4 files with different seeds
        # change sample size to .8 * 4000
        # then run PC selection code


        self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
        self.fracSibs=np.array([0])
        self.h2s = [.1]
        self.var_hidden = [0]
        self.num_causal = np.array([10])
        self.tempdirectory=tempdirectory
        self.output_prefix=output_prefix
        self.resultpickle=self.output_prefix+"_res.pickle"
        self.short_fn = short_fn
        self.numIndividuals=numIndividuals
        self.num_folds=num_folds
        self.seed=seed
        self.penalty=penalty
        self.param_sweeps = list(itertools.product(self.FSTs, self.fracSibs, self.h2s, self.var_hidden, self.num_causal))
        print "Param sweep processed, this will result in %d jobs" % len(self.param_sweeps)
    
    def __str__(self):
        return self.tempdirectory

    def copyinputs(self, copier):
        pass

    def copyoutputs(self, copier):     
        copier.output(self.resultpickle)
        pass

    @property
    def work_count(self):
        return len(self.param_sweeps)

    def work_sequence(self):
        for i_p, p in enumerate(self.param_sweeps):
            yield lambda p=p:self.dowork(params=p,tasknum=i_p)

    def dowork(self, params,tasknum):
        fst, fracSibs, h2, hh2, causal = params
        options, args = sim.parseArgs()
        # here we just call parseargs because I'm not sure whether
        # the default params are set in the init or not. It's just for safety
        # we override the important stuff anyway
        options.h2 = h2
        options.fracSibs = fracSibs
        options.csnps = causal
        options.fst = fst
        options.var_hidden = hh2
        options.short_fn=self.short_fn
        options.numIndividuals=self.numIndividuals
        options.num_folds=self.num_folds
        options.randomseed=self.seed
        options.penalty=self.penalty

        snps,Y,simsnps,simphen,i_SNPs,nonchild_index_list = sim.generate_data(options,args)

        #import fastlmm.external.sklearn.externals.decomposition.mainpca as mainpca
        import simulation.mainpca as mainpca
        geno = {'snps':snps}
        bestPCCount = mainpca.FindBestPCCount(geno, predictSNPs=True,useChrom=False)
        return bestPCCount

    def reduce(self, result_sequence):
        results_array = []
        for result in result_sequence:
            results_array.append(result)
        #pickle the results_array so that copier can copy them
        file = open(self.resultpickle,"wb")
        cPickle.dump(results_array,file)
        file.close()
        # return final (pickleable) object
        return results_array

if __name__ == "__main__":
    short_fn=1
    seed=int(sys.argv[1])
    num_folds=int(sys.argv[2])
    penalty=float(sys.argv[3])

    from simulation.carlk_runner import *

    iidFraction = .8
    #iidFraction = .05
    job_env = ParamSweepX(tempdirectory="te2", short_fn=short_fn,output_prefix="out",numIndividuals=4000 * iidFraction, num_folds=num_folds,seed=seed,penalty=penalty)
    if "local" in sys.argv:
        runner = run.Local()
    elif "hadoop" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.Hadoop(workcount, mapmemory=32*1024,reducememory=8*1024,mkl_num_threads=8)
    elif "debug" in sys.argv:
        #short_fn=0
        #debug = True
        workcount = 2
        runner = run.LocalMultiProc(workcount)
    elif "hpc" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.HPC(workcount, 'RR1-N13-09-H44',r'\\msr-arrays\scratch\msr-pool\eScience3',priority='AboveNormal',unit='node')
    else:
        raise NotImplementedError("please specify where to run")
        #import pdb;pdb.set_trace()
    result = runner.run(job_env)
    
    #on Hadoop the result does not get returned
    file = open(job_env.resultpickle,'r')
    result = cPickle.load(file)
    file.close()
    print result
