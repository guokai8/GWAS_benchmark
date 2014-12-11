import fastlmm.util.runner as run
import numpy as np
import simulator as sim
import itertools
import numpy as np
import cPickle
import sys

class ParamSweep(object):
    """
    This is a super simple class that shows how to use the util.runner library.
    """

    def __init__(self,tempdirectory, debug=False,short_fn=0,output_prefix="out",numIndividuals=2000,num_folds=10, seed=5,penalty=0.0):

        if debug:#small fast local run
            self.FSTs = [0.01]
            self.fracSibs=[0.3]
            self.h2s = np.arange(0.1, 0.2, 0.1)
            self.var_hidden = np.arange(0.3, 0.6, 0.3)
            self.num_causal = [10,100]
        else:#big param sweep for cluster
            self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            self.fracSibs=np.array([0.0,0.05,0.1,0.2])
            self.h2s = np.arange(0.1, 0.7, 0.1)
            self.var_hidden = np.arange(0.0, 1.0, 0.3)
            self.num_causal = np.array([10, 50, 100, 500, 1000])

            #self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            #self.fracSibs=np.array([0.2])
            #h2s = np.array([.1])
            #var_hidden = np.array([0])
            ##self.num_causal = np.array([10])
            #self.num_causal = np.array([200])
            ##self.num_causal = np.array([10000])
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
        results_df = sim.sim_main(options, args)
        result = {
               "options":options,
               "res":results_df
               }
        return result

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
    debug = False
    short_fn=1
    seed=int(sys.argv[1])
    num_folds=int(sys.argv[2])
    penalty=float(sys.argv[3])

    from simulation.carlk4_runner import *

    job_env = ParamSweep(tempdirectory="te2", debug=debug,short_fn=short_fn,output_prefix="out",numIndividuals=4000, num_folds=num_folds,seed=seed,penalty=penalty)
    if "local" in sys.argv:
        runner = run.Local()
    elif "LocalInParts" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.LocalInParts(0,workcount,mkl_num_threads=None)
    elif "hadoop" in sys.argv:
        workcount = len(job_env.param_sweeps)
        if "default" in sys.argv:
            queue="default"
        else:
            queue="shared"
        runner = run.Hadoop(workcount, mapmemory=32*1024,reducememory=8*1024,mkl_num_threads=8,queue=queue)
    elif "debug" in sys.argv:
        short_fn=0
        debug = True
        workcount = 2
        runner = run.LocalMultiProc(workcount)
    elif "hpc" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.HPC(workcount, 'RR1-N13-09-H44',r'\\msr-arrays\scratch\msr-pool\eScience3',priority='Normal',unit='node')
    else:
        raise NotImplementedError("please specify where to run")
        #import pdb;pdb.set_trace()
    result = runner.run(job_env)
    
    #on Hadoop the result does not get returned
    file = open(job_env.resultpickle,'r')
    result = cPickle.load(file)
    file.close()
    print result
