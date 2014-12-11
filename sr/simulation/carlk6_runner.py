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
        #N4000S50000c1000h0.10s0.00p0.50F0.0100FH0.2000v0.00_9.meta.pickle.oos', 
        if 1:# FS and PS
            self.FSTs = np.array([0.0, 0.005, 0.01, 0.05, 0.1])
            self.fracSibs= np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9])
            self.h2s = np.arange(0.1, 0.7, 0.1) # it's a range!
            self.var_hidden = np.array([0.0, 0.3])
            self.num_causal = np.array([10, 50, 100, 500, 1000])
            print "AAA"
        elif 0:# FS and PS
            self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            self.fracSibs= np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9]) #np.array([0.0,0.05,0.1,0.2,0.4]) #
            self.h2s = np.arange(0.1, 0.7, 0.1) # it's a range!
            self.var_hidden = np.array([0.0]) #np.arange(0.0, 1.0, 0.3) # it's a range! # #
            self.num_causal = np.array([10, 50, 100, 500, 1000])
            print "AAA"
        elif 0:#small fast local run
            self.FSTs = np.array([0.0])
            self.fracSibs= np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9]) #np.array([0.0,0.05,0.1,0.2,0.4]) #
            self.h2s = np.arange(0.1, 0.7, 0.1) # it's a range!
            self.var_hidden = np.array([0.0, 0.3]) #np.arange(0.0, 1.0, 0.3) # it's a range! # #
            self.num_causal = np.array([10, 50, 100, 500, 1000])
            print "AAA"
            """
            self.FSTs = np.array([0., 0.005, 0.01, 0.05, 0.1])
            self.fracSibs= np.array([0.2]) #np.array([0.0,0.05,0.1,0.2])
            self.h2s = np.arange(0.1, 0.7, 0.1) # it's a range!
            self.var_hidden = np.arange(0.0, 1.0, 0.3) # it's a range!
            self.num_causal = np.array([50, 500]) #np.array([10, 50, 100, 500, 1000])
            """
        else:#big param sweep for cluster
            self.FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            self.fracSibs=np.array([0.0,0.05,0.1,0.2])
            self.h2s = np.arange(0.1, 0.7, 0.1) # it's a range!
            self.var_hidden = np.arange(0.0, 1.0, 0.3) # # it's a range!
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
        self.resultpickle=self.output_prefix+"_new_res.pickle"
        self.short_fn = short_fn
        self.numIndividuals=numIndividuals
        self.num_folds=num_folds
        self.seed=seed
        self.penalty=penalty
        self.param_sweeps = list(itertools.product(self.FSTs, self.fracSibs, self.h2s, self.var_hidden, self.num_causal)) 
        #self.param_sweeps = self.param_sweeps[::-1]
        print "Param sweep processed, this will result in %d jobs" % len(self.param_sweeps)
    
        #import pdb; pdb.set_trace(); print self.resultpickle
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
        options.make_binary = False #ck01082014
        options.vertex_cutoff = 0.1
        if fst == 0.0:
            options.fst_hidden = 0.0
            print "WARNING: SETTING fst_hidden to 0.0"
        #options.recompute = 1
        #print "WARNING: force recompute!!"
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

    #This is a little python code to create new binary phen.txt files from old real-valued phen.txt files.
    #import os
    #inputDir = "oldPheno"
    #outputDir = "."
    #for fileName in os.listdir(inputDir):
    #    with open(os.path.join(inputDir, fileName),"r") as infile:
    #        with open(fileName,"w") as outfile:
    #            for inLine in infile:
    #                fields = inLine.split('\t')
    #                fields[2] = "0" if float(fields[2]) < 0 else "1"
    #                outfile.write('\t'.join(fields)+"\n")


    debug = False
    #print "WARNING: debugging mode!!!"
    short_fn=1
    seed=int(sys.argv[1])
    num_folds=int(sys.argv[2])
    penalty=float(sys.argv[3])

    from simulation.carlk6_runner import *

    #!!ck01082014 Change "out" and tempdirectory
    job_env = ParamSweep(tempdirectory="te6", debug=debug,short_fn=short_fn,output_prefix="out6",numIndividuals=4000, num_folds=num_folds,seed=seed,penalty=penalty)
    if "local" in sys.argv:
        runner = run.Local()
    elif "LocalInParts" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.LocalInParts(1803,workcount,mkl_num_threads=None)
        ##parts = [42,405,413,415,419,429,432,433,434,442,445,447,450,465,472,489,493,498-1,500-1,506-1,515-1,524-1,528-1,529-1,531-1,548-1]
        #parts = [548-1]
        #runners = []
        #for part in parts:
        #    runner = run.LocalInParts(part,workcount,mkl_num_threads=None)
        #    runner.run(job_env)
        #    runner = 0
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
        workcount = len(job_env.param_sweeps) #RR1-N13-09-H44 #RR1-N13-16-H44
        runner = run.HPC(workcount, 'RR1-N13-16-H44',r'\\msr-arrays\scratch\msr-pool\eScience3',priority='highest',unit='node',min=150)
    elif "hpc2" in sys.argv:
        workcount = len(job_env.param_sweeps)
        runner = run.HPC(workcount, 'RR1-N13-09-H44',r'\\msr-arrays\scratch\msr-pool\eScience3',priority='abovenormal',unit='node',min=25)
    else:
        raise NotImplementedError("please specify where to run")
        #import pdb;pdb.set_trace()
    result = runner.run(job_env)
    
    #on Hadoop the result does not get returned
    file = open(job_env.resultpickle,'r')
    result = cPickle.load(file)
    file.close()
    print result
