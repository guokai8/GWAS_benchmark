import fastlmm.util.runner as run
import numpy as np


class SimpleExample(object):
    """
    This is a super simple class that shows how to use the util.runner library.
    """

    def __init__(self, number_of_jobs, print_intermediate_results=False):
        # where number_of_jobs can for example be: cv folds, param_sweeps etc.
        # you can pass whatever parameter or argument you want, but it probably 
        # won't accept functions as input, because it must be pickleable.
        # self.jobs is an array of input values (here vals 0 to 4)
        self.jobs=np.arange(number_of_jobs)
        self.print_intermediate = print_intermediate_results

    def copyinputs(self, copier):
        # we don't exactly know what copier is. Carl: can you please add a docstring 
        # here or in the parent class to explain a bit? we couldn't find "copier" 
        pass

    def copyoutputs(self, copier):
        # same as above
        pass    

    @property
    def work_count(self):
        # this returns the effective number of jobs
        return len(self.jobs)

    def work_sequence(self):
    # this is a python generator. Cycles through the job list and calls the dowork()
    # or whatever other name. This is a mapper.
        for job in self.jobs:
            yield lambda job=job:self.dowork(job) # this is a generator

    def dowork(self, input):
        r = np.sum(np.arange(input))
        return r

    def reduce(self, result_sequence):
        results_array = []
        for result in result_sequence:
            results_array.append(result)
            if self.print_intermediate:
                print result
        
        # return final (pickleable) object
        return results_array


runner = run.Local()
job_env = SimpleExample(5)
result = runner.run(job_env)
print result