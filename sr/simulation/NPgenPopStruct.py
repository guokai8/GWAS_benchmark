import random
import math
import os
import pdb
import subprocess
random.seed()
import time
from optparse import OptionParser

try:
    import numpy as np
    usenumpy = True
except:
    usenumpy = False

def runCommand(command, arguments, showOutput=True):
	# if (getOutput): p.StartInfo.RedirectStandardOutput = True
	#print 'Running:', command, arguments
	# p.StartInfo.FileName = command
	# p.StartInfo.Arguments  = arguments
	# p.Start()
	# if (getOutput): return p.StandardOutput.ReadToEnd()
	# p.WaitForExit()

	#print '\n'+command + " " + arguments+'\n'
	if showOutput: subprocess.call([command]+arguments.split())
	#else:
	#	with open(os.devnull, "w") as fnull:
	#		subprocess.call([command]+arguments.split(), stdout = fnull, stderr = fnull)


	
def createFrqFile(freqs, fileName, numCausativeSNPs, numHiddenSnps):
	f = open(fileName, 'w')	
	numCausativeSnpsCreated=0
	numSnpsCreated=0
	for i in range(len(freqs[0])):
		snpName = "snp"+str(i+1)
		if (numCausativeSnpsCreated < numCausativeSNPs-numHiddenSnps):
			snpName = 'c'+snpName
			numCausativeSnpsCreated += 1
			isCausative = True
		if (i >= numCausativeSNPs-numHiddenSnps and i < numCausativeSNPs):continue
		if (numSnpsCreated >= numSnps+numCausativeSNPs-numHiddenSnps): snpName = 'd'+snpName	
		f.write(' '.join([str(x) for x in [snpName, freqs[0][i], freqs[1][i], abs(freqs[0][i]-freqs[1][i])]])+'\n')
		numSnpsCreated+=1
	f.close()


def generatePhenotype(options,mafCounts,numCausativeSNPs,sampleSize,snpEffectSizes,epsilon,freqs,covars):
    #pdb.set_trace()
    liabilities = np.zeros(sampleSize)
    sumBeta2 = 0.0
    
    #SNPs*snpEffectSizes
    liabilities += ((mafCounts-2.0*freqs)/np.sqrt(2.0*freqs*(1.0-freqs))).dot(snpEffectSizes)
    sumBeta2+=np.sum(snpEffectSizes*snpEffectSizes)
    
    #Add covariate terms (covars*options.covarCoeff)
    liabilities += covars.sum(1)*options.covarCoeff
    sumBeta2 += covars.shape[1]*options.covarCoeff*options.covarCoeff

    #add interaction term
    if (options.interactionCoeff > 0.001 or options.interactionCoeff < -0.001):
        raise Exception("not implemented")
        pass
    
    #add noise
    if (sumBeta2>0.0): liabilities *= math.sqrt((1.0-epsilon*epsilon) / sumBeta2)
    residual = np.random.normal(0.0, epsilon,liabilities.shape)
    liabilities += residual
    return liabilities

def generateIndividualPhenotype(options,individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon, freqs, covars):	
    '''
    generate the continuous liability for a  single individual
    noise ~ N(0, epsilon)
    liability = SNPs*snpEffectSizes + covars*options.covarCoeff + term*SNPsoptions.interactionCoeff + noise
    ------------------------------------------------------------------------------- 
    freqs:      N-dim array of allele frequencies
    covars :    covariates
    options:    structure
                options.covarCoeff (constant)
                options.interactionCoeff (constant)
    -------------------------------------------------------------------------------
    return liability
    '''
    liability = 0.0
    sumBeta2 = 0.0
    #SNPs*snpEffectSizes
    for i in range(numCausativeSNPs):
        beta = snpEffectSizes[i]
        sumBeta2 += beta*beta
        liability += beta*(individualMafCounts[i]-2.0*freqs[i]) / math.sqrt(2.0*freqs[i]*(1.0-freqs[i]))
        
    #Add covariate terms (covars*options.covarCoeff)
    for covar in covars:
        beta = options.covarCoeff
        sumBeta2 += beta*beta
        liability += beta*covar #We assume all covariates already have zero mean and unit variance

    #add interaction term
    if (options.interactionCoeff > 0.001 or options.interactionCoeff < -0.001):
        beta = options.interactionCoeff
        sumBeta2 += beta*beta
        term = (individualMafCounts[0]-2.0*freqs[0]) * (individualMafCounts[1]-2.0*freqs[1])
        term /= math.sqrt(2.0*freqs[0]*(1.0-freqs[0]) * 2.0*freqs[1]*(1.0-freqs[1]))
        liability += beta*term

    if (sumBeta2>0.0): liability *= math.sqrt((1.0-epsilon*epsilon) / sumBeta2)
    residual = random.gauss(0, epsilon)
    liability += residual
    return liability

def generateSnps(numSnps, numPeople, mafs):
    randomNumbers = np.zeros((numPeople,len(mafs)),dtype = 'int8')
    randomNumbers[np.random.random((numPeople,len(mafs)))<mafs]+=1
    randomNumbers[np.random.random((numPeople,len(mafs)))<mafs]+=1
    return randomNumbers
    


def generateIndividualSnps(numSnps, mafs):
    if usenumpy:#faster version requires numpy
        randomNumbers = np.zeros(len(mafs),dtype = 'int8')
        randomNumbers[np.random.random(len(mafs))<mafs]+=1
        randomNumbers[np.random.random(len(mafs))<mafs]+=1
        return randomNumbers.tolist()
    else:#slower version can be without numpy
        randomNumbers = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in mafs]
        return randomNumbers

def generateCovars(freq, numAlleles,numPeople,covarCoeff,covarStd):
	if (abs(covarCoeff) < 1e-20): return np.zeros((numPeople,0))
	if (covarStd < 1e-20): covars = (np.random.random(numPeople,1)-0.5)*sqrt12
	else: covars = ((numAlleles-2*freq) / math.sqrt(2.0*freq*(1.0-freq)) + np.random.normal(loc = 0.0, scale = options.covarStd,size = (numPeople,1))) / (1.0+covarStd)
	return covars


def generateIndividualCovars(freq, numAlleles):
	if (abs(options.covarCoeff) < 1e-20): return []
	if (options.covarStd < 1e-20): covar = (random.random()-0.5)*sqrt12
	else: covar = ((numAlleles-2*freq) / math.sqrt(2.0*freq*(1.0-freq)) + random.gauss(0, options.covarStd)) / (1.0+options.covarStd)
	return [covar]
	
def createCovarsFile(covarsList, fileName):
	f = open(fileName, 'w')
	f.write('FID IID ')
	for i in range(len(covarsList[0])): f.write('covar'+str(i+1)+ ' ')
	f.write('\n')
	for personCounter in range(len(covarsList)):
		line = 'FAM1 person'+str(personCounter+1) + ' ' + ' '.join([str(x) for x in covarsList[personCounter]])
		f.write(line+'\n')
	f.close()
	
def createConfoundersFile(snpsList, liabilitiesList, fileName, numCausativeSNPs, numHiddenSnps, popList):
	f = open(fileName, 'w')
	f.write('FID IID ')
	for i in range(numHiddenSnps): f.write('confounder'+str(i+1)+ ' ')
	f.write('\n')
	
	for personCounter in range(len(liabilitiesList)):
		line = 'FAM1 person'+str(personCounter+1) + ' ' + str(popList[personCounter])
		i = -1
		for maf in snpsList[personCounter]:
			i+=1
			#if (i>= (numCausativeSNPs-numHiddenSnps) and i<numCausativeSNPs): line += ' ' + str(maf)
				
		f.write(line+'\n')
	f.close()
	

	
def createPedFile(liabilitiesList, snpsList, fileName, numCausativeSNPs, numHiddenSnps, numLD):
	f = open(fileName, "w")	
	for personCounter in range(len(liabilitiesList)):
		line = 'FAM1 person'+str(personCounter+1)+' 0 0 1 '+str(liabilitiesList[personCounter])
		i = -1
		randomDraws = [random.random()<0.02 for r in range(2*numLD*len(snpsList[personCounter]))]
		drawIndex=0
		for maf in snpsList[personCounter]:
			i+=1
			if (i>= (numCausativeSNPs-numHiddenSnps) and i<numCausativeSNPs): continue
			if (maf==0): line += ' 1 1'
			elif (maf==1): line += ' 1 2'
			elif (maf==2): line += ' 2 2'

			#create LD markers
			if (i >= numCausativeSNPs or True):				
				for LDMArkerInd in range(numLD):
					if (maf==0): snp1, snp2 = 0, 0
					elif (maf==2): snp1, snp2 = 1, 1
					else: snp1, snp2 = 1, 0
					p1 = randomDraws[drawIndex]
					drawIndex+=1
					p2 = randomDraws[drawIndex]
					drawIndex+=1
					if p1: snp1=1-snp1
					if p2: snp2=1-snp2
					line += ' ' + str(snp1+1) + ' ' +str(snp2+1)
		f.write(line+'\n')
		
	f.close()
		

def createMapFile(numSnps, numDiffSnps, numCausativeSNPs, fileName, numHiddenSnps, numLD):
	f = open(fileName, 'w')
	
	numCausativeSnpsCreated = 0
	numSnpsCreated = 0
	
	chromosome = '1'
	pos = 0.00
	base = 0
	for i in range(numSnps+numCausativeSNPs+numDiffSnps):
		isCausative = False
		snpName = "snp"+str(i+1)
		if (numCausativeSnpsCreated < numCausativeSNPs):
			snpName = 'c'+snpName
			numCausativeSnpsCreated += 1
			isCausative = True
		if (numSnpsCreated >= numSnps+numCausativeSNPs): snpName = 'd'+snpName		
		lineArr = [chromosome, snpName, str(pos), str(base)]
		if (i< (numCausativeSNPs-numHiddenSnps) or i>=numCausativeSNPs):
			f.write(' '.join(lineArr)+'\n')		

			#create LD markers
			if (not isCausative or True):
				for LDMArkerInd in range(numLD):
					ldsnpName = snpName + '_'+str(LDMArkerInd+1)
					lineArr = [chromosome, ldsnpName, str(pos+0.01*(LDMArkerInd+1)), str(base+1+LDMArkerInd)]
					f.write(' '.join(lineArr)+'\n')
		
		
		pos += 5
		base += (numLD+10) 		
		numSnpsCreated += 1
	f.close()
	
def writeGenesFile(fileName, numSnps, numDiffSnps, numCausativeSNPs, chrom, posDiff):
	f = open(fileName, 'w')
	totalNumSnps = numSnps+numDiffSnps+numCausativeSNPs
	maxDist = ((totalNumSnps-1)//options.geneSize)*posDiff
	pos = posDiff-10
	while True:
		f.write(str(chrom) + ' ' + str(pos) + '\n')
		pos += posDiff
		if (pos >= maxDist):
			if (pos-posDiff < maxDist): f.write(str(chrom) + ' ' + str(pos) + '\n')			
			break
	f.close()
	
def createPhenotypeFile(liabilitiesList, filename, binary=False, threshold=0.0, addColumns=False):
	f = open(filename, 'w')
	for personCounter in range(len(liabilitiesList)):	
		phenotype = liabilitiesList[personCounter]
		if (binary): phenotype = [1,2][phenotype>=threshold]
		if (addColumns): lineArr = ["FAM1", "person"+str(personCounter+1), '0', '0', '1', str(phenotype)]
		else: 			 lineArr = ["FAM1", "person"+str(personCounter+1), str(phenotype)]
		f.write(' '.join(lineArr)+'\n')
	f.close()
	
	
	
def findThreshold(diseasePrevalence, sampleSize, numCausativeSNPs, snpEffectSizes, epsilon, freqs, realFreqs,covarCoeff,covarStd):
    if 1:
        lastObservedSnpIndex = options.csnps-1 - options.numHiddenSnps
        mafCounts = generateSnps(numCausativeSNPs,sampleSize,freqs)
        if (lastObservedSnpIndex >= 0): 
            covars = generateCovars(freqs[lastObservedSnpIndex],  mafCounts[:,lastObservedSnpIndex], sampleSize,covarCoeff,covarStd)
        else:
            covars = np.zeros((sampleSize,0))
        liabilities = generatePhenotype(options,mafCounts,numCausativeSNPs,sampleSize,snpEffectSizes,epsilon,realFreqs,covars)
        prevIndex = int(round(sampleSize * (1-diseasePrevalence)))
        liabilities.sort()
        return liabilities[prevIndex]
    if 0:
        liabilitiesArr = []
        for i in range(sampleSize):
            if (i % 100000 == 0 and i>0): print 'generating individual', i, 'out of', sampleSize		
            individualMafCounts = generateIndividualSnps(numCausativeSNPs, freqs)
        
            #Add covariates	
            covars_list = []
            if (lastObservedSnpIndex >= 0): covars = generateIndividualCovars(freqs[lastObservedSnpIndex], individualMafCounts[lastObservedSnpIndex])
            
            individualLiability = generateIndividualPhenotype(options,individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon, realFreqs, covars_list)
            liabilitiesArr.append(individualLiability)	
        liabilitiesArr.sort()
        prevIndex = int(round(sampleSize * (1-diseasePrevalence)))
        return liabilitiesArr[prevIndex]

def generateRandomAlleles(fst, options,freqs):    
    numSnps = options.numSnps
    numCausativeSNPs = options.csnps
    numDiffSnps = options.numDiff
    N = options.numIndividuals
    casesFrac = options.caseFrac
    hidCoeff = options.hidCoeff
    numCases = int(N*casesFrac)
    numControls = N - numCases
    #Find the subpopulation with the required number of cases
    numCaseCaseSiblings = int(options.numSiblingsRatio * numCases)
    numCaseControlSiblings = numCaseCaseSiblings
    numControlControlSiblings = numCaseCaseSiblings
    siblingPop = -1
    maxPopCases = -1
    for popIndex in range(len(fst)):
        if (len(popCasesLiabilities[popIndex]) > maxPopCases):
            siblingPop = popIndex
            maxPopCases = len(popCasesLiabilities[popIndex])

    #Generate random alleles
    for popIndex in range(len(fst)):
        caseIndex = 0
        controlIndex = 0
        #Generate siblings
        if (popIndex == siblingPop or True):
            while (numCaseCaseSiblings > 0):
                (fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
                for i in range(options.numSibs):
                    popCasesSnps[popIndex][caseIndex] += sibSnps[i]
                    caseIndex+=1
                    numCaseCaseSiblings -= 1
            while (numControlControlSiblings > 0):
                (fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
                for i in range(options.numSibs):
                    popControlsSnps[popIndex][controlIndex] += sibSnps[i]
                    controlIndex+=1
                    numControlControlSiblings -= 1
            while (numCaseControlSiblings > 0 and controlIndex < len(popControlsSnps[popIndex])):
                (fatherSnps, motherSnps, sibSnps) = generateBrothers(options.numSibs, freqs, popIndex, numCausativeSNPs)
                for sibIndex in xrange(0, options.numSibs, 2):
                    popCasesSnps[popIndex][caseIndex] += sibSnps[sibIndex]
                    caseIndex+=1
                    popControlsSnps[popIndex][controlIndex] += sibSnps[sibIndex+1]
                    controlIndex+=1
                    numCaseControlSiblings -= 2

        #Generate non-siblings controls
        while (controlIndex < len(popControlsSnps[popIndex])):
            popControlsSnps[popIndex][controlIndex] += [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
            controlIndex+=1
        #Generate non-siblings cases
        while (caseIndex < len(popCasesSnps[popIndex])):
            popCasesSnps[popIndex][caseIndex] += [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
            caseIndex+=1

def generateBrothers(numBrothers, freqs, popIndex, numCausativeSNPs):
    if usenumpy:#faster version requires numpy
        fatherSnps = np.zeros(len(freqs[popIndex][numCausativeSNPs:]),dtype = 'int8')
        fatherSnps[np.random.random(len(freqs[popIndex][numCausativeSNPs:]))<freqs[popIndex][numCausativeSNPs:]]+=1
        fatherSnps[np.random.random(len(freqs[popIndex][numCausativeSNPs:]))<freqs[popIndex][numCausativeSNPs:]]+=1
        motherSnps = np.zeros(len(freqs[popIndex][numCausativeSNPs:]),dtype = 'int8')
        motherSnps[np.random.random(len(freqs[popIndex][numCausativeSNPs:]))<freqs[popIndex][numCausativeSNPs:]]+=1
        motherSnps[np.random.random(len(freqs[popIndex][numCausativeSNPs:]))<freqs[popIndex][numCausativeSNPs:]]+=1
        brotherSnps = np.zeros((numBrothers,len(fatherSnps)),dtype = 'int8')
        for i in range(numBrothers):
            brotherSnps[i][0.5*fatherSnps>=np.random.random(len(fatherSnps))]+=1
            brotherSnps[i][0.5*motherSnps>=np.random.random(len(motherSnps))]+=1
        return (fatherSnps.tolist(), motherSnps.tolist(), brotherSnps.tolist())
    else:#slower version can be without numpy
        fatherSnps = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]
        motherSnps = [[0,1][random.random()<m] + [0,1][random.random()<m] for m in freqs[popIndex][numCausativeSNPs:]]        
        brotherSnps = []
        for i in range(numBrothers):
            childFatherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in fatherSnps]
            childMotherSnps = [[[0,1][s==2], random.randint(0,1)][s==1] for s in motherSnps]
            brotherSnps.append(map(int.__add__, childFatherSnps, childMotherSnps))
        return (fatherSnps, motherSnps, brotherSnps)

def parseArgs():
    parser = OptionParser()
    parser.add_option('--bfile', metavar='bfile', help='output bfile name (default = out)', default='out')
    parser.add_option('--csnps', metavar='csnps', type=int, default=50, help='number of causal SNPs')
    parser.add_option('--fst', metavar='fst',  type=float, default=0.01, help='Fst distance between the two populations')
    parser.add_option('--her', metavar='her', type=float, default=0.5, help='trait heritability (in liability space)')
    parser.add_option('--prev', metavar='prev',  type=float, default = 0.5, help='trait prevalence')
    parser.add_option('--caseFrac', metavar='caseFrac', type=float, default=0.5, help='fraction of cases')
    parser.add_option('--snpCoeff', metavar='snpCoeff', type=float, default=1.0, help='ratio of first causal SNP coefficient to the others (if effstd option is not used)')
    parser.add_option('--hidCoeff', metavar='hidCoeff', type=float, default=5.0, help='ratio of hidden causal SNP coefficient to the others')
    parser.add_option('--numSnps', metavar='numSnps', type=int, default=10000, help='number of non-causal SNPs')
    parser.add_option('--numIndividuals', metavar='numIndividuals', type=int, default=2000, help='number of individuals')
    parser.add_option('--numDiff', metavar='numDiff', type=int, default=100, help='number of unusually differentiated SNPs')
    parser.add_option('--numSiblingsRatio', metavar='numSiblingsRatio', type=float, default=0.2, help='fraction of case-case siblings out of population 2 (e.g., 0.2 means that in population 2, 40% of individuals are unrelated, 20% are case-case sibs, 20% are case-control sibs and 20% are control-controls sibs)')
    parser.add_option('--numHiddenSnps', metavar='numHiddenSnps', type=int, default=1, help='number of hidden SNPs')
    parser.add_option('--numLD', metavar='numLD', type=int, default=0, help='num SNPs in LD (not in use)')
    parser.add_option('--minFreq', metavar='minFreq', type=float, default=0.1, help='minimum minor allele frequency')
    parser.add_option('--computeOdds', metavar='computeOdds', type=int, default=0, help='compute odds ratio (not in use)')
    parser.add_option('--interactionCoeff', metavar='interactionCoeff', type=float, default=0.0, help='interaction coefficient between two first causal SNPs')
    parser.add_option('--compress', metavar='compress', type=int, default=0, help='compress files (requires tar and bz2)')
    parser.add_option('--covarCoeff', metavar='covarCoeff', type=float, default=0.0, help='coefficient of clinical covariate (0 means no clinical covariate)')
    parser.add_option('--covarStd', metavar='covarStd', type=float, default=0.0, help='if greater than 0, clinical covar value is drawn from a guassian centered on the number of minor alleles of the second causal SNP')
    parser.add_option('--printCorr', metavar='printCorr', type=int, default=0, help='print maf correlations between the two populations')
    parser.add_option('--diffSize', metavar='diffSize', type=float, default=-1, help='minor allele frequency difference between populations for unusually differentiated SNPs (-1 means 0.6, with maf for pop1 randomly drawn between 0 and 0.4)')
    parser.add_option('--numRare', metavar='numRare', type=int, default=0, help='num of rare causal SNPs')
    parser.add_option('--rareFreq', metavar='rareFreq', type=float, default=0.01, help='frequency of rare causal SNPs')
    parser.add_option('--geneSize', metavar='geneSize', type=int, default=10, help='number of SNPs in each gene (only useful when using rare causal SNPs)')
    parser.add_option('--effStd', metavar='effStd', type=float, default=0.0, help='if greater than zero, effect sizes are drawn from a zero-mean gaussian with this std')
    parser.add_option('--numSibs', metavar='numSibs', type=int, default=2, help='number of siblings in each sib-group')
    parser.add_option('--diffCause', metavar='diffCause', type=int, default=1, help='1 if causal SNPs are differentiated, 0 otherwise')
    (options, args) = parser.parse_args()
    return (options,args)

if __name__ == "__main__":
    t0 = time.time()
    print ("starting (%.2fs)"%(t0-time.time()))
    (options, args) = parseArgs()

    sqrt12 = math.sqrt(12.0)

    bfile = options.bfile
    numSnps = options.numSnps
    numCausativeSNPs = options.csnps
    numDiffSnps = options.numDiff
    N = options.numIndividuals
    casesFrac = options.caseFrac
    hidCoeff = options.hidCoeff
    numCases = int(N*casesFrac)
    numControls = N - numCases
    if (options.numSibs == 0): options.numSiblingsRatio = 0

    fst = [options.fst, options.fst]
    epsilon = [math.sqrt(1.0-options.her), math.sqrt(1.0-options.her)]
    prevalence = options.prev

    #Generate SNP effect sizes
    print ("Generate SNP effect sizes (%.2fs)"%(t0-time.time()))
    if (options.effStd == 0.0): #fixed effects
        snpEffectSizes = np.ones(numCausativeSNPs)
        try:
            snpEffectSizes[0] = options.snpCoeff
            snpEffectSizes[(numCausativeSNPs-options.numHiddenSnps):numCausativeSNPs]*=hidCoeff
        except: pass
    else:                       #random effects
        snpEffectSizes = math.sqrt(options.effStd)*np.random.randn(numCausativeSNPs)
    lastObservedSnpIndex = options.csnps-1 - options.numHiddenSnps
    if (options.covarStd>1e-20): snpEffectSizes[lastObservedSnpIndex] = 0.0 #that SNP corresponds to the covariate for some reason

    #Generate allele frequencies for the old and new populations
    print ("Generate allele frequencies for the old and new populations (%.2fs)"%(t0-time.time()))
    ancestralFreqs = np.random.random(numSnps+numCausativeSNPs)*(0.5-options.minFreq) + options.minFreq
    ancestralFreqs[1:min(options.numRare+1, options.csnps-1)] = options.rareFreq
    
    freqs = [np.zeros(numDiffSnps+numSnps+numCausativeSNPs) for x in range(len(fst))]
    fstCoeffs = np.array([(1.0-fst[popIndex]) / fst[popIndex] for popIndex in range(len(fst)) if fst[popIndex]>0])
    for popIndex in range(len(fst)):
        if fst[popIndex]<1e-20:
            freqs[popIndex] = ancestralFreqs
        else:
            #sample froam a beta distribution
            freqs[popIndex][0:numSnps+numCausativeSNPs] = np.random.beta(ancestralFreqs*fstCoeffs[popIndex],(1.0-ancestralFreqs)*fstCoeffs[popIndex])
            iNotOK = np.logical_or(freqs[popIndex][0:numSnps+numCausativeSNPs]<options.minFreq,freqs[popIndex][0:numSnps+numCausativeSNPs]>1.0-options.minFreq)
            while iNotOK.any():#redo the frequencies that failed the cutoff untill all pass
                freqs[popIndex][0:numSnps+numCausativeSNPs][iNotOK] = np.random.beta(ancestralFreqs[iNotOK]*fstCoeffs[popIndex],(1.0-ancestralFreqs[iNotOK])*fstCoeffs[popIndex])
                iNotOK = np.logical_or(freqs[popIndex][0:numSnps+numCausativeSNPs]<options.minFreq,freqs[popIndex][0:numSnps+numCausativeSNPs]>1.0-options.minFreq)


    #Generate allele frequencies for unusually differentiated SNPs
    print ("Generate allele frequencies for unusually differentiated SNPs (%.2fs)"%(t0-time.time()))    
    if (options.diffSize >= 0):
        freqs[0][numSnps+numCausativeSNPs::] = 0.5 - options.diffSize
        freqs[1][numSnps+numCausativeSNPs::] = 0.5 + options.diffSize
    else:
        freqs[0][numSnps+numCausativeSNPs::] = np.random.random(numDiffSnps)*0.4
        freqs[1][numSnps+numCausativeSNPs::] = freqs[0][numSnps+numCausativeSNPs::]+0.6
    print ("Generate aaffection tresholds (%.2fs)"%(t0-time.time()))    
    
    #Find the affection threshold
    threshold = 0.0
    for popIndex in range(len(fst)):	
        popThreshold = findThreshold(prevalence, int(3000.0/prevalence), numCausativeSNPs, snpEffectSizes, epsilon[popIndex], freqs[popIndex][:numCausativeSNPs], ancestralFreqs[:numCausativeSNPs],options.covarCoeff,options.covarStd)
        threshold = max(threshold, popThreshold)

    #Generate individuals
    print ("Generate individuals (%.2fs)"%(t0-time.time()))
    
    casesCounter = 0
    controlsCounter = 0
    casesLiabilities = np.zeros((numCases))
    controlLiabilities = np.zeros((numControls))
    casesSnps = np.zeros((numCases,numCausativeSnps),dtype = 'int8')
    controlSnps = np.zeros((numCases,numControlSnps),dtype = 'int8')
    if (lastObservedSnpIndex >= 0): 
        casesCovars = np.zeros((numCases,1))
        cotrolCovars =np.zeros((numControls,1))
    else:
        casesCovars = np.zeros((numCases,0))
        cotrolCovars =np.zeros((numControls,0))
    #this also subsamples individuals: (they are being thrown out if too many/little cases)
    while (casesCounter < numCases or controlsCounter < numControls):
        for popIndex in range(len(fst)):

            mafCounts = generateSnps(numCausativeSnps, numPeople, freqs[popIndex][:numCausativeSNPs])
            if (lastObservedSnpIndex >= 0): 
                covars = generateIndividualCovars(freqs[popIndex][lastObservedSnpIndex], individualMafCounts[lastObservedSnpIndex])
            else:
                covars = np.zeros((numPeople,0))
            liabilities = generatePhenotype()
            i_control = liabilities<treshold
            n_control = i_control.sum()
            n_case = numPeople - n_control
            if (controlsCounter < numControls):
                numControlsAdd = min(n_control,numControls-controlsCounter)
                casesSnps[controlsCounter:controlsCounter+numControlsAdd,:] = mafCounts[i_control[controlsCounter:controlsCounter+numControlsAdd]]
                controlsCounter += numControlsAdd
                pass
            if (casesCounter < numCases):
                numCasesAdd = min(n_case,numCases-casesCounter)
                casesSnps[casesCounter:casesCounter+numCasesAdd,:] = mafCounts[i_control[casesCounter:casesCounter+numCasesAdd]]
                controlsCounter += numControlsAdd

    casesCounter = 0
    controlsCounter = 0
    popCasesSnps = []
    popControlsSnps = []
    popCasesLiabilities = []
    popControlsLiabilities = []
    popControlsCovars = []
    popCasesCovars = []

    print ("Find the subpopulation with the required number of cases (%.2fs)"%(t0-time.time()))
    
    for popIndex in range(len(fst)):
        popCasesSnps.append([])
        popControlsSnps.append([])
        popCasesLiabilities.append([])
        popControlsLiabilities.append([])
        popControlsCovars.append([])
        popCasesCovars.append([])


    #this also subsamples individuals: (they are being thrown out if too many/little cases)
    while (casesCounter < numCases or controlsCounter < numControls):
        for popIndex in range(len(fst)):

            individualMafCounts = generateIndividualSnps(numCausativeSNPs, freqs[popIndex][:numCausativeSNPs])

            #Add covariates
            covars = []
            if (lastObservedSnpIndex >= 0): covars = generateIndividualCovars(freqs[popIndex][lastObservedSnpIndex], individualMafCounts[lastObservedSnpIndex])
        
            individualLiability = generateIndividualPhenotype(options,individualMafCounts, numCausativeSNPs, snpEffectSizes, epsilon[popIndex], ancestralFreqs[:numCausativeSNPs], covars)
            if (controlsCounter < numControls and individualLiability < threshold):
                controlsCounter += 1
                popControlsSnps[popIndex].append(individualMafCounts)
                popControlsLiabilities[popIndex].append(individualLiability)
                popControlsCovars[popIndex].append(covars)
            elif (casesCounter < numCases and individualLiability >= threshold):
                casesCounter += 1
                popCasesSnps[popIndex].append(individualMafCounts)
                popCasesLiabilities[popIndex].append(individualLiability)
                popCasesCovars[popIndex].append(covars)

    print [len(popControlsLiabilities[0]), len(popControlsLiabilities[1])], [len(popCasesLiabilities[0]), len(popCasesLiabilities[1])]


    print ("Generate random alleles (%.2fs)"%(t0-time.time()))
    t1 = time.time()
    generateRandomAlleles(fst, options,freqs)
    print ("done Generate random alleles (%.2fs)"%(t1-time.time()))
    pdb.set_trace()
    #Create combined lists for all individuals
    allLiabilities = []
    allSnps = []
    allCovars = []
    popList = []
    for popIndex in range(len(fst)):
        allLiabilities += popControlsLiabilities[popIndex]
        allLiabilities += popCasesLiabilities[popIndex]
        popList += [popIndex+1 for i in range(len(popControlsLiabilities[popIndex]) + len(popCasesLiabilities[popIndex]))]
        allSnps += popControlsSnps[popIndex]
        allSnps += popCasesSnps[popIndex]	
        allCovars += popControlsCovars[popIndex]
        allCovars += popCasesCovars[popIndex]
    
    #Create Plink Files
    print ("Creating plink files (%.2fs)"%(t0-time.time()))
    createMapFile(numSnps, numDiffSnps, numCausativeSNPs, bfile+'.map', options.numHiddenSnps, options.numLD)
    writeGenesFile(bfile+'.genes', numSnps, numDiffSnps, numCausativeSNPs, 1, 10*options.geneSize)
    createPedFile(allLiabilities, allSnps, bfile+'.ped', numCausativeSNPs, options.numHiddenSnps, options.numLD)
    createPhenotypeFile(allLiabilities, bfile+'.phe', binary=True, threshold=threshold)
    createPhenotypeFile(allLiabilities, bfile+'.phe.liab', binary=False)
    #createPhenotypeFile(allLiabilities, bfile+'_pca.ped', binary=True, threshold=threshold, addColumns = True)
    createCovarsFile(allCovars, bfile+'.covar')
    createFrqFile(freqs, bfile+'.frq', options.csnps, options.numHiddenSnps)
    createConfoundersFile(allSnps, allLiabilities, bfile+'.confounders', numCausativeSNPs, options.numHiddenSnps, popList)
    print ("Running plink (converting to bed) (%.2fs)"%(t0-time.time()))
    runCommand('plink', '--noweb --silent --file ' + bfile + ' --make-bed --out ' + bfile)

    #Compress the files
    if (options.compress == 1):
        print 'Compressing the files...'
        runCommand('tar', '-cjf ' +  bfile+'.bz2 ' + ' '.join([bfile+'.phe', bfile+'.phe.liab', bfile+'.bed', bfile+'.bim', bfile+'.fam', bfile+'.covar', bfile+'.confounders']))

    if (options.printCorr == 1):
        x1, x2 = [], []
        for c in allSnps: x1.append(c[0]); x2.append(c[1])
        import numpy
        print 'corr:', numpy.corrcoef(x1, x2)[0,1]

    if (options.computeOdds == 1):
        #compute odd ratios for causal SNPs
        import scipy.stats
        import numpy
        from numpy import *
        sumBeta2=0
        for i in xrange(numCausativeSNPs): sumBeta2 += snpEffectSizes[i]**2
        snpVec = array([0.0 for i in xrange(numCausativeSNPs)])
        freqsArr = array(ancestralFreqs[:numCausativeSNPs])
        effVec = array(snpEffectSizes)
        currSnp = numCausativeSNPs-1
        snpAffectProbs = [[0.0,0.0,0.0] for i in xrange(numCausativeSNPs)]
        recipSqrtEff = numpy.sqrt(1.0/(2.0*freqsArr*(1.0-freqsArr)))
        sumBeta2 = numpy.sum(dot(freqsArr,freqsArr))
        liabCoeff = numpy.sqrt((1.0-epsilon[0]*epsilon[0]) / sumBeta2)

        while True:	
            geneticLiab = dot((snpVec-2*freqsArr)*recipSqrtEff, effVec) * liabCoeff
            probAff     = scipy.stats.norm.sf(threshold-geneticLiab, 0, epsilon[0])		
            prob = numpy.prod(freqsArr**snpVec * (1.-freqsArr)**(2.0-snpVec) * array([[1.0,2.0][s==1.0] for s in snpVec]))

            for i in xrange(len(snpVec)):		
                probWithoutSnp = prob / (freqsArr[i]**snpVec[i]  * (1-freqsArr[i])**(2-snpVec[i]) * [1,2][snpVec[i]==1.0])		
                snpAffectProbs[i][int(snpVec[i])] += probAff*probWithoutSnp

            snpVec[currSnp]+=1
            if (snpVec[currSnp] <= 2): continue

            while (snpVec[currSnp]>=2 and currSnp>=0):
                snpVec[currSnp] = 0.0
                currSnp-=1
            if (currSnp<0):break
            snpVec[currSnp]+=1
            currSnp=numCausativeSNPs-1

        f = open(bfile+'.coeff.txt', 'w')
        for i in xrange(numCausativeSNPs - options.numHiddenSnps):		
            odds = [math.log(p/(1-p)) for p in snpAffectProbs[i]]		
            beta = math.exp((odds[2] - odds[0]) / 2.0)
            f.write('csnp'+str(i+1) + ' ' + str(beta) + '\n')
        f.close()

