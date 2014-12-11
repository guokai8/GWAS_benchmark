import pylab as pl
import cPickle
import scipy as sp
import numpy as np
import scipy.stats as st
import sys
import pandas
import matplotlib.lines as Line2D
import logging
import os

markers = ['x','v','o','+','*','<','>','d']
lines = ['-', '-', '--', ':','-.']


#filenames_in=["./../simres3/out_res.pickle","./../simres4/out_res.pickle","./../simres5/out_res.pickle","./../simres_6_mult/out_res.pickle"]

#filenames_in=["./../simres_25folds/out_res_3.pickle","./../simres_25folds/out_res_high3.pickle","./../simres_25folds/out_res_10k3.pickle","./../simres_25folds/out_res_4.pickle"]
#filenames_in=["z:/heckerman4-c/simres_25folds/out_res_3.pickle","z:/heckerman4-c/simres_25folds/out_res_high3.pickle","z:/heckerman4-c/simres_25folds/out_res_4.pickle"]

#!!ck01102014
#1
#filenames_in=[r"\\carlk3\cachebio\genetics\synthetic\carlk\out_res.12102013.pickle",r"\\carlk3\cachebio\genetics\synthetic\carlk\out_res_3.pickle"]


#filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed9_fst0_rank_lr.pickle"] #, r"C:\Users\chwidmer\Documents\results\synthetic\seed9_fixed_delta.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\seed2_fix.pickle"] #, r"C:\Users\chwidmer\Documents\results\synthetic\seed6_out6_new_res.pickle"] 
#,
#
# stage2:
#filenames_in=[r"C:\Users\chwidmer\Documents\results\synthetic\stage2\out_res.01082014.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\stage2\out4_res.noauto01102014.pickle"]
#filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\original_out6_res.pickle"] 
#                r"C:\Users\chwidmer\Documents\results\synthetic\out9_res.pickle"]

# stage3:
#filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed6_fst0.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\seed2_fst0.fh0.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\seed9_fst0_lmm_rank.pickle"] #  r"C:\Users\chwidmer\Documents\results\synthetic\seed6_out6_new_res.pickle", 

# stage4:
#filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed3_fixed.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\seed9_fixed_delta.pickle", r"C:\Users\chwidmer\Documents\results\synthetic\seed2_fix.pickle"]
#filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed5_fam10_lots_fs.pickle"]
#seed1_fst0_hiddenvar0.pickle
"""
filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed1_fst0_hiddenvar0.pickle", 
                r"C:\Users\chwidmer\Documents\results\synthetic\seed2_fst0_hiddenvar0.pickle", 
                r"C:\Users\chwidmer\Documents\results\synthetic\seed4_fst0_hiddenvar0.pickle", 
                r"C:\Users\chwidmer\Documents\results\synthetic\seed6_fst0_hiddenvar0.pickle", 
                r"C:\Users\chwidmer\Documents\results\synthetic\seed9_fst0_hiddenvar0.pickle", ]
"""

filenames_in = [r"C:\Users\chwidmer\Documents\results\synthetic\seed7_two_val_hv_true_causals.pickle"]

do_merge = True

xvalbar=[1e-7,1e-6,1e-5,1e-4,1e-3]

## These are for .90 intervals
#barlow=[4.16667E-08,
#8.22917E-07,
#9.45833E-06,
#9.83125E-05,
#0.000994688]
#barhigh=[1.66667E-07,
#1.17708E-06,
#1.05417E-05,
#0.000101688,
#0.001005313]

##These are for .99 intervals
#barlow=[2.08333E-08,
#7.39583E-07,
#9.16667E-06,
#0.000097375,
#0.000991688]
#barhigh=[1.97917E-07,
#1.28125E-06,
#1.08542E-05,
#0.000102646,
#0.001008333]

#For the .95 intervals
barlow=[3.125E-08,
7.91667E-07,
9.36458E-06,
9.79896E-05,
0.000993677]
barhigh=[1.77083E-07,
1.21875E-06,
1.06458E-05,
0.00010201,
0.001006333]


methodmap={
 "linreg_0":"Linreg",
 "linreg_auto": r"Linreg + PC_geno",
 "linreg": r"Linreg + PC_pheno",
 "lmm_ins_0":"LMM(all)",
 "lmm_ins_auto": r"LMM(all) + PC_geno", 
 "lmm_oos_auto": r"Standard PCg oos",
 "lmm_ins": r"LMM(all) + PC_pheno",
 "insample_0":"LMM(select)",
 "insample_auto": r"LMM(select) + PC_geno",
 "insample": r"LMM(select) + PC_pheno",
 "oosample_0":"SVoos no PC",
 "oosample_auto": r"SVoos PCg",
 "oosample": r"SVoos PCp",
 "lmm_fs_cond_full_auto": r"LMM(select+all)",
 "lmm_fs_cond_full_auto_2": r"LMM(select+all) sub2",
 "lmm_fs_cond_full_auto_8": r"LMM(select+all) sub8",
 "lmm_true_causals_auto": "LMM(true_causals)"
 }

datamap = {
           'h2':"Signal-to-noise ratio",
           'confounding_var':"Frac of noise variance due to pop structure",
           'fraction_Trios':"Fraction of individuals belonging to a trio",
           'number_causal':"Number of causal SNPs",
           'Fst': "$F_{st}$"
           }


#filenames_in=[r"\\carlk3\cachebio\genetics\synthetic\carlk01b\out_res.01082014.pickle"]
def getKey(key="lambda",method = None, mstring=None,index=None):
    if method is None:
        ret =  np.array([res[i]["res"][key][beststr[:,i]] for i in range(len(res))],dtype="float")
    elif mstring is not None:
        ret =  np.array([res[i]["res"][key][mstring] for i in range(len(res))],dtype="float")
    else:
        index_methods=[methods.index(m) for m in method]
        ret =  np.array([res[i]["res"][key][beststr[index_methods,i]] for i in range(len(res))],dtype="float")
    if index is not None:
        return ret[index]
    else:
        return ret

def computeTest(key="Power corr @ 1.0e-07",test="wilcox",index=None): #,zero_method="pratt"
    vals = getKey(key=key,index=index)
    n = vals.shape[1]

    stats = np.zeros((n,n))
    pv = np.zeros((n,n))
    for i in xrange(n):
        for j in range(i+1,n):
            if test=="wilcox":
                stats[i,j],pv[i,j]=st.wilcoxon(vals[:,i],vals[:,j]) #,zero_method="pratt"
            elif test=="ttest":
                stats[i,j],pv[i,j]=st.ttest_rel(vals[:,i],vals[:,j])
                #raise NotImplementedError(test)
                #stats[i,j],pv[i,j]=st.wilcoxon(vals[:,i],vals[:,j],zero_method=zero_method)
            else:
                raise NotImplementedError(test)
    return stats,pv

def myplot(x,y,styles):
    for i in range(len(styles)):
        pl.plot(x,y[:,i],styles[i])

def plot_ind(ind,val,name,styles):

        #logging.info("A")
        if name=="fraction_Trios":
            val*=3.0
        f=val
        i=ind
        if datamap.has_key(name):
            mapped_name = datamap[name]
        else:
            mapped_name = name
        #logging.info("B")
        pl.figure(figsize=(14,5))
        pl.subplot(1,2,1)
        pl.title("%s=%s"%(mapped_name,str(float(f))),fontsize="large")
        pl.ylabel("type I error",fontsize="large")
        pl.xlabel(r"$\alpha$",fontsize="large")
        myplot(pv_thresholds,best_t1[:,i,:].mean(1).T,styles=styles)
        pl.xscale("log")
        pl.yscale("log")
        legends = [methodmap[methods[ii]] for ii in xrange(len(methods))]
        leg = pl.legend(legends,loc="lower right")
        labels = leg.get_texts()
        for label in labels:
            label.set_fontsize("large")
        
        import scipy.stats as stats
        rt = pv_thresholds[::-1]

        num_trials = 50000 * sum(ind)
        lower = [max(1e-7,(stats.distributions.binom.ppf(0.025, num_trials, t)-1)/float(num_trials)) for t in rt]
        #lower = [stats.distributions.binom.ppf(0.025, num_trials, t)/float(num_trials) for t in rt]
        
        #upper = [max(1e-7,(stats.distributions.binom.ppf(0.975, num_trials, t)/float(num_trials))) for t in rt]
        upper = [stats.distributions.binom.ppf(0.975, num_trials, t)/float(num_trials) for t in rt]

        
        pl.fill_between(rt, lower, upper, alpha=0.7, facecolor='#DDDDDD')

        pl.plot([1e-6,1e-3],[1e-6,1e-3],"k")#"k" is abreviation for black
        pl.xlim([1e-6,1e-3])

        #import pdb
        #pdb.set_trace()
        #logging.info("C")
        pl.subplot(1,2,2)
        pl.title("%s=%s"%(mapped_name,str(float(f))),fontsize="large")
        pl.ylabel("power",fontsize="large")
        y=best_power[:,i,:].mean(1).T
        #myplot(pv_thresholds,y,styles=styles)
        tmp_t1 = best_t1[:,i,:].mean(1).T

        for k in xrange(y.shape[1]):
            pl.plot(tmp_t1[:,k],y[:,k],styles[k])
        #pl.ylim(y.min(),y.max())
        pl.ylim(0.0,1.0)
        pl.xscale("log")
        #pl.yscale("log")
        #pl.legend(methods,loc="lower right")
        pl.xlabel("type I error",fontsize="large")
        pl.xlim([1e-6,1e-3])
        #logging.info("D")
    

        """
        pl.subplot(3,1,3)
        pl.title("for %s=%s"%(mapped_name,str(float(f))),fontsize="large")
        pl.ylabel("power - GC corrected",fontsize="large")
        y=best_power_c[:,i,:].mean(1).T
        #myplot(pv_thresholds,y,styles=styles)
        for k in xrange(y.shape[1]):
            pl.plot(tmp_t1[:,k],y[:,k],styles[k])
        #pl.ylim(y.min(),y.max())
        pl.ylim(.1,1.0)
        pl.xscale("log")
        pl.yscale("log")
        #pl.legend(methods,loc="lower right")
        pl.xlabel(r"$\alpha$",fontsize="large")
        pl.xlim([1e-7,1e-3])
        """
        #pl.subplot(2,2,4)
        #pl.title("global average GC corrected for %s=%.4f"%(name,f))
        #pl.ylabel("Type-1 error")
        #pl.xlabel(r"$\alpha$")
        #myplot(pv_thresholds,best_t1_c[:,i,:].mean(1).T,styles=styles)
        #pl.xscale("log")
        #pl.yscale("log")
        #pl.legend(methods,loc="upper left")
        #pl.plot([1e-8,1e-3],[1e-8,1e-3],"k")
        pl.savefig("%s_%.4f.pdf"%(name,f))


def PlotForOneMethodsList():
    logging.info("methods: {0}".format(",".join(methods)))
    
    logging.info("gather counts and create numPCs.txt file")
    with open("numPCs.txt", "w") as numPCsFile:
        for i_r, r in enumerate(res):

            results_df = pandas.DataFrame(index=methods, columns=cols)
            for i_m, method in enumerate(methods):
                #print "method", method
                #import pdb; pdb.set_trace();

                if method=='insample' or method=='oosample' or method=='linreg' or method=='lmm_ins' or method=='lmm_oos':
                    minval = 9999999999.9
                    for n in npc:
                        fullstr = "%s_%i"%(method,n)
                        if minval>r["res"]["obj"][fullstr]:
                            minval = r["res"]["obj"][fullstr]
                            bestpc[i_m,i_r]=n
                            beststr[i_m,i_r]=fullstr
                            best_k[i_m,i_r]=r["res"]["k"][fullstr]
                elif method[-3::]=="_GC":
                    minval = 9999999999.9
                    for n in npc:
                        fullstr = "%s_%i"%(method[:-3],n)
                        if minval>max(1.0,r["res"]["lambda"][fullstr]):
                            minval = max(1.0,r["res"]["lambda"][fullstr])
                            bestpc[i_m,i_r]=n
                            beststr[i_m,i_r]=fullstr
                            best_k[i_m,i_r]=r["res"]["k"][fullstr]
                elif method[-2::]=="_0" or method[-2::]=="_1" or method[-2::]=="_2":
                    beststr[i_m,i_r]=method
                    best_k[i_m,i_r]=r["res"]["k"][method]
                    bestpc[i_m,i_r]= int(method[-1])          
                elif "_auto" in  method:
                    beststr[i_m,i_r]=method
                    try:
                        best_k[i_m,i_r]=r["res"]["k"][method]
                    except:
                        print "hello"
                    #import pdb; pdb.set_trace()
                    bestpc[i_m,i_r]= r["res"]["numPCs"][method] # wasint(method[-1])          
                elif method[-4::]=='_lin':
                    minval = 9999999999.9
                    for n in npc:
                        fullstr = "%s_%i"%(method[:-4],n)
                        selectstr_pc="linreg_%i"%n
                        if minval>r["res"]["obj"][selectstr_pc]:
                            minval = r["res"]["obj"][selectstr_pc]
                            bestpc[i_m,i_r]=n
                            beststr[i_m,i_r]=fullstr
                            best_k[i_m,i_r]=r["res"]["k"][fullstr]
                else:
                    print "unknown method", method
                    raise NotImplementedError(method)
                #import pdb
                #pdb.set_trace()
                best_t1[i_m,i_r,:]=r["res"][l_t1].T[beststr[i_m,i_r]]    
                best_power[i_m,i_r,:]=r["res"][l_power].T[beststr[i_m,i_r]]    
                best_t1_c[i_m,i_r,:]=r["res"][l_t1_c].T[beststr[i_m,i_r]]    
                best_power_c[i_m,i_r,:]=r["res"][l_power_c].T[beststr[i_m,i_r]]
                numPCsFile.write("{0}\t{1}\t{2}\n".format(i_m, i_r, bestpc[i_m,i_r]))
    if 1:
        logging.info("creating the i_'s")
        i_fst = np.zeros(len(res),dtype = int)
        i_sib = np.zeros(len(res),dtype = int)
        i_h2 = np.zeros(len(res),dtype = int)
        i_varh = np.zeros(len(res),dtype = int)        
        i_ncausal = np.zeros(len(res),dtype = int)
        for i_r, r in enumerate(res):
            i_fst[i_r] = np.nonzero(FSTs==r['options'].fst)[0][0]
            i_sib[i_r] = np.nonzero(fracSibs==r['options'].fracSibs)[0][0]
            i_h2[i_r] =  np.nonzero(h2s==r['options'].h2)[0][0]
            i_varh[i_r] = np.nonzero(var_hidden==r['options'].var_hidden)[0][0]
            i_ncausal[i_r]=np.nonzero(num_causal==r['options'].csnps)[0][0]
    
    if 1: #plot global averages    
        logging.info("creating all.pdf")
        #pl.ion()
    
        pl.figure(figsize=(7,21))
        pl.subplot(3,1,2)
        pl.title("",fontsize="large") #"all data sets"
        pl.ylabel("power",fontsize="large")
        y=best_power.mean(1).T
     
        myplot(pv_thresholds,y,styles=styles)
        #pl.ylim(y.min(),y.max())
        pl.ylim(.1,.4)
        try:
            pl.xscale("log")
            #pl.yscale("log")
        except:
            pass
        #pl.legend(methods,loc="upper left")
        #pl.xlabel(r"$\alpha$")
        pl.xlabel(r"$\alpha$",fontsize="large")
        pl.xlim([1e-7,1e-3])
        pl.subplot(3,1,1)
        pl.title("",fontsize="large")#"all data sets"
        pl.ylabel("type I error",fontsize="large")
        pl.xlabel(r"$\alpha$",fontsize="large")
        myplot(pv_thresholds,best_t1.mean(1).T,styles=styles)
        pl.fill_between(xvalbar,barlow,barhigh,color = "grey",alpha=0.5)
        try:
            pl.xscale("log")
            pl.yscale("log")
        except:
            pass
        pl.xlim([1e-7,1e-3])
        legends = [methodmap[methods[ii]] for ii in xrange(len(methods))]
        leg = pl.legend(legends,loc="lower right")
        labels = leg.get_texts()
        for label in labels:
            label.set_fontsize("large")
        pl.plot([1e-8,1e-3],[1e-8,1e-3],"k")
    
        #pl.subplot(1,2,1)
        #pl.title("global average GC corrected")
        #pl.ylabel("Type-1 error")
        #pl.xlabel(r"$\alpha$")
        #myplot(pv_thresholds,best_t1_c.mean(1).T,styles=styles)
        #try:
        #    pl.xscale("log")
        #    pl.yscale("log")
        #except:
        #    pass
        ##pl.legend(methods,loc="upper left")
        #pl.plot([1e-7,1e-3],[1e-7,1e-3],"k")
    
        pl.subplot(3,1,3)
        pl.title("",fontsize="large")#"all data sets"
        pl.ylabel("power - GC corrected",fontsize="large")
        y=best_power_c.mean(1).T
        myplot(pv_thresholds,y,styles=styles)
        #pl.ylim(y.min(),y.max())
        pl.ylim(.1,.4)
        pl.xlim([1e-7,1e-3])
        try:
            pl.xscale("log")
            #pl.yscale("log")
        except:
            pass
        #pl.legend(methods,loc="upper left")
        pl.xlabel(r"$\alpha$",fontsize="large")
    
        pl.savefig("all.pdf")
    

 
    if 1:#new all plot
        thr=pv_thresholds[4]
        
        counter = 0
        logging.info("creating grid plots")
        k = np.ones(len(i_fst), dtype=np.bool)

        try:
            tmp_name = "all_new"
                    
            plot_ind(ind=k,val=0,name=tmp_name,styles=styles)
        except Exception, detail:
            print "failed to generate plot", tmp_name
            print detail
            pass
 
    if 1:#new all plot
        thr=pv_thresholds[4]
        
        counter = 0
        logging.info("creating grid plots")
        k = np.ones(len(i_fst), dtype=np.bool)
        #k = np.bitwise_and(k, i_varh==0)#############################################!!!!!!!!!!!!!!!!!!!
        k = np.bitwise_and(k, i_varh==0)
        k = np.bitwise_and(k, i_sib==0)
        #k = np.bitwise_and(k, i_ncausal==2)
        #import pdb; pdb.set_trace()
        num_set = sum(k)
        try:
            tmp_name = "no_fs_var_hidden=0_num%i" % num_set
            
            plot_ind(ind=k,val=0,name=tmp_name,styles=styles)
        except Exception, detail:
            print "failed to generate plot", tmp_name
            print detail
            pass

    if 1:#single threshold only
        thr=pv_thresholds[4]
        
        counter = 0
        logging.info("creating grid plots")
        for i_f,f in enumerate(fracSibs):
            for j_f,f2 in enumerate(num_causal):
                
                #ii=i_fst==i_f
                ii=i_sib==i_f
                jj=i_ncausal==j_f
                k = np.bitwise_and(ii,jj)
                k = np.bitwise_and(k, i_varh==0)#############################################!!!!!!!!!!!!!!!!!!!
                #k = np.bitwise_and(k, i_fst==0)
                #k = np.bitwise_and(k, i_h2==3)
                #import pdb; pdb.set_trace()
                num_set = sum(k)
                try:
                    tmp_name = "grid_fracSibs=%f_numCausal=%i_numSet=%i" % (f,f2,num_set)
                    
                    plot_ind(ind=k,val=f2*f,name=tmp_name,styles=styles)
                except Exception, detail:
                    print "failed to generate plot", tmp_name
                    print detail
                    pass
    
        #population structure (Fst)
        logging.info("creating FSTs plots")
        for i_f,f in enumerate(FSTs):
            i=i_fst==i_f
            try:
                plot_ind(ind=i,val=f,name="Fst",styles=styles)
            except:
                pass

        #family structure (sibs)
        logging.info("creating Trios plots")
        for i_f,f in enumerate(fracSibs):
            i=i_sib==i_f
            try:
                plot_ind(ind=i,val=f,name="fraction_Trios",styles=styles)
            except:
                pass
        #heritability
        logging.info("creating h2 plots")
        for i_f,f in enumerate(h2s):
            i=i_h2==i_f
            try:
                plot_ind(ind=i,val=f,name="h2",styles=styles)
            except:
                pass
        #var_hidden 
        logging.info("creating confounding var plots")
        for i_f,f in enumerate(var_hidden):
            i=i_varh==i_f
            try:
                plot_ind(ind=i,val=f,name="confounding_var",styles=styles)
            except:
                pass
        #number causal
        logging.info("creating number causal plots")
        for i_f,f in enumerate(num_causal):
            i=i_ncausal==i_f
            i = np.bitwise_and(i, i_varh==0)
            i = np.bitwise_and(i, i_sib>0)
            try:
                plot_ind(ind=i,val=f,name="number_causal",styles=styles)
            except:
                pass
    if 1:
        best_k_mean=sp.zeros((len(methods),len(num_causal)))
        frac_k_all=sp.zeros((len(methods),len(num_causal)))
        logging.info("creating k_ plots")
        for i_f,f in enumerate(num_causal):
            i=i_ncausal==i_f
            for i_m, method in enumerate(methods):
                best_k_mean[i_m,i_f]=best_k[i_m,i].mean()
                frac_k_all[i_m,i_f]=(best_k[i_m,i]>20000).mean()
                if 0:#method=='insample' or method=='oosample':
                    pl.figure()
                    pl.hist(best_k[i_m,i][best_k[i_m,i]<20000],50)
                    pl.title("all SNPs selected in %i out of %i cases"%((best_k[i_m,i]>=20000).sum(),best_k[i_m,i].shape[0]))
                    name = 'k_nC%i_%s.pdf'%(f,method)
                    pl.savefig(name)
    
    pl.close("all")
    
    if 1:
        pl.figure(figsize=(7,7));
        for iplot in [0]:#will plot method[0], so set to what you want
            pl.plot(num_causal,frac_k_all[iplot],styles[iplot])
        pl.xlabel("Number of causal SNPs")
        pl.ylabel("Fraction of data sets where all SNPs are selected")
        pl.savefig("frac_all_SNPs.pdf")
    return best_k, best_power, best_power_c, best_t1, best_t1_c, bestpc, beststr, frac_k_all, i_fst, i_h2, i_ncausal, i_sib, i_varh, ii, npc, styles

def create_stats():
    out = np.zeros((6,10),dtype='|S30')
    out[0,:] = ["key","diff of means","ttest 2-sided","ttest 1-sided","ranksum","N_reduced","wilcox 2-sided","z_score_wiki","wilcox 2-sided wiki","wilcox 1-sided"]
    for ii,thres in enumerate(["1.0e-03","1.0e-04","1.0e-05","1.0e-06","1.0e-07"]):
        i = ii+1
        key="Type I @ %s"%thres
        gk = getKey(key=key,method = None, mstring=None,index=None)
        ttest = computeTest(key=key,test="ttest",index=None)
        wilcox = computeTest(key=key,test="wilcox",index=None)
        out[i,0] = key
        out[i,1] = str(gk[:,0].mean()-gk[:,1].mean())
        out[i,2] = str(ttest[1][0,1])
    
        if float(out[i,1]) < 0:
            out[i,3] = str(1-.5*float(out[i,2]))
        else:
            out[i,3] = str(.5*float(out[i,2]))
    
        diff = gk[:,0]-gk[:,1]
        diff = diff[diff != 0] #remove any zero's
        N_r=diff.shape[0]
        sign = (diff < 0) * 2 - 1 #make -1,1
        absdiff = np.abs(diff)
        rank = st.mstats.rankdata(absdiff)
        ranksum = -(sign * rank).sum()
        out[i,4] = str(ranksum)
        out[i,5] = str(N_r)
        out[i,6] = str(wilcox[1][0,1])
        sigma_w=np.sqrt(N_r*(N_r+1)*(2.0*N_r+1)/6.0)
        z_wiki = (np.abs(ranksum)-0.7)/sigma_w
        wilcox_wiki = st.norm.sf(z_wiki)
        out[i,7]=str(z_wiki)
        out[i,8]= str(wilcox_wiki)
        if ranksum < 0:
            out[i,9] = str(1-.5*float(out[i,8]))
        else:
            out[i,9] = str(.5*float(out[i,8]))
    np.savetxt("stats.txt",out,fmt="%s",delimiter="\t")
    return ii, out, thres

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    if 1:
        res=[]
        if not do_merge:
            for fn in filenames_in:
                logging.info("reading result pickle file: {0}".format(fn))
                file = open(fn,"r")
                res+=(cPickle.load(file))
                file.close()
        else:
            res=None
            for fn in filenames_in:
                logging.info("reading result pickle file: {0}".format(fn))
                with open(fn) as file:
                    items = (cPickle.load(file))
                if res is None:
                    res = items
                else:
                    if len(res) != len(items) : raise Exception("Expect all files to have the same number of items")
                    for i_r, r in enumerate(res):
                        item = items[i_r]
                        #import pdb; pdb.set_trace()
                        if not r['options'].__dict__ == item['options'].__dict__: raise Exception("The options differ for result #{0}".format(i_r))
                        res[i_r]['res'] = r['res'].append(item['res'])


        #import pdb;pdb.set_trace()
        logging.info("setting run parameters")
        if 0:
            FSTs = np.array([0.01, 0.05, 0.1])
            fracSibs=np.array([0.0, 0.15,0.3])
            h2s = np.arange(0.1, 0.5, 0.1)
            var_hidden = np.arange(0.0, 0.9, 0.3)
            num_causal = np.array([10, 100, 500, 1000])
        elif 0:
            FSTs = np.array([0.01, 0.05, 0.1])
            fracSibs=np.array([0.0, 0.15,0.3])
            h2s = np.arange(0.1, 0.7, 0.1)
            var_hidden = np.arange(0.0, 1.0, 0.3)
            num_causal = np.array([10, 100, 500, 1000])
        elif 0:
            FSTs = np.array([0.0])
            fracSibs=np.array([0.0,0.05,0.1,0.2])
            h2s = np.arange(0.1, 0.7, 0.1) #! this is a range!  6 entries
            #var_hidden = np.array([0.0]) #
            var_hidden = np.arange(0.0, 1.0, 0.3) # ! this is a range! 4 entries
            num_causal = np.array([10, 50, 100, 500, 1000])
        elif 0:
            FSTs = np.array([0.0])
            fracSibs=np.array([0.0,0.05,0.1,0.2,0.4])
            h2s = np.arange(0.1, 0.7, 0.1) #! this is a range!  6 entries
            var_hidden = np.array([0.0]) #
            #var_hidden = np.arange(0.0, 1.0, 0.3) # ! this is a range! 4 entries
            num_causal = np.array([10, 50, 100, 500, 1000])
        elif 1:
            FSTs = np.array([0.0])
            fracSibs=np.array([0.0, 0.5, 0.6, 0.7, 0.8, 0.9])
            h2s = np.arange(0.1, 0.7, 0.1) #! this is a range!  6 entries
            var_hidden = np.array([0.0, 0.3]) #
            #var_hidden = np.arange(0.0, 1.0, 0.3) # ! this is a range! 4 entries
            num_causal = np.array([10, 50, 100, 500, 1000])
        else:
            #FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            #fracSibs=np.array([0.2])
            #h2s = np.array([.1])
            #var_hidden = np.array([0])
            #num_causal = np.array([200])
            FSTs = np.array([0.005, 0.01, 0.05, 0.1])
            fracSibs=np.array([0.0,0.05,0.1,0.2])
            h2s = np.arange(0.1, 0.7, 0.1) #! this is a range!  6 entries
            var_hidden = np.arange(0.0, 1.0, 0.3) # ! this is a range! 4 entries
            num_causal = np.array([10, 50, 100, 500, 1000])


        pv_thresholds = [1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8]
        sim_params = ['h2', 'h2 hidden', '% trios', '#causal', 'Fst']
        l_power  = ['Power @ %.1e' % i for i in pv_thresholds]
        l_t1     = ['Type I @ %.1e' % i for i in pv_thresholds] 
        l_power_c= ['Power corr @ %.1e' % i for i in pv_thresholds]
        l_t1_c   = ['Type I corr @ %.1e' % i for i in pv_thresholds]
        cols = sum([sim_params, l_power, l_t1, ['lambda', 'k', 'obj'],l_power_c,l_t1_c], [])
#        methods = ["linreg_0","lmm_0","linreg_1","lmm_1","linreg_2","lmm_2","insample","oosample"]
#        methods = ["linreg_1","lmm_ins_0","lmm_ins_1","lmm_oos_0","lmm_oos_1","insample_0","oosample_0","insample_1","oosample_1","linreg","lmm_ins","lmm_oos","insample","oosample"]
#        methods = ["linreg","linreg_1","lmm_ins_0","lmm_ins","lmm_oos","insample","oosample"]
#        methods = ["linreg","linreg_1","lmm_ins_0","lmm_ins","lmm_oos","insample","oosample","insample_1","oosample_1"]
#        methods = ["linreg","linreg_1","lmm_ins_0","lmm_ins","lmm_oos","insample","oosample","insample_1","oosample_1","insample_GC","oosample_GC"]
#        methods = ["linreg","insample","oosample","insample_GC","oosample_GC","insample_lin","oosample_lin","insample_2","oosample_2"]
#        if "2meth" in sys.argv:
#            methods = ["insample","oosample"]
#        else:
#            methods = ["lmm_ins_0","lmm_oos","insample","oosample","insample_0","oosample_0","insample_1","oosample_1"]
#        #methods = ["linreg_1","lmm_ins_0","lmm_ins","insample","oosample","insample_1","oosample_1","linreg"]
#        methods = ["oosample_2","lmm_ins"]
##import pdb
#        #methods = ["oosample_auto"]k
#        #pdb.set_trace()
#        #methods = ["lmm_0","lmm_1","insample_0","oosample_0","insample_1","oosample_1","insample","oosample"]
#        methods = ["oosample_auto","insample_auto", "lmm_ins_auto", "linreg_auto"] # "lmm_oss_auto", 
#        methods = ["oosample_auto","insample_auto", "lmm_ins_auto", "linreg_auto","lmm_oos_auto"] # "lmm_oss_auto", 

#        #methods = ["lmm_oos", "lmm_oos_0", "lmm_oos_auto"]
#        #methods = ["oosample", "oosample_0", "oosample_auto"]
#        #methods = ["linreg", "linreg_0", "linreg_auto"]
#        #methods = ["oosample_auto","lmm_ins_auto", "linreg_auto"]

        stylesall=[':x', '-v', '--o', ':+','-.*',':>', '-<', '--d', ':x','-.v',':o', '-+', '--*', ':d','-.x']

        for folder_name, methods in [
                                                #["stage1", ["linreg_auto", "lmm_ins_auto", "insample_auto"]],
                                                #["stage2_linreg", ["linreg_0", "linreg", "linreg_auto"]],
                                                #["stage2_lmm", ["lmm_ins_0", "lmm_ins", "lmm_ins_auto"]],
                                                #["stage2_select", ["insample_0", "insample", "insample_auto"]]
                                                #["stage3", ["linreg_auto", "lmm_ins_auto", "insample_auto", "lmm_fs_cond_full_auto"]]
                                                #["stage4", ["linreg_auto", "lmm_ins_auto", "insample_auto", "lmm_fs_cond_full_auto"]]
                                                #["auto", ["linreg_auto", "lmm_ins_auto", "insample_auto", "lmm_fs_cond_full_auto"]]
                                                ["auto", ["linreg_auto", "lmm_ins_auto", "insample_auto", "lmm_fs_cond_full_auto", "lmm_true_causals_auto"]]
                                                #["auto", ["lmm_fs_cond_full_auto_2", "lmm_fs_cond_full_auto_8", "lmm_fs_cond_full_auto", "linreg_auto", "lmm_ins_auto", "insample_auto"]],

                                                #["01172014d", ["insample_auto", "oosample_auto"]],
                                                #["01172014e", ["oosample_0", "oosample_auto"]],
                                                ]:
            os.mkdir(folder_name)
            os.chdir(folder_name)

            styles = stylesall[:len(methods)]
            npc = [0,1,2]
            bestpc = np.zeros((len(methods),len(res)),dtype='int')
            beststr = np.zeros((len(methods),len(res)),dtype='|S30')
            best_power = np.zeros((len(methods),len(res),len(l_power)))
            best_t1 = np.zeros((len(methods),len(res),len(l_t1)))
            best_power_c = np.zeros((len(methods),len(res),len(l_power)))
            best_t1_c = np.zeros((len(methods),len(res),len(l_t1)))
            best_k = np.zeros((len(methods),len(res)),dtype='int')

            PlotForOneMethodsList()

            create_stats()

            os.chdir("..")

