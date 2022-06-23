# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 17:01:19 2022

@author: samiu
"""

# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""

import numpy as np
import matplotlib.pyplot as plt
import stomtools as stm
import scipy as sp
import pandas as pd
from scipy.optimize import curve_fit
from tqdm import tqdm

from scipy.stats import chi2

def exp_func(x, A, lamb):
    ''' 
    Return a set of expectation values for the background distribution for the 
    passed in x values. 
    '''
    return A*np.exp(-x/lamb)

def chisq(x,y,A,L,s):
    diff=((y-stm.get_B_expectation(x, A, L))/s)**2
    chi_sum=np.sum(diff)
    return chi_sum
   #Divide by Ndof for reduced chi squared
    
nbins = 45
ndof = nbins-2
    
data = stm.generate_data()


#%%
bin_heights, bin_edges = np.histogram(data, bins = nbins, range = [104,155])
x = bin_edges[0:-1]+(((155-104)/nbins)/2)
# above was a test to find the interval of A and Lambda to iterate over

#%% sampling the data before 120 GeV to find the A and Lambda
x_samp=np.array([])
y_samp=np.array([])
for val in range(0,len(x)):
    if x[val]<120:
        x_samp=np.append(x_samp,x[val])
        y_samp=np.append(y_samp,bin_heights[val])
    else:
        continue
#%% iterating over a range of lambda and A values to find the optimized background fit

h_err = np.sqrt(bin_heights)
iterA = np.linspace(33300,34000,1000)
iterLamb = np.linspace(29,33,40)
chi=chisq(x_samp,y_samp, iterA[0], iterLamb[0],h_err[0:len(x_samp)])
for i in range(0,len(iterA)):
    for j in range(0,len(iterLamb)):
        test_chi = chisq(x_samp, y_samp, iterA[i], iterLamb[j],h_err[0:len(x_samp)])
        if test_chi<chi:
            chi=chisq(x_samp,y_samp, iterA[i], iterLamb[j],h_err[0:len(x_samp)])
            minA=iterA[i]
            minLamb=iterLamb[j]
       
chi_min=chi
print(chi_min)
print(minA)
print(minLamb)
#%% finding chi min with signal
iterA2 = np.linspace(41200,41300,100)
iterLamb2 = np.linspace(29,31,100)

chi=chisq(x,bin_heights, iterA2[0], iterLamb2[j],h_err)
for i in range(0,len(iterA2)):
    for j in range(0,len(iterLamb2)):
        test_chi = chisq(x, bin_heights, iterA2[i], iterLamb2[j],h_err)
        if test_chi<chi:
            chi=chisq(x,bin_heights, iterA2[i], iterLamb2[j],h_err)
            minA2=iterA2[i]
            minLamb2=iterLamb2[j]

all_chi_min=chi
print(all_chi_min)
print(minA2)
print(minLamb2)

#%%

for k in tqdm(range (0,100)):
    results = pd.DataFrame()
    chi_dist = []
    prob_dist = []
    average_probability=[]
    for l in tqdm(range(0,100)):     
        data_chi_t=stm.generate_data(k*4)
        bin_heights2 = np.histogram(data_chi_t, bins = nbins, range = [104,155])[0]
        
        h_err2 = np.sqrt(bin_heights2)
        iterA3 = np.linspace(3000,6000,1000)
        iterLamb3 = np.linspace(29,33,100)
        chi=chisq(x,bin_heights2, iterA3[0], iterLamb3[0],h_err2)
        for i in range(0,len(iterA3)):
            for j in range(0,len(iterLamb3)):
                test_chi = chisq(x, bin_heights2, iterA3[i], iterLamb3[j],h_err2)
                if test_chi<chi:
                    chi=chisq(x,bin_heights2, iterA3[i], iterLamb3[j],h_err2)
                    minA3=iterA3[i]
                    minLamb3=iterLamb3[j]
                    
        bg_chi_min = chi
        # prob = 1-chi2.cdf(bg_chi_min,ndof)
        chi_dist.append(bg_chi_min)
        # prob_dist.append(prob)
    # avg_p = np.mean(np.array(prob_dist))
    # average_probability.append(avg_p)
    results['Chi min'] = chi_dist
    # results['Probability'] = prob_dist
    # results['Average Probability'] = average_probability
    results.to_csv("C:\\Users\\samiu\\Desktop\\Chi_pdf_Amplitude{}.csv".format(k*4), index=True)
    

#%%

for k in tqdm(range (51,89)):
    results = pd.DataFrame()
    chi_dist = []
    prob_dist = []
    average_probability=[]
    for l in tqdm(range(0,10000)):     
        data_chi_t=stm.generate_data(k*4)
        bin_heights2 = np.histogram(data_chi_t, bins = nbins, range = [104,155])[0]
        
        h_err2 = np.sqrt(bin_heights2)
        
        fits, cov = curve_fit(exp_func, x, bin_heights2, [4000,31], maxfev=50000)
        #plt.hist(data_chi_t, bins = nbins, range = [104,155])
        #opt_fit=exp_func(x,fits[0],fits[1])
        #plt.plot(x,opt_fit)
        #plt.show()
        chi=chisq(x,bin_heights2, fits[0], fits[1],h_err2)
        # prob = 1-chi2.cdf(chi,ndof)
        chi_dist.append(chi)
        # prob_dist.append(prob)
    # avg_p = np.mean(np.array(prob_dist))
    # average_probability.append(avg_p)
    results['Chi min'] = chi_dist
    # results['Probability'] = prob_dist
    results.to_csv("C:\\Users\\samiu\\Desktop\\Results\\Chi_pdf_Amplitude{}.csv".format(k*4), index=True)
    
#%%

plt.hist(chi_dist, bins=20)
plt.show()
#%%
print(1-chi2.cdf(chi_min,ndof))
print(1-chi2.cdf(all_chi_min,ndof))
#%%
# 1.235618750418813
# 0.1384031934777934
# 2.3107538834522416
# 2.3591800105693395e-06
# [124.74057212   1.44872224 405.99120767]
# 0.0004994792380509996 0.0004998468087683226 0.12124901519410582
# 0.0008680087789532109
# 0.9999995164905917

