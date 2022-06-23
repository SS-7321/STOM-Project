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
from scipy.stats import chi2
from scipy.signal import find_peaks


def chisq(x,y,A,L,s):
    diff=((y-stm.get_B_expectation(x, A, L))/s)**2
    chi_sum=np.sum(diff)
    return chi_sum


data = stm.generate_data()
nbins = 45
Ndof = nbins-2

#%%
bin_heights, bin_edges, patches = plt.hist(data, range = [104,155], bins = nbins)
plt.show()
# #%%
x = bin_edges[0:-1]+(((155-104)/nbins)/2)

# plt.scatter(x[:-1], bin_heights, marker='x')
# plt.show()


test_data = stm.get_B_expectation(x, 37000, 30) # exp dist with arbitrary vals for A and Lamb
plt.plot(x,test_data)
plt.show()
# above was a test to find the interval of A and Lambda to iterate over

#%% finding freq uncertainty for each bin

w_err=np.array([])
for i in range(0,len(bin_edges)-1):
    temp=np.array([])
    for val in data:
        if (val<=bin_edges[i+1] and val>=bin_edges[i]):
            temp=np.append(temp,val)
            
    unc=(max(temp)-min(temp))/2
    w_err=np.append(w_err,unc)  

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
iterA = np.linspace(33000,35000,2000)
iterLamb = np.linspace(29,33,1000)
minA=0
minLamb=0
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


#%% plot background exp graph from optimised parameters
bin_heights, bin_edges, patches = plt.hist(data, range = [104,155], bins = nbins)
test_data = stm.get_B_expectation(x, minA, minLamb) # exp dist with arbitrary vals for A and Lamb
plt.plot(x,test_data)
plt.show()


#%% finding chi min with signal
iterA2 = np.linspace(40000,45000,2000)
iterLamb2 = np.linspace(29,31,1000)

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

#%%
print(1-chi2.cdf(all_chi_min,Ndof))
print(1-chi2.cdf(chi_min,Ndof))

