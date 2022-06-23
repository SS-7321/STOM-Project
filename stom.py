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
import matplotlib.image as im


def chisq(x,y,A,L,s):
    diff=((y-stm.get_B_expectation(x, A, L))/s)**2
    chi_sum=np.sum(diff)
    return chi_sum

def chisqSB(x,y,s,A,L, mu, sig, amp):
    diff=((y-stm.get_SB_expectation(x, A, L, mu, sig, amp))/s)**2
    chi_sum=np.sum(diff)
    return chi_sum

def signal_gaus(x, mu, sig, signal_amp):
    return signal_amp/(np.sqrt(2.*np.pi)*sig)*np.exp(-np.power((x - mu)/sig, 2.)/2)

def get_SB_expectation(xs, A, lamb, mu, sig, signal_amp):
    ys = []
    for x in xs:
        ys.append(A*np.exp(-x/lamb) + signal_gaus(x, mu, sig, signal_amp))
    return ys

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
#%%]
data = stm.generate_data()


bin_heights, bin_edges = np.histogram(data, range = [104,155], bins = nbins)

fits, cov = curve_fit(get_SB_expectation, x, bin_heights, [minA,minLamb,125,1.5,700], maxfev=50000)
plt.plot(x, get_SB_expectation(x, fits[0], fits[1], fits[2], fits[3], fits[4]))
plt.show()
plt.scatter(x,bin_heights,label="Data")
plt.plot(x,stm.get_B_expectation(x, fits[0], fits[1]),label="Best Fit for Background",color="red")
plt.legend()
plt.title("Background Fit against Data")
plt.xlabel("Energy (GeV)")
plt.ylabel("Frequency of Occurences")
plt.grid()
plt.savefig("C:\\Users\\samiu\\Documents\\Edu\\Codes\\Spyder\\bgVdata.png", dpi=1200)
plt.show()
print(chisq(x,bin_heights,fits[0], fits[1],h_err))
print(1-chi2.cdf(chi_min,Ndof))
#%%
chi=chisqSB(x,bin_heights,h_err, fits[0], fits[1], fits[2], fits[3], fits[4])
print(chi/Ndof)
print(1-chi2.cdf(chi,Ndof))


pure_bg=bin_heights-stm.signal_gaus(x,fits[2],fits[3],fits[4])
bh, be = np.histogram(pure_bg, range = [104,155], bins = nbins)
fits2, cov = curve_fit(get_SB_expectation, (be[0:-1]+(be[1]-be[0])/2), pure_bg, [minA,minLamb,125,1.5,700], maxfev=50000)
xnew = be[0:-1]+(be[1]-be[0])/2
opt_fit=get_SB_expectation(xnew,fits2[0], fits2[1], fits[2], fits[3], fits[4])
# plt.plot(xnew,opt_fit,label="Data")
# plt.scatter(xnew,pure_bg,label="Estimated Background")
plt.show()
chi=chisqSB(xnew,pure_bg,np.sqrt(pure_bg), fits2[0], fits2[1], fits[2], fits[3], fits[4])
print(chi/Ndof)
print(1-chi2.cdf(chi,Ndof))
std = 5.93


plt.scatter(x,bin_heights)
plt.plot(xnew,get_SB_expectation(xnew,fits[0],fits[1],fits[2],fits[3],fits[4]))
plt.show()
opt_SB=get_SB_expectation(xnew,fits[0],fits[1],fits[2],fits[3],fits[4])
pure_sig=bin_heights-pure_bg
plt.scatter(xnew,pure_sig,label="Data",marker='o')
fit3, cov3 = curve_fit(stm.signal_gaus, xnew, pure_sig, [125, 1.5, 160],sigma=1/pd.Series(np.sqrt(pure_bg)),absolute_sigma=True)
plt.plot(xnew, stm.signal_gaus(xnew, fit3[0], fit3[1], fit3[2]), label="Best Fit")
plt.xlabel("Energy (GeV)")
plt.ylabel("Frequency of Occurences")
plt.legend()
plt.grid()
plt.title("Isolated Signal")
plt.savefig("C:\\Users\\samiu\\Documents\\Edu\\Codes\\Spyder\\isoSig.png", dpi=1200)

plt.show()
print(fit3)
print(np.sqrt(cov3[0][0]),np.sqrt(cov3[1][1]),np.sqrt(cov3[2][2]))
print(1-chi2.cdf(all_chi_min,Ndof))
print(1-chi2.cdf(chi_min,Ndof))
#%%
plt.scatter(x,bin_heights)
min_chi=2000
m_energy=0
for m in np.linspace(0,160,1000):
    opt_SB=get_SB_expectation(x,fits[0],fits[1],m,fits[3],fits[4])
    plt.plot(x,get_SB_expectation(x,fits[0],fits[1],m,fits[3],fits[4]))
    get_chi=chisqSB(x,bin_heights,h_err,fits[0],fits[1],m,fits[3],fits[4])
    if get_chi<min_chi:
        min_chi=get_chi
        m_energy=m
        
    
        