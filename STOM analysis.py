# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:39:43 2022

@author: Freddie
"""

import pandas as pd
import numpy as np 
from scipy.optimize import curve_fit
import stomtools as stm
import os
import matplotlib.pyplot as plt
import re
from scipy.integrate import quad
#Define a gaussian for fitting the chi PDFs
def gauss(x,m,s):
    factor=1/(s*np.sqrt(2*np.pi))
    exp=np.exp((-(x-m)**2)/(2*(s**2)))
    return factor*exp

def gauss_fact(x,m,s,A):
    factor=1/(s*np.sqrt(2*np.pi))
    exp=np.exp((-(x-m)**2)/(2*(s**2)))
    return A*factor*exp
def quadratic(x,a,b,c,):
    x=pd.Series(x)
    return a*x**2+b*x+c
#Iterate over all filepaths and assigns the directories to a list
filepaths=[]
for root, dirs, files in os.walk("C:\\Users\\Freddie\\OneDrive\\Desktop\\Results",topdown=True):    
    for file in files:                              
        filepaths.append(os.path.join(root,file)) 
print(filepaths)
#%%
#Plot a gaussian best fit for each chi square distribution
g_means=[]
g_stds=[]
signals=[]
p_values=[]
for dr in filepaths:
    temp=pd.read_csv(dr)
    chi_sq=temp.iloc[:,1]
    bin_heights,bin_edges=np.histogram(chi_sq,bins=50,density=True)
    x=bin_edges[0:-1]+(bin_edges[1]-bin_edges[0])/2
    plt.scatter(x,bin_heights)
    fit,cov=curve_fit(gauss,x,bin_heights,p0=[45,10])
    opt_fit=gauss(x,fit[0],fit[1])
    p_val=quad(gauss,fit[0],500,args=(45.69966252997833,9.972451533314945))
    print(fit[0],fit[1])
    p_values.append(p_val[0])
    plt.plot(x,opt_fit)
    plt.show()
    g_means.append(fit[0])
    g_stds.append(fit[1])
    sig_temp=re.findall(r'\d+',dr)
    for val in sig_temp:
        signals.append(int(val))
    print(sig_temp)
    print(p_values)
#%%
df=pd.DataFrame([signals,g_means,g_stds,p_values])
df=df.transpose()
df["Signal"]=df.iloc[:,0]
df=df.sort_values(by = "Signal", ascending=True)
df_sigs=df["Signal"]
df_means=df.iloc[:,1]
df_stds=df.iloc[:,2]
df_pvals=df.iloc[:,3]
print(df)
plt.scatter(signals,g_means)
fit1,cov=curve_fit(quadratic,df_sigs,df_means,p0=[0.1,1,0])
dummy_sig=np.linspace(0,np.max(df_sigs),int(max(df_sigs)))
print(dummy_sig)
opt_fit=quadratic(dummy_sig,fit1[0],fit1[1],fit1[2])
plt.plot(dummy_sig,opt_fit)
exact_pvals=[]
for sig in dummy_sig:
    l_b=quadratic(sig,fit1[0],fit1[1],fit1[2])
    exact_p=quad(gauss,float(l_b),500,args=(45.69966252997833,9.972451533314945))
    exact_pvals.append(exact_p)
    if 0.05>=exact_p[0] and 0.048<=exact_p[0]:
        print(sig)
print(sig)
plt.xlabel("Signal Amplitude")
plt.ylabel("Chi_Squared")
plt.grid()
plt.show()
plt.scatter(df_sigs,df_pvals)
plt.plot(dummy_sig,exact_pvals)
plt.xlabel("Signal Amplitude")
plt.ylabel("P-value")
plt.grid()

