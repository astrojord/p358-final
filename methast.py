#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 19:06:12 2018

@author: cameronkass
"""

import numpy as np

#---------------------------------------------------------------
def rnorm(R):
    # returns normal random deviates following the Box-Mueller method.
    x = []
    r = 0
    while r < R:
        #obtian two random variables on [0,1]
        U = np.random.random_sample((2,))
        #compute the array of desired deviates
        xr = ((-2.*np.log(U))**0.5)*np.array([np.cos(2.*np.pi*U[1]),np.sin(2.*np.pi*U[1])]) 
        x.append(xr)
        r+=1
    x = np.array(x)
    return x

#--------------------------------------------------------------
def methast(fTAR,scl,R,x0,delta):
    # Returns an array of random variables sampled according to a
    # target distribution fTAR.
    # input: fTAR      : function pointer to the target distribution.
    #                    The function must take arguments fTAR(x,bounds),
    #                    and must return the value of fTAR at x.
    #        scl       : scaling factor for fTAR
    #        R         : number of samples to be generated
    #        x0        : initial guess
    #        delta     : step size, or scale (width of local proposal distribution Q)
    # output: xr      : random variables drawn from fTAR.
    xr = []
    A = rnorm(R) #array of normal deviates
    for r in range(R):
        P0 = fTAR(x0,scl) #evaluate the initial guess
        xt = x0+(delta*A[r,np.random.choice([0,1])]) #pick a trial state
        Pt = fTAR(xt,scl) #evaluate the trial state
        if Pt/P0 >= 1.:
            xr.append(xt)
            x0 = xt
        else:
            u = np.random.random()
            if u <= Pt/P0:
                xr.append(xt)
                x0 = xt
            else:
                xr.append(x0)
    xr = np.array(xr)
    return xr
