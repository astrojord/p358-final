#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 18:40:14 2018

@author: cameronkass
"""

import numpy as np
import matplotlib.pyplot as plt
import make_halos
import methast

rs       = 17.  #kpc



#overlay dens_NFW plot on dens_NFW hist
def make_hist_NFW():
    delta = [0.05,0.01,0.005,.001]
    for i in range(len(delta)):
        r_NFW         = methast.methast(make_halos.NFW_dist,rs,1000000,1.,delta[i])
        hist_NFW,edges= np.histogram(r_NFW,1000,normed=False)
        x             = 0.5*(edges[0:edges.size-1]+edges[1:edges.size])
        #tothist       = np.sum(hist_NFW.astype(float))
        #hist_NFW      = hist_NFW/tothist
        shellvols     = np.zeros(x.size)
        shellvols[0]  = (4./3.)*np.pi*edges[1]**3
        shellvols[1:] = ((4./3.)*np.pi*edges[2:]**3)-((4./3.)*np.pi*edges[1:edges.size-1]**3)
        hist_NFW      = np.divide(hist_NFW,shellvols)
        plt.figure(i)
        
        plt.plot(x,make_halos.NFW_dist(x,rs),linestyle='-',color='red',linewidth=1.,label='NFW_dist')
        plt.plot(x,hist_NFW,'.',label='methast')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel('particle density')
        plt.legend()
        plt.savefig('NFW_hist_%s.png' % str(delta[i]))
        plt.show()

#overlay dens_cusp plot on dens_cusp hist


#new delta list is advisable for this profile
def make_hist_cusp():
    delta = [0.0025,0.001,.00075,.0005]
    for i in range(len(delta)):
        r_cusp = methast.methast(make_halos.cusp_dist,rs,1000000,1.,delta[i])
        hist_cusp,edges= np.histogram(r_cusp,1000,normed=False)
        x             = 0.5*(edges[0:edges.size-1]+edges[1:edges.size])
        tothist       = np.sum(hist_cusp.astype(float))
        hist_cusp     = hist_cusp/tothist
        shellvols     = np.zeros(x.size)
        shellvols[0]  = (4./3.)*np.pi*edges[1]**3
        shellvols[1:] = ((4./3.)*np.pi*edges[2:]**3)-((4./3.)*np.pi*edges[1:edges.size-1]**3)
        hist_cusp     = np.divide(hist_cusp,shellvols)
        plt.figure(i)
        
        plt.plot(x,make_halos.cusp_dist(x,rs),linestyle='-',color='red',linewidth=1.,label='cusp_dist')
        plt.plot(x,hist_cusp,'.',label='methast')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('r')
        plt.ylabel('normalized density')
        plt.legend()
        plt.savefig('cusp_hist_%s.png' % str(delta[i]))
        plt.show()

        
def main():
    make_hist_NFW()
    make_hist_cusp()
    
main()