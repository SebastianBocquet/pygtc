from matplotlib import pyplot as plt
import numpy as np
import scipy.ndimage
import ConfigParser
import glob
import os
import sys
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm

def pyGTC(chains, param_names=None, priors=None, have_weight=False):
    
    if param_names is not None:
        have_names = True
    else:
        have_names = False
        
    ##########
    # Setup figure and color scheme

    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['text.usetex']=True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',r'\sansmath']

    colors = [['#7fa5e3','#4c72b0'], ['#88db9b','#55a868'], ['#ffc982','#f5964f'], ['#f78185','#c44e52'], ['#b4a5e5','#8172b2'], ['#333333','#000000'], ['#7fa5e3','#4c72b0'], ['#88db9b','#55a868'], ['#ffc982','#f5964f'], ['#f78185','#c44e52'], ['#b4a5e5','#8172b2'], ['#333333','#000000']]
    maincolors = ['#4c72b0', '#55a868', '#f5964f', '#c44e52', '#8172b2', '#000000', '#4c72b0', '#55a868', '#f5964f', '#c44e52', '#8172b2', '#000000']
    lightblack = '#333333'


    ##########
    Nbins = 30.
    xaxis = np.linspace(0., Nbins**2, Nbins**2)

    Nchains = len(chains)

    # Number of dimensions
    if not have_weight:
	ndim = len(chains[0][0,:])
	for i in range(Nchains):
	    Nsamples = len(chains[i][:,0])
	    chains[i] = np.insert(chains[i], 0, np.ones(Nsamples), axis=1)
    else:
	ndim = len(chains[0][0,1:])

    size_pt = ndim*70.
    fig_width = size_pt / 72.27
    figu = plt.figure(figsize=(fig_width,fig_width))

    xmin=np.empty(ndim)
    xmax=np.empty(ndim)

    ########## 2D contour plots
    for i in range(ndim): # row
        for j in range(ndim): # column
            if j<i:
                fig = figu.add_subplot(ndim,ndim,(i*ndim)+j+1)
            
                data, extent = [], []
                sigma1, sigma2 = [], []
                for k in reversed(range(Nchains)):
                    plotdata, xedges, yedges = np.histogram2d(chains[k][:,j+1], chains[k][:,i+1], weights = chains[k][:,0],bins=Nbins)
                    extent.append([xedges[0],xedges[-1],yedges[0],yedges[-1]])
                    plotdata = plotdata/np.sum(plotdata)
                    histogram_ordered = np.sort(plotdata.flat)
                    histogram_cumulative = np.cumsum(histogram_ordered)
                    temp = np.interp(.32, histogram_cumulative, xaxis)
                    sigma1temp = np.interp(temp, xaxis, histogram_ordered)
                    temp = np.interp(.05, histogram_cumulative, xaxis)
                    sigma2temp = np.interp(temp, xaxis, histogram_ordered)
                    tempdata = scipy.ndimage.gaussian_filter(plotdata.T, sigma=1)
                    xbins = np.delete(xedges+.5*(xedges[1]-xedges[0]), -1)
                    ybins = np.delete(yedges+.5*(yedges[1]-yedges[0]), -1)
                
                    fig.contourf(xbins,ybins,tempdata, levels=[sigma2temp,sigma1temp,1], colors=colors[k])
                
                    data.append(tempdata)
                    sigma1.append(sigma1temp)
                    sigma2.append(sigma2temp)
                for k in range(Nchains):
                    fig.contour(data[Nchains-1-k], [sigma1[Nchains-1-k]], extent=extent[Nchains-1-k], origin='lower', colors=colors[k][1])
                    fig.contour(data[Nchains-1-k], [sigma2[Nchains-1-k]], extent=extent[Nchains-1-k], origin='lower', colors=colors[k][0])
                
            
                # Labels and ticklables
                if i!=ndim-1:
                    fig.get_xaxis().set_ticklabels([])
                else:
                    if have_names:
                        fig.set_xlabel('$'+param_names[j]+'$')
                if j!=0:
                    fig.get_yaxis().set_ticklabels([])
                else:
                    if have_names:
                        fig.set_ylabel('$'+param_names[i]+'$')
            
                my_locator = MaxNLocator(5)
                fig.xaxis.set_major_locator(my_locator)
                my_locator = MaxNLocator(5)
                fig.yaxis.set_major_locator(my_locator)
            
                xmin[j], xmax[j] = fig.get_xlim()
            


    ########## 1D histograms
    for i in range(ndim):
        fig = figu.add_subplot(ndim,ndim,(i*ndim)+i+1)
    
        # 1D histogram
        for k in reversed(range(Nchains)):
            hist1d, edges = np.histogram(chains[k][:,i+1], weights = chains[k][:,0], normed=True, bins=Nbins)
            centers = np.delete(edges+.5*(edges[1]-edges[0]), -1)
            data = scipy.ndimage.gaussian_filter1d((centers,hist1d), sigma=1)
            plt.fill_between(data[0],data[1],0, color=colors[k][0])
            plt.plot(data[0],data[1], ls=':', color=colors[k][0])
    
        if priors is not None:
            if i < len(width) and width[i]>0:
                if i==ndim-1:
                    arr = np.linspace(fig.get_xlim()[0],fig.get_xlim()[1],20)
                    plt.plot(arr,norm.pdf(arr,mean[i],width[i]), color='k')
                arr = np.linspace(xmin[i],xmax[i],20)
                plt.plot(arr,norm.pdf(arr,mean[i],width[i]), color='k')

        ########## Ticks, labels, range
        fig.get_yaxis().set_ticklabels([])
    
        fig.set_ylim(bottom=0)
    
        if i!=ndim-1:
            fig.set_xlim(xmin[i],xmax[i])
            fig.get_xaxis().set_ticklabels([])
        else:
            if have_names:
                fig.set_xlabel('$'+param_names[i]+'$')
        if i==0:
            if have_names:
                fig.set_ylabel('$'+param_names[i]+'$')
                
        my_locator = MaxNLocator(5)
        fig.xaxis.set_major_locator(my_locator)
        my_locator = MaxNLocator(5)
        fig.yaxis.set_major_locator(my_locator)
    

    #######

    figu.subplots_adjust(hspace=0)
    figu.subplots_adjust(wspace=0)

    plt.savefig('GTC.pdf', bbox_inches='tight')

