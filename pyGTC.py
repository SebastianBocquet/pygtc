from matplotlib import pyplot as plt
import numpy as np
import scipy.ndimage
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm

def plotGTC(chains, **kwargs):

    assert np.shape(chains) in [2,3], "chains shape unexpected"

    #increase dimensionality by 1 if user only supplies one chain
    if len(np.shape(chains)) == 2:
        chains = [chains]

    #set defaults
    param_names=None
    truths=None
    priors=None
    have_weight=False
    plotname=None
    figuresize=None
    chainlabels=None
    confidencelevels=2
    smoothing_kernel=1

    #Read in column names from Pandas DataFrame if exists
    if hasattr(chains[0], 'columns'):
        param_names = list(chains[0].columns.values)

    #parse kwargs
    if kwargs is not None:
        for key, val in kwargs.iteritems():
            if key == 'param_names':
                if all(isinstance(s, basestring) for s in val):
                    if len(val) == chains[0][0,:]
                        param_names = list(val)
                    else:
                        raise ValueError("param_names length must match number of parameters in chains")
                else:
                    raise TypeError("param_names must be a list of strings")
            if key == 'truths':
                truths = val
            if key == 'priors':
                priors = val
            if key == 'have_weight':
                have_weight = val
            if key == 'plotname':
                plotname = val
            if key == 'figuresize':
                figuresize = val
            if key == 'chainlabels':
                chainlabels = val
            if key == 'confidencelevels':
                assert confidencelevels in [1,2,3], "ERROR, confidencelevels must be 1, 2, or 3"
                confidencelevels = val
            if key == 'smoothing_kernel':
                smoothing_kernel = val


    ##########
    #Magic numbers defining 2D Gaussian confidence levels
    conflevels = (.3173, .0455, .0027)

    ##########
    # Setup figure and colors

    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['text.usetex']=True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath']

    colors = [['#4c72b0','#7fa5e3','#b2d8ff'],
        ['#55a868','#88db9b','#bbffce'],
        ['#f5964f','#ffc982','#fffcb5'],
        ['#c44e52','#f78185','#ffb4b8'],
        ['#8172b2','#b4a5e5','#37d8ff'],
        ['#000000','#333333','#666666']]

    lightblack = '#333333'

    #TODO: Add arbitrary number of truths
    truthcolor = '#666666'


    ########## Get number of chains and dimensions and take care of sample weights (optional)

    # Number of chains
    Nchains = len(chains)
    assert Nchains <= 6, "Currently only supports up to 6 chains"

    # Number of dimensions
    if not have_weight:
        ndim = len(chains[0][0,:])
        for i in range(Nchains):
            Nsamples = len(chains[i][:,0])
            chains[i] = np.insert(chains[i], 0, np.ones(Nsamples), axis=1)
    else:
        ndim = len(chains[0][0,1:])


    ##########
    # These are needed to compute the confidence levels
    Nbins = 30.
    xaxis = np.linspace(0., Nbins**2, Nbins**2)

    # Left and right panel boundaries
    xmin, xmax = np.empty(ndim), np.empty(ndim)


    ########## Define the figure size
    # If no size if given, choose something nice
    if figuresize is None:
        fig_width = ndim*70. / 72.27
    else:
        if not type(figuresize=='str'):
            fig_width = figuresize
        else:
            if figuresize=='APJ_column':
                fig_width = 245.26653 / 72.27
            elif figuresize=='APT_page':
                fig_width = 513.11743 / 72.27 # TO DO!
            elif figuresize=='MNRAS_column':
                fig_width = 240. / 72.27
            elif figuresize=='MNRAS_page':
                fig_width = 504. / 72.27
            else:
                print "ERROR: Figuresize",figuresize,"not known!"
                return

    fig = plt.figure(figsize=(fig_width,fig_width))




    ########## 2D contour plots
    for i in range(ndim): # row
        for j in range(ndim): # column
            if j<i:
                ax = fig.add_subplot(ndim,ndim,(i*ndim)+j+1)

                ##### The actual contour plots
                data, extent, chainlevels = [], [], []

                # Draw filled contours in reversed order to have first chain in list on top
                for k in reversed(range(Nchains)):
                    # Create 2d histogram
                    plotdata, xedges, yedges = np.histogram2d(chains[k][:,j+1], chains[k][:,i+1], weights=chains[k][:,0], bins=Nbins)
                    # image extent
                    extent.append([xedges[0], xedges[-1], yedges[0], yedges[-1]])
                    # Normalize
                    plotdata = plotdata/np.sum(plotdata)
                    # Computed cumulative 1d distribution
                    histogram_ordered = np.sort(plotdata.flat)
                    histogram_cumulative = np.cumsum(histogram_ordered)

                    # Compute confidence levels (from low to high for technical reasons)
                    templist = []
                    for l in reversed(range(confidencelevels)):
                        # Find location of confidence level in 1d histogram_cumulative
                        temp = np.interp(conflevels[l], histogram_cumulative, xaxis)
                        # Find "height" of confidence level
                        templist.append( np.interp(temp, xaxis, histogram_ordered) )
                    templist.append(1)

                    # Get bin center of histogram edges
                    xbins = np.delete(xedges+.5*(xedges[1]-xedges[0]), -1)
                    ybins = np.delete(yedges+.5*(yedges[1]-yedges[0]), -1)

                    # Apply Gaussian smoothing and plot
                    tempdata = scipy.ndimage.gaussian_filter(plotdata.T, sigma=smoothing_kernel)
                    ax.contourf(xbins, ybins, tempdata, levels=templist, colors=colors[k][:confidencelevels][::-1])

                    # Store for next step (contour line)
                    data.append(tempdata)
                    chainlevels.append(templist[:-1])


                # Draw contour lines in order to see contours lying on top of each other
                for k in range(Nchains):
                    for l in range(confidencelevels):
                        ax.contour(data[Nchains-1-k], [chainlevels[k][l]], extent=extent[Nchains-1-k], origin='lower', colors=colors[k][confidencelevels-1-l])


                # Truth lines
                if truths is not None:
                    if i < len(truths):
                        if truths[i] is not None:
                            ax.axhline(truths[i], color=truthcolor, ls='--')
                    if j < len(truths):
                        if truths[j] is not None:
                            ax.axvline(truths[j], color=truthcolor, ls='--')
                            lo, hi = ax.get_xlim()
                            if lo>truths[j]: lo = truths[j]-.05*(hi-lo)
                            if hi<truths[j]: hi = truths[j]+.05*(hi-lo)
                            ax.set_xlim(lo, hi)


                # Labels and ticklables
                if i!=ndim-1:
                    ax.get_xaxis().set_ticklabels([])
                else:
                    if param_names is not None:
                        ax.set_xlabel(param_names[j])
                if j!=0:
                    ax.get_yaxis().set_ticklabels([])
                else:
                    if param_names is not None:
                        ax.set_ylabel(param_names[i])

                # No more than 5 ticks per panel
                my_locator = MaxNLocator(5)
                ax.xaxis.set_major_locator(my_locator)
                my_locator = MaxNLocator(5)
                ax.yaxis.set_major_locator(my_locator)

                # Limits to be applied to 1d histograms
                xmin[j], xmax[j] = ax.get_xlim()



    ########## 1D histograms
    for i in range(ndim):
        ax = fig.add_subplot(ndim,ndim,(i*ndim)+i+1)

        # 1D histogram
        for k in reversed(range(Nchains)):
            # create 1d histogram
            hist1d, edges = np.histogram(chains[k][:,i+1], weights = chains[k][:,0], normed=True, bins=Nbins)
            # Get bin center of histogram edges
            centers = np.delete(edges+.5*(edges[1]-edges[0]), -1)
            # Gaussian smoothing
            data = scipy.ndimage.gaussian_filter1d((centers,hist1d), sigma=smoothing_kernel)
            # Filled histogram
            plt.fill_between(data[0],data[1],0, color=colors[k][1])
            # Dotted line for hidden histogram
            plt.plot(data[0],data[1], ls=':', color=colors[k][1])


        # Truth line
        if truths is not None:
            if i < len(truths):
                if truths[i] is not None:
                    ax.axvline(truths[i], color=truthcolor, ls='--')


        # Gaussian prior
        if priors is not None:
            if i < len(priors):
                if priors[i]:
                    if priors[i][1]>0:
                        if i==ndim-1:
                            arr = np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],40)
                            plt.plot(arr,norm.pdf(arr,priors[i][0],priors[i][1]), color=lightblack)
                        else:
                            arr = np.linspace(xmin[i],xmax[i],40)
                            plt.plot(arr,norm.pdf(arr,priors[i][0],priors[i][1]), color=lightblack)


        ########## Ticks, labels, range
        # No labels on y-axes, lower limit 0
        ax.get_yaxis().set_ticklabels([])
        ax.set_ylim(bottom=0)

        # Panels in middle of triangle have no x label
        if i!=ndim-1:
            ax.set_xlim(xmin[i],xmax[i])
            ax.get_xaxis().set_ticklabels([])
        # x label for last panel
        else:
            if param_names is not None:
                ax.set_xlabel(param_names[i])
        # y label for first panel
        if i==0:
            if param_names is not None:
                ax.set_ylabel(param_names[i])

        # No more than 5 ticks per panel
        my_locator = MaxNLocator(5)
        ax.xaxis.set_major_locator(my_locator)
        my_locator = MaxNLocator(5)
        ax.yaxis.set_major_locator(my_locator)



    ########## Legend
    if chainlabels is not None:
        # Dummy plot for label line color
        labelcolors = []
        ax = fig.add_subplot(ndim,ndim,ndim)

        # Label for each chain
        for k in range(Nchains):
            ax.plot(0,0, color=colors[k][0], label=chainlabels[k])
            labelcolors.append(colors[k][0])
        ax.axis('off')

        # Legend and label colors according to plot
        leg = plt.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.)
        for color,text in zip(labelcolors,leg.get_texts()):
            text.set_color(color)



    ##########

    # No space between panels
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # Save figure
    if plotname is not None:
        plt.savefig(plotname, bbox_inches='tight')

    return fig
