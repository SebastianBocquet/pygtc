from matplotlib import pyplot as plt
import numpy as np
import scipy.ndimage
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm

def plotGTC(chains, **kwargs):
    """Create a beautiful triangle plot - aka Giant Triangle Confusogram (GTC).
    
    Parameters:
    -----------
    chains: 2d array or list of 2d arrays
        Sample points (Npoints x ndim) drawn from a distribution
    kwargs:
        See below for all arguments.

    Returns:
    --------
    fig: matplotlib.figure
        the GTC in all its glory
    """
    assert len(np.shape(chains)) in [2,3], "Your chains' all fucked up, brah!"

    #increase dimensionality by 1 if user only supplies one chain
    if len(np.shape(chains)) == 2:
        chains = [chains]

    #set defaults
    ParamNames=None
    truths=None
    priors=None
    UseWeights=False
    PlotName=None
    FigureSize=None
    ChainLabels=None
    TruthLabels=None
    NConfidenceLevels=2
    SmoothingKernel=1

    #Read in column names from Pandas DataFrame if exists
    if hasattr(chains[0], 'columns'):
        ParamNames = chains[0].columns.values

    #parse kwargs
    if kwargs is not None:
        for key, val in kwargs.iteritems():
            if key == 'ParamNames':
                ParamNames = val
            if key == 'truths':
                truths = val
            if key == 'priors':
                priors = val
            if key == 'UseWeights':
                UseWeights = val
            if key == 'PlotName':
                PlotName = val
            if key == 'FigureSize':
                FigureSize = val
            if key == 'ChainLabels':
                ChainLabels = val
            if key == 'TruthLabels':
                TruthLabels = val
            if key == 'NConfidenceLevels':
                assert NConfidenceLevels in [1,2,3], "ERROR, NConfidenceLevels must be 1, 2, or 3"
                NConfidenceLevels = val
            if key == 'SmoothingKernel':
                SmoothingKernel = val
        
    

    ##########
    #Magic numbers defining Gaussian confidence levels
    GaussConfLevels = [.3173, .0455, .0027]


    ##########
    # Setup figure and colors

    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath']

    colors = [['#4c72b0','#7fa5e3','#b2d8ff'],
        ['#55a868','#88db9b','#bbffce'],
        ['#f5964f','#ffc982','#fffcb5'],
        ['#c44e52','#f78185','#ffb4b8'],
        ['#8172b2','#b4a5e5','#37d8ff'],
        ['#000000','#333333','#666666']]

    LightBlack = '#333333'

    TruthColor = ['r','c','g','b','m']


    ##########
    # Increase dimensionality of truth list by one if single list
    if truths is not None:
        if not isinstance(truths[0], list):
            truths = [truths]
        else:
            assert len(truths)<len(TruthColor), "ERROR: Currently only supported up to "+str(len(TruthColor))+" truths"
            
    
    ########## Get number of chains and dimensions and take care of sample weights (optional)

    # Number of chains
    Nchains = len(chains)
    assert Nchains<len(colors), "ERROR: Currently only supported up to "+str(len(colors))+" chains"

    # Number of dimensions
    if not UseWeights:
        ndim = len(chains[0][0,:])
        for i in range(Nchains):
            Nsamples = len(chains[i][:,0])
            chains[i] = np.insert(chains[i], 0, np.ones(Nsamples), axis=1)
    else:
        ndim = len(chains[0][0,1:])


    ##########
    # These are needed to compute the confidence levels
    Nbins = 30.
    NbinsFlat = np.linspace(0., Nbins**2, Nbins**2)

    # Left and right panel boundaries
    xmin, xmax = np.empty(ndim), np.empty(ndim)


    ########## Define the figure size
    # If no size if given, choose something nice
    if FigureSize is None:
        fig_width = ndim*70. / 72.27
    else:
        if not type(FigureSize=='str'):
            fig_width = FigureSize
        else:
            if FigureSize=='APJ_column':
                fig_width = 245.26653 / 72.27
            elif FigureSize=='APT_page':
                fig_width = 513.11743 / 72.27 # TO DO!
            elif FigureSize=='MNRAS_column':
                fig_width = 240. / 72.27
            elif FigureSize=='MNRAS_page':
                fig_width = 504. / 72.27
            else:
                print "ERROR: Figuresize",FigureSize,"not known!"
                return

    fig = plt.figure(figsize=(fig_width,fig_width))




    ########## 2D contour plots
    ChainLevels = np.ones((Nchains,NConfidenceLevels+1))
    extents = np.empty((Nchains,4))
    
    for i in range(ndim): # row
        for j in range(ndim): # column
            if j<i:
                ax = fig.add_subplot(ndim,ndim,(i*ndim)+j+1)

                ##### The actual contour plots
                SmoothData = []
                
                # Draw filled contours in reversed order to have first chain in list on top
                for k in reversed(range(Nchains)):
                    # Create 2d histogram
                    hist2d, xedges, yedges = np.histogram2d(chains[k][:,j+1], chains[k][:,i+1], weights=chains[k][:,0], bins=Nbins)
                    # image extent, needed below for contour lines
                    extents[k] = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    # Normalize
                    hist2d = hist2d/np.sum(hist2d)
                    # Cumulative 1d distribution
                    HistogramOrdered = np.sort(hist2d.flat)
                    HistogramCumulative = np.cumsum(HistogramOrdered)
                    
                    # Compute confidence levels (from low to high for technical reasons)
                    for l in range(NConfidenceLevels):
                        # Find location of confidence level in 1d HistogramCumulative
                        temp = np.interp(GaussConfLevels[l], HistogramCumulative, NbinsFlat)
                        # Find "height" of confidence level
                        ChainLevels[k][NConfidenceLevels-1-l] = np.interp(temp, NbinsFlat, HistogramOrdered)
                    
                    # Get bin center of histogram edges
                    xbins = np.delete(xedges+.5*(xedges[1]-xedges[0]), -1)
                    ybins = np.delete(yedges+.5*(yedges[1]-yedges[0]), -1)
                    
                    # Apply Gaussian smoothing and plot
                    SmoothData.append( scipy.ndimage.gaussian_filter(hist2d.T, sigma=SmoothingKernel) )
                    ax.contourf(xbins, ybins, SmoothData[-1], levels=ChainLevels[k], colors=colors[k][:NConfidenceLevels][::-1])


                # Draw contour lines in order to see contours lying on top of each other
                for k in range(Nchains):
                    for l in range(NConfidenceLevels):
                        ax.contour(SmoothData[Nchains-1-k], [ChainLevels[k][NConfidenceLevels-1-l]], extent=extents[k], origin='lower', colors=colors[k][l])


                # Truth lines
                if truths is not None:
                    for k in range(len(truths)):
                        if i < len(truths[k]):
                            if truths[k][i] is not None:
                                ax.axhline(truths[k][i], color=TruthColor[k])
                        if j < len(truths[k]):
                            if truths[k][j] is not None:
                                ax.axvline(truths[k][j], color=TruthColor[k])
                                lo, hi = ax.get_xlim()
                                if lo>truths[k][j]: lo = truths[k][j]-.05*(hi-lo)
                                if hi<truths[k][j]: hi = truths[k][j]+.05*(hi-lo)
                                ax.set_xlim(lo, hi)


                # Labels and ticklables
                if i!=ndim-1:
                    ax.get_xaxis().set_ticklabels([])
                else:
                    if ParamNames is not None:
                        ax.set_xlabel(ParamNames[j])
                if j!=0:
                    ax.get_yaxis().set_ticklabels([])
                else:
                    if ParamNames is not None:
                        ax.set_ylabel(ParamNames[i])

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
            # Bin center between histogram edges
            centers = np.delete(edges+.5*(edges[1]-edges[0]), -1)
            # Gaussian smoothing
            PlotData = scipy.ndimage.gaussian_filter1d((centers,hist1d), sigma=SmoothingKernel)
            # Filled histogram
            plt.fill_between(PlotData[0], PlotData[1], 0, color=colors[k][1])
            # Dotted line for hidden histogram
            plt.plot(PlotData[0], PlotData[1], ls=':', color=colors[k][1])


        # Truth line
        if truths is not None:
            for k in range(len(truths)):
                if i < len(truths[k]):
                    if truths[k][i] is not None:
                        ax.axvline(truths[k][i], color=TruthColor[k])


        # Gaussian prior
        if priors is not None:
            if i < len(priors):
                if priors[i]:
                    if priors[i][1]>0:
                        if i==ndim-1:
                            arr = np.linspace(ax.get_xlim()[0],ax.get_xlim()[1],40)
                            plt.plot(arr,norm.pdf(arr,priors[i][0],priors[i][1]), color=LightBlack)
                        else:
                            arr = np.linspace(xmin[i],xmax[i],40)
                            plt.plot(arr,norm.pdf(arr,priors[i][0],priors[i][1]), color=LightBlack)


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
            if ParamNames is not None:
                ax.set_xlabel(ParamNames[i])
        # y label for first panel
        if i==0:
            if ParamNames is not None:
                ax.set_ylabel(ParamNames[i])

        # No more than 5 ticks per panel
        my_locator = MaxNLocator(5)
        ax.xaxis.set_major_locator(my_locator)
        my_locator = MaxNLocator(5)
        ax.yaxis.set_major_locator(my_locator)



    ########## Legend
    if (ChainLabels is not None) or (TruthLabels is not None):
        # Dummy plot for label line color
        LabelColors = []
        ax = fig.add_subplot(ndim,ndim,ndim)
        ax.axis('off')
    
        # Label the data sets
        if ChainLabels is not None:
            # Label for each chain
            for k in range(Nchains):
                ax.plot(0,0, color=colors[k][0], label=ChainLabels[k])
                LabelColors.append(colors[k][0])
        
        # Label the truth lines
        if TruthLabels is not None:
            # Label for each truth
            for k in range(len(TruthLabels)):
                ax.plot(0,0, color=TruthColor[k], label=TruthLabels[k])
                LabelColors.append(TruthColor[k])

        # Legend and label colors according to plot
        leg = plt.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.)
        for color,text in zip(LabelColors,leg.get_texts()):
            text.set_color(color)



    ##########

    # No space between panels
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # Save figure
    if PlotName is not None:
        plt.savefig(PlotName, bbox_inches='tight')

    return fig
