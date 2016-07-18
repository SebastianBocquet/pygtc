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
        Sample points (length x Ndim) or multiple sets of samples points
        Note: If you are using emcee (http://dan.iel.fm/emcee/current/) - and you should! - you need to pass the EnsembleSampler.flatchain object.
    kwargs:
        See below for all arguments.

    Returns:
    --------
    fig: matplotlib.figure
        the GTC in all its glory
    """
    # Set defaults
    #TODO: make default labels be 1, 2, 3, etc...

    ParamNames=None # label the x and y axes, supports latex
    truths=None # Highlight a point (or several) in parameter space by lines
    priors=None # Draw a Gaussian distribution (or several) in the 1d panels
    weights=None # Provide weight factors for your sample points
    PlotName=None # Save plot as PlotName, else return matplotlib.figure object
    FigureSize=None # Width=height of figure in inches
    ChainLabels=None # Create legend with names for the plotted data
    TruthLabels=None # Label the truth lines in the legend
    NConfidenceLevels=2 # Draw 2d contours out to the NConfidenceLevels level, support 1, 2, 3
    SmoothingKernel=1 # Gaussian smoothing kernel (in pixels)
    TruthColors = ['r','c','g','b','m'] #Default colors for plotting truths

    #Numpy really doesn't like lists of Pandas DataFrame objects
    #so if it gets one, extract array vals and throw away the rest
    try: #Not a list of DFs, but might be a single DF
        shape_length = len(np.shape(chains))
        assert shape_length in [2,3], "unexpected chains shape"
        if shape_length == 2:
            chains = [chains]

        # Read in column names from Pandas DataFrame if exists
        #Also convert DataFrame to simple numpy array to avoid later conflicts
        if hasattr(chains[0], 'columns'):
            ParamNames = list(chains[0].columns.values)
            chains = [df.values for df in chains]

    except ValueError: #Probably a list of pandas DFs
        if hasattr(chains[0], 'columns') and hasattr(chains[0], 'values'):
            ParamNames = list(chains[0].columns.values)
            chains = [df.values for df in chains]

        #If not a DF but something else, this will fail again, which is good
        shape_length = len(np.shape(chains))
        assert shape_length in [2,3], "unexpected chains shape"

    # Parse kwargs
    if kwargs is not None:
        for key, val in kwargs.iteritems():
            if key == 'ParamNames':
                if all(isinstance(s, basestring) for s in val):
                    if len(val) == len(chains[0][0,:]):
                        ParamNames = list(val)
                    else:
                        raise ValueError("ParamNames length must match number of parameters in chains")
                else:
                    raise TypeError("ParamNames must be a list of strings")
            elif key == 'truths':
                truths = val
            elif key == 'priors':
                priors = val
            elif key == 'weights':
                weights = val
            elif key == 'PlotName':
                PlotName = val
            elif key == 'FigureSize':
                FigureSize = val
            elif key == 'ChainLabels':
                ChainLabels = val
            elif key == 'TruthLabels':
                TruthLabels = val
            elif key == 'NConfidenceLevels':
                assert NConfidenceLevels in [1,2,3], "NConfidenceLevels must be 1, 2, or 3"
                NConfidenceLevels = val
            elif key == 'SmoothingKernel':
                SmoothingKernel = val
            elif key == 'TruthColors':
                TruthColors = val
            else:
                raise NameError("illegal keyword argument: " + key)

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

    # Number of chains
    Nchains = len(chains)

    assert Nchains<len(colors), "currently only supports up to "+str(len(colors))+" chains"

    # Check that chains are 2d arrays
    for i in range(Nchains):
        assert len(np.shape(chains[i]))==2, "chain "+str(i)+" has unexpected shape"

    # Number of dimensions
    ndim = len(chains[0][0,:])

    # Manage the sample point weights
    if weights==None:
        # Set unit weights if no weights are provided
        weights = []
        for i in range(Nchains):
            weights.append( np.ones(len(chains[i])) )
    else:
        if len(np.shape(weights))==len(chains[0]):
            weights = [weights]
        for i in range(Nchains):
            if len(weights[i])!=len(chains[i]):
                raise ValueError("missmatch in chain/weights #%d: len(chain) %d, len(weights) %d"%(i,len(chains[i]),len(weights[i])))

    # Use the 68%, 95%, 99% confidence levels
    GaussConfLevels = [.3173, .0455, .0027]

    # Increase dimensionality of truth list by one if single list
    try: #calling len(scalar) will raise a TypeError
        for ts in truths:
            if len(ts)<len(TruthColors):
                raise ValueError("More truths than available colors. Set colors with TruthColors = [colors...]")
    except TypeError: #Probably a single list, so raise dimensionality
        truths = [truths]

    #Now if this raises a TypeError it's for a good reason
    for ts in truths:
        if len(ts)<len(TruthColors):
            raise ValueError("More truths than available colors. Set colors with TruthColors = [colors...]")

    # These are needed to compute the confidence levels
    Nbins = 30.
    NbinsFlat = np.linspace(0., Nbins**2, Nbins**2)

    # Left and right panel boundaries
    xmin, xmax = np.empty(ndim), np.empty(ndim)

    # Figure size (or resolution)
    if FigureSize is None:
        # If no figure size is given, use resolution of 70 ppp (pixel per panel)
        fig_width = ndim*70. / 72.27
    else:
        # User-defined width=height in inches
        if not isinstance(FigureSize, basestring):
            fig_width = FigureSize
        else:
            # Choose from a couple of presets to fit your publication
            if FigureSize=='APJ_column':
                fig_width = 245.26653 / 72.27
            elif FigureSize=='APJ_page':
                fig_width = 513.11743 / 72.27
            elif FigureSize=='MNRAS_column':
                fig_width = 240. / 72.27
            elif FigureSize=='MNRAS_page':
                fig_width = 504. / 72.27
            else:
                raise ValueError("Figuresize %s unknown!"%FigureSize)
    fig = plt.figure(figsize=(fig_width,fig_width))



    ########## 2D contour plots

    ChainLevels = np.ones((Nchains,NConfidenceLevels+1))
    extents = np.empty((Nchains,4))

    for i in range(ndim): # row
        for j in range(ndim): # column
            if j<i:
                ax = fig.add_subplot(ndim,ndim,(i*ndim)+j+1)


                ##### The filled contour plots
                SmoothData = []
                # Draw filled contours in reversed order to have first chain in list on top
                for k in reversed(range(Nchains)):
                    # Create 2d histogram
                    hist2d, xedges, yedges = np.histogram2d(chains[k][:,j], chains[k][:,i], weights=weights[k], bins=Nbins)
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


                ###### Draw contour lines in order to see contours lying on top of each other
                for k in range(Nchains):
                    for l in range(NConfidenceLevels):
                        ax.contour(SmoothData[Nchains-1-k], [ChainLevels[k][NConfidenceLevels-1-l]], extent=extents[k], origin='lower', colors=colors[k][l])


                ##### Truth lines
                if truths is not None:
                    for k in range(len(truths)):
                        # horizontal line
                        if i < len(truths[k]):
                            if truths[k][i] is not None:
                                ax.axhline(truths[k][i], color=TruthColors[k])
                        # vertical line
                        if j < len(truths[k]):
                            if truths[k][j] is not None:
                                ax.axvline(truths[k][j], color=TruthColors[k])
                                # If needed, readjust limits of x_axis
                                lo, hi = ax.get_xlim()
                                if lo>truths[k][j]: lo = truths[k][j]-.05*(hi-lo)
                                if hi<truths[k][j]: hi = truths[k][j]+.05*(hi-lo)
                                ax.set_xlim(lo, hi)


                ##### Ticks & labels

                # x-labels at bottom of plot only
                if i==ndim-1:
                    if ParamNames is not None:
                        ax.set_xlabel(ParamNames[j])
                else:
                    ax.get_xaxis().set_ticklabels([])

                # y-labels for left-most panels only
                if j==0:
                    if ParamNames is not None:
                        ax.set_ylabel(ParamNames[i])
                else:
                    ax.get_yaxis().set_ticklabels([])

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

        ##### 1D histogram
        for k in reversed(range(Nchains)):
            # create 1d histogram
            hist1d, edges = np.histogram(chains[k][:,i], weights = weights[k], normed=True, bins=Nbins)
            # Bin center between histogram edges
            centers = np.delete(edges+.5*(edges[1]-edges[0]), -1)
            # Gaussian smoothing
            PlotData = scipy.ndimage.gaussian_filter1d((centers,hist1d), sigma=SmoothingKernel)
            # Filled histogram
            plt.fill_between(PlotData[0], PlotData[1], 0, color=colors[k][1])
            # Dotted line for hidden histogram
            plt.plot(PlotData[0], PlotData[1], ls=':', color=colors[k][1])


        ##### Truth line
        if truths is not None:
            for k in range(len(truths)):
                if i < len(truths[k]):
                    if truths[k][i] is not None:
                        ax.axvline(truths[k][i], color=TruthColors[k])


        ##### Gaussian prior
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


        ##### Ticks, labels, range

        # No labels on y-axes, lower limit 0
        ax.get_yaxis().set_ticklabels([])
        ax.set_ylim(bottom=0)

        # x-label for bottom-right panel only
        if i==ndim-1:
            if ParamNames is not None:
                ax.set_xlabel(ParamNames[i])
        else:
            ax.set_xlim(xmin[i],xmax[i])
            ax.get_xaxis().set_ticklabels([])

        # y label for top-left panel
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
        ##### Dummy plot for label line color
        LabelColors = []
        ax = fig.add_subplot(ndim,ndim,ndim)
        ax.axis('off')

        ##### Label the data sets
        if ChainLabels is not None:
            # Label for each chain
            for k in range(Nchains):
                ax.plot(0,0, color=colors[k][0], label=ChainLabels[k])
                LabelColors.append(colors[k][0])

        ##### Label the truth lines
        if TruthLabels is not None:
            # Label for each truth
            for k in range(len(TruthLabels)):
                ax.plot(0,0, color=TruthColors[k], label=TruthLabels[k])
                LabelColors.append(TruthColors[k])

        ##### Legend and label colors according to plot
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
