from matplotlib import pyplot as plt
import numpy as np
import scipy.ndimage
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm


#################### Create a full GTC

def plotGTC(chains, **kwargs):
    """Create a beautiful triangle plot - aka Giant Triangle Confusogram (GTC).

    Parameters:
    -----------
    chains: 2d array or list of 2d arrays
        Sample points (length x Ndim) or multiple sets of samples points
        Note: If you are using emcee (http://dan.iel.fm/emcee/current/) - and you should! - you need to pass the EnsembleSampler.flatchain object.
    kwargs:
        weights
        chainLabels
        paramNames
        truthColors
        truths
        truthLabels
        priors
        plotName
        nConfidenceLevels
        nBins
        smoothingKernel
        figureSize
        paramRanges
        colorsOrder
        do1dplots
        doonly1dplot

    Returns:
    --------
    fig: matplotlib.figure
        the GTC in all its glory
    """
    # Setup matplotlb rcParams TODO: make sure this list is exhaustive
    ipyNotebook = kwargs.pop('paramRanges', None)
    
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath']

    colorsDict = { 'blues' : ('#4c72b0','#7fa5e3','#b2d8ff'),
                    'greens' : ('#55a868','#88db9b','#bbffce'),
                    'yellows' : ('#f5964f','#ffc982','#fffcb5'),
                    'reds' : ('#c44e52','#f78185','#ffb4b8'),
                    'purples' : ('#8172b2','#b4a5e5','#37d8ff')}

    colorsOrder = ['blues', 'greens', 'yellows', 'reds', 'purples']
    
    colors = [colorsDict[cs] for cs in colorsOrder]
    
    lightBlack = '#333333'

    tickAngle = 45 #Angle of tick labels


    #Dictionary of size types or whatever:
    mplPPI = plt.rcParams['figure.dpi'] #Matplotlib dots per inch
    figSizeDict = { 'APJ_column' : 245.26653 / mplPPI,
                    'APJ_page' : 513.11743 / mplPPI,
                    'MNRAS_column' : 240. / mplPPI,
                    'MNRAS_page' : 504. / mplPPI}


    #Check the validity of the chains argument:

    #Numpy really doesn't like lists of Pandas DataFrame objects
    #so if it gets one, extract array vals and throw away the rest
    dfColNames = None
    try: #Not a list of DFs, but might be a single DF
        #shapeLength = len(np.shape(chains))
        #assert shapeLength in [2,3], "unexpected chains shape"
        try:
            if chains.ndim == 2:
                chains = [chains]
        except:
            for i in range(len(chains)):
                assert len(chains[i].shape)==2, "unexpected shape of chain %d"%(chains[i])

        # Read in column names from Pandas DataFrame if exists
        #Also convert DataFrame to simple numpy array to avoid later conflicts
        if hasattr(chains[0], 'columns'):
            #Set param names from DataFrame column names, can be overridden later
            dfColNames = list(chains[0].columns.values)
            chains = [df.values for df in chains]

    except ValueError: #Probably a list of pandas DFs
        if hasattr(chains[0], 'columns') and hasattr(chains[0], 'values'):
            dfColNames = list(chains[0].columns.values)
            chains = [df.values for df in chains]

        #If not a DF but something else, this will fail again, which is good
        assert len(np.shape(chains)) in [2,3], "unexpected chains shape"

    #Get number of chains
    nChains = len(chains)
    assert nChains<=len(colorsOrder), "currently only supports up to "+str(len(colorsOrder))+" chains"
    
    # Double-check that chains are 2d arrays
    for i in range(nChains):
        assert len(chains[i].shape)==2, "chain "+str(i)+" has unexpected shape"

    # Number of dimensions (parameters)
    nDim = len(chains[0][0,:])

    #Process kwargs and set defaults
    chainLabels = kwargs.pop('chainLabels', None) #Labels for multiple chains, goes in plot legend
    if chainLabels is not None:
        if not (isinstance(chainLabels, tuple) or isinstance(chainLabels, list)):
            chainLabels = [chainLabels]
        assert len(chainLabels) == nChains, "chainLabels mismatch with number of chains"

    # Label the x and y axes, supports latex
    paramNames = kwargs.pop('paramNames', None)
    if paramNames is not None:
        if not (isinstance(paramNames, tuple) or isinstance(paramNames, list)):
            paramNames = [paramNames]
        #if all(isinstance(s, basestring) for s in paramNames):
            #if len(paramNames) == len(chains[0][0,:]):
            #    paramNames = list(val)
            #else:
        if len(paramNames) != nDim:
            raise ValueError("paramNames length (%d) must match number of parameters in chains (%d)"%(len(paramNames),nDim))
        #else:
        #    raise TypeError("paramNames must be a list of strings")
    elif dfColNames is not None:
        paramNames = dfColNames

    # Custom parameter range
    paramRanges = kwargs.pop('paramRanges', None)
    
    # User-defined color ordering
    customColorsOrder = kwargs.pop('colorsOrder', None) #Labels for multiple chains, goes in plot legend
    if customColorsOrder is not None:
        if not all(color in colorsDict.keys() for color in customColorsOrder):
            raise ValueError("Bad color name in colorsOrder=%s, pick from %s"%(customColorsOrder,colorsDict.keys()))
        if not (isinstance(customColorsOrder, tuple) or isinstance(customColorsOrder, list)):
            customColorsOrder = [customColorsOrder]
        lencustomColorsOrder = len(customColorsOrder)
        assert lencustomColorsOrder<=nChains, "colorsOrder mismatch with number of chains"
        colorsOrder[:lencustomColorsOrder] = customColorsOrder[:lencustomColorsOrder]    
        colors = [colorsDict[cs] for cs in colorsOrder]
        
    # Colors of truth lines
    truthColors = kwargs.pop('truthColors', ['r','c','g','b','m']) #Default supports up to five truths TODO: prettier colors
    truths = kwargs.pop('truths', None) # Highlight a point (or several) in parameter space by lines
    if truths is not None:
        if not (isinstance(truths[0], tuple) or isinstance(truths[0], list)):
            truths = [truths]
        else:
            if len(truths)>len(truthColors):
                raise ValueError("More truths than available colors. Set colors with truthColors = [colors...]")
    
    # Fill up truths lists with None for missing entries
    if truths is not None:
        truthsTemp = []
        for k in range(len(truths)):
            tempList = []
            for i in range(nDim):
                if i<len(truths[k]):
                    temp = truths[k][i] if truths[k][i] is not None else None
                else:
                    temp = None
                tempList.append( temp )
            truthsTemp.append(tempList)
        truths = np.array(truthsTemp)
    
    # Labels for the different truth lines
    truthLabels = kwargs.pop('truthLabels', None) #Labels for multiple truths, goes in plot legend
    if truthLabels is not None:
        if not (isinstance(truthLabels, tuple) or isinstance(truthLabels, list)):
            truthLabels = [truthLabels]
        assert len(truthLabels) == len(truths), "truthLabels mismatch with number of truths"

    #Show Gaussian priors on plots (assuming flat priors)
    priors = kwargs.pop('priors', None)

    # Manage the sample point weights
    weights = kwargs.pop('weights', None)
    if weights==None:
        # Set unit weights if no weights are provided
        weights = []
        for i in range(nChains):
            weights.append( np.ones(len(chains[i])) )
    else:
        if len(weights)==len(chains[0]):
            weights = [weights]
        else:
            try:
                for i in range(nChains):
                    if len(weights[i])!=len(chains[i]):
                        raise ValueError("missmatch in chain/weights #%d: len(chain) %d, len(weights) %d"%(i,len(chains[i]),len(weights[i])))
            except: raise ValueError("your weights don't make sense")

    # Set plotName to save the plot to plotName
    plotName = kwargs.pop('plotName', None) #Um... the name of the plot?!
    if plotName is not None:
        assert isinstance(plotName, basestring), "plotName must be a string type"

    # Define which confidence levels to show
    nConfidenceLevels = kwargs.pop('nConfidenceLevels', 2) #How many of the above confidence levels to show
    assert nConfidenceLevels in [1,2,3], "nConfidenceLevels must be 1, 2, or 3"

    # Data binning and smoothing
    nBins = kwargs.pop('nBins', 30) # Number of bins for 1d and 2d histograms. 30 works...
    smoothingKernel = kwargs.pop('smoothingKernel', 1) #Don't you like smooth data?

    # Figure size: choose size to fit journal, use reasonable default, or provide your own
    figureSize = kwargs.pop('figureSize', None) #Figure size descriptor or figure width=height in inches
    if figureSize is None:
        # If no figure size is given, use resolution of 70 ppp (pixel per panel)
        figureWidth = nDim*70. / mplPPI
    else:
        # User-defined width=height in inches
        if not isinstance(figureSize, basestring):
            figureWidth = figureSize
        else:
            # Choose from a couple of presets to fit your publication
            if figureSize in figSizeDict.keys():
                figureWidth = figSizeDict[figureSize]
            else:
                raise ValueError("figureSize %s unknown!"%figureSize)
                
    # Plot 1d histograms
    do1dplots = kwargs.pop('do1dplots', True)
    
    # Plot ONLY 1d histograms
    doonly1dplot = kwargs.pop('doonly1dplot', False)
    if doonly1dplot:
        for i in range(nChains):
            assert chains[i].shape[1]==1, "Provide chains of shape(Npoints,1) if you only want the 1d histogram"
        do1dplots = True

    #Check to see if there are any remaining keyword arguments
    keys = ''
    for key in kwargs.iterkeys():
        keys = keys + key + ' '
        raise NameError("illegal keyword arguments: " + keys)

    # These are needed to compute the confidence levels
    nBinsFlat = np.linspace(0., nBins**2, nBins**2)

    # Left and right panel boundaries
    xmin, xmax = np.empty(nDim), np.empty(nDim)

    #Create the figure
    fig = plt.figure(figsize=(figureWidth,figureWidth))


    
    
    ########## 2D contour plots
    if not doonly1dplot:
        for i in range(nDim): # row
            for j in range(nDim): # column
                if j<i:
                    ##### Create subplot
                    if do1dplots:
                        ax = fig.add_subplot(nDim,nDim,(i*nDim)+j+1)
                    else:
                        ax = fig.add_subplot(nDim-1,nDim-1,((i-1)*nDim)+j+1)
                    

                    ##### Draw contours and truths
                    # Extract 2d chains
                    chainsForPlot2D = [[chains[k][:,j], chains[k][:,i]] for k in range(nChains)]

                    # Extract 2d truths
                    if truths is not None:
                        truthsForPlot2D = [[truths[k,i], truths[k,j]] for k in range(len(truths))]
                    else:
                        truthsForPlot2D = None
                
                    # Plot!
                    ax = __plot2d(ax, nChains, chainsForPlot2D, weights, nBins, nBinsFlat, smoothingKernel, colors, nConfidenceLevels, truthsForPlot2D, truthColors)
                
                
                    ##### Range
                    if paramRanges is not None:
                        if j<len(paramRanges):
                            if paramRanges[j]:
                                ax.set_xlim(paramRanges[j][0],paramRanges[j][1])
                        if i<len(paramRanges):
                            if paramRanges[i]:
                                ax.set_ylim(paramRanges[i][0],paramRanges[i][1])

                    ##### Ticks & labels
                    ax.get_xaxis().get_major_formatter().set_useOffset(False)
                    ax.get_xaxis().get_major_formatter().set_scientific(False)

                    ax.get_yaxis().get_major_formatter().set_useOffset(False)
                    ax.get_yaxis().get_major_formatter().set_scientific(False)

                    # x-labels at bottom of plot only
                    if i==nDim-1:
                        if paramNames is not None:
                            ax.set_xlabel(paramNames[j])
                    else:
                        ax.get_xaxis().set_ticklabels([])

                    for xLabel in ax.get_xticklabels():
                        xLabel.set_rotation(tickAngle)

                    # y-labels for left-most panels only
                    if j==0:
                        if paramNames is not None:
                            ax.set_ylabel(paramNames[i])
                    else:
                        ax.get_yaxis().set_ticklabels([])

                    for yLabel in ax.get_yticklabels():
                        yLabel.set_rotation(tickAngle)

                    # No more than 5 ticks per panel
                    myLocator = MaxNLocator(5)
                    ax.xaxis.set_major_locator(myLocator)
                    myLocator = MaxNLocator(5)
                    ax.yaxis.set_major_locator(myLocator)
                    
                    # Remove first and last tick location
                    ax.xaxis.set_ticks(ax.xaxis.get_ticklocs()[1:-1])
                    ax.yaxis.set_ticks(ax.yaxis.get_ticklocs()[1:-1])

                    # Limits to be applied to 1d histograms
                    xmin[j], xmax[j] = ax.get_xlim()
                


    if do1dplots:
        ########## 1D histograms
        for i in range(nDim):
            ##### Create subplot
            ax = fig.add_subplot(nDim,nDim,(i*nDim)+i+1)


            ##### Plot histograms, truths, Gaussians
            # Extract 1d chains
            chainsForPlot1D = [chains[k][:,i] for k in range(nChains)]

            # Extract 1d truths
            if truths is not None:
                truthsForPlot1D = [truths[k,i] for k in range(len(truths))]
            else:
                truthsForPlot1D = None
            
            # Extract 1d prior
            if priors is not None:
                if i<len(priors):
                    if priors[i] and priors[i][1]>0:
                        prior1d = priors[i]
            else:
                prior1d = None
                    
            # Plot!
            ax = __plot1d(ax, nChains, chainsForPlot1D, weights, nBins, smoothingKernel, colors, truthsForPlot1D, truthColors, prior1d, lightBlack)
        

            ##### Ticks, labels, range
            ax.get_xaxis().get_major_formatter().set_useOffset(False)
            ax.get_xaxis().get_major_formatter().set_scientific(False)

            # No ticks or labels on y-axes, lower limit 0
            ax.get_yaxis().set_ticklabels([])
            ax.yaxis.set_ticks_position('none')
            ax.set_ylim(bottom=0)
            ax.xaxis.set_ticks_position('bottom')

            # x-label for bottom-right panel only
            if i==nDim-1:
                if paramNames is not None:
                    ax.set_xlabel(paramNames[i])
            else:
                ax.set_xlim(xmin[i],xmax[i])
                ax.get_xaxis().set_ticklabels([])

            for xLabel in ax.get_xticklabels():
                xLabel.set_rotation(tickAngle)

            # y label for top-left panel
            if i==0:
                if doonly1dplot:
                    ax.set_ylabel('Probability')
                elif paramNames is not None:
                    ax.set_ylabel(paramNames[i])

            # No more than 5 ticks per panel
            myLocator = MaxNLocator(5)
            ax.xaxis.set_major_locator(myLocator)

            # Remove first and last tick location
            ax.xaxis.set_ticks(ax.xaxis.get_ticklocs()[1:-1])
            


    ########## Legend
    if (chainLabels is not None) or (truthLabels is not None):
        ##### Dummy plot for label line color
        labelColors = []
        if not doonly1dplot:
            ax = fig.add_subplot(nDim,nDim,nDim)
            ax.axis('off')
        else:
            xmin, xmax = ax.get_xlim()

        ##### Label the data sets
        if chainLabels is not None:
            # Label for each chain
            for k in range(nChains):
                ax.plot(0,0, color=colors[k][0], label=chainLabels[k])
                labelColors.append(colors[k][0])

        ##### Label the truth lines
        if truthLabels is not None:
            # Label for each truth
            for k in range(len(truthLabels)):
                ax.plot(0,0, color=truthColors[k], label=truthLabels[k])
                labelColors.append(truthColors[k])
        
        # Set xlim back to what the data wanted
        if doonly1dplot:
            ax.set_xlim(xmin, xmax)
            

        ##### Legend and label colors according to plot
        leg = plt.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.)
        for color,text in zip(labelColors,leg.get_texts()):
            text.set_color(color)
    
    



    ##########

    # No space between panels
    fig.subplots_adjust(hspace=0)
    fig.subplots_adjust(wspace=0)

    # Save figure
    if plotName is not None:
        plt.savefig(plotName, bbox_inches='tight')

    return fig
    


#################### Create single 1d panel

def __plot1d(ax, nChains, chains1d, weights, nBins, smoothingKernel, colors, truths1d, truthColors, prior1d, lightBlack):

    ##### 1D histogram
    for k in reversed(range(nChains)):
        # create 1d histogram
        hist1d, edges = np.histogram(chains1d[k], weights = weights[k], normed=True, bins=nBins)
        # Bin center between histogram edges
        centers = np.delete(edges+.5*(edges[1]-edges[0]), -1)
        # Gaussian smoothing
        plotData = scipy.ndimage.gaussian_filter1d((centers,hist1d), sigma=smoothingKernel)
        # Filled histogram
        plt.fill_between(plotData[0], plotData[1], 0, color=colors[k][1])
        # Dotted line for hidden histogram
        plt.plot(plotData[0], plotData[1], ls=':', color=colors[k][1])


    ##### Truth line
    if truths1d is not None:
        for k in range(len(truths1d)):
            if truths1d[k] is not None:
                ax.axvline(truths1d[k], color=truthColors[k])


    ##### Gaussian prior
    if prior1d is not None:
        # Plot prior in -4 to +4 sigma range
        arr = np.linspace(prior1d[0]-4*prior1d[1], prior1d[0]+4*prior1d[1], 40)
        plt.plot(arr,norm.pdf(arr,prior1d[0],prior1d[1]), color=lightBlack)
    
    return ax
    
    
    
#################### Create single 2d panel

def __plot2d(ax, nChains, chains2d, weights, nBins, nBinsFlat, smoothingKernel, colors, nConfidenceLevels, truths2d, truthColors):
    
    #generateLabels defaults to True for standalone plotting
    #pltGTC sets generateLabels=False and handles its own
    
    # Use the 68%, 95%, and 99% confidence levels, which look different in 2D
    gaussConfLevels = [.3173, .0455, .0027]
    
    # 
    chainLevels = np.ones((nChains,nConfidenceLevels+1))
    extents = np.empty((nChains,4))
    
    ##### The filled contour plots
    smoothData = []
    # Draw filled contours in reversed order to have first chain in list on top
    for k in reversed(range(nChains)):
        # Create 2d histogram
        hist2d, xedges, yedges = np.histogram2d(chains2d[k][0], chains2d[k][1], weights=weights[k], bins=nBins)
        # image extent, needed below for contour lines
        extents[k] = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        # Normalize
        hist2d = hist2d/np.sum(hist2d)
        # Cumulative 1d distribution
        histOrdered = np.sort(hist2d.flat)
        histCumulative = np.cumsum(histOrdered)

        # Compute confidence levels (from low to high for technical reasons)
        for l in range(nConfidenceLevels):
            # Find location of confidence level in 1d histCumulative
            temp = np.interp(gaussConfLevels[l], histCumulative, nBinsFlat)
            # Find "height" of confidence level
            chainLevels[k][nConfidenceLevels-1-l] = np.interp(temp, nBinsFlat, histOrdered)

        # Get bin center of histogram edges
        xbins = np.delete(xedges+.5*(xedges[1]-xedges[0]), -1)
        ybins = np.delete(yedges+.5*(yedges[1]-yedges[0]), -1)

        # Apply Gaussian smoothing and plot
        smoothData.append( scipy.ndimage.gaussian_filter(hist2d.T, sigma=smoothingKernel) )
        ax.contourf(xbins, ybins, smoothData[-1], levels=chainLevels[k], colors=colors[k][:nConfidenceLevels][::-1])


    ###### Draw contour lines in order to see contours lying on top of each other
    for k in range(nChains):
        for l in range(nConfidenceLevels):
            ax.contour(smoothData[nChains-1-k], [chainLevels[k][nConfidenceLevels-1-l]], extent=extents[k], origin='lower', colors=colors[k][l])


    ##### Truth lines
    if truths2d is not None:
        for k in range(len(truths2d)):
            # horizontal line
            if truths2d[k][0] is not None:
                ax.axhline(truths2d[k][0], color=truthColors[k])
            # vertical line
            if truths2d[k][1] is not None:
                ax.axvline(truths2d[k][1], color=truthColors[k])

    return ax
