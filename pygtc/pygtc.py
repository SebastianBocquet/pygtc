from matplotlib import pyplot as plt
import numpy as np
import scipy.ndimage
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm


#################### Create a full GTC

def plotGTC(chains, **kwargs):
    r"""Make a great looking Giant Triangle Confusogram (GTC) with one line of code!
    A GTC is a lot like a triangle (or corner) plot, but you get to put as many
    sets of data, and overlay as many truths as you like. That's what can make
    it so *confusing*!

    Parameters
    ----------
    chains : array-like[nDims , samples] or a list[[nDims, samples1], ...]
        All chains (where a chain is [nDims, samples]) in the list must have
        the same number of dimensions. Note: If you are using ``emcee``
        (http://dan.iel.fm/emcee/current/) - and you should! - each element
        of chains is an ``EnsembleSampler.flatchain`` object.

    weights : array-like[nSamples] or a list[[nSamples1], ...]
        Weights for the sample points. If a 1d array is passed, the same
        weights are used for all dimension in chains. If a list of 1d arrays
        is passed, there must be a weights array for each dimension of
        `chains`. Default weight is 1.

    chainLabels : array-like[nChains]
        A list of text labels describing each chain passed to chains.
        len(chainLabels) must equal len(chains). chainLabels supports LaTex
        commands enclosed in $..$. Additionally, you can pass None as a
        label. Default is ``None``.

    paramNames : list-like[nDims]
        A list of text labels describing each dimension of chains.
        len(paramNames) must equal len(chains[0]). paramNames supports LaTex
        commands enclosed in $..$. Additionally, you can pass None as a
        label. Default is None, however if you pass a ``pandas.DataFrame``
        object, `paramNames` defaults to the ``DataFrame`` column names.

    truths : list-like[nDims] or [[nDims], ...]
        A list of parameter values, one for each parameter in `chains` to
        hilite in the GTC parameter space, or a list of lists of values to
        hilite in the parameter space. For each set of truths passed to
        `truths`, there must be a value corresponding to every dimension
        in `chains`, although any value may be `None`. Default is
        ``None``.

    truthLabels : list-like[nTruths]
        A list of labels, one for each list passed to truths. truthLabels
        supports LaTex commands enclosed in $..$. Additionally, you can pass
        ``None`` as a label. Default is ``None``.

    truthColors : list-like[nTruths]
        User-defined colors for the truth lines, must be one per set of
        truths passed to `truths`. Default colors are
        ``['r','c','g','b','m']``.

    priors : list of tuples [(mu1, sigma1), ...]
        Each tuple describes a Gaussian to be plotted over that parameter's
        histogram. The number of priors must equal the number of dimensions
        in `chains`. Default is ``None``.

    plotName : string
        A path to save the GTC to in pdf form. Default is ``None``.

    nConfidenceLevels : int
        The number of contours to plot in the 2d histograms. Each contour
        corresponds to a sigma. May be 1, 2, or 3. Default is 2.

    nBins : int
        An integer describing the number of bins used to compute the
        histograms. Default is 30.

    smoothingKernel : int
        An integer describing the size of the Gaussian smoothing kernel in
        bins. Default is 1.

    figureSize : float or string
        A number in inches describing the length = width of the GTC, or a
        string indicating a predefined journal setting. Default is 70/dpi
        where ``dpi = plt.rcParams['figure.dpi']``.

    panelSpacing : string
        Options are ``'loose'`` or ``'tight'``. Determines whether there is
        some space between the subplots of the GTC or not. Default is
        ``'tight'``.

    paramRanges : list of tuples [nDim]
        Set the boundaries of each paramter range. Must provide a tuples for
        each dimension of `chains`. If ``None`` is provided for a
        parameter, the range defaults to 4 sigma bounds.

    colorsOrder : list-like[nDims]
        The color order for chains passed to `chains`. Default is
        ``['blues', 'greens', 'yellows', 'reds', 'purples']``. Currently,
        ``pygtc`` is limited to these color values, so you can reorder them,
        but can't yet define your own colors.

    do1dPlots : bool
        Whether or not 1d histrograms are plotted on the diagonal. Default
        is ``True``.

    doOnly1dPlot : bool
        Pot only ONE 1d histogram. If this is True, then chains must have
        shape ``(samples,1)``. Default is ``False``.

    Returns
    -------
    fig : ``matplotlib.figure`` object
        You can do all sorts of fun things with this in terms of
        customization after it gets returned. If you are using a ``JuPyter``
        notebook with inline plotting enabled, you should assign a variable
        to catch the return or else the figure will plot twice.

    """

    ##### Matplotlib and figure setting
    # Mtplotlb rcParams TODO: make sure this list is exhaustive
    plt.rcParams['legend.fontsize'] = 9
    plt.rcParams['axes.labelsize'] = 9
    plt.rcParams['xtick.labelsize'] = 6
    plt.rcParams['ytick.labelsize'] = 6
    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}', r'\sansmath']

    #Set up some colors
    truthsDefaultColors = ['r','c','g','b','m']
    colorsDict = { 'blues' : ('#4c72b0','#7fa5e3','#b2d8ff'),
                    'greens' : ('#55a868','#88db9b','#bbffce'),
                    'yellows' : ('#f5964f','#ffc982','#fffcb5'),
                    'reds' : ('#c44e52','#f78185','#ffb4b8'),
                    'purples' : ('#8172b2','#b4a5e5','#37d8ff')}
    colorsOrder = ['blues', 'greens', 'yellows', 'reds', 'purples']
    colors = [colorsDict[cs] for cs in colorsOrder]
    lightBlack = '#333333'

    #Angle of tick labels
    tickAngle = 45

    #Dictionary of size types or whatever:
    mplPPI = plt.rcParams['figure.dpi'] #Matplotlib dots per inch
    figSizeDict = { 'APJ_column' : 245.26653 / mplPPI,
                    'APJ_page' : 513.11743 / mplPPI,
                    'MNRAS_column' : 240. / mplPPI,
                    'MNRAS_page' : 504. / mplPPI}


    ##### Check the validity of the chains argument:

    #Numpy really doesn't like lists of Pandas DataFrame objects
    #so if it gets one, extract array vals and throw away the rest
    dfColNames = None
    try: #Not a list of DFs, but might be a single DF
        try:
            # Check if single numpy 2d chain
            if chains.ndim == 2:
                chains = [chains]
        except:
            pass

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

    #Get number of chains
    nChains = len(chains)
    assert nChains<=len(colorsOrder), "currently only supports up to "+str(len(colorsOrder))+" chains"

    # Check that each chain looks reasonable (2d shape)
    for i in range(nChains):
        assert len(chains[i].shape)==2, "unexpected shape of chain %d"%(chains[i])

    # Number of dimensions (parameters), check all chains have same nDim
    nDim = len(chains[0][0,:])
    for i in range(nChains):
        nDimi = len(chains[i][0,:])
        assert nDimi==nDim, "chain %d has unexpected number of dimensions %d"%(i,nDimi)

    # Labels for multiple chains, goes in plot legend
    chainLabels = kwargs.pop('chainLabels', None)
    if chainLabels is not None:
        # Convert to list if only one label
        if isstr(chainLabels):
            chainLabels = [chainLabels]
        # Check that number of labels equals number of chains
        assert len(chainLabels) == nChains, "chainLabels mismatch with number of chains"
        # Check that it's a list of strings
        assert all(isstr(s) for s in chainLabels), "chainLabels must be list of strings"

    # Label the x and y axes, supports latex
    paramNames = kwargs.pop('paramNames', None)
    if paramNames is not None:
        # Convert to list if only one name
        if isstr(paramNames):
            paramNames = [paramNames]
        # Check that number of paramNames equals nDim
        assert len(paramNames) == nDim, "paramNames mismatch with number of dimensions"
        # Check that it's a list of strings
        assert all(isstr(s) for s in paramNames), "paramNames must be list of strings"
    elif dfColNames is not None:
        paramNames = dfColNames

    # Custom parameter range
    paramRanges = kwargs.pop('paramRanges', None)
    if paramRanges is not None:
        assert len(paramRanges)==nDim, "paramRanges must match number of parameters"

    # User-defined color ordering
    customColorsOrder = kwargs.pop('colorsOrder', None) #Labels for multiple chains, goes in plot legend
    if customColorsOrder is not None:
        # Convert to list if only one entry
        if isstr(customColorsOrder):
            customColorsOrder = [customColorsOrder]
        lencustomColorsOrder = len(customColorsOrder)
        if not all(color in colorsDict.keys() for color in customColorsOrder):
            raise ValueError("Bad color name in colorsOrder=%s, pick from %s"%(customColorsOrder,colorsDict.keys()))
        colorsOrder[:lencustomColorsOrder] = customColorsOrder[:lencustomColorsOrder]
        colors = [colorsDict[cs] for cs in colorsOrder]

    # Highlight a point (or several) in parameter space by lines
    # Colors of truth lines
    truthColors = kwargs.pop('truthColors', truthsDefaultColors) #Default supports up to five truths TODO: prettier colors
    truths = kwargs.pop('truths', None)
    if truths is not None:
        # Convert to list if needed
        if len(np.shape(truths))==1:
            truths = [truths]
        truths = np.array(truths)
        assert np.shape(truths)[0]<=len(truthColors), "More truths than available colors. Set colors with truthColors = [colors...]"
        assert np.shape(truths)[1]==nDim, "Each list of truths must match number of parameters"

    # Labels for the different truth lines
    truthLabels = kwargs.pop('truthLabels', None) #Labels for multiple truths, goes in plot legend
    if truthLabels is not None:
        # Convert to list if only one label
        if isstr(truthLabels):
            truthLabels = [truthLabels]
        # Check that it's a list of strings
        assert all(isstr(s) for s in truthLabels), "truthLabels must be list of strings"
        assert len(truthLabels) == len(truths), "truthLabels mismatch with number of truths"

    #Show Gaussian priors on plots (assuming flat priors)
    priors = kwargs.pop('priors', None)
    if priors is not None:
        assert len(priors)==nDim, "List of priors must match number of parameters"
        for i in range(nDim):
            if priors[i]:
                assert priors[i][1]>0, "Prior width must be positive"

    # Manage the sample point weights
    weights = kwargs.pop('weights', None)
    if weights==None:
        # Set unit weights if no weights are provided
        weights = [np.ones(len(chains[i])) for i in range(nChains)]
    else:
        if len(weights)==len(chains[0]):
            weights = [weights]
        for i in range(nChains):
            assert len(weights[i])==len(chains[i]), "missmatch in chain/weights #%d: len(chain) %d, len(weights) %d"%(i,len(chains[i]),len(weights[i]))

    # Set plotName to save the plot to plotName
    plotName = kwargs.pop('plotName', None) #Um... the name of the plot?!
    if plotName is not None:
        assert isstr(plotName), "plotName must be a string type"

    # Define which confidence levels to show
    nConfidenceLevels = kwargs.pop('nConfidenceLevels', 2) #How many of the above confidence levels to show
    assert nConfidenceLevels in [1,2,3], "nConfidenceLevels must be 1, 2, or 3"

    # Data binning and smoothing
    nBins = kwargs.pop('nBins', 30) # Number of bins for 1d and 2d histograms. 30 works...
    smoothingKernel = kwargs.pop('smoothingKernel', 1) #Don't you like smooth data?
    if smoothingKernel>=nBins/10:
        print("Wow, that's a huge smoothing kernel! You sure you want its scale to be %.1f percent of the plot?!"%(100.*float(smoothingKernel)/float(nBins)))

    # Figure size: choose size to fit journal, use reasonable default, or provide your own
    figureSize = kwargs.pop('figureSize', None) #Figure size descriptor or figure width=height in inches
    if figureSize is None:
        # If no figure size is given, use resolution of 70 ppp (pixel per panel)
        figureWidth = nDim*70. / mplPPI
    else:
        # User-defined width=height in inches
        if not isstr(figureSize):
            figureWidth = figureSize
        else:
            # Choose from a couple of presets to fit your publication
            if figureSize in figSizeDict.keys():
                figureWidth = figSizeDict[figureSize]
            else:
                raise ValueError("figureSize %s unknown"%figureSize)

    # Space between panels
    panelSpacing = kwargs.pop('panelSpacing', 'tight')

    # Plot 1d histograms
    do1dPlots = kwargs.pop('do1dPlots', True)

    # Plot ONLY 1d histograms
    doOnly1dPlot = kwargs.pop('doOnly1dPlot', False)
    if doOnly1dPlot:
        for i in range(nChains):
            assert chains[i].shape[1]==1, "Provide chains of shape(Npoints,1) if you only want the 1d histogram"
        do1dPlots = True

    #Check to see if there are any remaining keyword arguments
    keys = ''
    for key in iter(kwargs.keys()):
        keys = keys + key + ' '
        raise NameError("illegal keyword arguments: " + keys)

    # These are needed to compute the confidence levels
    nBinsFlat = np.linspace(0., nBins**2, nBins**2)

    # Left and right panel boundaries
    xmin, xmax = np.empty(nDim), np.empty(nDim)

    #Create the figure
    fig = plt.figure(figsize=(figureWidth,figureWidth))




    ########## 2D contour plots
    if not doOnly1dPlot:
        for i in range(nDim): # row
            for j in range(nDim): # column
                if j<i:
                    ##### Create subplot
                    if do1dPlots:
                        ax = fig.add_subplot(nDim,nDim,(i*nDim)+j+1)
                    else:
                        ax = fig.add_subplot(nDim-1,nDim-1,((i-1)*(nDim-1))+j+1)


                    ##### Draw contours and truths
                    # Extract 2d chains
                    chainsForPlot2D = [[chains[k][:,j], chains[k][:,i]] for k in range(nChains)]

                    # Extract 2d truths
                    truthsForPlot2D = None
                    if truths is not None:
                        truthsForPlot2D = [[truths[k,i], truths[k,j]] for k in range(len(truths))]

                    # Plot!
                    ax = __plot2d(ax, nChains, chainsForPlot2D, weights, nBins, nBinsFlat, smoothingKernel, colors, nConfidenceLevels, truthsForPlot2D, truthColors)


                    ##### Range
                    if paramRanges is not None:
                        if paramRanges[j]:
                            ax.set_xlim(paramRanges[j][0],paramRanges[j][1])
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



    if do1dPlots:
        ########## 1D histograms
        for i in range(nDim):
            ##### Create subplot
            ax = fig.add_subplot(nDim,nDim,(i*nDim)+i+1)


            ##### Plot histograms, truths, Gaussians
            # Extract 1d chains
            chainsForPlot1D = [chains[k][:,i] for k in range(nChains)]

            # Extract 1d truths
            truthsForPlot1D = None
            if truths is not None:
                truthsForPlot1D = [truths[k,i] for k in range(len(truths))]

            # Extract 1d prior
            prior1d = None
            if priors is not None:
                if priors[i] and priors[i][1]>0:
                    prior1d = priors[i]

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
                if doOnly1dPlot:
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
        if not doOnly1dPlot:
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
        if doOnly1dPlot:
            ax.set_xlim(xmin, xmax)


        ##### Legend and label colors according to plot
        leg = plt.legend(loc='upper right', fancybox=True)
        leg.get_frame().set_alpha(0.)
        for color,text in zip(labelColors,leg.get_texts()):
            text.set_color(color)





    ##########

    # No space between panels
    space = 0
    if panelSpacing=='loose':
        space = .05
    fig.subplots_adjust(hspace=space)
    fig.subplots_adjust(wspace=space)

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



#################### Check for basestring in python 2/3 compatible way
def isstr(s):
    try:
        isinstance("", basestring)
        return isinstance(s, basestring)
    except NameError:
        return isinstance(s, str)
