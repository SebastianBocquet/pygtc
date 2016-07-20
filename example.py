from matplotlib import pyplot as plt
import numpy as np
import pyGTC

""" Create beautiful triangle plots - aka Giant Triangle Confusograms (GTCs)

    Basic usage:
    ----------
    GTC = pyGTC.plotGTC(chains, kwargs)

    Parameters:
    ----------
    chains: 2d array or list of 2d arrays
        Sample points (length x Ndim) or multiple sets of samples points
        Note: If you are using emcee (http://dan.iel.fm/emcee/current/) - and you should! - you need to pass the EnsembleSampler.flatchain object.
    kwargs:
        weights: weights for the sample points. Either 1d array or list of 1d arrays
        chainLabels: one label per chain, supports latex
        paramNames: one label per parameter, supports latex
        truths: show points in parameter space by lines, multiple sets possible
        truthLabels: one label per truth line, supports latex
        truthColors: user-defined colors for the truth lines
        priors: plot Gaussian priors in 1d panels
        plotName: provide filename for direct output
        nConfidenceLevels: in 2d panels, plot nConfidenceLevels confidence levels (1,2, or 3)
        nBins: change number of bins for histograms
        smoothingKernel: size of Gaussian smoothing kernel (in bins)
        figureSize: provide either figure width, or choose from predefined journal setting (recommended)
        panelSpacing: "loose" / "tight" (default)
        paramRanges: provide ranges for subset or all parameters 
        colorsOrder: change the default order of colors for the contours
        do1dPlots: set to False if you don't want the 1d panels
        doOnly1dPlot: only plot ONE 1d histogram. Provide chain(s) of shape (Npoints,1)

    returns:
        fig: matplotlib.figure
            the GTC in all its glory
"""



def create_random_samples(ndim, Npoints):
    means = np.random.rand(ndim)
    cov = .5 - np.random.rand(ndim**2).reshape((ndim,ndim))
    cov = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov = np.dot(cov,cov)
    samples =  np.random.multivariate_normal(means, cov, Npoints)
    return samples


# Create two sets of sample points with 8 parameters and 50000 points
samples1 = create_random_samples(8, 50000)
samples2 = 1+create_random_samples(8, 60000)

# List of parameter names, supports latex
names = ['param name', '$B_\mathrm{\lambda}$', '$C$', '$\\lambda$', 'C', 'D', 'M', '$\\gamma$']

# Labels for the different chains
#chainLabels = "data1 $\lambda$"
chainLabels = ["data1 $\lambda$", "data 2"]

# List of Gaussian curves to plot (to represent priors): mean, width
# List can be shorter than number of parameters
# Empty () or None if no prior to plot
priors = (None, (2, 1), (.5, 2), (), (0, .4), None, None, None)

# List of truth values, to mark best-fit or input values
#truths = (None, .5, None, .1, 0, None, None, 0)
truths = ((4, .5, None), (None, None, .3))
#truths = ((4, .5, None, .1, 0, None, None, 0), (None, None, .3, 1, None, None, None, None))

# Labels for the different truths
#truthLabels = 'the truth'
truthLabels = ( 'the truth', 'alternative truth')

# List of parameter ranges to show, empty () or None to let pyGTC decide
paramRanges = ((-3,5),None,(-2,4))

########## Do the magic
##### Full GTC
# inputarr must either be:
# -list of length nChains, each being an array of shape (Npoints, nDim)
# -single array of shape (Npoints, nDim)
#GTC = pyGTC.plotGTC(chains=[samples1[:,:2]], do1dPlots=True, paramNames=names[:2], truths=truths, priors=priors, chainLabels=chainLabels[0], truthLabels=truthLabels, paramRanges=paramRanges, figureSize='APJ_column')
GTC = pyGTC.plotGTC(chains=[samples1[:,:3],samples2[:,:3]], do1dPlots=True, paramNames=names[:3], truths=truths[:3], priors=priors[:3], chainLabels=chainLabels, truthLabels=truthLabels, paramRanges=paramRanges, figureSize='APJ_column')


##### Only one 1d histogram
# inputarr must be list of length nChains, each being array of shape (Npoints, 1)
#inputarr = [np.array([samples1[:,0]]).T, np.array([samples2[:,0]]).T]
#GTC = pyGTC.plotGTC(chains=inputarr, paramNames=[names[0]],  chainLabels=chainLabels,figureSize='APJ_column', doOnly1dPlot=True)


#plt.show()
plt.savefig('GTC.pdf', bbox_inches='tight')

#plt.clf()


#fig = plt.figure()
#fig = pyGTC.plot1d(fig, 1, [samples1[:,0]], [np.ones(len(samples1))], 30, 1, [('#4c72b0','#7fa5e3','#b2d8ff')], None, None, None, None)
#plt.show()