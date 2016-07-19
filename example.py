from matplotlib import pyplot as plt
import numpy as np
import pyGTC


def create_random_samples(ndim, Npoints):
    """Return samples drawn from random multivariate Gaussian.

    Parameters:
    -----------
    ndim: int
        number of dimensions (parameters)
    Npoints: int
        number of random points to draw

    Returns:
    --------
    samples: array
        samples of points drawn from the distribution
    """
    means = np.random.rand(ndim)
    cov = .5 - np.random.rand(ndim**2).reshape((ndim,ndim))
    cov = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov = np.dot(cov,cov)
    samples =  np.random.multivariate_normal(means, cov, Npoints)
    return samples


# Create two sets of sample points with 8 parameters and 50000 points
samples1 = create_random_samples(8, 50000)
samples2 = 1+create_random_samples(8, 50000)

# List of parameter names, supports latex
names = ['param name', '$B_\mathrm{\lambda}$', '$C$', '$\\lambda$', 'C', 'D', 'M', '$\\gamma$']

# Labels for the different chains
#chainLabels = "data1 $\lambda$"
chainLabels = ["data1 $\lambda$", "data 2"]

# List of Gaussian curves to plot (to represent priors): mean, width
# List can be shorter than number of parameters
# Empty () or None if no prior to plot
priors = ((2, 1), (.5, 2), (), (0, .4), None, ())

# List of truth values, to mark best-fit or input values
# NOT a python array because of different lengths
#truths = (4, .5, None, .1, 0, None, None, 0)
truths = ((4, .5, None, .1, 0, None, None, 0), (None, None, .3, 1))

# List of parameter ranges to show, empty () or None to let pyGTC decide
paramRanges = ((-3,5),None,(-2,4),())

# Labels for the different truths
#truthLabels = 'the truth'
truthLabels = ( 'the truth', 'alternative truth')

# New order of colors
colorsOrder = ['blues', 'yellows']

########## Do the magic
##### Full GTC
# inputarr must either be:
# -list of length nChains, each being an array of shape (Npoints, nDim)
# -single array of shape (Npoints, nDim)
#inputarr=[samples1,samples2]
#GTC = pyGTC.plotGTC(chains=inputarr, paramNames=names, truths=truths, priors=priors, chainLabels=chainLabels, truthLabels=truthLabels, paramRanges=paramRanges, colorsOrder=colorsOrder)

##### Only 2 parameters
# inputarr must be list of length nChains, each being an array of shape (Npoints, 2)
inputarr = [np.array([samples1[:,0],samples1[:,1]]).T]
# Covariance and 1d histograms
#GTC = pyGTC.plotGTC(chains=inputarr, paramNames=names[:2], chainLabels='data', figureSize='APJ_column')
# Only covariance
#GTC = pyGTC.plotGTC(chains=inputarr, paramNames=names[:2], chainLabels='data', figureSize='APJ_column', do1dplots=False)

##### Only one 1d histogram
# inputarr must be list of length nChains, each being array of shape (Npoints, 1)
inputarr = [np.array([samples1[:,0]]).T, np.array([samples2[:,0]]).T]
print np.shape(inputarr)
GTC = pyGTC.plotGTC(chains=inputarr, paramNames=[names[0]],  chainLabels=chainLabels,figureSize='APJ_column', doonly1dplot=True)


#plt.show()
plt.savefig('GTC.pdf', bbox_inches='tight')

#plt.clf()


#fig = plt.figure()
#fig = pyGTC.plot1d(fig, 1, [samples1[:,0]], [np.ones(len(samples1))], 30, 1, [('#4c72b0','#7fa5e3','#b2d8ff')], None, None, None, None)
#plt.show()