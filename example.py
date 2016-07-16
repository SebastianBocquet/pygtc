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


# List of parameter names, supports latex
names = ['param name', '$B_\mathrm{\lambda}$', '$C$', '$\\lambda$', 'C', 'D', 'M', '$\\gamma$']


# List of priors: mean, width
# List can be shorter than number of parameters
priors = [[2, 1], [.5, 2], [], [0, .4], [], []]


# List of truth value
# NOT a python array because of different lengths
#truths = [4, .5, None, .1, None, None, None, None, 0]
truths = [[4, .5, None, .1, None, None, None, None, 0], [None, None, .3, 1]]


# Create two sets of sample points with 8 parameters and 50000 points
samples1 = create_random_samples(8, 50000)
samples2 = 1+create_random_samples(8, 50000)


# Labels for the different chains
ChainLabels = ["data1 $\lambda$", "data 2"]


# Labels for the different truths
TruthLabels = ['the truth', 'alternative truth']

# Do the magic
# Unused arguments: NConfidenceLevels, figuresize
GTC = pyGTC.plotGTC(chains=[samples1,samples2], ParamNames=names, truths=truths, priors=priors, ChainLabels=ChainLabels, TruthLabels=TruthLabels)


#plt.show()
plt.savefig('GTC.pdf', bbox_inches='tight')
