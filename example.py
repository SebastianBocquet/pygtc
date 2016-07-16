from matplotlib import pyplot as plt
import numpy as np
import pyGTC


# List of parameter names, supports latex
names = ['param name', '$B_\mathrm{\lambda}$', '$C$', '$\\lambda$', 'C', 'D', 'M', '$\\gamma$']


# List of priors: mean, width
# List can be shorter than number of parameters
priors = [[2, 1], [.5, 2], [], [0, .4], [], []]


# List of truth value
# NOT a python array because of different lengths
#truths = [4, .5, None, .1, None, None, None, None, 0]
truths = [[4, .5, None, .1, None, None, None, None, 0], [None, None, .3, 1]]


# Create two sets of sample points with 8 parameters
ndim = 8

means = np.random.rand(ndim)
cov = .5 - np.random.rand(ndim**2).reshape((ndim,ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov,cov)
samples1 = np.random.multivariate_normal(means, cov, 50000)

means = np.random.rand(ndim) + 1
cov = .5 - np.random.rand(ndim**2).reshape((ndim,ndim))
cov = np.triu(cov)
cov += cov.T - np.diag(cov.diagonal())
cov = np.dot(cov,cov)
samples2 = np.random.multivariate_normal(means, cov, 50000)


# Labels for the different chains
ChainLabels = ["data1 $\lambda$", "data 2"]

# Do the magic
# Unused arguments: NConfidenceLevels, figuresize
GTC = pyGTC.plotGTC(chains=[samples1,samples2], ParamNames=names, truths=truths, priors=priors, ChainLabels=ChainLabels)

#plt.show()
plt.savefig('GTC.pdf', bbox_inches='tight')
