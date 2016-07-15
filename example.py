from matplotlib import pyplot as plt
import numpy as np
import pyGTC

# List of parameter names, latex-enabled
names = ['normalization', '$B_\mathrm{X}$', '$C_\mathrm{X}$', '$D_\mathrm{X}$', '$A_\mathrm{SZ}$', '$B_\mathrm{SZ}$', '$C_\mathrm{SZ}$', '$D_\mathrm{SZ}$', '$\\rho_\mathrm{SZ-X}$', '$\mathsf{\Omega}_\mathrm{m}$']

# List of priors: mean, width
# List can be shorter than number of parameters
priors = [[6,38, .61], [.57, .03], [], [.12, .08], [], []]

# List of truth value
truths = [6, .57, None, .12, None, None, None, None, 0]

# Load chain
chain1 = np.loadtxt("samples.txt")
chain2 = np.loadtxt("samples2.txt")

# Labels for the different chains
chainlabels = ["data1 $\lambda$", "chain 2"]

# Do the magic
# Unused arguments: confidencelevels, figuresize
GTC = pyGTC.plotGTC(chains=[chain1, chain2], param_names=names, truths=truths, priors=priors, chainlabels=chainlabels)

#plt.show()
plt.savefig('GTC.pdf', bbox_inches='tight')
