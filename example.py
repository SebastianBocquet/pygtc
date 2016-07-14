import numpy as np
from pyGTC import pyGTC

names = "A_\mathrm{X} B_\mathrm{X} C_\mathrm{X} D_\mathrm{X} A_\mathrm{SZ} B_\mathrm{SZ} C_\mathrm{SZ} D_\mathrm{SZ} \rho_\mathrm{SZ-X} \mathsf{\Omega}_\mathrm{m} \ln(10^{10}A_s) h \mathsf{\Omega}_\mathrm{b} n_s \sigma_8".split()

chain = np.loadtxt("samples.txt")
pyGTC(chains=[chain], param_names=names)
