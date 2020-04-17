# Make sure we always use the same backend for image comparison tests
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest
import warnings

import pygtc

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

HAS_SCIPY = pygtc.haveScipy

MPLVER = int(matplotlib.__version__.split('.')[0])
if MPLVER < 2:
    warnings.warn("Several tests are known to fail under matplotlib versions " +
                  "less than 2.0. The plots should still look good!",
                  UserWarning)

image_comp = pytest.mark.mpl_image_compare


# Set up some global variables for testing
def _make_random_chain(ndim=4, Npoints=10000, seed=0):
    np.random.seed(seed)
    means = np.random.rand(ndim)
    cov = .5 - np.random.rand(ndim**2).reshape((ndim, ndim))
    cov = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov = np.dot(cov, cov)
    samples = np.random.multivariate_normal(means, cov, Npoints)
    return samples


# Create two sets of sample points with 4 parameters and 10000 points
SAMPLES_1 = 2*_make_random_chain(seed=1)
SAMPLES_2 = 1+_make_random_chain(seed=2)
SAMPLES_1[:, 3] += 1e8
SAMPLES_2[:, 3] += 1e8

# Specify kwargs for savefig. We change two things:

# 1: bbox tight ensures that the labels don't get cut off.

# 2: Set a dpi that won't suck on retina displays and will look fine on anything
# else too. This is only really an issue for raster graphics. Sane people will
# use a vector format, but testing is faster with raster.

SFKWARGS = {'bbox_inches': 'tight',
            'dpi': 300}


# If this one fails, something is really wrong with matplotlib
@image_comp(filename='img.png', savefig_kwargs=SFKWARGS)
def test_img():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(1, 1)
    return fig


# A test for (almost) every keyword argument
@image_comp(filename='bare.png', savefig_kwargs=SFKWARGS)
def test_GTC_bare():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         smoothingKernel=0)


@image_comp(filename='pandas.png', savefig_kwargs=SFKWARGS)
def test_GTC_pandas():
    namesNoTex = ['param name', 'B_labmda', 'C', 'lambda']

    if HAS_PANDAS:
        samples1_pd = pd.DataFrame(SAMPLES_1, columns=namesNoTex)
        samples2_pd = pd.DataFrame(SAMPLES_2, columns=namesNoTex)
    else:
        pytest.skip("Can't test pandas auto-name without pandas.")

    return pygtc.plotGTC(chains=[samples1_pd, samples2_pd],
                         smoothingKernel=0)


@image_comp(filename='paramNames_noTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_paramNames_noTex():
    namesNoTex = ['param name', 'B_labmda', 'C', 'lambda']
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         paramNames=namesNoTex,
                         smoothingKernel=0)


@image_comp(filename='paramNames_withTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_paramNames_withTex():
    namesWithTex = ['param name', '$B_\\mathrm{\\lambda}$',
                    '$Q^a$', '$\\lambda$']
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         paramNames=namesWithTex,
                         smoothingKernel=0)


@image_comp(filename='chainLabels_noTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_chainLabels_noTex():
    chainLabelsNoTex = ['data1', 'data 2']
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         chainLabels=chainLabelsNoTex,
                         smoothingKernel=0)


@image_comp(filename='chainLabels_withTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_chainLabels_withTex():
    chainLabelsWithTex = ['data1 $\\lambda$', 'data 2']
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         chainLabels=chainLabelsWithTex,
                         smoothingKernel=0)


@image_comp(filename='truthLabels_noTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_truthLabels_noTex():
    truths = ((4, .5, None, .1),
              (None, None, .3, 1))
    truthLabelsNoTex = ('the truth', 'alternative truth')
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         truths=truths,
                         truthLabels=truthLabelsNoTex,
                         smoothingKernel=0)


@image_comp(filename='truthLabels_withTex.png', savefig_kwargs=SFKWARGS)
def test_GTC_truthLabels_withTex():
    truths = ((4, .5, None, .1),
              (None, None, .3, 1))
    truthLabelsWithTex = ('the truth $f_0$', 'alternative truth $\\lambda$')
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         truths=truths,
                         truthLabels=truthLabelsWithTex,
                         smoothingKernel=0)

# TODO: Add a test for truthColors


@image_comp(filename='truthLineStyles.png', savefig_kwargs=SFKWARGS)
def test_GTC_truthLineStyles():
    truthLineStyles = ['-', '-']
    truths = ((4, .5, None, .1),
              (None, None, .3, 1))
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         truths=truths,
                         truthLineStyles=truthLineStyles,
                         smoothingKernel=0)


@image_comp(filename='priors.png', tol=5e-3, savefig_kwargs=SFKWARGS)
def test_GTC_priors():
    if not HAS_SCIPY:
        pytest.skip("Can't test priors without scipy installed.")

    priors = (None, (2, 1), (.5, 2), ())
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         priors=priors,
                         smoothingKernel=0)


# TODO: Think up a good way to test plotName

@image_comp(filename='nContourLevels.png', savefig_kwargs=SFKWARGS)
def test_GTC_nContourLevels():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         nContourLevels=3,
                         smoothingKernel=0)


@image_comp(filename='sigmaContourLevels.png', savefig_kwargs=SFKWARGS)
def test_GTC_sigmaContourLevels():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         sigmaContourLevels=True,
                         smoothingKernel=0)


@image_comp(filename='nBins.png', savefig_kwargs=SFKWARGS)
def test_GTC_nBins():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         nBins=20,
                         smoothingKernel=0)


@image_comp(filename='smoothingKernel.png', savefig_kwargs=SFKWARGS)
def test_GTC_smoothingKernel():
    if not HAS_SCIPY:
        pytest.skip("Can't test smoothing without scipy.")

    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         smoothingKernel=2)


@image_comp(filename='filledPlots.png', savefig_kwargs=SFKWARGS)
def test_GTC_filledPlots():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         filledPlots=False,
                         smoothingKernel=0)


@image_comp(filename='plotDensity.png', savefig_kwargs=SFKWARGS)
def test_GTC_plotDensity():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         plotDensity=True,
                         smoothingKernel=0)


@image_comp(filename='figureSize.png', savefig_kwargs=SFKWARGS)
def test_GTC_figureSize():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         figureSize='APJ_page',
                         smoothingKernel=0)


@image_comp(filename='panelSpacing.png', savefig_kwargs=SFKWARGS)
def test_GTC_panelSpacing():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         panelSpacing='loose',
                         smoothingKernel=0)


# TODO: Add a test for legendMarker

# TODO: Add a test for paramRanges

@image_comp(filename='labelRotation.png', savefig_kwargs=SFKWARGS)
def test_GTC_labelRotation():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         labelRotation=(False, False),
                         smoothingKernel=0)


@image_comp(filename='tickShifts.png', savefig_kwargs=SFKWARGS)
def test_GTC_tickShifts():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         tickShifts=(0.2, 0.2),
                         smoothingKernel=0)


@image_comp(filename='colorsOrder.png', savefig_kwargs=SFKWARGS)
def test_GTC_colorsOrder():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         colorsOrder=['purples', 'yellows'],
                         smoothingKernel=0)


@image_comp(filename='do1dPlots.png', savefig_kwargs=SFKWARGS)
def test_GTC_do1dPlots():
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         do1dPlots=False,
                         smoothingKernel=0)


@image_comp(filename='doOnly1dPlot.png', savefig_kwargs=SFKWARGS)
def test_GTC_doOnly1dPlot():
    input_chains = [np.array([SAMPLES_1[:, 0]]).T,
                    np.array([SAMPLES_2[:, 0]]).T]
    return pygtc.plotGTC(chains=input_chains,
                         doOnly1dPlot=True,
                         smoothingKernel=0)


@image_comp(filename='mathTextFontSet.png', savefig_kwargs=SFKWARGS)
def test_GTC_mathTextFontSet():
    namesWithTex = ['param name', '$B_\\mathrm{\\lambda}$',
                    '$Q^a$', '$\\lambda$']
    return pygtc.plotGTC(chains=[SAMPLES_1, SAMPLES_2],
                         paramNames=namesWithTex,
                         mathTextFontSet=None,
                         smoothingKernel=0)

# TODO: Could add a few more tests to deal with label font customization...
