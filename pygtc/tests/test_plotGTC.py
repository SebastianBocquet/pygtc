import numpy as np
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

import pygtc


#Set up some global variables for testing
def _make_random_chain(ndim=4, Npoints=10000, seed = 0):
    np.random.seed(seed)
    means = np.random.rand(ndim)
    cov = .5 - np.random.rand(ndim**2).reshape((ndim,ndim))
    cov = np.triu(cov)
    cov += cov.T - np.diag(cov.diagonal())
    cov = np.dot(cov,cov)
    samples =  np.random.multivariate_normal(means, cov, Npoints)
    return samples

# Create two sets of sample points with 4 parameters and 10000 points
SAMPLES_1 = 2*_make_random_chain(seed = 1)
SAMPLES_2 = 1+_make_random_chain(seed = 2)
SAMPLES_1[:,3]+=1e8
SAMPLES_2[:,3]+=1e8

#Specify kwargs for savefig. We change two things:

#1: bbox tight ensures that the labels don't get cut off.

#2: Set a dpi that won't suck on retina displays and will look fine on anything
#else too. This is only really an issue for raster graphics. Sane people will
#use a vector format, but testing is faster with raster.

SFKWARGS = {'bbox_inches':'tight',
            'dpi':300}


#If this one fails, something is really wrong with matplotlib
@image_comparison(baseline_images=['img'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_img():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(1,1)

#A test for (almost) every keyword argument
@image_comparison(baseline_images=['bare'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_bare():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2])

@image_comparison(baseline_images=['pandas'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_pandas():
    namesNoTex = ['param name', 'B_labmda', 'C', 'lambda']

    if HAS_PANDAS:
        samples1_pd = pd.DataFrame(SAMPLES_1, columns=namesNoTex)
        samples2_pd = pd.DataFrame(SAMPLES_2, columns=namesNoTex)
    else:
        samples1_pd = SAMPLES_1
        samples2_pd = SAMPLES_2

    pygtc.plotGTC(chains=[samples1_pd,samples2_pd])

@image_comparison(baseline_images=['paramNames_noTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_paramNames_noTex():
    namesNoTex = ['param name', 'B_labmda', 'C', 'lambda']
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    paramNames = namesNoTex)

@image_comparison(baseline_images=['paramNames_withTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_paramNames_withTex():
    namesWithTex = ['param name', '$B_\mathrm{\lambda}$', '$Q^a$', '$\\lambda$']
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    paramNames = namesWithTex)

@image_comparison(baseline_images=['chainLabels_noTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_chainLabels_noTex():
    chainLabelsNoTex = ['data1', 'data 2']
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    chainLabels = chainLabelsNoTex)

@image_comparison(baseline_images=['chainLabels_withTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_chainLabels_withTex():
    chainLabelsWithTex = ['data1 $\lambda$', 'data 2']
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    chainLabels = chainLabelsWithTex)

@image_comparison(baseline_images=['truthLabels_noTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_truthLabels_noTex():
    truths = ((4, .5, None, .1),
                (None, None, .3, 1))
    truthLabelsNoTex = ('the truth', 'alternative truth')
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    truths = truths,
                    truthLabels = truthLabelsNoTex)

@image_comparison(baseline_images=['truthLabels_withTex'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_truthLabels_withTex():
    truths = ((4, .5, None, .1),
                (None, None, .3, 1))
    truthLabelsWithTex = ('the truth $f_0$', 'alternative truth $\\lambda$')
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    truths = truths,
                    truthLabels = truthLabelsWithTex)

#TODO: Add a test for truthColors

@image_comparison(baseline_images=['truthLineStyles'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_truthLineStyles():
    truthLineStyles = ['-', '-']
    truths = ((4, .5, None, .1),
                (None, None, .3, 1))
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    truths = truths,
                    truthLineStyles = truthLineStyles)


@image_comparison(baseline_images=['priors'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_priors():
    priors = (None, (2, 1), (.5, 2), ())
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    priors = priors)

#TODO: Think up a good way to test plotName

@image_comparison(baseline_images=['nConfidenceLevels'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_nConfidenceLevels():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    nConfidenceLevels = 3)

@image_comparison(baseline_images=['gaussianConfLevels'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_gaussianConfLevels():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    gaussianConfLevels = True)

@image_comparison(baseline_images=['nBins'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_nBins():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    nBins = 20)

@image_comparison(baseline_images=['smoothingKernel'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_smoothingKernel():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    smoothingKernel = 2)

@image_comparison(baseline_images=['filledPlots'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_filledPlots():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    filledPlots = False)

@image_comparison(baseline_images=['plotDensity'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_plotDensity():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    plotDensity = True)

@image_comparison(baseline_images=['figureSize'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_figureSize():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    figureSize = 'APJ_page')

@image_comparison(baseline_images=['panelSpacing'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_panelSpacing():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    panelSpacing = 'loose')

#TODO: Add a test for legendMarker

#TODO: Add a test for paramRanges

@image_comparison(baseline_images=['labelRotation'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_labelRotation():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    labelRotation = (False, False))

@image_comparison(baseline_images=['tickShifts'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_tickShifts():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    tickShifts = (0.2, 0.2))

@image_comparison(baseline_images=['colorsOrder'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_colorsOrder():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    colorsOrder = ['purples', 'yellows'])

@image_comparison(baseline_images=['do1dPlots'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_do1dPlots():
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    do1dPlots = False)

@image_comparison(baseline_images=['doOnly1dPlot'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_doOnly1dPlot():
    input_chains = [np.array([SAMPLES_1[:,0]]).T, np.array([SAMPLES_2[:,0]]).T]
    pygtc.plotGTC(chains=input_chains,
                    doOnly1dPlot = True)

@image_comparison(baseline_images=['mathTextFontSet'], extensions=['png'], savefig_kwarg=SFKWARGS)
def test_GTC_mathTextFontSet():
    namesWithTex = ['param name', '$B_\mathrm{\lambda}$', '$Q^a$', '$\\lambda$']
    pygtc.plotGTC(chains=[SAMPLES_1,SAMPLES_2],
                    paramNames = namesWithTex,
                    mathTextFontSet = None)

#TODO: Could add a few more tests to deal with label font customization...

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=['-s', '--with-doctest'])
