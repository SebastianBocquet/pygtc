---
  title: 'pyGTC: beautiful Giant Triangle Confusograms aka. likelihood surface plots'
  tags:
    - Python
    - Parameter space exploration
    - data analysis
    - visualization
  authors:
   - name: Sebastian Bocquet
     orcid: 0000-0002-4900-805X
     affiliation: Argonne National Laboratory
   - name: Faustin W. Carter
     orcid: 0000-0002-1236-3517
     affiliation: Argonne National Laboratory
  date: 27 July 2016
  bibliography: paper.bib
---

# Summary

Multi-dimensional model parameter spaces are most efficiently sampled using Markov Chain Monte-Carlo methods and similar algorithms (see, for instance, `emcee` [@emcee]). A common way of visualizing the output of such a code is to show the one-dimensional parameter posteriors as well as the pairwise projections of the likelihood surface. The individual panels are most intuitively presented on a grid in which the diagonal shows the 1-dimensional posteriors and the lower-left half displays the pairwise covariances. Due to the triangular appearance, such plots are typically referred to as "triangle" or "corner" plots. If the parameter space is large, the resulting plot can be visually overwhelming; we refer to such a figure as a Giant Triangle Confusogram (GTC).

Several packages exist to produce these plots and range in functionality from the lightweight `corner` [@corner], which can display a single set of likelihood surfaces from a single analysis, to the feature-heavy `getdist` [@getdist], which allows arbitrary large data sets and contains nearly every conceivable option and customization (at the expense of a steep learning curve). `pyGTC` aims to fill the parameter space in between by producing beautiful, publication-ready figures with as little as a single line of code, while allowing for an arbitrary number of sets of likelihood surfaces (a crucial feature in our field of cosmology). Various user-friendly options allow for further adjustments and ensure that details such as font sizes are chosen to match the layout of several astrophysical journals.

`pyGTC` [@pygtc] uses `matplotlib` [@matplotlib] to produce an aesthetically appealing GTC, and is designed to accept the standard output of the `emcee` [@emcee] package.

Figure 1 shows an example figure (comparing three sets of data from different analyses) that was produced with the pyGTC package and formatted for the Astrophysical Journal, all with the following line of code:

```python
GTC = pygtc.plotGTC(chains=[dat1,dat2,dat3],
                    chainLabels=["data1", "data2", "data3"],
                    paramNames=['foo', 'bar', '$C$', '$\\nu$'],
                    truths=(None, .5, None, .1),
                    truthLabels='ref',
                    figureSize='APJ_page')
```

# References
