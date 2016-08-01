pygtc.py
===
Make a publication-ready giant triangle confusogram with just one line of code.


**What is a Giant Triangle Confusogram?**

A Giant Triangle Confusogram (GTC, aka triangle/corner plot) is a way of displaying the results of an MCMC
or similar analysis. The recovered parameter constraints are displayed on a grid in which the diagonal shows
the one-dimensional posteriors and the lower-left triangle shows the pairwise projections.


Installation
---
pygtc depends on ``numpy``, ``matplotlib``, and ``scipy``. Installation is as easy as ``pip install pygtc``.


Example Usage
---
We recommend to familiarize yourself with pygtc's main features by studying our tutorial
`demo.ipynp <https://github.com/SebastianBocquet/pygtc/blob/master/demo.ipynb>`_


Contributions, problems, questions, recommendations
--
Please report any problems on `GitHub <https://github.com/SebastianBocquet/pygtc/issues>`_ where pygtc is being developed.


API
===
.. autofunction :: pygtc.plotGTC
