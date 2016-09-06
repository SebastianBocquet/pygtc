pygtc.py
========
Make a publication-ready giant-triangle-confusogram (GTC) with just one line of code!


**What is a Giant Triangle Confusogram?**

A Giant-Triangle-Confusogram (GTC, aka triangle/corner plot) is a way of
displaying the results of a Monte-Carlo Markov Chain (MCMC) sampling or similar
analysis. The recovered parameter constraints are displayed on a grid in which
the diagonal shows the one-dimensional posteriors (and, optionally, priors) and
the lower-left triangle shows the pairwise projections. You might want to look
at a plot like this if you are fitting model to data and want to see the
parameter covariances along with the priors.

Although several other packages exists to make such a plot, we were unsatisfied
with the amount of extra work required to massage the result into something we
were happy to publish. With ``pygtc``, we hope to take that extra legwork out of
the equation by providing a package that gives a figure that is publication
ready on the first try!

Here's an example of a GTC:

.. image:: _static/demo_files/demo_8_0.png

Contents:
---------
.. toctree::
   :maxdepth: 2

   installation.rst
   demo.rst
   api.rst


Contributions, problems, questions, recommendations
---------------------------------------------------
Please `report any problems on GitHub <https://github.com/SebastianBocquet/pygtc/issues>`_ where pygtc is being developed.
