Version 0.2.3:
  * It is now possible to omit a parameter from a chain by using None as a
    placeholder. Demo notebook and documentation updated appropriately.
  * Added changelog to manifest.

Version 0.2.2:
  * No longer fails if you ask for priors and scipy is missing. Instead, it
    issues a warning and ignores the request.
  * Test suite now skips tests that required pandas or scipy if pandas or scipy
    is missing.

Version 0.2.1:
  * Refactored code to make contour choice options more clear.
  * Updated documentation with discussion about contour level choices.
  * Added tests for most Keyword Argument options. Tests require nose.
  * Updated documentation to include notes about testing.
  * Added a requirements.txt file to the package specifying version numbers and
    removed those specifications from setup.py.

Version 0.2.0:
  * Option to turn off tick label rotation
  * Long tick labels are now lined up correctly with their ticks
  * Option to fine tune position of label relative to tick
  * Can choose to use 1d or 2d Gaussian sigma contours
  * scipy is now an optional dependency
  * Can choose to fill contours and 1d histograms independently
  * Can plot the actual data-point density underneath the contours
  * Full control over fonts for labels, ticks, and legend
  * Use built-in matplotlib LaTex renderer for equations in labels so no LaTex
    installation is required locally.
  * Added a new example to the documentation

Version 0.1.1:
  * Prettier truth colors
  * Parameter labels line up even if parameter numbers are very different in size
  * Fixed bug where parameter labels near the edge could overlap


Version 0.1.0:
 * Initial version
