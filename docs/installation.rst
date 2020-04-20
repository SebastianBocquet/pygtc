============
Installation
============

Required packages
~~~~~~~~~~~~~~~~~

pygtc is compatible with Python 2.7 and 3.6 or greater and requires the following packages:

* numpy >= 1.5
* matplotlib >= 1.5.3 (preferably >= 2.0)
* scipy (optional)
* pandas (optional -- provides easy auto-naming of chains)
* pytest (optional -- only needed for running unit tests)
* pytest-mpl (optional -- only needed for running unit tests)


Downloading and installing
~~~~~~~~~~~~~~~~~~~~~~~~~~

The latest stable version of pygtc is hosted at `PyPi
<http://pypi.python.org/pypi/pygtc/>`_.

If you like pip, you can install with::

  pip install pygtc

Or, alternatively, download and extract the tarball from the above link and then
install from source with::

  cd /path/to/pygtc
  python setup.py install


Development
~~~~~~~~~~~

Development happens at `github <https://github.com/SebastianBocquet/pygtc>`_. You can
clone the repo with::

  git clone https://github.com/SebastianBocquet/pygtc

And you can install it in developer mode with pip::

  pip install -e /path/to/pygtc

or from source::

  cd /path/to/pygtc
  python setup.py develop

Running tests
~~~~~~~~~~~~~
For tests to *for sure* run properly, you'll want to have matplotlib >= 3
installed and be working on python >= 3.6. Additionally you'll need to ensure
that you are using freetype 2.9.1. You'll need ``pytest`` and ``pytest-mpl``
installed to run the tests, although pygtc functions fine without them. You
also should have the Arial font installed, as that is pygtc's default font and
tests will "fail" if matplotlib falls back on Deja Vu Sans (even though the
images produced might look perfectly fine). Test base images were produced on
Windows 10 using the ``Agg`` backend and if you are on another system there is
no guarantee that you will get a pixel-perfect copy of what the Agg backend
produces. However, the images produced by the tests should still look great!
They are saved to a folder called ``result_images`` in whatever directory you
ran the tests from.

To run the test suite, use the pytest utility from within the pygtc package
directory::

  cd /path/to/pygtc
  pytest --mpl --mpl-results-path=results

If you are running tests on OSX, you may have better results comparing against
the legacy_baseline images like so::

  pytest --mpl --mpl-baseline-path=pygtc/tests/legacy_baseline --mpl-results-path=result_images

There are 25 tests to run, and it should take between 20-30 seconds to run them
all. If the first test fails, there may be something wrong with your matplotlib
install in general (or maybe something weird in your rcParams). If you are
missing pandas or scipy, a few tests will be skipped. If you get a ton of
errors make sure you read the first paragraph in this section and you have all
the prerequisites installed. If matplotlib can't find Arial and you recently
installed it, delete your matplotlib font cache and try again. If you are
running on matplotlib 2.x, then expect the chain label tests to fail because
mpl 2 and 3 handle margins slightly differently. If errors persist, let us know
at GitHub.

Contribution and/or bug reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you like this project and want to contribute, send a message over at the
gitHub repo and get involved. If you find a bug, please open an issue at gitHub.
Better yet, you can try and fix the bug and submit a pull request!
