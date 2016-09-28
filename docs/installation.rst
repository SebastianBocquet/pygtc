============
Installation
============

Required packages
~~~~~~~~~~~~~~~~~

pygtc is written in Python 2.7 and requires the following packages:

* numpy >= 1.5
* matplotlib >= 1.5 and < 2.0
* scipy (optional)


Downloading and installing
~~~~~~~~~~~~~~~~~~~~~~~~~~

The latest stable version of pygtc is hosted at `PyPi
<http://pypi.python.org/pypi/pygtc/>`_.

If you like pip, you can install with::

  pip install pygtc

Or, alternatively, download and extract the tarball from the above link and then
install from source with::

  cd directory_you_unpacked_ptgtc_into
  python setup.py install


Development
~~~~~~~~~~~

Development happens at `github <https://github.com/SebastianBocquet/pygtc>`_. You can
clone the repo with::

  git clone https://github.com/SebastianBocquet/pygtc

And you can install it in developer mode with pip::

  pip install -e directory_where_you_cloned_pygtc

or from source::

  cd directory_cloned_pygtc_into
  python setup.py develop

Running tests
~~~~~~~~~~~~~

For tests to *for sure* run properly, you'll want to have matplotlib v1.5.3
installed, as they fixed a bug in their ``image_comparison`` decorator. In any
case, you'll need ``nose`` installed to run the tests. There are two ways to do
the test suite. You can use the nosetests utility::

  nosetests directory_cloned_pygtc_into

Or, you can run the tests as a script::

  cd directory_cloned_pygtc_into/pygtc/tests
  python test_plotGTC.py

In either case, there are 25 tests to run, and it should take between 20-30
seconds to run them all. If the first test fails, there is likely something
wrong with matplotlib in general.

Contribution and/or bug reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you like this project and want to contribute, send a message over at the
gitHub repo and get involved. If you find a bug, please open an issue at gitHub.
Better yet, you can try and fix the bug and submit a pull request!
