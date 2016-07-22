from setuptools import setup
import pygtc

setup(
    name = 'pyGTC',
    description = 'Make an awesome giant triangle confusogram (gtc)!',
    version = '0.1.0',
    author = 'Sebastian Bocquet and Faustin Carter',
    author_email = 'sebastian.bocquet@gmail.com',
    license = 'MIT',
    url = 'http://github.com/sebastianbocquet/pygtc',
    packages = ['pygtc'],
    long_description = open('README.rst').read(),
    install_requires = ['numpy', 'matplotlib', 'scipy'],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Visualization'
    ]

)
