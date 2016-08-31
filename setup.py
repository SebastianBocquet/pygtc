from setuptools import setup
import pygtc

setup(
    name = 'pygtc',
    description = 'Make an awesome giant triangle confusogram (gtc)!',
    version = '0.1.1',
    author = 'Sebastian Bocquet and Faustin Carter',
    author_email = 'sebastian.bocquet@gmail.com',
    license = 'MIT',
    url = 'http://github.com/sebastianbocquet/pygtc',
    packages = ['pygtc'],
    long_description = open('README.rst').read(),
    install_requires = [
        'numpy>=1.5',
        'matplotlib>=1.5',
        'scipy>=0.14'
    ],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Visualization'
    ]

)
