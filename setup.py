from setuptools import setup
import os
import pygtc

version_file = open(os.path.join('.', 'VERSION.txt'))
version_number = version_file.read().strip()
version_file.close()

setup(
    name = 'pygtc',
    description = 'Make an awesome giant triangle confusogram (gtc)!',
    version = version_number,
    author = 'Sebastian Bocquet and Faustin Carter',
    author_email = 'sebastian.bocquet@gmail.com',
    license = 'MIT',
    url = 'http://github.com/sebastianbocquet/pygtc',
    packages = ['pygtc'],
    long_description = open('README.rst').read(),
    install_requires = [
        'numpy>=1.5',
        'matplotlib>=1.5, <2.0'
    ],
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Visualization'
    ]

)
