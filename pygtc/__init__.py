import inspect
import os

__author__ = "Sebastian Bocquet and Faustin W. Carter"
__contributors__ = [
    'Samuel Flender',
    'Ben Hoyle'
]
__all__ = ['plotGTC']
version_filename = os.path.join(os.path.dirname(os.path.abspath(inspect.stack()[0][1])), 'VERSION.txt')
with open(version_filename, 'r') as version_file:
    __version__ = version_file.read().strip()

from .pygtc import plotGTC
