from setuptools import setup
import os

with open(os.path.join("pygtc", "VERSION.txt")) as version_file:
    version_number = version_file.read().strip()

setup(
    name="pyGTC",
    description="Make an awesome giant triangle confusogram (gtc)!",
    version=version_number,
    author="Sebastian Bocquet and Faustin Carter",
    author_email="sebastian.bocquet@gmail.com",
    license="MIT",
    url="http://github.com/sebastianbocquet/pygtc",
    packages=["pygtc"],
    package_data={
        "pygtc": [
            "VERSION.txt",
            "tests/*.py",
            "tests/baseline/*.png",
            "tests/legacy_baseline/*.png",
        ]
    },
    long_description=open("README.rst").read(),
    install_requires=[
        "packaging",
        "numpy>=1.7.1",
        "matplotlib>=2.2",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
