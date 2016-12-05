PyPSD.py
========
By [Sean Anderson](https://github.com/roguephysicist/PyPSD)

usage: python pypsd.py [-h] inputfile binsfile

Requirements: 
* python3
* numpy
* matplotlib
* argparse

A Python script for calculating the particle size distribution (PSD) of any
sample.

The script requires two inputs. The first, a text file with your particle areas
determined with, for example, the [ImageJ](https://fiji.sc) scientific imaging
program. The second input is a text file with the 'bins' used to classify your
particles by diameter. Sample input (`sample.txt`) and bins (`bins.txt`) files
are included.

This script has been tested with [Anaconda Python]
(https://www.continuum.io/downloads) on macOS/OS X, Windows 10, and Linux.
