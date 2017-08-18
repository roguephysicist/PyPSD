PyPSD
========

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

PyPSD is a Python program for calculating the particle size distribution (PSD) of any sample. Measuring particle sizes and calculating their distribution in a sample is crucial for understanding diverse physical and chemical processes. Particle sizes can be determined with [ImageJ](https://fiji.sc), a free and open-source scientific imaging program (or any other software with this functionality). That data is then used by PyPSD in order to produce the relevant distributions.

![](http://i.imgur.com/eM8dUPk.png?1)


Installation and use
------------------------------------

PyPSD has been tested with Python 3 and [Anaconda Python 4+](https://www.continuum.io/downloads) on macOS, Linux, and Windows. It should work on any system with the required Python packages installed.

Python requirements:
`numpy`, `matplotlib`, `argparse`, `os`

Usage:
`python pypsd.py [--help] -b BINSFILE -i INPUTFILE [-o OUTPUTDIR]`

For maximum compatibility, PyPSD is designed to be executed directly from the command line which is available on all operating systems. The `--help` option will display a useful help message. Use the `-b` and `-i` flags to specify your bin and input files. The optional `-o` argument allows you to specify a directory to output your results.

Your input data should be in a text (`.txt`) file with a single column containing the measured particle areas. Your bins should also be in a text file with a single column of bin limits. These are used to classify your particles by diameter. Your bins should obviously be in units compatible with the areas listed in your input file. Sample input (`sample_input.txt`) and bin (`bins.txt`) files are included.


License
------------------------------------

Copyright 2017 [Sean M. Anderson](mailto:sma@cio.mx) and Liliana Villafaña-López.

PyPSD is free software made available under the BSD-3-Clause License. For details please see the LICENSE file.
