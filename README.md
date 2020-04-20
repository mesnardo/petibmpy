# PetibmPy: Python processing tools for PetIBM

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/mesnardo/petibmpy/raw/master/LICENSE)
[![Travis](https://img.shields.io/travis/mesnardo/petibmpy/master.svg?style=flat-square&logo=travis)](https://travis-ci.org/mesnardo/petibmpy)
[![Coverage Status](https://coveralls.io/repos/github/mesnardo/petibmpy/badge.svg?branch=master)](https://coveralls.io/github/mesnardo/petibmpy?branch=master)

Small Python package to perform pre- and post-processing steps for [PetIBM](https://github.com/barbagroup/PetIBM).

## Dependencies

* PetIBM (last tested: `0.5.1`)
* Python (`3.6.10`)
* H5Py (`2.10.0`)
* lxml (`4.5.0`)
* Matplotlib (`3.1.1`)
* NumPy (`1.18.1`)
* Pyyaml (`5.3`)
* SciPy (`1.3.2`)

## Installation

With Anaconda:

```shell
conda env create --name=py36-petibmpy --file=environment.yaml
conda activate py36-petibmpy
python setup.py install
```

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) if you have any questions, suggestions, or feedback.

To report bugs, please use the GitHub issue tracking system.
We also welcome pull-requests.
