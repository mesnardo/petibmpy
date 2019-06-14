# petibmpy: Python processing tools for PetIBM

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/mesnardo/petibmpy/raw/master/LICENSE)
[![Travis](https://img.shields.io/travis/org/mesnardo/petibmpy/master.svg?style=flat-square&logo=travis)](https://travis-ci.org/mesnardo/petibmpy)

Small Python package to perform pre- and post-processing steps for [PetIBM](https://github.com/barbagroup/PetIBM).

## Dependencies

* PetIBM (last tested: `0.4.1`)
* Python (`3.6.8`)
* H5Py (`2.9.0`)
* lxml (`4.3.1`)
* Matplotlib (`3.0.2`)
* NumPy (`1.16.2`)
* Pyyaml (`3.13`)
* SciPy (`1.2.1`)

## Installation

With Anaconda:

```shell
conda env create --name=py36-petibmpy --file=environment.yaml
conda activate py36-petibmpy
python setup.py develop
```

## Contact

Please e-mail [Olivier Mesnard](mailto:mesnardo@gwu.edu) if you have any questions, suggestions, or feedback.

To report bugs, please use the GitHub issue tracking system.
We also welcome pull-requests.
