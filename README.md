# petibmpy

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/mesnardo/petibmpy/raw/master/LICENSE)

Small Python package to perform pre- and post-processing steps for a [PetIBM](https://github.com/barbagroup/PetIBM) run.

## Dependencies

* PetIBM (last tested: `0.4`)
* Python (`3.6.8`)
* H5Py (`2.9.0`)
* lxml (`4.3.1`)
* Matplotlib (`3.0.2`)
* NumPy (`1.16.2`)
* Pyyaml (`3.13`)

## Installation

With Anaconda:

```shell
conda create env --name=py36-petibmpy --file=environment.yaml
conda activate py36-petibmpy
python setup.py develop
```