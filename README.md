# gmx-top4py 

<p align="left">
<a href="https://github.com/psf/black"><img alt="Code style: black" src="https://img.shields.io/badge/code%20style-black-000000.svg"></a>
</p>

GROMACS topology files for python

## Description
The package gmx-top4py provides an python interface to
* read and write toplogy and force field information from GROMACS-type top-files
* alter force field parameters 

## Tutorials
The following tutorials are avaiable as Google Colab notebooks and hence do not require a local installation:
* [Basic tutorial](https://colab.research.google.com/drive/19zjIw2H5O6InLDQWLXEvrwrZ5w0AL7i1?usp=sharing) - Basic accessing of topology information with gmx-top4py
* [Advanced tutorial](https://colab.research.google.com/drive/1uWgji2O1N4PnxnYpGc1yV68Qzm-PfVX1?usp=sharing) - Advanded accessing of topology information with gmx-top4py
* [Force field parameter tutorial](https://colab.research.google.com/drive/1YOPPkCDudqSwd2QV0P7qEcd-1bogrvtY?usp=sharing) - Adapting force field parameters with gmx-top4py

## Installation
### From PyPi
```bash
pip install gmx-top4py
```

### From TestPyPi
```bash
pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple gmx-top4py
```
### From source via uv
Clone repository and move into
```bash
git clone git@github.com:graeter-group/gmx-top4py.git
cd gmx-top4py
```
Install repository
```bash
uv sync
```
Activate virtual environment
```bash
source .venv/bin/activate
```
Verify install by running the tests
```bash
pytest tests
```
