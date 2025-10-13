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
* [ff_params.py](examples/ff_params.py) - Adapting force field parameter with gmx-top4py

## Installation
### From pipy
Coming soon
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
