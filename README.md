# gmx-top4py
## Scope
Python interface to
* Read and write GROMACS-type top-files
* Alter force field parameters 

## Installation
### From source via uv
Clone repository and move into
```bash
git clone git@github.com:graeter-group/gmx-top4py.git
cd gmx-top4py
```
Install repository
a) uv.lock
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
