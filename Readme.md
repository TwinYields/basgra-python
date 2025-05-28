
# Python interface to BASGRA_N model

This is a python interface to BASGRA_N model from: https://github.com/MarcelVanOijen/BASGRA_N . The main aim is to be able to control the model from Python each day for e.g. data assimilation purposes, this is implemented in `src/bglib.f90`.

## Installation note

The newest tested build system works with python 3.10 and numpy 1.26.0. Update to new build system is needed:

```
pip install setuptools==69.5.1 numpy==1.26.0
python setup.py install
```