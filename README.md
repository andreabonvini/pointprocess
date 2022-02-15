[![Documentation Status](https://readthedocs.org/projects/pointprocess/badge/?version=latest)](https://pointprocess.readthedocs.io/en/latest/?badge=latest)
[![GitHub Super-Linter](https://github.com/andreabonvini/pointprocess/workflows/Lint%20Code%20Base/badge.svg)](https://github.com/marketplace/actions/super-linter)
![Build Status](https://img.shields.io/github/workflow/status/andreabonvini/pointprocess/CI?event=push&label=Build&logo=Github-Actions)
[![codecov.io](https://codecov.io/github/andreabonvini/pointprocess/coverage.svg?branch=master)](https://codecov.io/github/andreabonvini/pointprocess?branch=master)

![ppbig](docs/images/ppbig.png)

ALERT: I'm working on building a Python package in order to ease the use the software in an academic environment, write me in case you need a temporary working solution.

This repository contains a `C++` implementation of the `MATLAB` software provided by Riccardo Barbieri and Luca
Citi [here](http://users.neurostat.mit.edu/barbieri/pphrv).

*Refer to the following papers for more details:*

- [*A point-process model of human heartbeat intervals: new definitions of heart rate and heart rate
  variability*](https://pubmed.ncbi.nlm.nih.gov/15374824/)
- [*The time-rescaling theorem and its application to neural spike train data
  analysis*](https://pubmed.ncbi.nlm.nih.gov/11802915/)

# Requirements

The project can be built with `CMake`, the only dependencies are [Eigen](https://eigen.tuxfamily.org)
and [Boost](https://www.boost.org).

# Documentation

The technical and scientific *documentation* for this repository is available [here](https://pointprocess.readthedocs.io/en/latest/). (WIP)
