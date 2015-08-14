automatic-equilibration-detection
=================================

Automatic detection of equilibrated regions of molecular simulations

A preprint manuscript of this work is available [on bioRxiv](http://biorxiv.org/content/early/2015/07/04/021659).

## Manifest
* `examples/` - code to generate simulation data and figures used in the manuscript
* `manuscript/` - the manuscript source and figures
* `references/` - some references cited in the paper

## Reproducing the results

* Check out a copy of this repository, or [download a release](https://github.com/choderalab/automatic-equilibration-detection/releases) corresponding to the version of the paper of interest.
* Go to `examples/liquid-argon-density/` and check the `README.md` for instructions on running the `./reproduce.sh` script.

## Automatically computing equilibration times

To detect equilibration start times in your own code, you can use the [`detectEquilibration()` method](http://pymbar.readthedocs.org/en/latest/timeseries.html#pymbar.timeseries.detectEquilibration) in the [`timeseries` module](http://pymbar.readthedocs.org/en/latest/timeseries.html) of [`pymbar`](http://pymbar.readthedocs.org/en/latest/).

You can install pymbar via the [`conda` package manager](http://conda.pydata.org):
```bash
conda install -c omnia pymbar
```
and then use the [`detectEquilibration()` method](http://pymbar.readthedocs.org/en/latest/timeseries.html#pymbar.timeseries.detectEquilibration):
```python
from pymbar import timeseries
[t0, g, Neff_max] = timeseries.detectEquilibration(A_t)
A_t_equlibrated = A_t[t0:]
```
where `A_t` is the observable timeseries and `t_0` is the beginning of the equilibrated region.

The [`equilibration.py` module](https://github.com/choderalab/automatic-equilibration-detection/blob/master/examples/liquid-argon-density/equilibration.py) provided in this repo also provides a stand-alone implementation of this feature with a broader choice of methods for computing integrated autocorrelation times.
