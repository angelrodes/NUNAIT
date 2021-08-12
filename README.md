# NUNAIT
Supporting information for the manuscript *"The NUNAtak Ice Thinning (NUNAIT) calculator for cosmonuclide elevation profiles."* submitted to *Geosciences*.

The NUNAIT calculator is an easy-to-use MATLAB/Octave tool that constrains parameters describing the geological history of a nunatak from a set of surface exposure ages.

*Ángel Rodés, 2021*

## Graphical description of the model

The NUNAIT model generates vertical profiles of apparent exposure ages based on given parameter values (erosion rates, uplift rate, and current and maximum ice elevations) assuming that the ice surface has changed *parallel* to the climate proxy data.

![NUNAIT_cartoon](https://user-images.githubusercontent.com/53089531/129214127-8459cb71-7675-4239-bc63-5f1075c46a7f.png)

## Fitting types

There are four ways of fitting the model to the data: 
* ```0```: normal fitting.
* ```1```: minimum apparent ages. E.g. for nunataks with inhomogeneous erosion rates.
* ```2```: maximum apparent ages. E.g. for datasets containing boulders with variable inherited signals.
* ```3```: no fitting, used for testing purposes. E.g. generating expected profiles before sampling.

![fit types](https://user-images.githubusercontent.com/53089531/129227177-fae66375-82c1-48d5-9581-db72d63c778e.png)

## Preprint

A description of the calculator is available at:

Rodés, Á. [The NUNAtak Ice Thinning (NUNAIT) Calculator for Cosmonuclide Elevation Profiles.](https://www.preprints.org/manuscript/202107.0411) Preprints 2021, 2021070411 (doi: 10.20944/preprints202107.0411.v1)
