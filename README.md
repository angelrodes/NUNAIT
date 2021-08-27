# NUNAIT

![NUNAIT_logo](https://user-images.githubusercontent.com/53089531/129394469-bd579e0d-14f3-4eeb-98ef-528a92301586.png)

>Supporting information for the manuscript *"The NUNAtak Ice Thinning (NUNAIT) calculator for cosmonuclide elevation profiles."* [published in *Geosciences*](https://www.mdpi.com/2076-3263/11/9/362).

The NUNAIT calculator is an easy-to-use MATLAB/Octave tool that constrains parameters describing the geological history of a nunatak from a set of surface exposure ages.

*Ángel Rodés, 2021*\
[angelrodes.com](http://www.angelrodes.com/)

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

## Cite

A description of the calculator is available at:

**Rodés, Á.** (2021) "The NUNAtak Ice Thinning (NUNAIT) Calculator for Cosmonuclide Elevation Profiles" *Geosciences* 11, no. 9: 362. doi:[10.3390/geosciences11090362](https://doi.org/10.3390/geosciences11090362 )
