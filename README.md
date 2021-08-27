# NUNAIT

![NUNAIT_logo](https://user-images.githubusercontent.com/53089531/129394469-bd579e0d-14f3-4eeb-98ef-528a92301586.png)

>Supplementary material of \
>Rodés, Á. (2021) The NUNAtak Ice Thinning (NUNAIT) Calculator for Cosmonuclide Elevation Profiles. *Geosciences* 11, no. 9: 362. doi:[10.3390/geosciences11090362](https://doi.org/10.3390/geosciences11090362 )

The NUNAIT calculator is an easy-to-use MATLAB/Octave tool that constrains parameters describing the geological history of a nunatak from a set of surface exposure ages.

*Ángel Rodés, 2021*\
[**angelrodes.com**](http://www.angelrodes.com/)

## Graphical description of the model

The NUNAIT model generates vertical profiles of apparent exposure ages based on given parameter values (erosion rates, uplift rate, and current and maximum ice elevations) assuming that the ice surface has changed *parallel* to the climate proxy data.

![NUNAIT_cartoon](https://user-images.githubusercontent.com/53089531/129214127-8459cb71-7675-4239-bc63-5f1075c46a7f.png)

## Input data

Save your data in comma separated files ```.csv```, one for each nunatak or group of nunataks. Some examples of input files are included in the folder "Examples". The input file contains the following headers (first line, ignored in NUNAIT): 

* ```name```: Sample name without spaces or symbols.
* ```lat```: Latitude used to calculate the muon contributions ({decimal degrees}).
* ```site_elv```: Elevation of the sample above sea level (m).
* ```isotope```: Mass of the cosmogenic isotope. Currently accepting 3, 10, 14, 21, 26, and 36 for 3He, 10Be, 14C, 21Ne, 26Al, and 36Cl, respectively.
* ```base_level```: Current elevation of the glacier surface above sea level at the sampling site (m). This is used to calculate the ice position through time.
* ```apparent_years```: Apparent surface exposure age calculated with any cosmogenic calculator, any scaling scheme, and any production rate reference.
* ```dapparent_years```: External uncertainty of the previous age.

![image](https://user-images.githubusercontent.com/53089531/131122210-d86e1034-834f-433a-b599-4ffc4081a194.png)

## Fitting parameters

When running ```START.m``` in MATLAB or Octave (e.g. ```octave --persist START.m```) and choosing ```Run simulation```, the program asks the user to set maximum and minimum values for the parameters to be fitted:

![image](https://user-images.githubusercontent.com/53089531/131124678-7670b0f6-deef-4c93-9729-c8ee10168486.png)


## Fitting types

There are four ways of fitting the model to the data: 
* ```0```: normal fitting.
* ```1```: minimum apparent ages. E.g. for nunataks with inhomogeneous erosion rates.
* ```2```: maximum apparent ages. E.g. for datasets containing boulders with variable inherited signals.
* ```3```: no fitting, used for testing purposes. E.g. generating expected profiles before sampling.

![fit types](https://user-images.githubusercontent.com/53089531/129227177-fae66375-82c1-48d5-9581-db72d63c778e.png)

## Output

After running the models, the NUNAIT calculator produces text and graphical outputs, including probability plots for each parameter, glacier elevation model, and apparent age profiles:

![image](https://user-images.githubusercontent.com/53089531/131139758-52c595c0-8483-45e5-862f-532c752dd963.png)

The text output in the command window is tabulated, so it can be easily copied to any spreadsheet.


## Cite

A description of the calculator is available at:

**Rodés, Á.** (2021) "The NUNAtak Ice Thinning (NUNAIT) Calculator for Cosmonuclide Elevation Profiles" *Geosciences* 11, no. 9: 362. doi:[10.3390/geosciences11090362](https://doi.org/10.3390/geosciences11090362 )
