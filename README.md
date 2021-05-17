# Simulation-in-CASA
Construct interesting data sets to be used for feature demonstrations, test scripts, bug reports and scientific characterization of numerical features in CASA. Each folder in this repository contains its own README.md with information about how to set up and run each set of scripts and notebooks. 

*** 
__Simulations scripts and demo notebooks in this repository__ :

 - __Basic Simulations__ : A demo of how to construct a Measurement Set, populate it with visibilities predicted from a sky model, corrupt it with noise and/or antenna gains, and then calibrate and image it. 
 
 - __Heterogeneous Array Simulations__ : A demo of how to construct a heterogeneous array dataset, including cross-baselines, for a mosaic observation, including visibility prediction and imaging for the case of two different antenna sizes (using CASA's mosaic gridder). These simulations are scaled-down demos for use cases of ALMA 12m + 7m with cross baselines, the NGVLA with 18m and 6m antennas used together, and the specification of arbitrary primary beam models as images. Numerical results are compared against theoretically calculated and predicted estimates for an absolute measure of accuracy. 
 
 - __Multiscale Wideband Simulations__ : A demo and comparison of the Hogbom, MS-Clean, MT-MFS, (and soon) Asp and MT-Asp deconvolution algorithms. Simulation datasets include compact and extended emission with varying degrees of complexity in both spatial and spectral structure. The purpose is to perform a direct comparison of deconvolution algorithms and cube vs mtmfs for the recovery of spatial and spectral characteristics. 

 - __Widefield Wideband Simulations__ : Demos to quantify the accuracy at which wideband mosaic reconstructions are currently possible using CASA imaging options. Currently, a notebook exists to demo all currently available options. Additional demos will be added for thorough numerical characterization, as well as include options such as linear mosaics for multi-term wideband imaging. Imaging options that will be evaluated against simulated truth are the 'mosaic' and 'awproject' gridders, along with hogbom deconvolution of point sources. 

*** 
*** 
*** 
*** 
 
 
 - *To be added later* 
      - More detailed Wideband Mosaic Simulations : Scripts to construct a wideband mosaic simulated dataset using a sky model derived from SKADS and wide-field effects such as primary beam squint, rotation, frequency-dependence and sidelobe structure (using CASA's awproject gridder). These simulations were used for an [analysis of wideband mosaic imaging accuracy](http://iopscience.iop.org/article/10.3847/0004-6256/152/5/124/meta). 
      - Single Dish and Interferometer Combination : Scripts to simulate single dish images and wideband mosaic interferometer data to use for testing joint single dish and interferometer combination algorithms. These simulations were used to develop the [SDINT joint reconstruction algorithm](https://iopscience.iop.org/article/10.3847/1538-3881/ab1aa7/meta).
