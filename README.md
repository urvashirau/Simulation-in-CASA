# Simulation-in-CASA
Construct interesting data sets to be used for feature demonstrations, test scripts, bug reports and scientific characterization of numerical features in CASA. Each folder in this repository contains its own README.md with information about how to set up and run each set of scripts and notebooks. 

*** 
Simulations scripts and demo notebooks in this repository :

 - Basic Simulations : A demo of how to construct a Measurement Set, populate it with visibilities predicted from a sky model, corrupt it with noise and/or antenna gains, and then calibrate and image it. 
 
 - Heterogeneous Array Simulations : A demo of how to construct a heterogeneous array dataset, including cross-baselines, for a mosaic observation, including visibility prediction and imaging for the case of two different antenna sizes (using CASA's mosaic gridder). These simulations are scaled-down demos for use cases of ALMA 12m + 7m with cross baselines, the NGVLA with 18m and 6m antennas used together, and the specification of arbitrary primary beam models as images. 
 
*** 
 
 - *TO BE ADDED LATER* 
      - Wideband Mosaic Simulations : Scripts to construct a wideband mosaic simulated dataset using a sky model derived from SKADS and wide-field effects such as primary beam squint, rotation, frequency-dependence and sidelobe structure (using CASA's awproject gridder). These simulations were used for an [analysis of wideband mosaic imaging accuracy](http://iopscience.iop.org/article/10.3847/0004-6256/152/5/124/meta). 
      - Single Dish and Interferometer Combination : Scripts to simulate single dish images and wideband mosaic interferometer data to use for testing joint single dish and interferometer combination algorithms. These simulations were used to develop the [SDINT joint reconstruction algorithm](https://iopscience.iop.org/article/10.3847/1538-3881/ab1aa7/meta).
