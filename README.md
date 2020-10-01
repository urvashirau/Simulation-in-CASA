# Simulation-in-CASA
Construct interesting data sets to be used for feature demonstrations, test scripts, bug reports and scientific characterization of numerical features in CASA. Each folder in this repository contains its own README.md with information about how to set up and run each set of scripts and notebooks. 

*** 
Simulations scripts and demo notebooks in this repository :

 - Basic Simulations : A demo of how to construct a Measurement Set, populate it with visibilities predicted from a sky model, corrupt it with noise and/or antenna gains, and then calibrate and image it. 
 
 - Heterogeneous Array Simulations : A demo of how to construct a heterogeneous array dataset, including cross-baselines, for a mosaic observation, including visibility prediction and imaging for the case of two different antenna sizes and arbitrary primary beam models specified per antenna (using CASA's mosaic gridder)
 
*** 
 
 - *TO BE ADDED LATER* 
      - Wideband Mosaic Simulations : Scripts to construct a wideband mosaic simulated dataset using a sky model derived from SKADS and wide-field effects such as primary beam squint, rotation, frequency-dependence and sidelobe structure (using CASA's awproject gridder). 
      - Single Dish and Interferometer Combination : Scripts to simulate single dish images and wideband mosaic interferometer data to use for testing joint single dish and interferometer combination algorithms. 
