# Heterogeneous Array Mosaic Simulations and Imaging

CASR-524 : ALMA. CASR-301 : ngVLA : ALMA (CASR-514) and NGVLA (CASR-301) need the ability to simulate and image data from an array containing different dish sizes, arbitrarily specifiable names for the telescope and antenna types, support for cross baseline visibilities, and realistic noise predictions for each baseline type. 

__Topics Covered in this notebook__

    - Construct an MS with multiple sets of antenna sizes and cross-baselines
    - Predict visibilities with correct PB application during simulation
    - Image the data with correct PB application during imaging
    - Introduce noise per baseline type and track it through to imaging
    
__Use Cases__

    - ALMA single pointing + mosaic  cube
    - ngVLA single pointing + mosaic  cube
   
__Tests__
   
(1) Test het-array model prediction and imaging using the Mosaic gridder.
   
   For a 1 Jy source, 
   - For baselines within each antenna type, image and PB should match.
   - For cross baselines, PB(cross) = sqrt( PB(type1) * PB(type2) )  and similarly for the image flux. Image and PB should match.
   - For all baselines together, image and PB should match
   
(2) Test adding noise and measuring noise in cleaned images
   
   - Add noise to all visibilities. 
       - sm.setnoise(mode='tsys-atm') will do dish sizes automatically (but not yet accurate for cross-baselines)
       - custom function to compute noise based only on dish diameter
   - Image different subsets and measure noise level.
   - Calculated expected noise levels for cross baselines and joint imaging, using appropriate weighted and quadrature averages, and measure what we get. They should match. 

(3) Test setting arbitrary antenna and telescope names (and PB models) 
   - Known and Unknown observatory name
   - Arbitrary antenna names
   - For unknown PB models,
       - Default : Calculate an Airy Dish model from the dish diameter information in the ANTENNA subtable
       - Custom : Use the vptable to specify custom Voltage Patterns as images

   
# Code dependencies

To run the Demo notebook as is, we need 
  - __CASA6__ (version 6.1 or later, for ALMA,  version 6.2 (or tarball from CAS-13010) for NGVLA) : All simulation, visibility prediction, and imaging code uses casa6 tasks and tools. 
  - __astropy__ : Image displays use WCS from astropy for world-coordinates. (This is a very simple option for imaging viewing in a notebook)
  - __cngi-prototype__ : Used for data plots. Measurement Sets are converted to Zarr and read in as an XArray Dataset and plotted using XArray selection and plot methods. (This is a very simple option for plotting visibilities, uv-coverage, antenna names/locations, etc in a notebook)
  
## A minimal python environment for this notebook may be generated as follows 
```
  export PPY=`which python3`
  virtualenv -p $PPY --setuptools ./local_python3
  ./local_python3/bin/pip install --upgrade pip
  ./local_python3/bin/pip install ipython numpy
  ./local_python3/bin/pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatools
  ./local_python3/bin/pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatasks
  ./local_python3/bin/pip3 install jupyter
  ./local_python3/bin/pip install cngi-prototype
  ./local_python3/bin/pip install matplotlib astropy
```
### Installing casa6 from local pip wheels
   Download the latest pip wheels for casatools and casatasks from the JIRA ticket (in this case, CAS-13010 under development for casa 6.2). 
```
  ./local_python3/bin/pip install casatools_xxxxxx.whl
  ./local_python3/bin/pip install casatasks_xxxxxx.whl 
```
### Running without astropy and cngi-prototype
The Sim_Heterogeneous_Array_Scripts.py contains only CASA6 python code, with no dependencies on astropy or cngi-prototype.  If required, methods from this file alone will be able to run all the simulations and imaging steps shown in the demo notebook.  Display and analysis methods, however, will require scripts from the Display_Experiments_Script.py.  Alternatively, the casaviewer or casaplotms from casa5/casa6 may be used for image and MS viewing. 

***

Please refer to [CASA 6 installation instructions](https://casadocs.readthedocs.io/en/latest/notebooks/usingcasa.html) and [CNGI-Prototype Docs](https://cngi-prototype.readthedocs.io/en/latest/) for more details.
