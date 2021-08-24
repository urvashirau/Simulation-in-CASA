# Multi-Scale and Wideband Simulations and Image Reconstruction.

__Purpose__

To test and characterize Cube and MTMFS deconvolution with extended emission. Algorithms compared are hogbom, multiscale (and asp) for cubes, and mtmfs and mt-asp for wideband continuum imaging. 

__Simulations__
Three types of source structures are simulated into an MS containing 5 channels spanning 1-2 GHz. These sources have multi-scale structure along with spectral structure. Spectral index maps are constructed from the true input images, as well as from cube imaging results and the multi-term imaging outputs. The goal is to easily compare these outputs to assess comparitive algorithm quality. 

 
# Code dependencies

To run the Demo notebook as is, we need 
  - __CASA6__ (version 6.2) All simulation, visibility prediction, and imaging code uses casa6 tasks and tools. 
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
   Download the latest pip wheels for casatools and casatasks from the JIRA ticket (Currently wheels have been obtained from build "dev19" on CAS-940). 
```
  ./local_python3/bin/pip install casatools_xxxxxx.whl
  ./local_python3/bin/pip install casatasks_xxxxxx.whl 
```
### Running without astropy and cngi-prototype
The Sim_Heterogeneous_Array_Scripts.py contains only CASA6 python code, with no dependencies on astropy or cngi-prototype.  If required, methods from this file alone will be able to run all the simulations and imaging steps shown in the demo notebook.  Display and analysis methods, however, will require scripts from the Display_Experiments_Script.py.  Alternatively, the casaviewer or casaplotms from casa5/casa6 may be used for image and MS viewing. 

***

Please refer to [CASA 6 installation instructions](https://casadocs.readthedocs.io/en/latest/notebooks/usingcasa.html) and [CNGI-Prototype Docs](https://cngi-prototype.readthedocs.io/en/latest/) for more details.
