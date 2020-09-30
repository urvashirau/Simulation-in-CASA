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

To run these examples, we need 
  - CASA6 (version 6.1 or later, for ALMA,  version 6.2 (or tarball from CAS-13010) for NGVLA)
  - astropy
  - cngi-prototype
  
