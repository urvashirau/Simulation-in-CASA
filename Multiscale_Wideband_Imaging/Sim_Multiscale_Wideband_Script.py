#!/usr/bin/env python
# coding: utf-8

# # Install and Import
# 
# __Option 1 : Install local python3__
# 
#     export PPY=`which python3`
#     virtualenv -p $PPY --setuptools ./local_python3
#     ./local_python3/bin/pip install --upgrade pip
#     ./local_python3/bin/pip install --upgrade numpy matplotlib ipython astropy
#     ./local_python3/bin/pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatools
#     ./local_python3/bin/pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatasks
#     ./local_python3/bin/pip3 install jupyter
# 
# ***

from __future__ import print_function

import os
import pylab as pl
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


# Import required tools/tasks
from casatools import simulator, image, table, coordsys, measures, componentlist, quanta, ctsys, ms, vpmanager
from casatasks import tclean, ft, imhead, listobs, exportfits, flagdata, bandpass, applycal,imstat, visstat,mstransform, imsmooth,widebandpbcor
from casatasks.private import simutil

from IPython.display import Markdown as md


# Instantiate all the required tools
sm = simulator()
ia = image()
tb = table()
cs = coordsys()
me = measures()
qa = quanta()
cl = componentlist()
mysu = simutil.simutil()
myms = ms()
vp = vpmanager()

def makeMSFrame(msn = 'sim_data',tel='ALMA'):
    """ 
    Construct an empty Measurement Set that has the desired observation setup
    with uvw/scan/field/ddid setup
 
    Construct an empty Measurement Set that has the desired observation setup. 
    This includes antenna configuration, phase center direction, spectral windows, 
    date and timerange of the observation, structure of scans/spws/obsidd/fieldids 
    (and all other MS metadata). Evaluate UVW coordinates for the entire observation 
    and initialize the DATA column to zero.
    """
        
    msname = msn + '_' + tel + '.ms'
        
    print("Making an MS named : "+msname)
        
    os.system('rm -rf '+msname)
        
    ## Open the simulator
    sm.open(ms=msname);

    ## Read/create an antenna configuration. 
    ## Canned antenna config text files are located here : /home/casa/data/trunk/alma/simmos/*cfg
    ## Fictitious telescopes can be simulated by specifying x, y, z, d, an, telname, antpos.
    ##     x,y,z are locations in meters in ITRF (Earth centered) coordinates. 
    ##     d, an are lists of antenna diameter and name.
    ##     telname and obspos are the name and coordinates of the observatory. 

    
    if tel=='ALMA':
        antennalist = os.path.join( ctsys.resolve("alma/simmos") ,"alma.all.cfg")   
        (x,y,z,d,an,an2, telname, obspos) = mysu.readantenna(antennalist)
        ## Pick a subset : Example : four 7m dishes and eight 12m dishes is to
        antlist = ['N601','N606','J505','J510', 'A001', 'A012','A025', 'A033','A045', 'A051','A065', 'A078']
        obsposname='ALMA'
        dir_pointing = me.direction(rf='J2000', v0='19h59m28.5s',v1='-40d44m01.5s')

    if tel=='VLA':
        antennalist = os.path.join( ctsys.resolve("alma/simmos") ,"vla.d.cfg")   
        (x,y,z,d,an,an2, telname, obspos) = mysu.readantenna(antennalist)
        ## Pick a subset : Example : four 7m dishes and eight 12m dishes is to
        antlist = ['W01','W02','W03','W05','W07','W09','E01','E02','E03','E05','E07','E09','N01','N02','N03','N05','N07','N09']
        #antlist = ['W02','W04','W06','W08','W010','W12','E02','E04','E06','E16','E18','N02','N03','N12','N14','N16','N18']
        obsposname='VLA'
        dir_pointing = me.direction(rf='J2000', v0='19h59m28.5s',v1='+40d40m00.0s')

    if (tel=='NGVLA'):
        antennalist = "ngvla-demo-revC.cfg"  ## Local copy
        ##### Or concat two of them programmatically and then read..... 
        #alist1 = os.path.join( ctsys.resolve("alma/simmos") ,"ngvla-core-revC.cfg")   
        #alist2 = os.path.join( ctsys.resolve("alma/simmos") ,"ngvla-sba-revC.cfg")   
        #os.system('cat '+alist1 + ' > ngvla-demo-revC.cfg')
        #os.system('cat '+alist2 + ' >> ngvla-demo-revC.cfg')
        (x,y,z,d,an,an2,telname, obspos) = mysu.readantenna(antennalist)
        ## Pick a subset : Example : seven 18m dishes and five 6m dishes is to
        antlist = ['m153','m155','m140','m142', 'm130', 'm122','m151', 's012','s006', 's013','s008', 's009']
        obsposname='VLA'
        dir_pointing = me.direction(rf='J2000', v0='19h59m28.5s',v1='+40d44m01.5s')
        telname='NGVLA1'

    antmask = np.isin(np.array(an),antlist)
    xsel=np.array(x)[antmask]
    ysel=np.array(y)[antmask]
    zsel=np.array(z)[antmask]
    dsel=np.array(d)[antmask]
    antsel=np.array(an)[antmask]
    newnamelist = np.array(an)[antmask]

    # To edit antenna names, edit this list.
    #        for ii in range(0,len(newnamelist)):
    #            newnamelist[ii] = newnamelist[ii]+'new'

    ## Set the antenna configuration
    sm.setconfig(telescopename=telname,
                     x=xsel,
                     y=ysel,
                     z=zsel,
                     dishdiameter=dsel,
                     mount=['alt-az'], 
                     antname=list(antsel),
                     #antname=list(newnamelist),
                     coordsystem='global',
                     referencelocation=me.observatory(obsposname));

    ## Set the polarization mode (this goes to the FEED subtable)
    sm.setfeed(mode='perfect X Y', pol=['']);

    ## Set the spectral window and polarization (one data-description-id). 
    ## Call multiple times with different names for multiple SPWs or pol setups.
    if tel=='ALMA':
        sm.setspwindow(spwname="Band3",
                       freq='90GHz',
                       deltafreq='2GHz',
                       freqresolution='1.0MHz',
                       nchannels=5,
                       stokes='XX YY');
        useband='Band3'
    if tel=='VLA':
        sm.setspwindow(spwname="LBand",
                       freq='1.0GHz',
                       deltafreq='0.2GHz',
                       freqresolution='1.0MHz',
                       nchannels=5,
                       stokes='RR LL');
        useband='LBand'
    if (tel=='NGVLA'):
        sm.setspwindow(spwname="Band5",
                       freq='40GHz',
                       deltafreq='2GHz',
                       freqresolution='1.0MHz',
                       nchannels=5,
                       stokes='XX YY');
        useband='Band5'

    ## Setup source/field information (i.e. where the observation phase center is)
    ## Call multiple times for different pointings or source locations.
    sm.setfield( sourcename="fake1", sourcedirection=dir_pointing);
    
    ## Set shadow/elevation limits (if you care). These set flags.
    sm.setlimits(shadowlimit=0.01, elevationlimit='1deg');

    ## Leave autocorrelations out of the MS.
    sm.setauto(autocorrwt=0.0);  

    ## Set the integration time, and the convention to use for timerange specification
    ## Note : It is convenient to pick the hourangle mode as all times specified in sm.observe()
    ##        will be relative to when the source transits.
    sm.settimes(integrationtime='2000s', 
                usehourangle=True,
                referencetime=me.epoch('UTC','2020/10/4/00:00:00'));

    ## Construct MS metadata and UVW values for one scan and ddid 
    ## Call multiple times for multiple scans.
    ## Call this with different sourcenames (fields) and spw/pol settings as defined above.
    ## Timesteps will be defined in intervals of 'integrationtime', between starttime and stoptime.
    sm.observe(sourcename="fake1",
               spwname=useband, 
               starttime='-5.0h', 
               stoptime='+5.0h');

    ## Close the simulator
    sm.close()
    
    ## Unflag everything (unless you care about elevation/shadow flags)
    flagdata(vis=msname,mode='unflag')
    

def makeCompList(clname_true='sim_onepoint.cl',tel='VLA',stype='basic'):
    """
    Make a component list.
    Evaluate this onto an image. ( Later, predict using the Mosaic gridder )
    """
    # Make sure the cl doesn't already exist. The tool will complain otherwise.
    os.system('rm -rf '+clname_true)
    cl.done()
    
    # Add sources, one at a time. 
    # Call multiple times to add multiple sources. ( Change the 'dir', obviously )
    if tel=='ALMA':
        cl.addcomponent(
            #            dir='J2000 19h59m28.5s -40d44m01.5s',  ## Center of PB
            dir='J2000 19h59m28.5s -40d44m21.5s', ## Offset by 100 pixels = 20 arcsec.
            flux=1.0,            # For a gaussian, this is the integrated area.
            fluxunit='Jy', 
            freq='95.0GHz', 
            shape='point',       ## Point source
            #                    shape='gaussian',   ## Gaussian
            #                    majoraxis="5.0arcmin", 
            #                    minoraxis='2.0arcmin', 
            spectrumtype="spectral index",
            index=0.0)

    if (tel=='NGVLA'):
        cl.addcomponent(
            #            dir='J2000 19h59m28.5s -40d44m01.5s',  ## Center of PB
            dir='J2000 19h59m28.5s +40d44m21.5s', ## Offset by 40 pixels = 20 arcsec.
            flux=1.0,            # For a gaussian, this is the integrated area.
            fluxunit='Jy', 
            freq='45.0GHz', 
            shape='point',       ## Point source
            #                    shape='gaussian',   ## Gaussian
            #                    majoraxis="5.0arcmin", 
            #                    minoraxis='2.0arcmin', 
            spectrumtype="spectral index",
            index=0.0)



    if tel=='VLA':
        if stype=='basic':
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h59m28.5s +40d45m0.0s', ## Offset by 150 pixels = 150*10.0arcec= 25arcmin.
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
                  spectrumtype="spectral index",
                index=0.0)
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h58m05.5s +40d55m0.0s', ## Offset by 150 pixels = 150*10.0arcec= 25arcmin.
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
                spectrumtype="spectral index",
                index=-0.5)
            cl.addcomponent(
                dir='J2000 19h59m28.5s +40d35m00.0s',
                flux=50.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="7.0arcmin", 
                minoraxis="7.0arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-1.0)

        if stype=='simple':
            ## One point source with spectral index = 0.0, offset from the center in one direction
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h59m28.5s +40d30m00.0s', 
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
            spectrumtype="spectral index",
                index=0.0)
            
            ## One point source with spectral index = -1.0, offset from the center in the other direction
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h59m28.5s +40d50m0.0s', ## Offset by 150 pixels = 150*10.0arcec= 25arcmin.
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
                #                    shape='gaussian',   ## Gaussian
                #                    majoraxis="5.0arcmin", 
                #                    minoraxis='2.0arcmin', 
                spectrumtype="spectral index",
                index=-1.0)

            cl.addcomponent(
                dir='J2000 19h59m38.0s +40d40m00.0s',
                flux=100.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="5.0arcmin", 
                minoraxis="5.0arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=+1.0)
            cl.addcomponent(
                dir='J2000 19h59m05.0s +40d40m00.0s',
                flux=100.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="5.0arcmin", 
                minoraxis="5.0arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-1.0)
            


        if stype=='jet':
            ## One point source with spectral index = 0.0, offset from the center in one direction
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h59m28.5s +40d30m30.0s', ## NO ! Offset by 75 pixels = 75*10.0arcec= 12arcmin,30asec
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
            spectrumtype="spectral index",
                index=0.0)
            
            ## One point source with spectral index = -1.0, offset from the center in the other direction
            cl.addcomponent(
                #dir='J2000 19h59m28.5s -40d40m00.0s',  ## Center of PB
                dir='J2000 19h59m20.0s +40d55m0.0s', ## Offset by 150 pixels = 150*10.0arcec= 25arcmin.
                flux=1.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='point',       ## Point source
                #                    shape='gaussian',   ## Gaussian
                #                    majoraxis="5.0arcmin", 
                #                    minoraxis='2.0arcmin', 
                spectrumtype="spectral index",
                index=-1.0)
            
            
            ## Add a few extended sources to construct a jet and lobe structure. 
            cl.addcomponent(
                dir='J2000 19h59m28.5s +40d38m00.0s',
                flux=6.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="4.0arcmin", 
            minoraxis="1.0arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-0.2)
            cl.addcomponent(
                dir='J2000 19h59m28.5s +40d43m00.0s',
                flux=10.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="5.5arcmin", 
                minoraxis="1.5arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-0.3)
            cl.addcomponent(
                dir='J2000 19h59m28.5s +40d48m00.0s',
                flux=20.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="6.0arcmin", 
                minoraxis="2.5arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-0.5)
            cl.addcomponent(
                dir='J2000 19h59m28.5s +40d54m00.0s',
                flux=50.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="7.0arcmin", 
                minoraxis="5.0arcmin", 
                positionangle="0.0deg",
                spectrumtype="spectral index",
                index=-0.6)
            cl.addcomponent(
                dir='J2000 19h59m05.0s +40d55m00.0s',
                flux=40.0,            # For a gaussian, this is the integrated area.
                fluxunit='Jy', 
                freq='1.5GHz', 
                shape='Gaussian',   ## Gaussian
                majoraxis="6.0arcmin", 
                minoraxis="4.0arcmin", 
                positionangle="15.0deg",
                spectrumtype="spectral index",
                index=-0.8)
            
            
    # Print out the contents of the componentlist
    #print('Contents of the component list')
    #print(cl.torecord())
    
    # Save the file
    cl.rename(filename=clname_true)
    cl.done()


def makeEmptyImage(imname_true='sim_onepoint_true.im', tel='ALMA'):
    ## Define the center of the image
    radir = '19h59m28.5s'

    if tel=='ALMA':
        decdir = '-40d44m21.5s'   ## Place the image center at the source location
        reffreq = '90.0GHz'
        freqinc = '1.0GHz'
        cell='0.2arcsec'
        imsize=2048
    if tel=='VLA':
        decdir = '+40d40m00.0s'   ## Place the image center at the source location
        reffreq = '1.0GHz'
        freqinc = '0.2GHz'
        cell='12.0arcsec'
        imsize=512
    if (tel=='NGVLA'):
        decdir = '+40d44m21.5s'   ## Place the image center at the source location
        reffreq = '40.0GHz'
        freqinc = '1.0GHz'
        cell='0.5arcsec'
        imsize=2048
    
    ## Make the image from a shape
    ia.close()
    ia.fromshape(imname_true,[imsize,imsize,1,5],overwrite=True)
    
    ## Make a coordinate system
    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad=qa.convert(qa.quantity(cell),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([qa.convert(radir,'rad')['value'],qa.convert(decdir,'rad')['value']],type="direction")
    cs.setreferencevalue(reffreq,'spectral')
    cs.setreferencepixel([0],'spectral')
    cs.setincrement(freqinc,'spectral')
    
    ## Set the coordinate system in the image
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.set(0.0)
    ia.close() 
    
### Note : If there is an error in this step, subsequent steps will give errors of " Invalid Table Operation : SetupNewTable.... imagename is already opened (is in the table cache)"
## The only way out of this is to restart the kernel (equivalent to exit and restart CASA).
## Any other way ? 



# In[10]:


def evalCompList(clname='sim_onepoint.cl', imname='sim_onepoint_true.im'):
    """
    # __Evaluate the component list onto the image cube__
    ##  Evaluate a component list
    """
    cl.open(clname)
    ia.open(imname)
    ia.modify(cl.torecord(),subtract=False)
    ia.close()
    cl.done()


def smoothModel(imname='sim_jet_vla.im', rbeam='50.0arcsec'):
    imsmooth(imagename=imname,
             outfile=imname+'.smoothed',
             kernel='gauss',
             major=rbeam,
             minor=rbeam,
             pa='0.0deg',
             targetres=True,
             overwrite=True)




# get_ipython().run_line_magic('pinfo', 'pl.contour')


# # Simulate visibilities from the sky model
# 
# Simulate visibilities for the true sky model, using the Mosaic gridder.
# __Use imager (or ft)__
# 
# Visibilities are predicted and saved in the MODEL_DATA column of the MS. The values must then be copied to the DATA column. Use this approach when non-standard gridders are required, typically when instrument-dependent effects are included, or when Taylor-coefficient wideband image models are to be used for visibility prediction.
# 
# __Step 1__ : Simulate visibilities into the MODEL column using tclean
# 
# tclean can be used for model prediction with all gridders ('standard', 'wproject', 'mosaic', 'awproject'). Wide-field and full-beam effects along with parallactic angle rotation may be included with appropriate settings. tclean can predict model visibilities only from input images and not component lists.

## Use an input model sky image - widefield gridders
def predictImager(msname='sim_data_ALMA.ms',
                  imname_true='sim_onepoint_true.im',
                  gridder='mosaic',
                  tel='ALMA',
                  vptable=''):

    if tel=='ALMA':
        cell='0.2arcsec'
        imsize=2048
    if tel=='VLA':
        cell='12.0arcsec'
        imsize=512
    if (tel=='NGVLA'):
        cell='0.5arcsec'
        imsize=2048

    
    os.system('rm -rf sim_predict.*')
    
    # Run tclean in predictModel mode. 
    tclean(vis=msname,
       startmodel=imname_true,
       vptable = vptable,
       imagename='sim_predict',
       savemodel='modelcolumn',
       imsize=imsize,
       cell=cell,
       specmode='cube',
       interpolation='nearest',
       nchan=-1,
       #start='90.0GHz',
       #width='1.0GHz',
       #nchan=10,
       #reffreq='95.0GHz',
       gridder=gridder,
       normtype='flatsky',  # sky model is flat-sky
       #cfcache='sim_predict.cfcache',
       wbawp=True,      # ensure that gridders='mosaic' and 'awproject' do freq-dep PBs
       pblimit=0.05,    
       conjbeams=False, 
       calcres=False, 
       calcpsf=True, 
       niter=0, 
       wprojplanes=1)


def copyModelToData(msname='sim_data_ALMA.ms'):
    """
    ### Copy visibilities from the MODEL column to the data columns
    ### This is required when predicting using tclean or ft as they will only write to the MODEL column
    Also, initialize the WEIGHT_SPECTRUM column so that the "visstat" task may be used later.....
    """
    tb.open(msname,nomodify=False);
    moddata = tb.getcol(columnname='MODEL_DATA');
    tb.putcol(columnname='DATA',value=moddata);
    #tb.putcol(columnname='CORRECTED_DATA',value=moddata);
    moddata.fill(0.0);
    tb.putcol(columnname='MODEL_DATA',value=moddata);
    tb.close();
    
    ## Add the WEIGHT_SPECTRUM column. This is needed only if you intend to use task_visstat on this dataset. 
    os.system('rm -rf tmp_addedcol.ms')
    mstransform(vis=msname,outputvis='tmp_addedcol.ms',datacolumn='DATA',usewtspectrum=True)
    os.system('rm -rf '+msname)
    os.system('cp -r tmp_addedcol.ms '+msname)




def addNoiseSim(msname='sim_data_ALMA.ms'):
    """
    ## Add Dish-size dependent noise - Use built in method
    """
    sm.openfromms(msname);
    sm.setnoise(mode='tsys-atm');
    sm.corrupt();   
    sm.close();


def getFreqList(imname=''):

      ia.open(imname)
      csys =ia.coordsys()
      shp = ia.shape()
      ia.close()

      if(csys.axiscoordinatetypes()[3] == 'Spectral'):
           restfreq = csys.referencevalue()['numeric'][3]#/1.0e+09; # convert more generally..
           freqincrement = csys.increment()['numeric'][3]# /1.0e+09;
           freqlist = [];
           for chan in range(0,shp[3]):
                 freqlist.append(restfreq + chan * freqincrement);
      elif(csys.axiscoordinatetypes()[3] == 'Tabular'):
           freqlist = (csys.torecord()['tabular2']['worldvalues']) # /1.0e+09;
      else:
           casalog.post('Unknown frequency axis. Exiting.','SEVERE');
           return False;

      csys.done()
      return freqlist

from scipy.optimize import curve_fit
import numpy as np

def func_powerlaw(nus, inu0, alpha, beta):
    nu = nus['nu']
    nu0 = nus['nu0']
    return inu0 * np.power((nu)/nu0 , (alpha+ beta*(np.log(nu/nu0))))

def fit_spectrum(cubename='', intensity='', alpha='', beta='', pixt=0.1):
    ia.open(cubename)
    datpix = ia.getchunk()
    csys = ia.coordsys()
    shp = ia.shape()
    ia.close()
    #print('Cube shape : ',shp)

    nu = np.array(getFreqList(cubename))/1e+9
    print('Frequencies (Ghz) : ',nu)
    nu0 = (nu[0] + nu[len(nu)-1])/2.0

    alpha_arr = np.zeros((shp[0],shp[1],1,1),'float')
    beta_arr = np.zeros((shp[0],shp[1],1,1),'float')
    int_arr = np.zeros((shp[0],shp[1],1,1),'float')

    target_func = func_powerlaw

    print("Start fitting...")
    for ii in range(0,shp[0]):
        for jj in range(0,shp[1]):
            if datpix[ii,jj,0,0] > pixt:
                y = datpix[ii,jj,0,:]
                popt, pcov = curve_fit(f=target_func, xdata={'nu':nu,'nu0':nu0}, ydata=y,p0=[1.0,0.0,0.0])
                alpha_arr[ii,jj,0,0] = popt[1]
                beta_arr[ii,jj,0,0] = popt[2]
                int_arr[ii,jj,0,0] = popt[0]
            
    pl.figure(figsize=(12,7))
    pl.clf()
    pl.subplot(131)
    pl.imshow(int_arr[:,:,0,0].transpose(),origin='lower')     
    pl.subplot(132)
    pl.imshow(alpha_arr[:,:,0,0].transpose(),origin='lower',cmap='jet',vmin=-1.0, vmax=0.2)     
    pl.subplot(133)
    pl.imshow(beta_arr[:,:,0,0].transpose(),origin='lower',cmap='jet',vmin=-0.2, vmax=0.5)     

    print("Save outputs")
    ia.fromarray(outfile=intensity, pixels=int_arr, csys=csys.torecord(),overwrite=True)
    ia.fromarray(outfile=alpha, pixels=alpha_arr, csys=csys.torecord(),overwrite=True)
    ia.fromarray(outfile=beta, pixels=beta_arr, csys=csys.torecord(),overwrite=True)

###################################################################################


def calcImageAccuracy( truthimage='', outputimage=''):
    '''
    Calculate image fidelity metrics to compare truth and reconstructed images

    truthimage : Model image convolved with a beam that represents the desired target resolution.
    outputimage : Restored image, convolved to the same target resoution as truthimage.
    '''

    ia.open(truthimage)
    shp = ia.shape()
    ia.close()
    #print("Shape:"+str(shp))

    for chan in range(0,1): ## shp[3]):
        ia.open(truthimage)
        imtrue = ia.getchunk(blc=[0,0,0,chan],trc=[shp[0],shp[1],0,chan])
        ia.close()
        ia.open(outputimage)
        imout = ia.getchunk(blc=[0,0,0,chan],trc=[shp[0],shp[1],0,chan])
        ia.close()

        tmax = np.max(imtrue)

        ## Relative Error. This is also the inverse of the formal definition of "imaging fidelity"
        imdiff = (imtrue - imout)/tmax

        ## Chi-square
        chisq = np.sqrt(np.mean((imtrue - imout)**2))
        

        pl.figure(figsize=(4,4))
        pl.clf()
        pl.subplot(111)
        pl.imshow(imdiff[:,:,0,0].transpose(),origin='lower',interpolation=None,cmap='jet')#,vmin=-0.005, vmax=+0.005)     
        pl.colorbar()
        pl.title('Relative Difference between\n true and output images')

        print("Standard deviation of diff (chan %d) = %3.6f."%(chan, chisq))
