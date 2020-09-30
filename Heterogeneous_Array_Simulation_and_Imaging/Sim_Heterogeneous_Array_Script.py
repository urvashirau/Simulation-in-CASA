#!/usr/bin/env python
# coding: utf-8

# # Heterogeneous Array Mosaic Simulations and Imaging

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
# __Option 2 : Install at runtime__ (for Google Colab)
# 
#     import os
#     print("installing pre-requisite packages...")
#     os.system("apt-get install libgfortran3")
#     print("installing casatasks...")
#     os.system("pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatools")
#     os.system("pip install --extra-index-url https://casa-pip.nrao.edu/repository/pypi-group/simple casatasks")
#     print("building config files...")
#     os.system("mkdir ~/.casa")
#     !echo home, datapath = \'/content/\', [\'/content/\'] > ~/.casa/toolrc.py
#     !more /root/.casa/toolrc.py
#     print('complete')
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
from casatasks import tclean, ft, imhead, listobs, exportfits, flagdata, bandpass, applycal,imstat, visstat,mstransform
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
        ## Two pointings, 20 arcsec away from the source, in either direction
        dir_pointing_1 = me.direction(rf='J2000', v0='19h59m28.5s',v1='-40d44m01.5s')
        dir_pointing_2 = me.direction(rf='J2000', v0='19h59m28.5s',v1='-40d44m51.5s')

 ### Add content for ngVLA too
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
        ## Two pointings, 40 arcsec away from the source, in either direction
        dir_pointing_1 = me.direction(rf='J2000', v0='19h59m28.5s',v1='+40d43m41.5s')
        dir_pointing_2 = me.direction(rf='J2000', v0='19h59m28.5s',v1='+40d45m11.5s')
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
    sm.setfield( sourcename="fake1",
                 sourcedirection=dir_pointing_1);
    sm.setfield( sourcename="fake2",
                 sourcedirection=dir_pointing_2);

    
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
               stoptime='+0.0h');
    sm.observe(sourcename="fake2",
               spwname=useband, 
               starttime='0.0h', 
               stoptime='+5.0h');

    ## Close the simulator
    sm.close()
    
    ## Unflag everything (unless you care about elevation/shadow flags)
    flagdata(vis=msname,mode='unflag')
    

def makeCompList(clname_true='sim_onepoint.cl',tel='ALMA'):
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
        cell='0.2arcsec'
        imsize=2048
    if (tel=='NGVLA'):
        decdir = '+40d44m21.5s'   ## Place the image center at the source location
        reffreq = '40.0GHz'
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
    cs.setincrement('1.0GHz','spectral')
    
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


def makeVPImage(imname='sim_vp_A',tel='ALMA',atype='A',dtype='float'):
    
    sample_vp = os.path.join( ctsys.resolve("alma/responses") ,"ALMA_0_DA__0_0_360_0_45_90_153.5_163_163_GHz_ticra2007_EFP.im")   
    #sample_vp = '/home/casa/data/trunk/alma/responses/ALMA_0_DA__0_0_360_0_45_90_153.5_163_163_GHz_ticra2007_EFP.im'

    ia.open(sample_vp)
    csys_template = ia.coordsys()
    ia.close()

    ## Make the image from a shape
    if tel=='ALMA':
        reffreq = '90.0GHz'
        freq = 90.0
        cell='0.2arcsec'
        cellval = 0.2
        imsize=2048
        if atype=='A':
            wid = calc_ang(freq, 12.0)  # arcmin
        else:
            wid = calc_ang(freq, 7.0)  # arcmin
        print("Calculated PB size for type A (dia=%2.2f) : %3.5f arcmin"%(D_A, calc_ang(freq,D_A)))
    if (tel=='NGVLA'):
        reffreq = '40GHz'
        freq = 40.0
        cell='0.5arcsec'
        cellval=0.5
        imsize=2048
        if atype=='A':
            wid= calc_ang(freq, 18.0)  #arcmin
        else:
            wid= calc_ang(freq, 6.0)   #arcmin
    if tel=='KSA':  ## Kitchen sink array.  Pick larger dishes for NGVLA1, so that the beams are smaller..... faster for verification tests. 
        reffreq = '40GHz'
        freq = 40.0
        cell='0.5arcsec'
        cellval=0.5
        imsize=512
        if atype=='A':
            wid= calc_ang(freq, 20.0)  #arcmin
        else:
            wid= calc_ang(freq, 10.0)   #arcmin

    
    ## Make the image from a shape
    ia.close()
    if dtype=='float':
        ia.fromshape(imname,[imsize,imsize,4,5],overwrite=True,type='f')
    else:
        ia.fromshape(imname,[imsize,imsize,4,5],overwrite=True,type='c')
        
    
    csys_template.setreferencevalue(reffreq,'spectral')
    csys_template.setreferencepixel([0],'spectral')
    csys_template.setincrement('1.0GHz','spectral')
    
    csys_template.setreferencepixel([imsize/2, imsize/2], 'direction')
    ## Set the coordinate system in the image
    ia.setcoordsys(csys_template.torecord())
    #    ia.setbrightnessunit("Jy/pixel")
    # ia.set(0.0)
    
    ## Now, fill the pixels with frequency dependent beams. 
    if dtype=='float':
        pix = np.zeros([imsize,imsize,4,5],float)   ## Real valued for now
    else:
        pix = np.zeros([imsize,imsize,4,5],complex)   ## Complex valued for now
    
    wid_arcsec = wid * 60.0
    wid_pix = wid_arcsec / cellval
    
    fov2 = imsize*cellval  # half field of view in arcsec.

    for chan in range(0,5):
        chanfreq = freq + 1.0*chan  # 1.0 GHz chanwidth specified above.
        x, y = np.meshgrid(np.linspace(-fov2,fov2,imsize), np.linspace(-fov2,fov2,imsize))
        d = np.sqrt(x*x+y*y)
        sigma, mu = 0.5*wid_arcsec * freq/chanfreq, 0.0
        g = np.exp(-( (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
        
        pix[:,:,0,chan] = g
        pix[:,:,3,chan] = g
        
    ia.putchunk(pix)
   
    ia.close() 
    

def makeVPTable(vptab='pbmod_ngvla.vp',pbA='', pbB='',dtype='float'):
    
    antlist = ['m153','m155','m140','m142', 'm130', 'm122','m151', 's012','s006', 's013','s008', 's009']

    os.system('rm -rf '+vptab)
    vp.reset()
    if dtype=='float':
        vp.setpbimage(telescope='NGVLA1', realimage=pbA, antnames=antlist[0:7])
        vp.setpbimage(telescope='NGVLA1', realimage=pbB, antnames=antlist[7:12])
    else:
        vp.setpbimage(telescope='NGVLA1', compleximage=pbA, antnames=antlist[0:7])
        vp.setpbimage(telescope='NGVLA1', compleximage=pbB, antnames=antlist[7:12])
    vp.summarizevps()
    vp.saveastable(vptab)


def addNoiseSim(msname='sim_data_ALMA.ms'):
    """
    ## Add Dish-size dependent noise - Use built in method
    """
    sm.openfromms(msname);
    sm.setnoise(mode='tsys-atm');
    sm.corrupt();   
    sm.close();

def addDishNoise(msname='sim_data_ALMA.ms'):
    """
    ## Add Dish-size dependent noise - Use custom method
    ### This cell reads SIGMA and adds noise to DATA accordingly. 
    ### The sm.observe() already sets up the SIGMA and WEIGHT correctly. 

    # Calculations for each baseline : 
    #  - Eff_Area =  ( Dia 1 * Dia 2 ) / (max_dia * max_dia)
    #  - SIGMA = 1/Eff_Area
    #  - WEIGHT = 1/SIGMA^2
    #  - Add Gaussian random noise to the real and imag parts of the visibilities, for each channel and pol. 
    """
    myms.open(msname,nomodify=False)
    myms.iterinit(interval=1000)  # Time chunk size in seconds
    myms.iterorigin()
    
    moretodo = True
    while moretodo:
        ## Extract the antenna indices, and the SIGMA and DATA columns
        dat = myms.getdata(items=['SIGMA','DATA'])
        
        ## Now calculate and add noise as per SIGMA to the DATA column
        shp = dat['data'].shape
        noise_re = np.random.normal(loc=0.0,scale=np.repeat(dat['sigma'][:,np.newaxis,:],shp[1],axis=1))
        noise_im = np.random.normal(loc=0.0,scale=np.repeat(dat['sigma'][:,np.newaxis,:],shp[1],axis=1))
        #print(dat['data'].shape)
        noise_re = noise_re/100.0
        noise_im = noise_im/100.0
        dat['data'] = dat['data'] + noise_re + 1j*noise_im

        ### code note : broadcasting doesn't work here as the random numbers need to be different for each channel.
        #dat['data'] = np.ones(dat['data'].shape) + np.random.normal(loc=0.0,scale=dat['sigma'].reshape(shp[0],1,shp[2]))
        
        ## Write DATA back into the MS
        myms.putdata(dat)
        
        moretodo = myms.iternext()

    myms.close()


def calcSigmaWeightAndAddDishNoise(msname='sim_data_noise_ALMA.ms'):
    """
    ### This cell computes SIGMA and WEIGHT according to the antenna diameter, and then does DATA noise
    ### This does not write anything into the MS
    ### The purpose is to check that the sm.observe calculation of weights and sigma is accurate. 
    """
    tb.open(msname+'/ANTENNA')
    ant_dia = tb.getcol('DISH_DIAMETER')
    tb.close()
    maxdia = np.max(ant_dia) # for some normalization
    
    myms.open(msname,nomodify=False)
    myms.iterinit(interval=1000)  # Time chunk size in seconds
    myms.iterorigin()
    
    moretodo = True
    while moretodo:
        ## Extract the antenna indices, and the SIGMA and WEIGHT and DATA columns
        antdat = myms.getdata(items=['ANTENNA1','ANTENNA2'])
        dat = myms.getdata(items=['SIGMA','WEIGHT','DATA'])
        
        ## Effective area per baseline, normalized to 1.0 for the largest diameter elements
        areas = (ant_dia[antdat['antenna1']] * ant_dia[antdat['antenna2']])/(maxdia**2)
        
        print('OLD')
        
        #print(dat['weight'].shape)
        print(dat['sigma'])
        
        ## SIGMA per visibility (baseline) is inversely proportional to Effective Area
        ## Here too, the SIGMA is normalized to 1.0 for the largest diameter baselines
        dat['sigma'] = np.repeat( (1.0/areas)[np.newaxis,:],2,axis=0 )
        
        ## WEIGHT is initialized to 1/sigma^2.  Again, this is 1.0 for the largest ants. 
        dat['weight'] = np.repeat( (1.0/(dat['sigma']**2))[np.newaxis,:],2,axis=0 )
        
        ## Now calculate and add noise of this level to the DATA column
        shp = dat['data'].shape
        noise = np.random.normal(loc=0.0,scale=np.repeat(dat['sigma'][:,np.newaxis,:],shp[1],axis=1))
        dat['data'] = np.ones(dat['data'].shape) + noise
        
        ### code note : broadcasting doesn't work here as the random numbers need to be different for each channel.
        #dat['data'] = np.ones(dat['data'].shape) + np.random.normal(loc=0.0,scale=dat['sigma'].reshape(shp[0],1,shp[2]))
        
        #print(dat['data'][:,:,0])
        print('NEW')
        print(dat['sigma'])
        
        ## Write SIGMA and WEIGHT back into the MS.
        #myms.putdata(dat)
        
        moretodo = myms.iternext()
        
    myms.close()




def getBaselineTypes(msname='',tel='ALMA'):
    """
    Define MSSelection strings for all 3 types of baselines, and all.
    Check them if an MS name is supplied. 
    """

    if tel=='ALMA':
        antsels = {'A':'A*&',
                   'B':'J*,N*&',
                   'cross':'A* && J*,N*',
                   'all':'*'}
    if (tel=='NGVLA' or tel=='KSA'):
        antsels = {'A':'m*&',
                   'B':'s*&',
                   'cross':'m* && s*',
                   'all':'*'}


    if msname != "":
        for atype in antsels.keys():
            asel = myms.msseltoindex(vis=msname, baseline=antsels[atype])
            print(antsels[atype])
            print(asel['antenna1'])
            print(asel['antenna2'])
            print(asel['baselines'].transpose())
   
    return antsels


def imageAnts(vis='sim_data_ALMA.ms',imname='try_ALMA', field='0', antsel='',tel='ALMA',vptable=''):
    """
    Run imager..... 
    TODO  set pblimit=-1 ... or set the image stats method to look inside the PB region only. 
    """
    antsels =  getBaselineTypes(tel=tel)

    if antsel not in antsels.keys():
        print('Pick an antenna selection from '+ str(antsels))
        return
   
    
    if tel=='ALMA':
        cell='0.4arcsec'
        imsize=1024
        phasecenter='J2000 +19h59m28.5s -40d44m21.5s'
    if (tel=='NGVLA'):
        cell='1.0arcsec'
        imsize=1024
        phasecenter='J2000 +19h59m28.5s +40d44m21.5s'

   
    if field=='0' or field=='1':
        ftype = 'single'
    else:
        ftype = 'mosaic'
    
    imname1 = imname+'_'+antsel+'_'+ftype
    
    os.system('rm -rf '+imname1+'.*')
     
    # Run tclean 
    tclean(vis=vis,
       antenna=antsels[antsel],
       field=field, 
       imagename=imname1,
       imsize=imsize,
       cell=cell,
       phasecenter=phasecenter,
       specmode='cube',
       interpolation='nearest',
       nchan=-1,
       gridder='mosaic',
       vptable=vptable,
       normtype='flatnoise', 
       wbawp=True,      
       pblimit=0.05,    
       conjbeams=False, 
       niter=1000,
       nsigma=3.0,
       datacolumn='data',
       weighting='natural')


def getVisStat(vis='sim_data_ALMA.ms',field='',spw='', antsels={}):
    """
    ### This cell examines the dataset and checks noise level per baseline type. 
    ## Use visstat with antenna selection ? 
    """

    vis_summ = {}
    
    for btype in antsels.keys():
        vrec = visstat(vis=vis,field=field,spw=spw,antenna=antsels[btype],axis='amp',datacolumn='data')['DATA_DESC_ID=0']
        wrec = visstat(vis=vis,field=field,spw=spw,antenna=antsels[btype],axis='weight_spectrum',datacolumn='data')['DATA_DESC_ID=0']
        #print(vrec)
        vis_summ[btype] = {'mean':vrec['mean'], 'stddev':vrec['stddev'], 'sumwt':vrec['sumOfWeights'], 'npts':vrec['npts'], 'meanwt':wrec['mean']}
    
 #   df = pd.DataFrame(vis_summ)
 #   print("Summary of stats from the input MS")
 #   with pd.option_context('display.float_format', '{:0.8f}'.format):
 #       print(df)
 #   print("\n")
    return vis_summ


def checkVals(vis='',field='',spw='',pb_A='', pb_B='', antsels={}, meas_peaks={}, meas_rms={},tel='ALMA'):
    """
    For the given (selected) MS and PB images, pick out the PB values at the source location for antenna types A and B
    Calculate expected Intensity for the A-A, B-B, A-B and All baseline imaging cases.
    Using sigmas (or weights) and input vis-rms values, calculate expected rms for all 4 imaging cases
    Compare expected intensity, PB and RMS values with measured values. 
    
    ALMA : A = 12m.   B = 7m
    ngVLA : A = 18m.  B = 6m 
   
    """

    if tel=='ALMA':
        D_A = 12.0
        D_B = 7.0
        label_A = '12m'
        label_B = '7m'
        freq = 90.0 # GHz
    if (tel=='NGVLA'):
        D_A = 18.0
        D_B = 6.0
        label_A = '18m'
        label_B = '6m'
        freq = 40.0 # GHz
    if (tel=='KSA'):
        D_A = 20.0
        D_B = 10.0
        label_A = '20m'
        label_B = '10m'
        freq = 40.0 # GHz

    antsels3 = getBaselineTypes(tel=tel)
    antsels3.pop('all')
    vis_summ = getVisStat(vis=vis,field=field,spw=spw, antsels=antsels3)
    
    calc_peak_A = vis_summ['A']['mean']
    calc_peak_B = vis_summ['B']['mean']
    calc_peak_cross = vis_summ['cross']['mean']
    
    
    calc_peak_all = ( vis_summ['A']['mean']* vis_summ['A']['meanwt'] *  vis_summ['A']['npts'] +                \
                              vis_summ['B']['mean']* vis_summ['B']['meanwt']* vis_summ['B']['npts'] +                        \
                              vis_summ['cross']['mean']* vis_summ['cross']['meanwt']* vis_summ['cross']['npts'] ) /               \
                                 (vis_summ['A']['npts']* vis_summ['A']['meanwt'] +                                                       \
                                  vis_summ['B']['npts']* vis_summ['B']['meanwt'] +                                                           \
                                  vis_summ['cross']['npts']* vis_summ['cross']['meanwt'])
    
    calc_rms_A = (vis_summ['A']['stddev']/np.sqrt(vis_summ['A']['npts']) )
    calc_rms_B = vis_summ['B']['stddev']/np.sqrt(vis_summ['B']['npts']) 
    calc_rms_cross = vis_summ['cross']['stddev']/np.sqrt(vis_summ['cross']['npts'])
    
    sumwt = (vis_summ['A']['npts']* vis_summ['A']['meanwt'] +                                                       \
                    vis_summ['B']['npts']* vis_summ['B']['meanwt'] +                                                           \
                    vis_summ['cross']['npts']* vis_summ['cross']['meanwt'])

    calc_rms_all = np.sqrt(( (calc_rms_A * vis_summ['A']['meanwt'] *  vis_summ['A']['npts']/sumwt)**2 +              \
                                         (calc_rms_B * vis_summ['B']['meanwt']* vis_summ['B']['npts']/sumwt)**2 +                      \
                                         (calc_rms_cross * vis_summ['cross']['meanwt']* vis_summ['cross']['npts']/sumwt)**2 )  )
  


    normn_A = vis_summ['A']['stddev']/vis_summ['A']['stddev']
    normn_B = vis_summ['B']['stddev']/vis_summ['A']['stddev']
    normn_cross = vis_summ['cross']['stddev']/vis_summ['A']['stddev']

    import pandas as pd
    from IPython.display import display, HTML
    #pd.set_option("display.precision", 3)
    results = { 'Baseline\nTypes' : [label_A+'-'+label_A,label_B+'-'+label_B , label_A+'-'+label_B, 'All'],
                'Vis Mean \n(V=I=P)' : [vis_summ['A']['mean'],vis_summ['B']['mean'], vis_summ['cross']['mean'], "-NA-" ],
           'Vis Std\n(sim)' : [vis_summ['A']['stddev'],vis_summ['B']['stddev'], vis_summ['cross']['stddev'], "-NA-" ],
                'Vis Std\n(norm sim)' : [normn_A, normn_B, normn_cross, "-NA-" ],
                'Vis Std\n(calc)' : [ 1.0, (D_A**2)/(D_B**2), (D_A**2)/(D_A*D_B)  , "-NA-" ],
           'Number of \n Data Pts' : [vis_summ['A']['npts'],vis_summ['B']['npts'],vis_summ['cross']['npts'],"-NA"],
           'Weight\n(calc)' :  [1.0, ((D_B**2)/(D_A**2))**2, ((D_A*D_B)/(D_A**2))**2 , "-NA-" ],
                'Weight\n(sim)' :[vis_summ['A']['meanwt'],vis_summ['B']['meanwt'], vis_summ['cross']['meanwt'], "-NA-" ],
           'Int Jy/bm\n(calc)' : [calc_peak_A, calc_peak_B, np.sqrt(calc_peak_A * calc_peak_B)  , calc_peak_all],
           'Int Jy/bm\n(meas)' : [meas_peaks['A'],meas_peaks['B'],meas_peaks['cross'],meas_peaks['all']],
           'RMS mJy\n(calc)' : [calc_rms_A*1e+3, calc_rms_B*1e+3, calc_rms_cross*1e+3, calc_rms_all*1e+3],
           'RMS mJy\n(meas)' : [meas_rms['A']*1e+3,meas_rms['B']*1e+3,meas_rms['cross']*1e+3,meas_rms['all']*1e+3]   }
    df = pd.DataFrame(results)
    with pd.option_context('display.float_format', '{:0.4f}'.format):
        display(HTML(df.to_html()))
        #print(pd.DataFrame(df).to_markdown())
    #print(df)
    

    print("Calculated PB size for type A (dia=%2.2f) : %3.5f arcmin"%(D_A, calc_ang(freq,D_A)))
    print("Calculated PB size for type B (dia=%2.2f) : %3.5f arcmin"%(D_B, calc_ang(freq,D_B)))
    

    return
    
    
def calc_ang(freq, dia):
    return  ( (3e+8/(freq*1e+9)) / dia ) * (180.0/3.14158) * 60.0


###################################################################################
###################################################################################
################  Interactive Displays... #####################################################
###################################################################################
###################################################################################


from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from IPython.display import display

#get_ipython().run_line_magic('matplotlib', 'inline')
from ipywidgets import interactive
import matplotlib.pyplot as plt
import numpy as np

def arrplot(iarr,pow,vran, pb, imspec, pbspec, yran):
    plt.figure(2,figsize=(10,5))
    plt.clf()
    plt.subplot(121)
    plt.imshow(np.sign(iarr) * (np.fabs(iarr)**(pow)), vmin=vran[0], vmax=vran[1], origin='lower')
    plt.contour(pb,[0.5],colors='magenta',origin='lower')
    plt.subplot(122)
    pl.plot(imspec,'bo-',label='Im')
    pl.plot(pbspec,'ro-',label='PB')
    pl.title('Spectrum at source peak')
    pl.xlabel('Channel')
    pl.ylim((yran[0],yran[1]))
    pl.legend()
    

def display_image_interact(imprefix=''):
    ia.open(imprefix+'.image')
    shp=ia.shape()
    impix = ia.getchunk()
    ia.close()
    ia.open(imprefix+'.pb')
    shp=ia.shape()
    impb = ia.getchunk()
    ia.close()

    ploc = np.where( impix == impix.max() )
    imspec = impix[ploc[0][0], ploc[1][0],0,:]
    pbspec = impb[ploc[0][0], ploc[1][0],0,:]
           
    istat = imstat(imprefix+'.residual')
    print('Residual RMS : %3.7f'%(istat['rms']))

    pow_slider = widgets.FloatSlider(
        value=1.0,
        min=0.1,
        max=2.0,
        step=0.05,
        description='Power Scale:',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
    )

    ran_slider = widgets.FloatRangeSlider(
        value=[-0.1,2.0],
        min=-0.1,
        max=1.0,
        step=0.01,
        description='Min-Max:',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
    )

    yran_slider = widgets.FloatRangeSlider(
        value=[0.6,1.1],
        min=0.0,
        max=1.5,
        step=0.01,
        description='yMin-yMax:',
        disabled=False,
        continuous_update=False,
        orientation='horizontal',
        readout=True,
        readout_format='.1f',
    )   
    

    interact(arrplot, iarr=fixed(impix[:,:,0,0].transpose()), 
                 pow=pow_slider, vran=ran_slider, 
                 pb=fixed(impb[:,:,0,0].transpose()),
                imspec=fixed(imspec), pbspec=fixed(pbspec),
                yran=yran_slider)


    #intplot
    #output = intplot.children[-1]
    #print(output.layout)
    #display(intplot)
    
    #interactive_plot = interactive(arrplot, iarr=fixed(pix), pow=1.0)
    #output = interactive_plot.children[-1]
    #output.layout.height = '350px'
    #interactive_plot
    







def plotData(msname='sim_data.ms', myplot='uv'):
    """
    Options : myplot='uv'
              myplot='data_spectrum'
              myplot='data_time'
    """
    from matplotlib.collections import LineCollection
    tb.open(msname)

    # UV coverage plot
    if myplot=='uv':
        pl.figure(figsize=(4,4))
        pl.clf()
        uvw = tb.getcol('UVW')
        pl.plot( uvw[0], uvw[1], '.')
        pl.plot( -uvw[0], -uvw[1], '.')
        pl.title('UV Coverage')
    
    # Spectrum of chosen column. Make a linecollection out of each row in the MS.
    if myplot=='data_spectrum' or myplot=='corr_spectrum' or myplot=='resdata_spectrum'  or myplot=='rescorr_spectrum' or myplot=='model_spectrum':
        dats=None
        if myplot=='data_spectrum':
            dats = tb.getcol('DATA')
        if myplot=='corr_spectrum':
            dats = tb.getcol('CORRECTED_DATA')
        if myplot=='resdata_spectrum':
            dats = tb.getcol('DATA') - tb.getcol('MODEL_DATA') 
        if myplot=='rescorr_spectrum':
            dats = tb.getcol('CORRECTED_DATA') - tb.getcol('MODEL_DATA') 
        if myplot=='model_spectrum':
            dats = tb.getcol('MODEL_DATA')
            
        xs = np.zeros((dats.shape[2],dats.shape[1]),'int')
        for chan in range(0,dats.shape[1]):
            xs[:,chan] = chan
    
        npl = dats.shape[0]
        fig, ax = pl.subplots(1,npl,figsize=(10,4))
        
        for pol in range(0,dats.shape[0]):
            x = xs
            y = np.abs(dats[pol,:,:]).T
            data = np.stack(( x,y ), axis=2)
            ax[pol].add_collection(LineCollection(data))
            ax[pol].set_title(myplot + ' \n pol '+str(pol))
            ax[pol].set_xlim(x.min(), x.max())
            ax[pol].set_ylim(y.min(), y.max())
            ax[pol].set_xlabel('Channel')
            ax[pol].set_ylabel('Vis Amp')
        pl.show()

    if myplot=="data_time":
        dats = tb.getcol('DATA')
        times = tb.getcol('TIME')
        times = times - times[0]
        pl.plot(times, np.abs(dats[0,0,:]),'.')


def msToZarr(vis='sim_data_ALMA.ms'):
    from cngi.conversion import convert_ms
    if os.path.exists(vis+'.zarr'):
        os.system('rm -rf '+vis+'.zarr')
    convert_ms(vis, vis+'.zarr')
    


def listobs_jupyter(vis='sim_data_ALMA.ms'):
    """
    Print out the contents of listobs.
    TODO : Convert the contents to a df (if possible) for a pretty display
    """

    listobs(vis=vis, listfile='obslist.txt', verbose=False, overwrite=True)
    ## print(os.popen('obslist.txt').read()) # ?permission denied?
    fp = open('obslist.txt')
    for aline in fp.readlines():
        print(aline.replace('\n',''))
    fp.close()

    tb.open(vis+'/ANTENNA')
    print("Dish diameter : " + str(tb.getcol('DISH_DIAMETER')))
    print("Antenna name : " + str(tb.getcol('NAME')))
    tb.close()



def dispAstropy(imname='sim_onepoint_true.im'):
    """
    Use working code to use AstroPy here....
    """
    exportfits(imagename=imname, fitsimage=imname+'.fits', overwrite=True)
    hdu = fits.open(imname+'.fits')[0]
    wcs = WCS(hdu.header,naxis=2)
    fig = pl.figure()
    fig.add_subplot(121, projection=wcs)
    pl.imshow(hdu.data[0,0,:,:], origin='lower', cmap=pl.cm.viridis)
    pl.xlabel('RA')
    pl.ylabel('Dec')
    


def dispImage(imname='',chan=0):
    ia.open(imname)
    pix = ia.getchunk()[:,:,0,chan]
    csys = ia.coordsys()
    ia.close()
    shp = pix.shape

    rad_to_deg =  180/np.pi
    w = WCS(naxis=2)
    w.wcs.crpix = csys.referencepixel()['numeric'][0:2]
    w.wcs.cdelt = csys.increment()['numeric'][0:2]*rad_to_deg
    w.wcs.crval = csys.referencevalue()['numeric'][0:2]*rad_to_deg
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']

    pl.subplot(projection=w)

    p1 = int(shp[0]*0.25)
    p2 = int(shp[0]*0.75)

    pl.imshow(pix[p1:p2,p1:p2].transpose(), origin='lower',  cmap=pl.cm.viridis)
    pl.xlabel('Right Ascension')
    pl.ylabel('Declination')
    

def display_image(imprefix='try_ALMA_A_single'):
    return displayImage(imname=imprefix+'.image', pbname=imprefix+'.pb', resname=imprefix+'.residual')

def displayImage(imname='try_ALMA_A_single.image', pbname='', resname=''):
    ia.open(imname)
    shp = ia.shape()
    csys = ia.coordsys()
    impix = ia.getchunk()
    ia.close()
    if pbname != '':
        ia.open(pbname)
        impb = ia.getchunk()
        ia.close()

    rad_to_deg =  180/np.pi
    w = WCS(naxis=2)
    w.wcs.crpix = csys.referencepixel()['numeric'][0:2]
    w.wcs.cdelt = csys.increment()['numeric'][0:2]*rad_to_deg
    w.wcs.crval = csys.referencevalue()['numeric'][0:2]*rad_to_deg
    w.wcs.ctype = ['RA---SIN','DEC--SIN']
    #w.wcs.ctype = ['RA','DEC']

    pl.figure(figsize=(12,5))
    pl.clf()
    pl.subplot(121, projection=w)

    p1 = int(shp[0]*0.25)
    p2 = int(shp[0]*0.75)

    pl.imshow(impix[p1:p2,p1:p2,0,0].transpose(), origin='lower')
    if pbname != '':
        pl.contour(impb[p1:p2,p1:p2,0,0].transpose(),[0.5],colors=['magenta'], origin='lower')
    pl.title('Image from channel 0')
    pl.xlabel('Right Ascension')
    pl.ylabel('Declination')
    
    
    pk = 0.0
    if shp[3]>1:
        pl.subplot(122)
        ploc = np.where( impix == impix.max() )
        pl.plot(impix[ploc[0][0], ploc[1][0],0,:],'bo-',label='Im')
        if pbname != '':
            pl.plot(impb[ploc[0][0], ploc[1][0],0,:],'ro-',label='PB')
        pl.title('Spectrum at source peak')
        pl.xlabel('Channel')
        pl.ylim((0.4,1.1))
        pl.legend()
        pk = impix[ploc[0][0], ploc[1][0],0,0]
        print('Peak Intensity (chan0) : %3.7f'%(pk))
        if pbname != '':
            pbk = impb[ploc[0][0], ploc[1][0],0,0]
            print('PB at location of Intensity peak (chan0) : %3.7f'%(pbk))

    else:
        ploc = np.where( impix == impix.max() )
        print("Image Peak : %3.4f"%(impix[ploc[0][0], ploc[1][0],0,0]))
        if pbname != '':
            print("PB Value : %3.4f"%(impb[ploc[0][0], ploc[1][0],0,0]))
        pk = impix[ploc[0][0], ploc[1][0],0,0]

    if resname !='':
        istat = imstat(resname)  ### Make this calc within the PB.
        rres = istat['rms'][0]
        print('Residual RMS : %3.7f'%(rres))
    else:
        rres = None
    
 
    return pk, rres   # Return peak intensity from channnel 0 and rms

    
###################################################


