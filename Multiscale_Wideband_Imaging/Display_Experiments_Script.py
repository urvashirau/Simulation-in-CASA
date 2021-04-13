#!/usr/bin/env python
# coding: utf-8


from __future__ import print_function

import os
import pylab as pl
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


# Import required tools/tasks
from casatools import simulator, image, table, coordsys, measures, componentlist, quanta, ctsys, ms
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

import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

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
    

def display_image(imprefix='try_ALMA_A_single',multiterm=False,scale=''):
    if multiterm==False:
        return displayImage(imname=imprefix+'.image', pbname=imprefix+'.pb', resname=imprefix+'.residual',scale=scale)
    else:
        displayImage(imname=imprefix+'.image.tt0',  resname=imprefix+'.residual.tt0',scale=scale)
        displayImage(imname=imprefix+'.image.tt1', scale=scale,alphaname=imprefix+'.alpha')


def displayImage(imname='try_ALMA_A_single.image', pbname='', resname='',scale='',alphaname=''):
    ia.open(imname)
    shp = ia.shape()
    csys = ia.coordsys()
    impix = ia.getchunk()
    ia.close()
#    if pbname != '':
#        ia.open(pbname)
#        impb = ia.getchunk()
#        ia.close()
    if resname != '':
        ia.open(resname)
        imres = ia.getchunk()
        ia.close()
        

    rad_to_deg =  180/np.pi
    w = WCS(naxis=2)
    w.wcs.crpix = csys.referencepixel()['numeric'][0:2]
    w.wcs.cdelt = csys.increment()['numeric'][0:2]*rad_to_deg
    w.wcs.crval = csys.referencevalue()['numeric'][0:2]*rad_to_deg
    w.wcs.ctype = ['RA---SIN','DEC--SIN']

    pl.figure(figsize=(9,4))
    pl.clf()

    p1 = int(shp[0]*0)#.25)
    p2 = int(shp[0]*1)#0.75)

    pl.subplot(121, projection=w)
    if scale=='sqrt':
        arr = impix[p1:p2,p1:p2,0,0].transpose()
        arr[arr<1e-07] = 1e-07
        pl.imshow(np.power(arr,0.1), origin='lower',interpolation='none') 
        #pl.imshow(np.sqrt(np.sqrt(impix[p1:p2,p1:p2,0,0])).transpose(), origin='lower')
    else:
        pl.imshow(impix[p1:p2,p1:p2,0,0].transpose(), origin='lower',interpolation='none')
#    if pbname != '':
#        pl.contour(impb[p1:p2,p1:p2,0,0].transpose(),[0.5],colors=['magenta'], origin='lower')
    pl.title('Image : '+imname)
    pl.xlabel('Right Ascension')
    pl.ylabel('Declination')
    
    if resname != "":
        pl.subplot(122, projection=w)
        pl.imshow(imres[p1:p2,p1:p2,0,0].transpose(), origin='lower',interpolation='none')
        pl.title('Residual : '+resname)
        pl.xlabel('Right Ascension')
        pl.ylabel('Declination')
        istat = imstat(resname)  ### Make this calc within the PB.
        rres = istat['rms'][0]
        print('Residual RMS : %3.7f'%(rres))

    if alphaname != '':
        pl.subplot(122, projection=w)
        ia.open(alphaname)
        imalpha = ia.getchunk()
        ia.close()
        pl.imshow(imalpha[p1:p2,p1:p2,0,0].transpose(), origin='lower',interpolation='none',vmin=-1.0,vmax=0.2,cmap='jet')
        pl.title('Spectral Index')
        pl.xlabel('Right Ascension')
        pl.ylabel('Declination')
    
 
    return

    


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
    




###############################################
import cngi.dio as dio


def msToZarr(vis='sim_data_ALMA.ms'):
    from cngi.conversion import convert_ms
    if os.path.exists(vis+'.zarr'):
        os.system('rm -rf '+vis+'.zarr')
    convert_ms(vis, vis+'.zarr')
    

def XPlot(vis='sim_data_ALMA.ms',ptype='amp-time',forceconvert=False):
    """
    Make a few types of plots
    Supported types : amp-time, amp-freq, uvcov, plotants
    forceconvert=True/False : Convert the input MS to a Zarr dataset and read it into an XArray for plotting. 
                                               If set to False, it will skip the conversion step (and all the output messages the conversion produces). 
    """
    zvis = vis+'.zarr'
    if not os.path.exists(zvis) or forceconvert==True:
        msToZarr(vis=vis)
        
    xdat = dio.read_vis(zvis,ddi=0)
    gxdat = dio.read_vis(zvis,ddi='global')
    xdat['DMAG'] = ((xdat['DATA'].real ** 2 + xdat['DATA'].imag ** 2) ** 0.5).mean(axis=3)
    xdat['U'] = xdat['UVW'][:,:,0]
    xdat['V'] = xdat['UVW'][:,:,1]
    xdat['-U'] = -xdat['UVW'][:,:,0]
    xdat['-V'] = -xdat['UVW'][:,:,1]
    xdat['UVDIST'] =  np.sqrt(xdat['UVW'][:,:,0]**2 + xdat['UVW'][:,:,1]**2)
    

    if ptype == 'amp-time':
        fig, axes = pl.subplots(ncols=1,figsize=(9,3))
        for fld in np.unique(xdat.field):
            xdat.where(xdat.field==fld).plot.scatter(x='time',y='DMAG',  marker='.', color='r',alpha=0.1)
            
        pl.title('Visibility ampllitude');

    if ptype == 'amp-freq':
        fig, axes = pl.subplots(ncols=2,figsize=(9,3))
        ax=0
        for fld in np.unique(xdat.field):
            xdat.where(xdat.field==fld).plot.scatter(x='chan',y='DMAG',  marker='.', color='r',alpha=0.1,ax=axes[ax])
            axes[ax].set_title('Visibility Spectrum for field : '+fld) #+ '\nRed (A-A), Blue (B-B), Purple (A-B)')
            ax = ax+1

    if ptype == 'amp-freq-uvdist':
        fig, axes = pl.subplots(ncols=2,figsize=(9,3))
        print(np.unique(xdat.field))
        for fld in np.unique(xdat.field):
            ax=0
            xdat.where(xdat.field==fld).plot.scatter(x='chan',y='DMAG',  marker='.', color='r',alpha=0.1,ax=axes[ax])
            axes[ax].set_title('Visibility Spectrum for field : '+fld) #+ '\nRed (A-A), Blue (B-B), Purple (A-B)')
            ax=1
            xdat.where(xdat.field==fld).plot.scatter(x='UVDIST',y='DMAG',  marker='.', color='b',alpha=0.1,ax=axes[ax])
            axes[ax].set_title('Amp vs UVdist for field : '+fld) #+ '\nRed (A-A), Blue (B-B), Purple (A-B)')



            
    if ptype == 'uvcov':
        fig, axes = pl.subplots(ncols=2,figsize=(9,4))
        ax=0
        for fld in np.unique(xdat.field):
            xdat.where(xdat.field==fld).plot.scatter(x='U',y='V',  marker='.', color='m',ax=axes[ax])
            xdat.where(xdat.field==fld).plot.scatter(x='-U',y='-V',  marker='.', color='m',ax=axes[ax])
            axes[ax].set_title('UV coverage for field : '+fld)
            ax=ax+1
 

    if ptype == 'plotants':
        if vis.count('ALMA')>0:
            typeA = 'A'
        else:
            typeA = 'm'
        fig, axes = pl.subplots(ncols=1,figsize=(6,5))
        gxdat['ANT_XPOS'] = gxdat['ANT_POSITION'][:,0] - gxdat['ANT_POSITION'][:,0].mean()
        gxdat['ANT_YPOS'] = gxdat['ANT_POSITION'][:,1] - gxdat['ANT_POSITION'][:,1].mean()
        gxdat.plot.scatter(x='ANT_XPOS', y='ANT_YPOS',color='k',marker="1",s=200,linewidth=3.0)
        for i, txt in enumerate(gxdat.ANT_NAME.values):
            col = ('b' if (txt.count(typeA)>0) else 'r')
            pl.annotate('   '+txt, (gxdat['ANT_XPOS'].values[i], gxdat['ANT_YPOS'].values[i]),fontsize=12,color=col)   
        pl.title('Antenna Positions');

