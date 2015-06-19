#!/usr/local/bin/python3
# stars - create a simulated starfield image
# Revision 2013-03-17
# Copyright (c) 2011 Paul Hickson  

usage = \
"""
  %prog [options] RA DEC 

Description
  This program creates an image containing a simulated star field generated from
  catalog data. The first five arguments are required. Any number of optional
  parameters may be entered, but none can be skipped over. The image will be
  written to the image directory, specified in the cfg.py file.
  
Parameters
  RA             right ascension of the image center (J2000)
  DEC            declination of the image center (J2000)
"""

import numpy as np
import healpy as hp
import math as mt
import matplotlib.pyplot as plt
import sys
from pylab import *
from array import *
import pyfits

def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)

def DeclRaToIndex(decl,RA,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))



#------------------------------------------------------------------------------
# main
#
def main():
    """
    This is the main routine.
    """
    
    '''
    # Parse the command line.
    parser = OptionParser(usage=usage)
    parser.disable_interspersed_args()
    parser.add_option('-q','--quiet',dest='quiet',action='store_true',
        default=False,help='quiet output')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',
        default=False,help='verbose output')
    (options,args) = parser.parse_args()
    nargs = len(args)
    '''
    
    # Read Galaxy Catalog Parameters
    
    filenameCat1 = '/Users/Elisa/c/EAntolini/Healpix/Catalogs/2MASS_Tully_NED_completed.txt'
    outfilename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSource.fits'

    #Name,Morphology,GAL_RA,GAL_DEC,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa,distance = np.loadtxt(filenameCat1,dtype=[('f0',str),('f1',str),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float),('f9',float)], unpack = True)
    
    GAL_RA,GAL_DEC,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa, vel, distance = np.loadtxt(filenameCat1,dtype=[('f0',float),('f1',float),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)], unpack = True)
    

    Hubble_Constant = 75.0 #[km/s/Mpc]
    
    # Generate my own Map
    '''
    subplot(221, projection="aitoff")
    title("Aitoff")
    grid(True)
    show()
    '''
    
    # Open LIGO map
    # By default, input maps are converted to RING ordering, if they are in NESTED ordering.
    
    map_dir  = '/Users/Elisa/Documents/Robotic_Telescope/LIGO/Healpix/'
    map_name = 'bayestar.fits'
    
    #wmap_map_Nested     = hp.read_map(map_dir+map_name, nest=True) #Remains NESTED
    wmap_map_Ring       = hp.read_map(map_dir+map_name,0)            #Change  to RING (Default), read the 0 columns of the file
    hp.mollview(np.log10(wmap_map_Ring),coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=-8,max=-6, xsize=4096)
    hp.graticule()
    plt.show()
    

    #Get RA and DEC values from LIGO map


    mypixels = np.asarray(np.log10(wmap_map_Ring))
    galpixels=0*mypixels
 
    
    print(len(mypixels))
    
    
    
    #2)
    #LIGO_RA=[]
    #LIGO_DEC=[]
    
    # Convert Pixels to RA and DEC in LIGO map

    '''
    for i in range(len(mypixels)):
        dec, ra = IndexToDeclRa(512,i)
        LIGO_RA.append(ra)
        LIGO_DEC.append(dec)
    

    with  open('/Users/Elisa/c/EAntolini/Healpix/LigoRADEC.txt', 'w') as fLigofile:
        for i in range(len(LIGO_RA)):
            fLigofile.write(str(LIGO_RA[i])+" "+str(LIGO_DEC[i])+"\n")
    '''
    
    #Load File with M31 and Andromeda (LIGO data set)
    
    LIGO_RA,LIGO_DEC = np.loadtxt('/Users/Elisa/c/EAntolini/Healpix/LigoRADEC.txt',dtype=[('f0',float),('f1',float)], unpack = True)


    pos    = (mt.pi/180.0)
    arcsec_to_radians = 4.84813681e-06
    cosdec_c= np.cos((LIGO_DEC)*pos) #cos(dec_c) radians
    sindec_c= np.sin((LIGO_DEC)*pos) #sin(dec_c) radians


    #2) Take RA and DEC from GALAXY Catalog and convert to Index
    
    #for r, d, radius in zip(GAL_RA[Name=='M31'],GAL_DEC[Name=='M31'],r_k20fe[Name=='M31']):
    
    #for r, d,radius, semi_mayor,polar_angle in zip(GAL_RA,GAL_DEC,r_k20fe,k_ba,k_pa):
    for r, d,semi_mayor,K_mag,ba,polar_angle,dist in zip(GAL_RA,GAL_DEC,r_k20fe,k_m_k20fe,k_ba,k_pa,distance):
    
        # Distance of the galaxy from the center [radians]
        dumy=np.arccos(np.cos(d*pos)*cosdec_c*np.cos((r-LIGO_RA)*pos)+np.sin(d*pos)*sindec_c)
        
        
        # Polar Angle (between North-South directions) [radians]
        dumx=np.arctan2(np.sin(d*pos)-np.cos(dumy)*sindec_c,np.cos(d*pos)*np.sin((r-LIGO_RA)*pos)*cosdec_c);
        
        dumx +=(polar_angle+90)*pos
        
        '''
        dumx -=(polar_angle+90)*pos
        dumx +=(90-polar_angle)*pos
        dumx -=(90-polar_angle)*pos
        '''
        
        semi_minor=ba*semi_mayor
        
        
        #Compute the semi-minor axes of the Glaxy from Catalog
  
    
        f_dumx = (semi_mayor * semi_minor)/np.sqrt(np.square(semi_minor*np.cos(dumx))+np.square(semi_mayor*np.sin(dumx)))
        
        
        LumK = np.power(10,(-0.4*(K_mag-5*np.log10(dist*1e5)-6.35))) # lUMINOSTY OF THE GALAXY IN SOLAR LUMINOSITY
        
        radius = f_dumx*pos/3600
    
        #Utilizzare dumy se non ho tutte le distanze delle galassie
        # Quando ho tutte le galassie posso usare Lumk
    
        #galpixels += np.exp(-dumy/radius)
    
        galpixels += (LumK/(semi_mayor * semi_minor))*np.exp(-dumy/radius)
        #galpixels += (1/(semi_mayor * semi_minor))*np.exp(-dumy/radius)
    



    hp.mollview(np.log10(galpixels),coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=4096)
    #hp.mollview(galpixels,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=4096)
    hp.graticule()
    plt.show()



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

