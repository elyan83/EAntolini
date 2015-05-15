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
    
    filename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSource.tbl'
    outfilename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSource.fits'
    
    Name,Morphology,GAL_RA,GAL_DEC,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa = np.loadtxt(filename,skiprows=174,dtype=[('f0',str),('f1',str),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)], unpack = True)
    
    
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

    '''
    for i in range(len(mypixels)):
        dec, ra = IndexToDeclRa(512,i)
        LIGO_RA.append(ra)
        LIGO_DEC.append(dec)
    

    with  open('/Users/Elisa/c/EAntolini/Healpix/LigoRADEC.txt', 'w') as fLigofile:
        for i in range(len(LIGO_RA)):
            fLigofile.write(str(LIGO_RA[i])+" "+str(LIGO_DEC[i])+"\n")
    '''
    
    #Postpone this when run again the for loop
    LIGO_RA,LIGO_DEC = np.loadtxt('/Users/Elisa/c/EAntolini/Healpix/LigoRADEC.txt',dtype=[('f0',float),('f1',float)], unpack = True)


    pos    = (mt.pi/180.0)

    #2) Take RA and DEC from GALAXY Catalog and convert to Index
    for r, d, radius in zip(GAL_RA[Name=='M31'],GAL_DEC[Name=='M31'],r_k20fe[Name=='M31']):
    #for r, d in zip(GAL_RA,GAL_DEC):
        #Add r_k20fe (arcsec)-> size of the galaxy -> is big if the galaxy is big -> 100,400 big numbers -> M31 600 arcsec -> for now leave in arcsec
        #dist = (r-LIGO_RA)**2+(d-LIGO_DEC)**2 #->
        #pos[0] = (mt.pi/180.0)*GAL_RA;
        #pos[1] =(mt.pi/180.0)*GAL_DEC;
        

        #np.cos(pos*r)
        #np.cos((r-LIGO_RA)*pos)
        dumy=np.arccos(np.cos(d*pos)*cosdec_c*np.cos((r-LIGO_RA)*pos)+np.sin(d*pos)*sindec_c)
        #sindec_c= np.sin((LIGO_DEC)*mt.pi/180.0)
        #dist = (((r-LIGO_RA)**2+(d-LIGO_DEC)**2)/(radius**2))*1e4
        #galpixels +=np.exp(-dist)
        galpixels +=np.exp(-dumy/radius)
    

    
    '''
    cosdec_c= np.cos((LIGO_DEC)*mt.pi/180.0)
    sindec_c= np.sin((LIGO_DEC)*mt.pi/180.0)
    cosposra_c = np.cos(pos - (LIGO_RA)*mt.pi/180.0)
    
    # distance in radiance of the galaxy from the center
    dumy=np.arccos(np.cos(pos)*cosdec_c*cosposra_c+np.sin(pos)*sindec_c);
    #dumx=np.atan2(sin(pos[1])-cos(dumy)*sin(dec_c),cos(pos[1])*sin(pos[0]-ra_c)*cos(dec_c));
        
    for r, d, radius in zip(GAL_RA[Name=='M31'],GAL_DEC[Name=='M31'],r_k20fe[Name=='M31']):
        
        #dist = (((r-LIGO_RA)**2+(d-LIGO_DEC)**2)/(radius**2))*1e4
        galpixels +=np.exp(-dumy/radius)
    

    #print(galpixels)
    '''

    hp.mollview(galpixels,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', xsize=4096)
    hp.graticule()
    plt.show()
    #savefig('/Users/Elisa/c/EAntolini/Healpix/GalMap.png')


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

