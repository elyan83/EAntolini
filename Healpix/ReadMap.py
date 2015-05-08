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
    
    '''
    #1) Take RA and DEC from LIGO and convert to Index
    for r, d in zip(Ra,Dec):
        i=DeclRaToIndex(d,r,512)
        galpixels[i]=galpixels[i]+1
    
    #1) Take only non zero values
    
    with  open('/Users/Elisa/c/EAntolini/Healpix/GalPixels1.txt', 'w') as fGalfile:
        for i in range(len(galpixels)):
            if galpixels[i] != 0 :
                fGalfile.write(str(i)+" "+str(galpixels[i])+"\n")
    
    '''
    
    '''
    #2)
    LIGO_RA=[]
    LIGO_DEC=[]

    
    for i in range(len(mypixels)):
        ra, dec = IndexToDeclRa(512,i)
        LIGO_RA.append(ra)
        LIGO_DEC.append(dec)
    

    
    #2) Take RA and DEC from GALAXY Catalog and convert to Index
    for r, d in zip(GAL_RA,GAL_DEC):
        i=DeclRaToIndex(d,r,512)
        dist = (r-LIGO_RA)**2+(d-LIGO_DEC)**2
        galpixels +=np.exp(-dist)


    '''
    
    
    #3)
    for i in range(len(mypixels)):
        ra,dec = IndexToDeclRa(512,i)
        dist=(GAL_RA-ra)**2+(GAL_DEC-dec)**2 #-> lenght of this is the number of galaxies I have
        galpixels[i]=np.sum(np.exp(-dist)) #-> draw a little circle
    
    '''
    with  open('/Users/Elisa/c/EAntolini/Healpix/GalPixels2.txt', 'w') as fGalfile:
        for g in galpixels:
            fGalfile.write(str(g)+"\n")
    
    '''

    hp.mollview(galpixels,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=0,max=2, xsize=4096)
    hp.graticule()
    plt.show()
    #savefig('/Users/Elisa/c/EAntolini/Healpix/GalMap.png')


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

