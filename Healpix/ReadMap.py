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

def IndexToDeclRa(NSIDE,index):
    
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)

def DeclRaToIndex(decl,RA):
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
    
    Name,Morphology,Ra,Dec,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa = np.loadtxt(filename,skiprows=174,dtype=[('f0',str),('f1',str),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)], unpack = True)
    
    
    # Generate my own Map
    
    subplot(221, projection="aitoff")
    title("Aitoff")
    grid(True)
    show()
    
    # Open LIGO map
    # By default, input maps are converted to RING ordering, if they are in NESTED ordering.
    
    map_dir  = '/Users/Elisa/Documents/Robotic_Telescope/LIGO/Healpix/'
    map_name = 'bayestar.fits'
    
    #wmap_map_Nested     = hp.read_map(map_dir+map_name, nest=True) #Remains NESTED
    wmap_map_Ring       = hp.read_map(map_dir+map_name,0)            #Change  to RING (Default), read the 0 columns of the file
    print(wmap_map_Ring[0])
    hp.mollview(np.log10(wmap_map_Ring),coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=-8,max=-6, xsize=4096)
    hp.graticule()
    plt.show()
    
    '''
    #View map with unseen pixels
    mask = hp.read_map(map_dir+map_name,0).astype(np.bool)
    wmap_map_Ring_masked = hp.ma(wmap_map_Ring)
    wmap_map_Ring_masked.mask = np.logical_not(mask)
    wmap_map_Ring_masked_filled = wmap_map_Ring_masked.filled()
    hp.mollview(wmap_map_Ring_masked.filled(),coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=1e-25,max=1e-6, xsize=4096,)
    hp.graticule()
    plt.hist(wmap_map_Ring_masked.compressed(), bins = 4096)

    plt.show()
    '''

    #Get RA and DEC values from LIGO map

    LIGO_RA  = array('i')
    LIGO_DEC = array('i')

    mypixels = np.asarray(np.log10(wmap_map_Ring))
    print(mypixels[0])
    
    print(len(mypixels))

    
    for i in range(len(mypixels)):
        ra,dec = IndexToDeclRa(512,int(mypixels[i])*-1)
        
        #LIGO_RA.append(int(ra))
        #LIGO_DEC.append(int(dec))

    #print(LIGO_RA," ",LIGO_DEC)

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

