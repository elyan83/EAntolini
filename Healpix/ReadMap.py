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
import matplotlib.pyplot as plt
import sys
from pylab import *



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
    
    # Generate my own Map
    
    subplot(221, projection="aitoff")
    title("Aitoff")
    grid(True)
    show()
    
    # Open LIGO map
    # By default, input maps are converted to RING ordering, if they are in NESTED ordering.
    
    map_dir  = '/Users/Elisa/Documents/Robotic_Telescope/LIGO/Healpix/'
    map_name = 'bayestar.fits'
    
    # wmap_map_Nested     = hp.read_map(map_dir+map_name, nest=True) #Remains NESTED
    wmap_map_Ring       = hp.read_map(map_dir+map_name,0)            #Change  to RING (Default), read the 0 columns of the file
    print(wmap_map_Ring)
    hp.mollview(wmap_map_Ring,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=1e-25,max=1e-6, xsize=4096)
    hp.graticule()
    plt.show()

    #View map with unseen pixels
    mask = hp.read_map(map_dir+map_name,0).astype(np.bool)
    wmap_map_Ring_masked = hp.ma(wmap_map_Ring)
    wmap_map_Ring_masked.mask = np.logical_not(mask)
    wmap_map_Ring_masked_filled = wmap_map_Ring_masked.filled()
    hp.mollview(wmap_map_Ring_masked.filled(),coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob', min=1e-25,max=1e-6, xsize=4096)
    hp.graticule()
    plt.hist(wmap_map_Ring_masked.compressed(), bins = 4096)

    plt.show()



    ''' This doesn't work'''
    '''
    nside = 512
    npix = hp.nside2npix(nside)
    pix = np.arange(npix)
    t,p = hp.pix2ang(nside,pix) #theta, phi

    r = hp.Rotator(deg=True, rot=[90, 30])

    map_rot = np.zeros(npix)

    for i in pix:
        trot, prot = r(t[i],p[i])
        tpix = int(trot*180./np.pi) #my data came in a theta, phi grid -- this finds its location there
        ppix = int(prot*180./np.pi)
        map_rot[i] = wmap_map_Ring_masked_filled[ppix,tpix] #this being the right way round may need double-checking


    map_rot = hp.mollview(wmap_map_Ring_masked.filled(),deg=True,rot=[90,30], return_projected_map=True)
    #hp.mollview( map_rot,coord='C', title='Histogram equalized Ecliptic', unit='prob', min=1e-25,max=1e-6, xsize=4096)

    '''

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

