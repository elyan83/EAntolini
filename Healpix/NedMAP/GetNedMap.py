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

import pyfits
import math as mt
import numpy as np
from array import *
from scipy.io import *
import healpy as hp
import matplotlib.pyplot as plt
import sys
from pylab import *
import binascii


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
    
   
    filenameCat = '/Users/Elisa/c/EAntolini/Healpix/NedMAP/NED_z_0.1_0.11.txt'
    
    RA,DEC = np.loadtxt(filenameCat,skiprows= 26,usecols = (2,3), delimiter = '|',dtype=[('f0',float),('f1',float)], unpack = True)


    pix = DeclRaToIndex(DEC,RA,512)
    galpixels= np.zeros(hp.nside2npix(512))
    galpixels[pix]=galpixels[pix]+1
    
    print(len(pix))

    print(hp.npix2nside(len(pix)))
    #print(hp.nside2npix(64))
    
    hp.mollview(galpixels,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob',xsize = 64)
    hp.graticule()
    plt.show()


    # Generate my own Map
    '''
    subplot(221, projection="aitoff")
    title("Aitoff")
    grid(True)
    show()
    '''



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

