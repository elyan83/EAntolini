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
    
    filename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSource.tbl'
    outfilename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSource.fits'
    
    Name,Morphology,Ra,Dec,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa = np.loadtxt(filename,skiprows=174,dtype=[('f0',str),('f1',str),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)], unpack = True)
    
    print(hp.pix2ang(16, 1440))
    
    print(IndexToDeclRa(16,1440))
    
    #print(Name) -> Doesn't work
    '''
    M=[Ra,Dec,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa]
    
    #Create Catalog in FITS format
    #c1=pyfits.Column(name='NAME',  format='E', array=Name)
    #c2=pyfits.Column(name='MORPHOLOGY', format='E', array=Morphology)
    c3=pyfits.Column(name='RA',  format='E', array=Ra)
    c4=pyfits.Column(name='DEC', format='E', array=Dec)
    c5=pyfits.Column(name='r_k20fe',  format='E', array=r_k20fe)
    c6=pyfits.Column(name='j_m_k20fe',  format='E', array=j_m_k20fe)
    c7=pyfits.Column(name='k_m_k20fe',  format='E', array=k_m_k20fe)
    c8=pyfits.Column(name='k_ba',  format='E', array=k_ba)
    c9=pyfits.Column(name='k_pa',  format='E', array=k_pa)

    cols = pyfits.ColDefs([c3, c4,c5,c6,c7,c8,c9])

    tbhdu = pyfits.new_table(cols)
    hdu = pyfits.PrimaryHDU(data=M)
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.writeto(outfilename)
    thdulist.close()
    '''

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

