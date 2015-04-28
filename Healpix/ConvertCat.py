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
import numpy as np
from array import *
from scipy.io import *




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
    
    filename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSourceRaDecNxNy.txt'
    outfilename = '/Users/Elisa/c/EAntolini/Healpix/IpacTableFromSourceRaDecNxNy.fits'

    data1 = []
    data2 = []
    data3 = []
    data4 = []
    
    #Open Catalog

    file = open(filename,'r')
    lines = file.readlines()



    for line in lines:
        p = line.split()
        data1.append(float(p[0]))
        data2.append(float(p[1]))
        data3.append(float(p[2]))
        data4.append(float(p[3]))
    
    
    M=[data1,data2,data3,data4]

    #Create Catalog in FITS format
    
    c1=pyfits.Column(name='RA',  format='E', array=data1)
    c2=pyfits.Column(name='DEC', format='E', array=data2)
    c3=pyfits.Column(name='nx',  format='E', array=data3)
    c4=pyfits.Column(name='ny',  format='E', array=data4)

    cols = pyfits.ColDefs([c1, c2, c3, c4])

    tbhdu = pyfits.new_table(cols)
    hdu = pyfits.PrimaryHDU(data=M)
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.writeto(outfilename)
    thdulist.close()

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

