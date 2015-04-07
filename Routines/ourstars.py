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

import math
import random
import sys
from optparse import OptionParser

sys.path.append('/ocs/modules')
import astro
import cfg
import image
import util
import pyfits

#------------------------------------------------------------------------------
# main
#
def main():
    """
    This is the main routine.
    """
    
    
    # Parse the command line.
    parser = OptionParser(usage=usage)
    parser.disable_interspersed_args()
    parser.add_option('-q','--quiet',dest='quiet',action='store_true',
        default=False,help='quiet output')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',
        default=False,help='verbose output')
    (options,args) = parser.parse_args()
    nargs = len(args)

    nx = cfg.ccd_nx         # number of columns in simulated image
    ny = cfg.ccd_ny         # number of rows in simulated image
    binf = 1                # binning factor
    fwhm = 1
    resolution = 5          # Number of images per FWHM for star trails
    band = 'V'              # Default band
    ra_rate = 0             # Default is untrailed image
    dec_rate = 0
    track_length = 0
    dx = 0
    dy = 0
    
    if nargs < 3:
        parser.print_help()

    else:
        ra = args[0]
        dec = args[1]
        filename = args[2]
        filedir  = args[3]
        binf=1
        nadd = int(resolution*track_length/(binf*fwhm))+1      # no of images in track
        nadd2 = nadd/2

        nx = int(nx/binf)
        print(nx)
        ny = int(ny/binf)
        print(ny)
        
        # Determine the ra and dec range of image
        ra0 = astro.deg(ra)
        #print(ra0*3600)
        dec0 = astro.deg(dec)
        #print(dec0*3600)
        cd = math.cos(math.radians(dec0))
        #ra_min = ra0-(cd*nx+dx)*binf*cfg.ccd_pixel_size/108000
        ra_min = ra0-(cd*nx+dx)*binf*cfg.ccd_pixel_size/58800
        #ra_max = ra0+(cd*nx+dx)*binf*cfg.ccd_pixel_size/108000
        ra_max = ra0+(cd*nx+dx)*binf*cfg.ccd_pixel_size/58800
        if (ra_min < 0): ra_min += 24
        if (ra_max >= 24): ra_max -= 24
        #print(ra_min*3600)
        #print(ra_max*3600)
        dec_min = dec0-(ny+dy)*binf*cfg.ccd_pixel_size/4050
        dec_max = dec0+(ny+dy)*binf*cfg.ccd_pixel_size/4050
        #print(dec_min*3600)
        #print(dec_max*3600)
        
        print("Catalog Size arcsec = ra : ",(ra_max-ra_min)*3600, "Dec : ",(dec_max-dec_min)*3600)

        # Determine x and y increments for each image
        dx /= nadd
        dy /= nadd

        # Add stars. The origin is in the SE corner of the image. Dec is
        # proportional to y and ra is proportional to -x.
        # First get a list of stars from the USNO catalog.
        stars = astro.usno(cfg.catalog_dir+'usno185r',ra_min,ra_max,dec_min,dec_max)
        
        
        # Print out the Catalog
        #fileout= open(cfg.catalog_dir+filename, 'a+')
        fileout= open(filedir+filename, 'a+')
        
        if not stars:
            print('no stars in range',ra_min,ra_max,dec_min,dec_max)
            exit(0)
        else:
            print(len(stars),'stars')
            fileout.write(str(len(stars))+"\n")

        for star in stars:
            if band.lower() == 'b':
                mag = star.b
            elif band.lower() == 'v':
                mag = 0.5*star.b+0.5*star.r
            elif band.lower() == 'r':
                mag = star.r
            flux = 10.0**(-0.4*(mag-20))
            x = nx/2-(star.ra-ra0)*cd*54000/(binf*cfg.ccd_pixel_size)
            y = ny/2+(star.dec-dec0)*3600/(binf*cfg.ccd_pixel_size)
            
            #print('{:4.1f} {:6.1f} {:6.1f} {:6.8f} {:6.8f}'.format(flux,x,y,star.ra,star.dec))
            fileout.write('{:4.1f} {:s} {:6.8f} {:s} {:6.8f} {:s}'.format(flux," ",(star.ra)*15," ",star.dec,"\n"))



    


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

