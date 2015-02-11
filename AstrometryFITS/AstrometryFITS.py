#!/opt/local/bin/python2.7

import os
import sys
from sys import argv
import pyfits
import time
from optparse import OptionParser

sys.path.append('/ocs/modules')
import cfg


#------------------------------------------------------------------------------
# Implementation of Functions
#

def CreateNewImage(Self,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Date):
    
    #EXTEND Generates an error
    del Self['EXTEND']
    
 
    # RA and DEC reference Points
 
    Self.append(('CRVAL1', cenVAL1, 'RA  of reference point'),end=True)
    Self.append(('CRVAL2', cenVAL2, 'DEC  of reference point'),end=True)
 
    # Pixels Reference Points
 
    Self.append(('CRPIX1', cenPIX1, 'X reference pixel'),end=True)
    Self.append(('CRPIX2', cenPIX2, 'Y reference pixel'),end=True)
 
    # Units
 
    Self.append(('CUNIT1', 'deg','X pixel scale units '),end=True)
    Self.append(('CUNIT2', 'deg','Y pixel scale units '),end=True)
 
    # Transformation Matrix
 
    Self.append(('CD1_1',a,'Transformation matrix'),end=True)
    Self.append(('CD1_2',b,' '),end=True)
    Self.append(('CD2_1',d,' '),end=True)
    Self.append(('CD2_2',e,' '),end=True)
 
    # Additional Constants of the Transformation Matrix
 
    Self.append(('CD1_3',c,'Additional Constants'),end=True)
    Self.append(('CD2_3',f,' '),end=True)
 
    # Date
 
    Self.append(('DATE', Date,'Date this file was created.'),end=True)


    return Self

def CreateWCS(Self,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Date):

    # RA and DEC reference Points

    Self.header['CRVAL1'] = cenVAL1,'RA  of reference point'
    Self.header['CRVAL2'] = cenVAL2,'DEC  of reference point'

    # Pixels Reference Points

    Self.header['CRPIX1'] = cenPIX1,'X reference pixel'
    Self.header['CRPIX2'] = cenPIX2,'Y reference pixel'

    # Units

    Self.header['CUNIT1'] = 'deg','X pixel scale units '
    Self.header['CUNIT2'] = 'deg','Y pixel scale units '

    # Transformation Matrix

    Self.header['CD1_1'] = a,'Transformation matrix'
    Self.header['CD1_2'] = b,' '
    Self.header['CD2_1'] = d,' '
    Self.header['CD2_2'] = e,' '

    # Additional Constants of the Transformation Matrix
    Self.header['CD1_3'] = c,'Additional Constants'
    Self.header['CD2_3'] = f,' '

    # Date
    Self.header['DATE'] = Date,'Date this file was created.'

    return Self

#------------------------------------------------------------------------------
# main
#

def main():
    
    """
    This is the main routine.
    """
    
    # Parse the command line.
    parser = OptionParser()
    parser.disable_interspersed_args()
    parser.add_option('-q','--quiet',dest='quiet',action='store_true',default=False,help='quiet output')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',default=False,help='verbose output')
    (options,args) = parser.parse_args()
    nargs = len(args)
    
    if nargs < 1:
        parser.print_help()
    
    else:

        # File Product Directory
        prod_dir = '/Users/Elisa/c/EAntolini/ProductFiles/'


        # Image to Analyze - Give the full path
        FileFitsImage  = argv[1]


        # Coordinates in degrees of Reference Point
        # Get transofrmation Matrix Parameters


        with open (prod_dir+'TriangleOut.txt',"r" ) as fileHandle:
            lineList = fileHandle.readlines()
        lastline = lineList[len(lineList)-1]
        lastline = lastline.split(" ")

        for i in range(1,7):
            lastline[i] = float(lastline[i])

        # Transform Reference Point from ICS to CCS

        CRVAL1 =  (cfg.ccd_field_centre[0]*lastline[1] + cfg.ccd_field_centre[1]*lastline[2] + lastline[3])/3600
        CRVAL2 =  (cfg.ccd_field_centre[0]*lastline[4] + cfg.ccd_field_centre[1]*lastline[5] + lastline[6])/3600



        # Access to the Header of the FITS Image
        
        
        hdulist = pyfits.open(FileFitsImage)
        prihdr  = hdulist[0].header
        

        # Get the Center of the Image from the Header and add the Reference Point (deg)


        CenterRA  = prihdr['TCSRA']
        CenterDEC = prihdr['TCSDEC']

        CenterRA = (CenterRA*15) + CRVAL1
        CenterDEC = CenterDEC + CRVAL2


        # Date

        tm = time.gmtime()
        date = str(tm[0])+"-"+str(tm[1])+"-"+str(tm[2])+"T"+str(tm[3])+":"+str(tm[4])+":"+str(tm[5])


        #Create NewImage File from the Header of the Image
        
        hdrNewImage = CreateNewImage(prihdr,CenterRA,CenterDEC,cfg.ccd_field_centre[0],cfg.ccd_field_centre[1],lastline[0],lastline[1],lastline[2],lastline[3],lastline[4],lastline[5],date)
        
        # Create WCS File from scratch
        # Create the standard Header of the fits file

        hdu = pyfits.PrimaryHDU()
        hdrWCS = CreateWCS(hdu,CenterRA,CenterDEC,cfg.ccd_field_centre[0],cfg.ccd_field_centre[1],lastline[0],lastline[1],lastline[2],lastline[3],lastline[4],lastline[5],date)


        #Write the FitsFiles
        
        
        if os.path.exists(prod_dir+'new-image.fits') == False :
            
            hdulist.writeto(prod_dir+'new-image.fits')
            print(prod_dir+'new-image.fits' + " has been created"+"\n")
        
        else :
        
            print(prod_dir+'new-image.fits' + " already exists"+"\n")
        
        
        if os.path.exists(prod_dir+'wcs.fits') == False :
            
            hdrWCS.writeto(prod_dir+'wcs.fits')
            print(prod_dir+'wcs.fits' + " has been created"+"\n")

        else :

            print(prod_dir+'wcs.fits' + " already exists"+"\n")


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

