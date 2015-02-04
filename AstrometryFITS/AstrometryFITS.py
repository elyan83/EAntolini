#!/opt/local/bin/python2.7

import sys
from sys import argv
import pyfits

sys.path.append('/ocs/modules')
import cfg





# File Product Directory
prod_dir = '/Users/Elisa/c/EAntolini/ProductFiles/'


# Image to Analyze - Give the full path
FileFitsImage  = argv[1]

# Create the standard Header of the fits file

hdu = pyfits.PrimaryHDU()

#Add all the informations in the Header

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


# Get the Center of the Image from the Header
hdulist = pyfits.open(FileFitsImage)
prihdr = hdulist[0].header
CenterRA  = prihdr['TCSRA']
CenterDEC = prihdr['TCSDEC']


CenterRA = (CenterRA*15) + CRVAL1
CenterDEC = CenterDEC + CRVAL2


# RA and DEC reference Points

hdu.header['CRVAL1'] = CenterRA,'RA  of reference point'
hdu.header['CRVAL2'] = CenterDEC,'DEC  of reference point'

# Pixels Reference Points

hdu.header['CRPIX1'] = cfg.ccd_field_centre[0],'X reference pixel'
hdu.header['CRPIX2'] = cfg.ccd_field_centre[1],'Y reference pixel'

# Units

hdu.header['CUNIT1'] = 'Deg','X pixel scale units '
hdu.header['CUNIT2'] = 'Deg','Y pixel scale units '

# Transformation Matrix

hdu.header['CD1_1'] = lastline[1],'Transformation matrix'
hdu.header['CD1_2'] = lastline[2]
hdu.header['CD2_1'] = lastline[4]
hdu.header['CD2_2'] = lastline[5]

# Additional Constants of the Transformation Matrix
hdu.header['CD1_3'] = lastline[3],'Additional Constants'
hdu.header['CD2_3'] = lastline[6]

#Write the FitsFile

hdulist = pyfits.HDUList([hdu])
hdulist.writeto('/Users/Elisa/c/EAntolini/AstrometryFITS/new.fits')
