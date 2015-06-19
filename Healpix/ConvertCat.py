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
    
    #filename = '/Users/Elisa/c/EAntolini/Healpix/Catalogs/2MASSComplete_Hmag.tbl'
    filenameCat1 = '/Users/Elisa/c/EAntolini/Healpix/Catalogs/2MASS_Tully_NED_completed.txt'
    outfilename = '/Users/Elisa/c/EAntolini/Healpix/Catalogs/2MASS_Tully_NED_completed.fits'
    
    #recordtype = dtype(['|S15', '|S15',np.float,np.float,np.float,np.float,np.float,np.float,np.float])
    
    #Name,Morphology,Ra,Dec,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa = np.loadtxt(filenameCat1,skiprows=174,dtype={'names':('s0', 's1', 'f2', 'f3','f4','f5','f6','f7','f8'),'formats':('|S15', '|S15',np.float,np.float,np.float,np.float,np.float,np.float,np.float)} , unpack = True)
    
    #Name,Morphology,Ra,Dec,r_k20fe,j_m_k20fe,h_m_k20fe,k_m_k20fe,k_ba,k_pa = np.loadtxt(filename,skiprows=174,dtype={'names':('s0', 's1', 'f2', 'f3','f4','f5','f6','f7','f8','f9'),'formats':('|S15', '|S15',np.float,np.float,np.float,np.float,np.float,np.float,np.float,np.float)} , unpack = True)
    
    
    RA,DEC,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa,vel,distance = np.loadtxt(filenameCat1,dtype=[('f0',float),('f1',float),('f2',float),('f3',float),('f4',float),('f5',float),('f6',float),('f7',float),('f8',float)], unpack = True)
    
    '''
    newfile = open(outfilename, 'w')
    
    c = 0
    
    
    for i in range(len(RA)):
        
        # if   Morphology[i] != b'null' and Morphology[i] != b'PN' and Morphology[i] != b'GC' and Morphology[i] != b'SNR' :
        if   H_mag[i] != 0.0 :
            #print(str(Name[i]))
            
            #name = str(Name[i])
            #morph = str(Morphology[i])
            newfile.write(str(RA[i])+" "+str(DEC[i])+" "+str(distance[i])+" "+str(H_mag[i])+"\n")
    
        else :
            c = c+1
 
 

    print(c)
    newfile.close();



    Prova = np.ndarray.tostring(Name)
    Prova2 = Prova.decode()
    str1 = ''.join(Prova2)
    str2 = [str1.replace('\x00\x00\x00\x00\x00\x00\x00\x00', ' ') for s in str1]
    str3 = [str1.replace('\x00\x00\x00', ' ') for s in str1]
    str4 = [str1.replace('\x00\x00', ' ') for s in str1]
    str5 = [str1.replace('\x00', ' ') for s in str1]
    #str2 = [str2.replace('\x00\x00\x00', ' ') for s in str2]
    
    #Prova3 = Prova2.split('\x00\x00\x00\x00\x00\x00\x00\x00')
    str6 = ''.join(str5)
    #Prova4 = Prova3.split('\x00\x00')
    str7 = str6.split(' ')

    #index = [ s.index('') for s in str7]
    
    str8 = [x for x in str7 if x != '']
    print(type(str8))
    
    
    '''
    
    #M=[Ra,Dec,r_k20fe,j_m_k20fe,h_m_k20fe,k_m_k20fe,k_ba,k_pa]
    M=[RA,DEC,r_k20fe,j_m_k20fe,k_m_k20fe,k_ba,k_pa,vel,distance]
    
    #Create Catalog in FITS format
    #c1=pyfits.Column(name='NAME',  format='A1000', array=str8)
    #c2=pyfits.Column(name='MORPHOLOGY', format='A1000', array=Morphology)
    
    c3=pyfits.Column(name='RA',  format='E', array=RA)
    c4=pyfits.Column(name='DEC', format='E', array=DEC)
    c5=pyfits.Column(name='r_k20fe',  format='E', array = r_k20fe)
    c6=pyfits.Column(name='j_m_k20fe',  format='E', array = j_m_k20fe)
    c7=pyfits.Column(name='k_m_k20fe',  format='E', array=k_m_k20fe)
    c8=pyfits.Column(name='k_ba',  format='E', array=k_ba)
    c9=pyfits.Column(name='k_pa',  format='E', array=k_pa)
    c10=pyfits.Column(name='vel',  format='E', array=vel)
    c11=pyfits.Column(name='distance',  format='E', array=distance)
   

    cols = pyfits.ColDefs([c3, c4,c5,c6,c7,c8,c9,c10,c11])
    
    
    tbhdu = pyfits.new_table(cols)
    print(M)
    hdu = pyfits.PrimaryHDU(data=M)
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.writeto(outfilename)
    thdulist.close()


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

