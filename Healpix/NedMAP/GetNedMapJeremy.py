#!/usr/local/bin/python3# stars - create a simulated starfield image# Revision 2013-03-17# Copyright (c) 2011 Paul Hickson  #usage = \"""  %prog [options] RA DEC Description  This program creates an image containing a simulated star field generated from  catalog data. The first five arguments are required. Any number of optional  parameters may be entered, but none can be skipped over. The image will be  written to the image directory, specified in the cfg.py file.  Parameters  RA             right ascension of the image center (J2000)  DEC            declination of the image center (J2000)"""import pyfitsimport math as mtimport numpy as npfrom array import *from scipy.io import *import healpy as hpimport matplotlib.pyplot as pltimport sysfrom pylab import *import binasciidef IndexToDeclRa(NSIDE,index):        theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)    return -np.degrees(theta-mt.pi/2.),np.degrees(mt.pi*2.-phi)def DeclRaToIndex(decl,RA,NSIDE):    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-decl+90.),np.radians(360.-RA))#------------------------------------------------------------------------------# main#def main():    """    This is the main routine.    """        '''    # Parse the command line.    parser = OptionParser(usage=usage)    parser.disable_interspersed_args()    parser.add_option('-q','--quiet',dest='quiet',action='store_true',        default=False,help='quiet output')    parser.add_option('-v','--verbose',dest='verbose',action='store_true',        default=False,help='verbose output')    (options,args) = parser.parse_args()    nargs = len(args)    '''           filenameCat = '/Users/Elisa/c/EAntolini/Healpix/NedMAP/NED_z_0.01_0.02.txt'        RA,DEC,MAG = np.loadtxt(filenameCat,skiprows= 26,usecols = (2,3,8), delimiter = '|',dtype=[('f0',float),('f1',float),('s2','str')], unpack = True)    print(len(MAG))    #Mag = np.array_str(MAG)        '''    Magnitude = MAG.replace("b'",'')    Magnitude = Magnitude.replace("'",'')    Magnitude = Magnitude.replace("g",'0')    Magnitude = Magnitude.replace("b",'0')    Magnitude = Magnitude.replace("R",'0')    Magnitude = Magnitude.replace("[",'')    Magnitude = Magnitude.replace("]",'')    MagStar = Magnitude.split(' ')    print(len(MagStar))    for i in MagStar:        print(float(i))    '''    pix = DeclRaToIndex(DEC,RA,512)    galpixels= np.zeros(hp.nside2npix(512))    galpixels[pix]=galpixels[pix]+1    '''    Lsun = 3.846e26 # Watt    MagSun = -26.832 # Bolometric        LumGal = Lsun*np.power(10,(-MagSun-Mag)/2.5) # lUMINOSTY OF THE GALAXY IN SOLAR LUMINOSITY    galpixels[pix]+=LumGal    print(len(pix))    #print(hp.nside2npix(512))    '''    hp.mollview(galpixels,coord='C',rot = [0,0.3], title='Histogram equalized Ecliptic', unit='prob',xsize = 512)    hp.graticule()    plt.show()    ## Add some amuunt proportional to the Luminosity but not with all the galaxies -> galpixels[pix]+=LumK    ## Select a list with a range 0.01 to 0.02 in redshift -> see if there is all sky    ## with size of 512 each pixels is that big 40000/3145728 = 0.01 sq deg -> 6 arc minutes by 6 arcminutes    ## also try with all redsifth in a range of magnitude 23-23.5 and see if the map is all sky#------------------------------------------------------------------------------# Start program execution.#if __name__ == '__main__':    main()