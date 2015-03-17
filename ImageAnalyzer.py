#!/opt/local/bin/python2.7

#usage = \



import os
import sys
import pyfits
from pylab import *
from array import *
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import math
# import fit
# import statistics as stat
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import time
from optparse import OptionParser
from subprocess import call
import subprocess
sys.path.append('/ocs/modules')
import cfg


#------------------------------------------------------------------------------
# Implementation of Functions
#


def bubble_sort_cat(DistX,DistY,Flux,DistRA,DistDEC):
    DistX, DistY, Flux, DistRA, DistDEC = zip(*sorted(zip(DistX, DistY, Flux, DistRA, DistDEC), reverse=True, key=lambda x: x[2]))
    
    return DistX, DistY, Flux, DistRA, DistDEC

def bubble_sort_image(DistRa,DistDec,Allpix,PeakVal,Size,Stat1,Stat2,Npixel):
    
    for i in range(len(Allpix)):
        for j in range(len(Allpix)-1-i):
            if Allpix[j] < Allpix[j+1]:
                Allpix[j], Allpix[j+1]   = Allpix[j+1], Allpix[j]
                PeakVal[j],PeakVal[j+1]  = PeakVal[j+1],PeakVal[j]
                DistRa[j], DistRa[j+1]   = DistRa[j+1], DistRa[j]
                DistDec[j], DistDec[j+1] = DistDec[j+1],DistDec[j]
                Size[j], Size[j+1] = Size[j+1],Size[j]
                Stat1[j], Stat1[j+1] = Stat1[j+1],Stat1[j]
                Stat2[j], Stat2[j+1] = Stat2[j+1],Stat2[j]
                Npixel[j], Npixel[j+1] = Npixel[j+1],Npixel[j]



def ExecuteCommand(Filename1,Filename2,command,out):
    
    if Filename2 == None :
        
        if os.path.exists(Filename1) == False :
        
            if out == False :
                call(command)
                Exec = True
                print(Filename1 + " has been created"+"\n")
        
            else :
                print(os.popen(command).read())
                Exec = True
                print(Filename1 + " has been created"+"\n")
        else:
        
            Exec = False
            print(Filename1 + " already exists"+"\n")
            print('\n')

    else :
         if os.path.exists(Filename1) == False or os.path.exists(Filename2) == False:
             
             if out == False :
                 call(command)
                 Exec = True
             
             else :
                 print(os.popen(command).read())
                 Exec = True
         else:
                         
            Exec = False
            print(Filename1+" and "+Filename2+ " already exists"+"\n")



    return Exec



def ExecuteCommandFileOut(Filename,command):
    
    
    if os.path.exists(Filename) == False :
        
        with open(Filename, 'w') as myoutfile:
            call(command, stdout = myoutfile)

        Exec = True
        print(Filename+" has been created "+"\n")

    
    else:
        
        Exec = False
        print(Filename + " already exists"+"\n")
    
    return Exec





def CreateWCS(Self,Dir,Filename,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Nx,Ny,Aorder,Z,Date):
    
    Self.header['WCSAXES'] = 2,' '
    
    Self.header['CTYPE1'] = 'RA---TAN-SIP','TAN (gnomic) projection + SIP distortions'
    Self.header['CTYPE2'] = 'DEC---TAN-SIP','TAN (gnomic) projection + SIP distortions'
    Self.header['EQUINOX'] = 2000.0,'Equatorial coordinates definition (yr)'
    
    Self.header['LONPOLE'] = 180.0,' '
    Self.header['LATPOLE'] = 0.0,' '
    
    
    
    
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
    #Self.header['CD1_3'] = c,'Additional Constants'
    #Self.header['CD2_3'] = f,' '
    
    Self.header['IMAGEW'] = Nx,'Image width,  in pixels'
    Self.header['IMAGEH'] = Ny,'Image height, in pixels.'
    
    # Polynomial Distortions
    
    #Self.header['A_ORDER'] = Aorder,'Polynomial order, axis 1'
    #Self.header['A_0']     = Z[0],' '
    #Self.header['A_1']     = Z[1],' '
    #Self.header['A_2']     = Z[2],' '
    
    # Date
    Self.header['DATE'] = Date,'Date this file was created.'
    
    if os.path.exists(Dir+Filename) == False :
        
        Self.writeto(Dir+Filename)
        
        print(Dir+Filename + " has been created"+"\n")
            
    else :
                
        print(Dir+Filename + " already exists"+"\n")
    
    return Self

def CreateNewImage(Self,hlist,Dir,Filename,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Nx,Ny,Aorder,Z,Date):
    
    #EXTEND Generates an error
    del Self['EXTEND']
    
    Self.append(('WCSAXES',2,' '),end=True)
    
    Self.append(('CTYPE1','RA---TAN-SIP','TAN (gnomic) projection + SIP distortions'),end=True)
    Self.append(('CTYPE2','DEC---TAN-SIP','TAN (gnomic) projection + SIP distortions'),end=True)
    Self.append(('EQUINOX', 2000.0,'Equatorial coordinates definition (yr)'),end=True)
    
    Self.append(('LONPOLE',180.0,' '),end=True)
    Self.append(('LATPOLE',0.0,' '),end=True)
    

    
    
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
    
    #Self.append(('CD1_3',c,'Additional Constants'),end=True)
    #Self.append(('CD2_3',f,' '),end=True)
                
    Self.append(('IMAGEW',Nx,'Image width,  in pixels'),end=True)
    Self.append(('IMAGEH',Ny,'Image height, in pixels.'),end=True)
    
    # Polynomial Distortions
    
    #Self.append(('A_ORDER',Aorder,'Polynomial order, axis 1'),end=True)
    #Self.append(('A_0',Z[0],' '),end=True)
    #Self.append(('A_1',Z[1],' '),end=True)
    #Self.append(('A_2',Z[2],' '),end=True)
    
    # Date
    
    Self.append(('DATE', Date,'Date this file was created.'),end=True)
    
                
    if os.path.exists(Dir+Filename) == False :
                
         hlist.writeto(Dir+Filename)
         
         print(Dir+Filename + " has been created"+"\n")
                
    else :
                
        print(Dir+Filename + " already exists"+"\n")
    
    return Self



#------------------------------------------------------------------------------
# main
#
def main():
    """
    This is the main routine.
    """
    time1 = time.time()
    
    # Parse the command line.
    parser = OptionParser()
    parser.disable_interspersed_args()
    parser.add_option('-q','--quiet',dest='quiet',action='store_true',
    default=False,help='quiet output')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',
    default=False,help='verbose output')
    (options,args) = parser.parse_args()
    nargs = len(args)
    
    if nargs < 1:
        parser.print_help()
    
    else:
    
    
    
        # Image to Analyze - Give the full path
        FileFitsImage  = args[0]
       
        #fileImage = args[1]
        
        
        # Catalog Path
        ctg = '/Users/Elisa/c/Catalogs/'
        
        # KD-3Programs Path
        kd3 = '/Users/Elisa/c/kd-match-0.3.0/'
        
        # Image Peak Path
        pso = '/Users/Elisa/c/Files/'
        
        # File Product Directory
        prod_dir = '/Users/Elisa/c/EAntolini/ProductFiles/'
        
        # Base Directory
        base_dir = '/Users/Elisa/c/EAntolini/'
        
        #FITS Dir
        
        fits_dir = '/Users/Elisa/c/EAntolini/AstrometryFITS/'
        
        # Routines Dir
        RoutDir = '/Users/Elisa/c/EAntolini/Routines/'
        
        # Transformed Catalog Coordinates and Fluxes
        X_Cat      = array('f')
        Y_Cat      = array('f')
        Flux_Cat   = array('f')
        RA_Cat     = array('f')
        DEC_Cat    = array('f')
        
        # Image Coordinates and Fluxes
        X_Imge       = array('f')
        Y_Image      = array('f')
        Flux_Image   = array('f')
        Peak_Image   = array('f')
        Size_Image   = array('f')
        Stat1_Image  = array('f')
        Stat2_Image  = array('f')
        Npixel_Image = array('f')
        
        # Image Transformed Coordinates
        X_Trans_Image     = array('f')
        Y_Trans_Image     = array('f')
        

        # Get the Center of the Image from the Header
        hdulist = pyfits.open(FileFitsImage)
        prihdr = hdulist[0].header
        CenterRA  = prihdr['TCSRA']
        CenterDEC = prihdr['TCSDEC']
        
        #Get the Size Image
        nx = prihdr['NAXIS1']
        ny = prihdr['NAXIS2']
        
        
        # Analyze the Image and create the statistical parameters
        #Create PeakStatObj.txt from Imageproc2f.c with the transformed and sorted full catalog
        
        FilePeakImage = prod_dir+'PeakStatObj.txt'
        
        
        if os.path.exists(FilePeakImage) == False :
        
            cmd = base_dir+"imageproc2f"
        
            execWait = True
            subproc = subprocess.Popen([cmd, FileFitsImage])
            print(FilePeakImage + " Has Been Created"+"\n")
 
        
        
        else :
            
            execWait = False
            print(FilePeakImage + " already exists"+"\n")



         # Generate the Catalog


        Remove = array('b')

        for i in range(0,len(FileFitsImage)):
            if FileFitsImage[i] == '/':
                Remove.append(i)

        stopRemove = Remove[len(Remove)-1]+1

        FileCatalog = FileFitsImage.replace(FileFitsImage[Remove[0]:stopRemove],'')

        FileCatalog = FileCatalog.replace('.fits','.cat')


        
        cmd = ['python3','/ocs/commands/ourstars',str(CenterRA),str(CenterDEC),FileCatalog]
        ExecuteCommand(cfg.catalog_dir+FileCatalog,None,cmd,False)




        # Generate a Copy of the Catalog removing the first row
        
        cmd = ['cp',cfg.catalog_dir+FileCatalog,prod_dir]
        ExecValue = ExecuteCommand(prod_dir+FileCatalog,None,cmd,False)


        # Read nstars from the original Catalog
        
        with open(cfg.catalog_dir+FileCatalog, 'r') as fin:
            data = fin.read().splitlines(True)
            nstars = data[0]
       
        print(" Number of Stars in the Catalog : "+ nstars+"\n")
        
        if ExecValue == True :
            with open(prod_dir+FileCatalog, 'w') as fout:
                fout.writelines(data[1:])




        #  Transform the RA and Dec columns of the catalog from Degrees to arcsec

        FileCatalogTrans = FileCatalog.replace('.cat','_Transformed.cat')


        cmd = [kd3+'transform',prod_dir+FileCatalog,'-x','2','-y','3','-c',str(CenterRA*15),str(CenterDEC)]

        ExecuteCommandFileOut(prod_dir+FileCatalogTrans,cmd)


        # Sort the transformed Ctalog from brighter to fainter stars


        X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat = np.loadtxt(prod_dir+FileCatalogTrans, unpack = True)


        X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat = bubble_sort_cat(X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat)
        
        SortedCat = prod_dir+FileCatalogTrans.replace('.cat','_sorted.cat')

        if os.path.exists(SortedCat) == False :
    
            with open( SortedCat, 'w') as fCat:
                for j in range(len(Flux_Cat)):
                    fCat.write(str(X_Cat[j])+" "+str(Y_Cat[j])+" "+str(Flux_Cat[j])+" "+str(RA_Cat[j])+" "+str(DEC_Cat[j])+"\n")
    
            print(SortedCat + " Has been created"+"\n")
            print('\n')

        else :
            print(SortedCat + " already exists"+"\n")
            print('\n')



        if execWait == True :
    
            subproc.wait()



        # Sort the PeakStatObj.txt from brightest to fainter fluxes

        X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image = np.loadtxt(FilePeakImage, unpack = True)
        bubble_sort_image(X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image)
        SortedImage = FilePeakImage.replace('.txt','Sorted.txt')

        if os.path.exists(SortedImage) == False :
    
            with open(SortedImage, 'w') as fImage:
                for j in range(len(Flux_Image)):
                    if Flux_Image[j] > (Peak_Image[j]*5):
                        fImage.write(str(X_Image[j])+" "+str(Y_Image[j])+" "+str(Flux_Image[j])+" "+str(Peak_Image[j])+" "+str(Size_Image[j])+" "+str(Stat1_Image[j])+" "+str(Stat2_Image[j])+" "+str(Npixel_Image[j])+"\n")


        # Generate Shortest Catalog for Triangle Process (first 200 Brightest Stars)

        SortedCat200 = SortedCat.replace('.cat','_200.cat')

        cmd = ['head','-200',SortedCat]

        ExecuteCommandFileOut(SortedCat200,cmd)



        # Generate shortest PeakStatObj for Triangle Process (first 100 Brightest Stars)

        SortedImage100 = SortedImage.replace('.txt','_100.txt')

        cmd = ['head','-100',SortedImage]

        ExecuteCommandFileOut(SortedImage100,cmd)



        # Generates the Trinagles and write the output to a file

        TriangleOut = 'TriangleOut.txt'

        cmd = [kd3+'triangle_kd',SortedImage100,SortedCat200,'-x1','1','-y1','2','-x2','1','-y2','2']
        ExecuteCommandFileOut(prod_dir+TriangleOut,cmd)


        #Take the transformation Parameters from ICS to CCS (-t parameter)

        with open(prod_dir+TriangleOut, 'r') as fin:
            data = fin.read().splitlines(True)
            Tasterism = data[len(data)-1]


         # Create newfile with the output of triangle_kd

        cmd = [kd3+'transform',SortedImage]+Tasterism.split()


        ExecuteCommandFileOut(prod_dir+'Newfile',cmd)



        #Load stars of the Image  in the same coordinate system of the catalog

        X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image = np.loadtxt(SortedImage, unpack = True)
        X_Trans_Image,Y_Trans_Image= np.loadtxt(prod_dir+'Newfile',usecols=[0,1],dtype=[('f0',float),('f1',float)], unpack = True)

        print len(X_Trans_Image)
        print len(X_Image)

        #Compute the Image coordinate transformation from the Image Coordinate System to the Catalog Coordinate Syatem

        with  open(prod_dir+'Newfile', 'w') as fNewfile:
            for j in range(len(X_Trans_Image)):
                fNewfile.write(str(X_Trans_Image[j])+" "+str(Y_Trans_Image[j])+" "+str(X_Image[j])+" "+str(Y_Image[j])+" "+str(Flux_Image[j])+" "+str(Peak_Image[j])+" "+str(Size_Image[j])+" "+str(Stat1_Image[j])+" "+str(Stat2_Image[j])+" "+str(Npixel_Image[j])+"\n")


        # Compute the matching between the objects in the image and the objects in the whole catalog using the Catalog Coordinate System (CCS)

        cmd = [kd3+'match_kd',prod_dir+'Newfile',SortedCat,'-x2','1','-y2','2']

        ExecuteCommandFileOut(prod_dir+'Match',cmd)


        # Compute Astrometry

        #cmd = ['python',RoutDir+'Astrometry.py',prod_dir,nx,ny,SortedCat,prod_dir+"PeakStatObj.txt"]
        cmd = ['python',RoutDir+'Astrometry.py',prod_dir,str(nx),str(ny),SortedCat,prod_dir+'PeakStatObj.txt',str(X_Cat),str(Y_Cat)]
        call(cmd)


        #Create FITS files with Astrometric Information


        # Create WCS File from scratch
        # Create the standard Header of the fits file



        # Transform Reference Point from ICS to CCS


        CRVAL1 =  (cfg.ccd_field_centre[0]*lastline[1] + cfg.ccd_field_centre[1]*lastline[2] + lastline[3])/3600
        CRVAL2 =  (cfg.ccd_field_centre[0]*lastline[4] + cfg.ccd_field_centre[1]*lastline[5] + lastline[6])/3600

        CentRA = (CenterRA*15) #+ CRVAL1

        CentDEC = CenterDEC #+ CRVAL2

        # Compute the Reference Pixels in arcseconds

        A = lastline[1]
        B = lastline[2]
        C = lastline[3]
        D = lastline[4]
        E = lastline[5]
        F = lastline[6]


        
        RefpixX = (B*F-C*E)/(A*E-B*D)
        RefpixY = (A*F - D*C)/(D*B - A*E)
   




        # Date

        tm = time.gmtime()
        date = str(tm[0])+"-"+str(tm[1])+"-"+str(tm[2])+"T"+str(tm[3])+":"+str(tm[4])+":"+str(tm[5])
        
        File1 = 'wcs.wcs'
        File2 = 'new-image.fits'

        hdu = pyfits.PrimaryHDU()
        hdrWCS = CreateWCS(hdu,prod_dir,File1,CentRA,CentDEC,RefpixX,RefpixY,A/3600,B/3600,C/3600,D/3600,E/3600,F/3600,nx,ny,2,z,date)
        
        hdrNewImage = CreateNewImage(prihdr,hdulist,prod_dir,File2,CentRA,CentDEC,RefpixX,RefpixY,A/3600,B/3600,C/3600,D/3600,E/3600,F/3600,nx,ny,2,z,date)



        time2 = time.time()
        duration = time2-time1
        print("The script runs from start to finish in "+str(duration)+" seconds "+"\n")

        '''
#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
