#!/usr/bin/python3

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

sys.path.append('/Users/Elisa/c/EAntolini/Routines')
import Astrometry


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


def create_dir(f):
    
    d = os.path.dirname(f)

    if  os.path.exists(d) == False:
    
        os.mkdir(d)
        print("The directory "+f+" has been created"+"\n")

    else:

        print("The directory "+f+" already exists"+"\n")


def RemovePath(Filename):
    
    Remove = array('b')
    
    for i in range(0,len(Filename)):
        if Filename[i] == '/':
            Remove.append(i)
    
    stopRemove = Remove[len(Remove)-1]+1

    return Filename.replace(Filename[Remove[0]:stopRemove],'')
    


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
    
    
        #### Input Parameters #####
        
        # Image to Analyze - Give the full path
        FileFitsImage  = args[0]
        
        #Fit Comparison
        CompareFit = (args[1] == "True")
        
        # Print Out Plots
        Show = (args[2] == "True")

        ###########################

        # Catalog Path
        ctg = '/Users/Elisa/c/Catalogs/'
        
        # KD-3Programs Path
        kd3 = '/Users/Elisa/c/kd-match-0.3.0/'
        
        # Image Peak Path
        pso = '/Users/Elisa/c/Files/'
        
        # Base Directory
        base_dir = '/Users/Elisa/c/EAntolini/'
        
        #FITS Dir
        
        fits_dir = '/Users/Elisa/c/EAntolini/AstrometryFITS/'
        
        # Routines Dir
        RoutDir = '/Users/Elisa/c/EAntolini/Routines/'
        
        #print FileFitsImage
        
        prod_dir = RemovePath(FileFitsImage)
        
        prod_dir = prod_dir.replace('.fits','/')
        
        # File Product Directory
        prod_dir = '/Users/Elisa/c/EAntolini/'+prod_dir
        
        create_dir(prod_dir)
        
        
        

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
            
            print(cmd, FileFitsImage,prod_dir+"PeakStatObj.txt")
        
            execWait = True
            subproc = subprocess.Popen([cmd, FileFitsImage,prod_dir+"PeakStatObj.txt"])
            print(FilePeakImage + " Has Been Created"+"\n")
 
        
        
        else :
            
            execWait = False
            print(FilePeakImage + " already exists"+"\n")



         # Generate the Catalog


        
        FileCatalog = RemovePath(FileFitsImage)
        
        FileCatalog = FileCatalog.replace('.fits','.cat')

        Astrometry.ourstars(CenterRA,CenterDEC,FileCatalog,prod_dir)

        #cmd = ['python3','/ocs/commands/ourstars',str(CenterRA),str(CenterDEC),FileCatalog,prod_dir]
        #ExecuteCommand(cfg.catalog_dir+FileCatalog,None,cmd,False)


        '''

        FileCatalog2 = FileCatalog.replace('_original.cat','.cat')


        # Generate a Copy of the Catalog removing the first row


        cmd = ['cp',prod_dir+FileCatalog,prod_dir]
        ExecValue = ExecuteCommand(prod_dir+FileCatalog2,None,cmd,False)


        # Read nstars from the original Catalog
        
        with open(prod_dir+FileCatalog, 'r') as fin:
            data = fin.read().splitlines(True)
            nstars = data[0]
       
        print(" Number of Stars in the Catalog : "+ nstars+"\n")



        if ExecValue == True :
            with open(prod_dir+FileCatalog2, 'w') as fout:
                fout.writelines(data[1:])


        '''

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



        # Sort the PeakStatObj.txt from brightest to fainter fluxes and write all information in Newfile

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



        #Compute the Image coordinate transformation from the Image Coordinate System to the Catalog Coordinate Syatem

        with  open(prod_dir+'Newfile', 'w') as fNewfile:
            for j in range(len(X_Trans_Image)):
                fNewfile.write(str(X_Trans_Image[j])+" "+str(Y_Trans_Image[j])+" "+str(X_Image[j])+" "+str(Y_Image[j])+" "+str(Flux_Image[j])+" "+str(Peak_Image[j])+" "+str(Size_Image[j])+" "+str(Stat1_Image[j])+" "+str(Stat2_Image[j])+" "+str(Npixel_Image[j])+"\n")



        # Compute the matching between the objects in the image and the objects in the whole catalog using the Catalog Coordinate System (CCS)

        cmd = [kd3+'match_kd',prod_dir+'Newfile',SortedCat,'-x2','1','-y2','2']

        ExecuteCommandFileOut(prod_dir+'Match',cmd)


        #Load matched infromations
        X_Image_Obj,Y_Image_Obj, distances,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat  = np.loadtxt(prod_dir+'Match',usecols=[2,3,10,11,12,13],unpack=True)

        ####################### Compute Astrometry #####################

         #Create the Box Image Size in the Catalog Coordinate Sytem
        lastline,vertices = Astrometry.GenBoxImge(prod_dir+'TriangleOut.txt',SortedCat,nx,ny,prod_dir)


         #Compute Mean and Sigma of the distances and Generate the Missing and Matching files and the statistical parameters

        Size_Goodness,Stat1_Goodness,Stat2_Goodness,Size_Good_Stars,Stat1_Good_Stars,Stat2_Good_Stars,Mag_Miss_Stars_Image,Mag_Missing_Image,X_Missing_Image,Y_Missing_Image = Astrometry.MatchMissObj(distances,prod_dir,vertices,X_Image_Obj,Y_Image_Obj,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat,X_Trans_Image,Y_Trans_Image,Flux_Image,Size_Image,Stat1_Image,Stat2_Image)



        # Plot the cumulative distribution of the distances
        Astrometry.MakePlot(1,1,"Cumulative distribution of the distances",prod_dir+"CatImageOverlap.pdf","Distances","Number of Times",False,sorted(distances),np.arange(len(distances)),0.0,0.0,0.0,0.0,Show)

        # Plot the cumulative distribution of the Matched Stars (Image and Catalog) and of the Stars in the Catalog inside the box image
        #Astrometry.MakePlot(2,3,"Image and Catalog Overlap",prod_dir+"MagnitudesDistribution.png","Magnitude","Number of Times",False,Mag_Box_Cat,np.arange(len(Mag_Box_Cat)),sorted(Mag_Matched_Cat), np.arange(len(Mag_Matched_Cat)),Mag_Matched_Image, np.arange(len(Mag_Matched_Image)),Show)


        # plot the Catalog and Image overlap in the same coordinate system
        Astrometry.MakePlot(3,3,"Image and Catalog Overlap",prod_dir+"CatImageOverlap.pdf","X coordinate","Y coordinate",False,X_Trans_Image,Y_Trans_Image,X_Cat,Y_Cat,X_Missing_Image,Y_Missing_Image,Show)

        # Plot Size_Goodness, Stat1_Goodness and Stat2_Goodness for all the Missing objects
        Astrometry.MakePlot(4,3,"Statistical Analysis of Possible Stars",prod_dir+"Goodness.png","Magnitude","Goodness",True,Mag_Missing_Image, Size_Goodness,Mag_Missing_Image, Stat1_Goodness,Mag_Missing_Image, Stat2_Goodness,Show)


        #  Plot Size_Goodness, Stat1_Goodness and Stat2_Goodness for the Star-like objects candidates belonging to the Missing Stars
        Astrometry.MakePlot(5,3,"Statistical Analysis of Possible Stars",prod_dir+"PossibleStarsGoodness.png","Magnitude","Goodness",True,Mag_Miss_Stars_Image,Size_Good_Stars,Mag_Miss_Stars_Image, Stat1_Good_Stars,Mag_Miss_Stars_Image, Stat2_Good_Stars,Show)




        #Get the Reference Point of the maps and the Transformation parameters
        RefpixX,RefpixY,A,B,C,D,E,F = Astrometry.GetRefPoint(lastline,cfg.ccd_field_centre,CenterRA*15,CenterDEC)

        date = Astrometry.GetDate()


         # Fit with least square methid  the difference in position between Matched Stars-objects (Catalog-Image) in Catalog Coordinate System

        A_11,A_12,A_13,B_11,B_12,B_13,DeltaU,DeltaV = Astrometry.LeastSquareFit(X_Image_Obj,Y_Image_Obj,X_Trans_Image,Y_Trans_Image,A,B,C,D,E,F)



        #Compare the previous Fit
        if CompareFit :

            p = [Astrometry.InitPar(str(A_11)),Astrometry.InitPar(str(A_12)),Astrometry.InitPar(str(A_13)),Astrometry.InitPar(str(B_11)),Astrometry.InitPar(str(B_12)),Astrometry.InitPar(str(B_13))]

            Astrometry.LeastSquareFitComparison(X_Image_Obj,Y_Image_Obj,DeltaU,DeltaV,p)

        Z = [A_11,A_12,A_13,B_11,B_12,B_13]

        #Create FITS files with Astrometric Information


        # Create WCS File from scratch
        # Create the standard Header of the fits file



        File1 = 'wcs.wcs'
        File2 = 'new-image.fits'

        hdu = pyfits.PrimaryHDU()
        hdrWCS = Astrometry.CreateWCS(hdu,prod_dir,File1,CenterRA*15,CenterDEC,RefpixX,RefpixY,A/3600,B/3600,C/3600,D/3600,E/3600,F/3600,nx,ny,2,2,Z,date)
        
        #hdrNewImage = Astrometry.CreateNewImage(prihdr,hdulist,prod_dir,File2,CenterRA*15,CenterDEC,RefpixX,RefpixY,A/3600,B/3600,C/3600,D/3600,E/3600,F/3600,nx,ny,2,2,Z,date)



        time2 = time.time()
        duration = time2-time1
        print("The script runs from start to finish in "+str(duration)+" seconds "+"\n")
        

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
