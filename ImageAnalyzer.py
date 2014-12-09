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
import statistics as stat

from optparse import OptionParser

from subprocess import call

import cfg


#------------------------------------------------------------------------------
# Implementation of Functions
#
'''
def bubble_sort_cat(DistRa,DistDec,Flux):
    
    for i in range(len(Flux)):
        for j in range(len(Flux)-1-i):
            if Flux[j] < Flux[j+1]:
                Flux[j], Flux[j+1] = Flux[j+1], Flux[j]
                DistRa[j], DistRa[j+1] = DistRa[j+1], DistRa[j]
                DistDec[j], DistDec[j+1] = DistDec[j+1], DistDec[j]
'''

def bubble_sort_cat(DistRa,DistDec,Flux):
    DistRa, DistDec, Flux = zip(*sorted(zip(DistRa, DistDec, Flux), reverse=True, key=lambda x: x[2]))
    
    return DistRa, DistDec, Flux

def bubble_sort_image(DistRa,DistDec,Allpix,PeakVal):
    
    for i in range(len(Allpix)):
        for j in range(len(Allpix)-1-i):
            if Allpix[j] < Allpix[j+1]:
                Allpix[j], Allpix[j+1]   = Allpix[j+1], Allpix[j]
                PeakVal[j],PeakVal[j+1]  = PeakVal[j+1],PeakVal[j]
                DistRa[j], DistRa[j+1]   = DistRa[j+1], DistRa[j]
                DistDec[j], DistDec[j+1] = DistDec[j+1],DistDec[j]

def ExecuteCommand(Filename,command,out):
    
    if os.path.exists(Filename) == False :
        
        if out == False :
            call(command)
            Exec = True
            print(Filename + " has been created")
        
        else :
            print(os.popen(command).read())
            Exec = True
            print(Filename + " has been created")
    else:
        
        Exec = False
        print(Filename + " already exists")

    return Exec

def ExecuteCommandFileOut(Filename,command):
    
    
    if os.path.exists(Filename) == False :
        
        myoutfile = open(Filename, 'w')
        
        call(command, stdout = myoutfile)
        Exec = True
        print(Filename+" has been created ")
    
    else:
        
        Exec = False
        print(Filename + " already exists")

    return Exec




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
        
        # Transformed Catalog Coordinates and Fluxes
        X_Cat     = array('f')
        Y_Cat     = array('f')
        Flux_Cat  = array('f')
        
        # Image Coordinates and Fluxes
        X_Imge      = array('f')
        Y_Image     = array('f')
        Flux_Image  = array('f')
        Peak_Image  = array('f')
        
        # Image Transformed Coordinates
        X_Trans_Imge      = array('f')
        Y_Trans_Image     = array('f')
        
        

        # Get the Center of the Image from the Header
        hdulist = pyfits.open(FileFitsImage)
        prihdr = hdulist[0].header
        CenterRA  = prihdr['TCSRA']
        CenterDEC = prihdr['TCSDEC']

         # Generate the Catalog
        FileCatalog = FileFitsImage.replace('/Users/Elisa/c/RealImage/','')
        FileCatalog = FileCatalog.replace('.fits','.cat')
        
        cmd = ['python3','/ocs/commands/ourstars',str(CenterRA),str(CenterDEC),FileCatalog]
        ExecuteCommand(cfg.catalog_dir+FileCatalog,cmd,False)

        # Generate a Copy of the Catalog removing the first row
        
        cmd = ['cp',cfg.catalog_dir+FileCatalog,ctg]
        ExecValue = ExecuteCommand(ctg+FileCatalog,cmd,False)
        
        # Read nstars from the original Catalog
        
        with open(cfg.catalog_dir+FileCatalog, 'r') as fin:
            data = fin.read().splitlines(True)
            nstars = data[0]
       
        print(" Number of Stars in the Catalog : "+ nstars)
        
        if ExecValue == True :
            with open(ctg+FileCatalog, 'w') as fout:
                fout.writelines(data[1:])


        #  Transform the RA and Dec columns of the catalog from Degrees to arcsec

        FileCatalogTrans = FileCatalog.replace('.cat','_Transformed.cat')


        cmd = [kd3+'transform',ctg+FileCatalog,'-x','2','-y','3','-c',str(CenterRA*15),str(CenterDEC)]

        ExecuteCommandFileOut(ctg+FileCatalogTrans,cmd)


        # Sort the transformed Ctalog from brighter to fainter stars


        X_Cat,Y_Cat,Flux_Cat = np.loadtxt(ctg+FileCatalogTrans,usecols=[0,1,2],dtype=[('f0',float),('f1',float),('f2',float)], unpack = True)

        X_Cat,Y_Cat,Flux_Cat = bubble_sort_cat(X_Cat,Y_Cat,Flux_Cat)
        
        SortedCat = ctg+FileCatalogTrans.replace('.cat','_sorted.cat')

        if os.path.exists(SortedCat) == False :
    
            fCat = open( SortedCat, 'w')
    
            for j in range(len(Flux_Cat)):
                fCat.write(str(X_Cat[j])+" "+str(Y_Cat[j])+" "+str(Flux_Cat[j])+"\n")
            
            print(SortedCat + " Has been created")

        else :
            print(SortedCat + " already exists")

         #  Add the number of stars in the transformed and sorted catalog (first row)


        f = open(SortedCat)
        text = f.read().splitlines(True)
        first_line = text[0]
        f.close()

        if first_line != nstars :
            #print("First Line Added")
            f = open(SortedCat, 'w')
            f.write(nstars)
            f.writelines(text[1:])
            f.close()
         
        #Create PeakStatObj.txt from Imageproc2f.c with the transformed and sorted full catalog

        FilePeakImage = 'PeakStatObj.txt'


        cmd = ["/Users/Elisa/c/imageproc2f "+FileFitsImage+" < "+SortedCat]

        ExecuteCommand(pso+FilePeakImage,cmd,True)
        print('\n')

        #  Remove the number of stars in the transformed and sorted catalog (first row)

        with open(SortedCat, 'r') as fin:
            data = fin.read().splitlines(True)
        first_line = data[0]

        if first_line == nstars :
            #print("First Line Removed")
            f = open(SortedCat, 'w')
            f.writelines(text[1:])
            f.close()


        # Sort the PeakStatObj.txt from brightest to fainter fluxes

        X_Image,Y_Image,Flux_Image,Peak_Image = np.loadtxt(pso+FilePeakImage,usecols=[0,1,2,3],dtype=[('f0',float),('f1',float),('f2',float),('f3',float)], unpack = True)
        bubble_sort_image(X_Image,Y_Image,Flux_Image,Peak_Image)
        SortedImage = pso+FilePeakImage.replace('.txt','Sorted.txt')

        if os.path.exists(SortedImage) == False :
    
            fImage = open(SortedImage, 'w')
    
            for j in range(len(Flux_Image)):
                if Flux_Image[j] > (Peak_Image[j]*5):
                    fImage.write(str(X_Image[j])+" "+str(Y_Image[j])+" "+str(Flux_Image[j])+"\n")

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

        cmd = ['/Users/Elisa/c/kd-match-0.3.0/triangle_kd',SortedImage100,SortedCat200,'-x1','1','-y1','2','-x2','1','-y2','2']
        ExecuteCommandFileOut(kd3+TriangleOut,cmd)


        #Take the results ftom Triangles (-t parameter)

        with open(kd3+TriangleOut, 'r') as fin:
            data = fin.read().splitlines(True)
            Tasterism = data[len(data)-1]


         # Create newfile with the output of triangle_kd

        #cmd = ['/Users/Elisa/c/kd-match-0.3.0/transform',SortedImage,TriangleOut]
        cmd = ['/Users/Elisa/c/kd-match-0.3.0/transform',SortedImage]+Tasterism.split()


        ExecuteCommandFileOut(kd3+'Newfile',cmd)

        #Load stars of the Image  in the same coordinate system od the catalog)

        X_Trans_Imge,Y_Trans_Image= np.loadtxt(kd3+'Newfile',usecols=[0,1],dtype=[('f0',float),('f1',float)], unpack = True)

        # Create match file with the whole catalog

        cmd = ['/Users/Elisa/c/kd-match-0.3.0/match_kd',kd3+'Newfile',SortedCat,'-x2','1','-y2','2']

        ExecuteCommandFileOut(kd3+'Match',cmd)


        # Sort objects from shorter to further distances (6th column of Match)

        dist=array('f')
        distances=array('f')

        with open(kd3+'Match') as fd:
            fd.readline()
            for line in fd:
                dist.append(float(line.split()[5]))

        col=[i for i in zip(dist)]
        distances = sorted(col)

        #for i in range(len(col)):
            #print(distances[i])

        # Generate the cumulative distribution of the distances


        #evaluate the cumulative
        cumulative = np.cumsum(distances)

        # plot the cumulative function
        plt.figure(1)
        #plt.plot(distances, cumulative, c='blue')
        plt.plot(distances, np.arange(len(distances)), c='blue')
        savefig("/Users/Elisa/c/Files/DistanceDistribution.png")
        plt.show()

        # plot the Catalog and Image overlap in the same coordinate system
        plt.figure(2)
        #Image
        plt.plot(X_Trans_Imge,Y_Trans_Image,linestyle = 'none',marker = 'o',c='red',markersize = 2.2)
        #Catalog
        plt.plot(X_Cat,Y_Cat,linestyle = 'none',marker = '^',c='blue', markersize = 2)
        savefig("/Users/Elisa/c/Files/CatImageOverlap.pdf")
        plt.show()


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
