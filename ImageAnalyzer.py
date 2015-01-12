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
import fit
import statistics as stat
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

from optparse import OptionParser

from subprocess import call

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



def ExecuteCommand(Filename,command,out):
    
    if os.path.exists(Filename) == False :
        
        if out == False :
            subproc = call(command)
            Exec = True
            #subproc.wait()
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


def gaus(r,sigma):
    return r*exp(-r**2/(2*sigma**2))


def point_in_poly(vertices,x,y):
    from matplotlib.path import Path as mpPath
    
    path = mpPath(vertices)
    point = (x,y)
    return path.contains_point(point)

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
        
        # File Product Directory
        
        prod_dir = '/Users/Elisa/c/EAntolini/ProductFiles/'
        
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
        X_Trans_Imge      = array('f')
        Y_Trans_Image     = array('f')
        

        # Get the Center of the Image from the Header
        hdulist = pyfits.open(FileFitsImage)
        prihdr = hdulist[0].header
        CenterRA  = prihdr['TCSRA']
        CenterDEC = prihdr['TCSDEC']
        
        #Get the Size Image
        nx = prihdr['NAXIS1']
        ny = prihdr['NAXIS2']

         # Generate the Catalog
        FileCatalog = FileFitsImage.replace('/Users/Elisa/c/RealImage/','')
        FileCatalog = FileCatalog.replace('.fits','.cat')
        
        cmd = ['python3','/ocs/commands/ourstars',str(CenterRA),str(CenterDEC),FileCatalog]
        ExecuteCommand(cfg.catalog_dir+FileCatalog,cmd,False)

        # Generate a Copy of the Catalog removing the first row
        
        cmd = ['cp',cfg.catalog_dir+FileCatalog,prod_dir]
        ExecValue = ExecuteCommand(prod_dir+FileCatalog,cmd,False)
        
        # Read nstars from the original Catalog
        
        with open(cfg.catalog_dir+FileCatalog, 'r') as fin:
            data = fin.read().splitlines(True)
            nstars = data[0]
       
        print(" Number of Stars in the Catalog : "+ nstars)
        
        if ExecValue == True :
            with open(prod_dir+FileCatalog, 'w') as fout:
                fout.writelines(data[1:])


        #  Transform the RA and Dec columns of the catalog from Degrees to arcsec

        FileCatalogTrans = FileCatalog.replace('.cat','_Transformed.cat')


        cmd = [kd3+'transform',prod_dir+FileCatalog,'-x','2','-y','3','-c',str(CenterRA*15),str(CenterDEC)]

        ExecuteCommandFileOut(prod_dir+FileCatalogTrans,cmd)


        # Sort the transformed Ctalog from brighter to fainter stars


        #X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat = np.loadtxt(ctg+FileCatalogTrans,usecols=[0,1,2],dtype=[('f0',float),('f1',float),('f2',float)], unpack = True)

        X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat = np.loadtxt(prod_dir+FileCatalogTrans, unpack = True)


        X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat = bubble_sort_cat(X_Cat,Y_Cat,Flux_Cat,RA_Cat,DEC_Cat)
        
        SortedCat = prod_dir+FileCatalogTrans.replace('.cat','_sorted.cat')

        if os.path.exists(SortedCat) == False :
    
            fCat = open( SortedCat, 'w')
    
            for j in range(len(Flux_Cat)):
                fCat.write(str(X_Cat[j])+" "+str(Y_Cat[j])+" "+str(Flux_Cat[j])+" "+str(RA_Cat[j])+" "+str(DEC_Cat[j])+"\n")
            
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


        cmd = "/Users/Elisa/c/EAntolini/imageproc2f "+FileFitsImage

        ExecuteCommand(prod_dir+FilePeakImage,cmd,True)
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

        X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image = np.loadtxt(prod_dir+FilePeakImage, unpack = True)
        bubble_sort_image(X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image)
        SortedImage = prod_dir+FilePeakImage.replace('.txt','Sorted.txt')

        if os.path.exists(SortedImage) == False :
    
            fImage = open(SortedImage, 'w')
    
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

        cmd = ['/Users/Elisa/c/kd-match-0.3.0/triangle_kd',SortedImage100,SortedCat200,'-x1','1','-y1','2','-x2','1','-y2','2']
        ExecuteCommandFileOut(prod_dir+TriangleOut,cmd)


        #Take the results ftom Triangles (-t parameter)

        with open(kd3+TriangleOut, 'r') as fin:
            data = fin.read().splitlines(True)
            Tasterism = data[len(data)-1]


         # Create newfile with the output of triangle_kd

        cmd = ['/Users/Elisa/c/kd-match-0.3.0/transform',SortedImage]+Tasterism.split()


        ExecuteCommandFileOut(prod_dir+'Newfile',cmd)

        #Load stars of the Image  in the same coordinate system of the catalog

        X_Image,Y_Image,Flux_Image,Peak_Image,Size_Image,Stat1_Image,Stat2_Image,Npixel_Image = np.loadtxt(SortedImage, unpack = True)
        X_Trans_Imge,Y_Trans_Image= np.loadtxt(prod_dir+'Newfile',usecols=[0,1],dtype=[('f0',float),('f1',float)], unpack = True)

        print len(X_Trans_Imge)
        print len(X_Image)

        #Write all Image Information in Newfile

        fNewfile = open(prod_dir+'Newfile', 'w')
        for j in range(len(Flux_Image)):
            fNewfile.write(str(X_Trans_Imge[j])+" "+str(Y_Trans_Image[j])+" "+str(X_Image[j])+" "+str(Y_Image[j])+" "+str(Flux_Image[j])+" "+str(Peak_Image[j])+" "+str(Size_Image[j])+" "+str(Stat1_Image[j])+" "+str(Stat2_Image[j])+" "+str(Npixel_Image[j])+"\n")

        fNewfile.close()

        # Create match file with the whole catalog

        cmd = ['/Users/Elisa/c/kd-match-0.3.0/match_kd',prod_dir+'Newfile',SortedCat,'-x2','1','-y2','2']

        ExecuteCommandFileOut(prod_dir+'Match',cmd)



        ##########   IMAGE AND CATALOG MATCH INFROMATION ######################


        #dist      = array('f')
        distances = array('f')
        x         = array('f')

        # Catalog Transformed Coordinates and Flux

        X_Trans_Cat        = array('f')
        Y_Trans_Cat        = array('f')
        F_Trans_Cat        = array('f')

        # Catalog Stars inside the box Image
        X_Box_Cat    = array('f')
        Y_Box_Cat    = array('f')
        F_Box_Cat    = array('f')
        Mag_Box_Cat  = array('f')


        # Stars in the Image which r > 5*mean
        X_Missing_Image    = array('f')
        Y_Missing_Image    = array('f')
        F_Missing_Image    = array('f')
        Mag_Missing_Image  = array('f')

        Size_Missing_Image     = array('f')
        Stat1_Missing_Image    = array('f')
        Stat2_Missing_Image    = array('f')

        # Stars in the Image which r <= 5*mean
        X_Matched_Image    = array('f')
        Y_Matched_Image    = array('f')
        F_Matched_Image    = array('f')
        Mag_Matched_Image  = array('f')

        Size_Matched_Image     = array('f')
        Stat1_Matched_Image    = array('f')
        Stat2_Matched_Image    = array('f')

        # Stars in the Catalog which r > 5*mean
        X_Missing_Cat    = array('f')
        Y_Missing_Cat    = array('f')
        F_Missing_Cat    = array('f')

        # Stars in the Catalog which r <= 5*mean
        X_Matched_Cat    = array('f')
        Y_Matched_Cat    = array('f')
        F_Matched_Cat    = array('f')
        Mag_Matched_Cat  = array('f')


        box_image_x        = np.empty(4)
        box_image_y        = np.empty(4)



        # Create the Box Image Size in the Catalog Coordinate Sytem


        # Obtain transformation parameters

        fileHandle = open (prod_dir+TriangleOut,"r" )
        lineList = fileHandle.readlines()
        fileHandle.close()
        lastline = lineList[len(lineList)-1]
        lastline = lastline.split(" ")

        for i in range(1,7):
            lastline[i] = float(lastline[i])

        #first Corner (0,0)

        box_image_x[0] = 0*lastline[1] + 0*lastline[2] + lastline[3]
        box_image_y[0] = 0*lastline[4] + 0*lastline[5] + lastline[6]

        #Second Corner (0,ny)
        box_image_x[1] = 0*lastline[1] + ny*lastline[2] + lastline[3]
        box_image_y[1] = 0*lastline[4] + ny*lastline[5] + lastline[6]

        #Third Corner  (nx,ny)
        box_image_x[2] = nx*lastline[1] + ny*lastline[2] + lastline[3]
        box_image_y[2] = nx*lastline[4] + ny*lastline[5] + lastline[6]

        #Fourth Corner (nx,0)
        box_image_x[3] = nx*lastline[1] + 0*lastline[2] + lastline[3]
        box_image_y[3] = nx*lastline[4] + 0*lastline[5] + lastline[6]


        # Fill the vertices of the Box
        vertices = [(box_image_x[0],box_image_y[0]),(box_image_x[1],box_image_y[1]),(box_image_x[2],box_image_y[2]),(box_image_x[3],box_image_y[3])]



        #Load matched infromations
        distances,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat  = np.loadtxt(kd3+'Match',usecols=[10,11,12,13],unpack=True)
    
        # Sort objects from shorter to further distances (11th column of Match)
        x = sorted(distances)
    
        mean  = np.median(np.array(distances))
        sigma = np.std(np.array(distances))
    
        print mean
        print sigma

        FCatStarsBoxImage = open(prod_dir+'FCatStarsBoxImage.txt', 'w')
        FMissingStars     = open(prod_dir+'FMissingStars.txt', 'w')
        FMatchedStars     = open(prod_dir+'FMatchedStars.txt', 'w')


        for i in range(len(distances)) :
            inside = point_in_poly(vertices,X_Trans_Cat[i],Y_Trans_Cat[i])
            
             # If the star Catalog is inside the Box Image
            if inside == 1 :
               
               #All the stars of the Catalog in the Box Image Size
               
                FCatStarsBoxImage.write(str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+"\n")
                X_Box_Cat.append(X_Trans_Cat[i])
                Y_Box_Cat.append(Y_Trans_Cat[i])
                F_Box_Cat.append(F_Trans_Cat[i])
               
                #  if the distance betweend the star in the catalog and the star in the image is  > 5 * mean
                if distances[i] > mean*5 :
                    
    
                    # Find Stars in the Image far away to the star Catalog inside the box Image size
           
                    X_Missing_Image.append(X_Trans_Imge[i])
                    Y_Missing_Image.append(Y_Trans_Image[i])
                    F_Missing_Image.append(Flux_Image[i])
                    Mag_Missing_Image.append(- 20 - ((log10(Flux_Image[i]))/0.4))
                    Size_Missing_Image.append(Size_Image[i])
                    Stat1_Missing_Image.append(Stat1_Image[i])
                    Stat2_Missing_Image.append(Stat2_Image[i])
                    
                    X_Missing_Cat.append(X_Trans_Cat[i])
                    Y_Missing_Cat.append(Y_Trans_Cat[i])
                    F_Missing_Cat.append(F_Trans_Cat[i])
                
                    FMissingStars.write(str(X_Trans_Imge[i])+" "+str(Y_Trans_Image[i])+" "+str(Flux_Image[i])+" "+str(- 20 - ((log10(Flux_Image[i]))/0.4))+" "+str(Size_Image[i])+" "+str(Stat1_Image[i])+" "+str(Stat2_Image[i])+" "+str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+" "+str(- 20 - ((log10(F_Trans_Cat[i]))/0.4))+"\n")
                    
                
                if distances[i] <= mean*5 :
                    
                 
                    # Find Stars in the Image inside the box closing to the star Catalog
                    
                    
                    X_Matched_Image.append(X_Trans_Imge[i])
                    Y_Matched_Image.append(Y_Trans_Image[i])
                    F_Matched_Image.append(Flux_Image[i])
                    Mag_Matched_Image.append(- 20 - ((log10(Flux_Image[i]))/0.4))
                    Size_Matched_Image.append(Size_Image[i])
                    Stat1_Matched_Image.append(Stat1_Image[i])
                    Stat2_Matched_Image.append(Stat2_Image[i])
                
                    X_Matched_Cat.append(X_Trans_Cat[i])
                    Y_Matched_Cat.append(Y_Trans_Cat[i])
                    F_Matched_Cat.append(F_Trans_Cat[i])
                    Mag_Matched_Cat.append(- 20 - ((log10(F_Trans_Cat[i]))/0.4))

                    FMatchedStars.write(str(X_Trans_Imge[i])+" "+str(Y_Trans_Image[i])+" "+str(Flux_Image[i])+" "+str(- 20 - ((log10(Flux_Image[i]))/0.4))+" "+str(Size_Image[i])+" "+str(Stat1_Image[i])+" "+str(Stat2_Image[i])+" "+str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+" "+str(- 20 - ((log10(F_Trans_Cat[i]))/0.4))+"\n")



        FMissingStars.close()
        FMatchedStars.close()
        FCatStarsBoxImage.close()




        if len(Size_Missing_Image) > len(Size_Matched_Image):
            lenDivide = int(len(Size_Missing_Image)/len(Size_Matched_Image))
        else :
            lenDivide = int(len(Size_Matched_Image)/len(Size_Missing_Image))



        lenPoint = 0
        count = 0


        Mag_Mean_Val = 0
        Size_Mean_Val = 0
        Stat1_Mean_Val = 0
        Stat2_Mean_Val = 0

        Mag_Mean   = array('f')
        Mag_Dev    = array('f')
        Size_Mean  = array('f')
        Size_Dev   = array('f')
        Stat1_Mean = array('f')
        Stat1_Dev  = array('f')
        Stat2_Mean = array('f')
        Stat2_Dev  = array('f')





        #while(lenPoint < (len(Mag_Matched_Image_Sorted))) :
        for j in range(0, len(Mag_Matched_Image_Sorted)):
            
                if(j < (lenPoint+lenDivide)) :
                    
                    count = count +1
                    
                    Mag_Mean_Val   += Mag_Matched_Image_Sorted[j]
                    Size_Mean_Val  += Size_Matched_Image[j]
                    Stat1_Mean_Val += Stat1_Matched_Image[j]
                    Stat2_Mean_Val += Stat2_Matched_Image[j]
                
                
                elif(j ==  (lenPoint+lenDivide)):
              
                    Mag_Mean.append(Mag_Mean_Val/count)
                  
                    Mag_Mean_Val = 0
                    Size_Mean.append(Size_Mean_Val/count)
                    Size_Mean_Val = 0
                    Stat1_Mean.append(Stat1_Mean_Val/count)
                    Stat1_Mean_Val = 0
                    Stat2_Mean.append(Stat2_Mean_Val/count)
                    Stat2_Mean_Val = 0
                                      
                    lenPoint = j
                    count = 0

        # Compute the Mean and standard deviation of Size, Stat1, Stat2 for the Missing Stars

        Mag_Dev.append(np.std(Mag_Mean))


        # Magnitude and Size Interpolation
        Mag_Size = np.interp(Mag_Mean, Mag_Mean,Size_Mean)
        Size_Dev.append(np.std(Size_Mean))
        #print Mag_Size
        
        # Magnitude and Stat1 Interpolation
        Mag_Stat1 = np.interp(Mag_Mean, Mag_Mean,Stat1_Mean)
        Stat1_Dev.append(np.std(Stat1_Mean))
        #print Mag_Stat1

        # Magnitude and Stat2 Interpolation
        Mag_Stat2 = np.interp(Mag_Mean, Mag_Mean,Stat1_Mean)
        Stat2_Dev.append(np.std(Stat2_Mean))
        #print Mag_Stat2


        Interp_Par = np.interp(Mag_Missing_Image,Mag_Mean,Mag_Size)

        # See if the No-matched stars are star-like objects


        Size_Goodness =  (Size_Missing_Image  - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Size))/np.interp(Mag_Missing_Image, Mag_Mean,Size_Dev)
        Stat1_Goodness = (Stat1_Missing_Image - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Stat1))/np.interp(Mag_Missing_Image, Mag_Mean,Stat1_Dev)
        Stat2_Goodness = (Stat2_Missing_Image - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Stat2))/np.interp(Mag_Missing_Image, Mag_Mean,Stat2_Dev)



        #Create a file

        # y_Cat_ok = np.arange(len(Cat_ok))

  
        # plot the cumulative function
        plt.figure(1)
        plt.plot(x, np.arange(len(distances)), c='blue')
        savefig(prod_dir+"DistanceDistribution.png")
        plt.show()



        # plot the Catalog and Image overlap in the same coordinate system
        plt.figure(2)
        #Image
        plt.plot(X_Trans_Imge,Y_Trans_Image,linestyle = 'none',marker = 'o',c='red',markersize = 2.2)
        #Catalog
        plt.plot(X_Cat,Y_Cat,linestyle = 'none',marker = '^',c='blue', markersize = 2)
        #Missing Stars from the Catalog
        plt.plot(X_Missing_Image,Y_Missing_Image,linestyle = 'none',marker = 'x',c='green', markersize = 2.4)

        savefig(prod_dir+"CatImageOverlap.pdf")
        plt.show()


        # Plot Size_Goodness, Stat1_Goodness and Stat2_Goodness
        plt.figure(9)
        host = host_subplot(111)
        host.set_xlabel("Magnitude")
        host.set_ylabel("Goodness")
        plt.gca().invert_xaxis()
        plt.scatter(Miss_Mag, Size_Goodness, marker = 'o', color = 'g')
        plt.scatter(Miss_Mag, Stat1_Goodness, marker = 'o', color = 'r')
        plt.scatter(Miss_Mag, Stat2_Goodness, marker = 'o', color = 'b')
        #plt.errorbar(Star_Size_Mean,Mag_Image_Sorted_Mean,xerr = Star_Size_StDev,yerr = Mag_Image_Sorted_StDev, marker = 'o', color = 'g', linestyle = "None")
        plt.title("Statistical Analysis")
        grid(False)
        savefig(prod_dir+"Goodness.png")
        plt.legend( loc='lower left')
        plt.show()

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
