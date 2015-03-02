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



def point_in_poly(vertices,x,y):
    from matplotlib.path import Path as mpPath
    
    path = mpPath(vertices)
    point = (x,y)
    return path.contains_point(point)



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
     

        ##########   IMAGE AND CATALOG MATCH INFROMATION ######################


        #dist      = array('f')
        distances = array('f')
        x         = array('f')

        # Catalog Transformed Coordinates and Flux

        X_Trans_Cat        = array('f')
        Y_Trans_Cat        = array('f')
        F_Trans_Cat        = array('f')
        X_Image_Obj        = array('f')
        Y_Image_Obj        = array('f')


        # Catalog Stars inside the box Image
        X_All_Cat     = array('f')
        Y_All_Cat     = array('f')
        F_All_Cat     = array('f')
        Ra_All_Cat    = array('f')
        Dec_All_Cat   = array('f')

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
        X_Matched_Image     = array('f')
        Y_Matched_Image     = array('f')
        X_Matched_Image_ICS = array('f')
        Y_Matched_Image_ICS = array('f')
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

        with open (prod_dir+TriangleOut,"r" ) as fileHandle:
            lineList = fileHandle.readlines()
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

         #Load the Transformed and Sorted Catalog informations

        X_All_Cat,Y_All_Cat,F_All_Cat,Ra_All_Cat,Dec_All_Cat  = np.loadtxt(SortedCat,unpack=True)

        #Load matched infromations
        X_Image_Obj,Y_Image_Obj, distances,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat  = np.loadtxt(prod_dir+'Match',usecols=[2,3,10,11,12,13],unpack=True)
    
        # Sort objects from shorter to further distances (11th column of Match)
        x = sorted(distances)


        #Compute Mean and Sigma of the distances

        mean  = np.median(np.array(distances))
        sigma = np.std(np.array(distances))
    
        print("The Mean of the Distances between the stars in the catalog and the object in the image is :"+str(mean)+" +/- "+str(sigma)+"\n")

        FCatStarsBoxImage = open(prod_dir+'FCatStarsBoxImage.txt', 'w')
        FMissingObjects     = open(prod_dir+'FMissingObjects.txt', 'w')
        FMatchedStars     = open(prod_dir+'FMatchedStars.txt', 'w')
        FGoodnessAll      = open(prod_dir+'FGoodnessAll.txt', 'w')
        FGoodnessStars    = open(prod_dir+'FGoodnessStars.txt', 'w')
 
        #Generate the file which will contain all the Stars in Catalog inside the size Image

        for i in range(len(X_All_Cat)) :
            inside = point_in_poly(vertices,X_All_Cat[i],Y_All_Cat[i])
            
            if inside == 1 :
               
                X_Box_Cat.append(X_All_Cat[i])
                Y_Box_Cat.append(Y_All_Cat[i])
                F_Box_Cat.append(F_All_Cat[i])
                Mag_Box_Cat.append(- 20 - ((log10(F_All_Cat[i]))/0.4))
                FCatStarsBoxImage.write(str(X_All_Cat[i])+" "+str(Y_All_Cat[i])+" "+str(F_All_Cat[i])+" "+str(- 20 - ((log10(F_All_Cat[i]))/0.4))+"\n")



        for i in range(len(distances)) :
            inside = point_in_poly(vertices,X_Trans_Cat[i],Y_Trans_Cat[i])
             # If the star Catalog is inside the Box Image
            if inside == 1 :
               
               
                #  if the distance betweend the star in the catalog and the star in the image is  > 5 * mean
                if distances[i] > mean*5 :
                    
    
                    # Find Stars in the Image far away from the star Catalog inside the box Image size
           
                    X_Missing_Image.append(X_Trans_Image[i])
                    Y_Missing_Image.append(Y_Trans_Image[i])
                    F_Missing_Image.append(Flux_Image[i])
                    Mag_Missing_Image.append(- 20 - ((log10(Flux_Image[i]))/0.4))
                    Size_Missing_Image.append(Size_Image[i])
                    Stat1_Missing_Image.append(Stat1_Image[i])
                    Stat2_Missing_Image.append(Stat2_Image[i])
                    
                    X_Missing_Cat.append(X_Trans_Cat[i])
                    Y_Missing_Cat.append(Y_Trans_Cat[i])
                    F_Missing_Cat.append(F_Trans_Cat[i])
                
                    FMissingObjects.write(str(X_Trans_Image[i])+" "+str(Y_Trans_Image[i])+" "+str(Flux_Image[i])+" "+str(- 20 - ((log10(Flux_Image[i]))/0.4))+" "+str(Size_Image[i])+" "+str(Stat1_Image[i])+" "+str(Stat2_Image[i])+" "+str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+" "+str(- 20 - ((log10(F_Trans_Cat[i]))/0.4))+"\n")
                    
                
                if distances[i] <= mean*5 :
                    
                 
                    # Find Stars in the Image inside the box closer to the stars in the Catalog
                    
                    
                    X_Matched_Image.append(X_Trans_Image[i])
                    Y_Matched_Image.append(Y_Trans_Image[i])
                    X_Matched_Image_ICS.append(X_Image_Obj[i])
                    Y_Matched_Image_ICS.append(Y_Image_Obj[i])
                    F_Matched_Image.append(Flux_Image[i])
                    Mag_Matched_Image.append(- 20 - ((log10(Flux_Image[i]))/0.4))
                    Size_Matched_Image.append(Size_Image[i])
                    Stat1_Matched_Image.append(Stat1_Image[i])
                    Stat2_Matched_Image.append(Stat2_Image[i])
                
                    X_Matched_Cat.append(X_Trans_Cat[i])
                    Y_Matched_Cat.append(Y_Trans_Cat[i])
                    F_Matched_Cat.append(F_Trans_Cat[i])
                    Mag_Matched_Cat.append(- 20 - ((log10(F_Trans_Cat[i]))/0.4))
                    
                 
                    FMatchedStars.write(str(X_Image_Obj[i])+" "+str(Y_Image_Obj[i])+" "+str(X_Trans_Image[i])+" "+str(Y_Trans_Image[i])+" "+str(Flux_Image[i])+" "+str(- 20 - ((log10(Flux_Image[i]))/0.4))+" "+str(Size_Image[i])+" "+str(Stat1_Image[i])+" "+str(Stat2_Image[i])+" "+str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+" "+str(- 20 - ((log10(F_Trans_Cat[i]))/0.4))+"\n")




        FMissingObjects.close()
        FMatchedStars.close()
        FCatStarsBoxImage.close()


        lenPoint = 0
        count = 0

          # Make binning with 50 stars for bin
        lenDivide = 50
        count2  = 0


        Mag_Mean_Val = 0
        Std_Mag_Mean_Val = 0
        Size_Mean_Val = 0
        Std_Size_Mean_Val = 0
        Stat1_Mean_Val = 0
        Std_Stat1_Mean_Val = 0
        Stat2_Mean_Val = 0
        Std_Stat2_Mean_Val = 0

        Mag_Mean   = array('f')
        Mag_Dev    = array('f')
        Size_Mean  = array('f')
        Size_Dev   = array('f')
        Stat1_Mean = array('f')
        Stat1_Dev  = array('f')
        Stat2_Mean = array('f')
        Stat2_Dev  = array('f')



        #Computing Mean and Standard Deviation for Magnitude, Size, Stat1 and Stat2 Matched Image stars

        for j in range(0, len(Mag_Matched_Image)):
            
                if(j < (lenPoint+lenDivide)) :
                    
                    count = count +1
                    
                    Mag_Mean_Val   += Mag_Matched_Image[j]
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
                    
                    
                    for k in range(lenPoint, (lenPoint+lenDivide)):
                        Std_Mag_Mean_Val   += math.pow((Mag_Matched_Image[k] -  Mag_Mean[count2]),2)
                        Std_Size_Mean_Val  += math.pow((Size_Matched_Image[k] - Size_Mean[count2]),2)
                        Std_Stat1_Mean_Val += math.pow((Stat1_Matched_Image[k] - Stat1_Mean[count2]),2)
                        Std_Stat2_Mean_Val += math.pow((Stat2_Matched_Image[k] - Stat2_Mean[count2]),2)
                    
                    
                    Mag_Dev.append(math.sqrt(Std_Mag_Mean_Val/count))
                    Size_Dev.append(math.sqrt(Std_Size_Mean_Val/count))
                    Stat1_Dev.append(math.sqrt(Std_Stat1_Mean_Val/count))
                    Stat2_Dev.append(math.sqrt(Std_Stat2_Mean_Val/count))

                    Std_Mag_Mean_Val   = 0
                    Std_Size_Mean_Val  = 0
                    Std_Stat1_Mean_Val = 0
                    Std_Stat2_Mean_Val = 0

                    count2 = count2+1
                                      
                    lenPoint = j
                    count = 0
                        
       

        # Compute the Mean and standard deviation of Size, Stat1, Stat2 for the Missing Objects


        # Magnitude and Size Interpolation
        Mag_Size = np.interp(Mag_Mean, Mag_Mean,Size_Mean)

        # Magnitude and Stat1 Interpolation
        Mag_Stat1 = np.interp(Mag_Mean, Mag_Mean,Stat1_Mean)
     
        # Magnitude and Stat2 Interpolation
        Mag_Stat2 = np.interp(Mag_Mean, Mag_Mean,Stat2_Mean)


        # Compute the Goodness of Size, Stat1, Stat2 for the Missing Objects

        Size_Goodness =  (Size_Missing_Image  - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Size))/np.interp(Mag_Missing_Image, Mag_Mean,Size_Dev)
        Stat1_Goodness = (Stat1_Missing_Image - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Stat1))/np.interp(Mag_Missing_Image, Mag_Mean,Stat1_Dev)
        Stat2_Goodness = (Stat2_Missing_Image - np.interp(Mag_Missing_Image, Mag_Mean,Mag_Stat2))/np.interp(Mag_Missing_Image, Mag_Mean,Stat2_Dev)

        Size_Good_Stars      = array('f')
        Stat1_Good_Stars     = array('f')
        Stat2_Good_Stars     = array('f')
        Mag_Miss_Stars_Image = array('f')




        for i in range(len(Size_Goodness)):
            FGoodnessAll.write(str(X_Missing_Image[i])+" "+str(Y_Missing_Image[i])+" "+str(F_Missing_Image[i])+" "+str(Mag_Missing_Image[i])+" "+str(Size_Goodness[i])+" "+str(Stat1_Goodness[i])+" "+str(Stat2_Goodness[i])+"\n")
            
            #Find potential Stars contained in th Missing Objects
            
            if (Size_Goodness[i] <= 1) and (Size_Goodness[i] >= -1):
                if (Stat1_Goodness[i] <= 1) and (Stat1_Goodness[i] >= -1):
                    if (Stat2_Goodness[i] <= 1) and (Stat2_Goodness[i] >= -1):
                        
                        Size_Good_Stars.append(Size_Goodness[i])
                        Stat1_Good_Stars.append(Stat1_Goodness[i])
                        Stat2_Good_Stars.append(Stat2_Goodness[i])
                        Mag_Miss_Stars_Image.append(Mag_Missing_Image[i])
                        
                        
                        FGoodnessStars.write(str(X_Missing_Image[i])+" "+str(Y_Missing_Image[i])+" "+str(F_Missing_Image[i])+" "+str(Mag_Missing_Image[i])+" "+str(Size_Goodness[i])+" "+str(Stat1_Goodness[i])+" "+str(Stat2_Goodness[i])+"\n")



        FGoodnessStars.close()
        FGoodnessAll.close()



        DeltaX      = array('f')
        DeltaY      = array('f')

        # Find the difference in position between Matched Stars (Catalog-Image) in Catalog Coordinate System

        for j in range(len(X_Matched_Image)):
            DeltaX.append(X_Matched_Image[j] - X_Matched_Cat[j])
            DeltaY.append(Y_Matched_Image[j] - Y_Matched_Cat[j])


        z = np.polyfit(DeltaX,DeltaY,2)
        #w = np.polyfit(Y_Matched_Image,Y_Matched_Cat,2)

        #z = np.polyfit(X_Matched_Image,Y_Matched_Image,2)
        #z = np.polynomial.polynomial.polyval2d(X_Matched_Image,X_Matched_Cat,[1,1])

        print z


        # Plot the cumulative distribution of the distances
        plt.figure(1)
        host = host_subplot(111)
        host.set_xlabel("Distances")
        host.set_ylabel("Number of Times")
        plt.plot(x, np.arange(len(distances)), c='blue')
        savefig(prod_dir+"DistanceDistribution.png")
        #plt.show()


        # Plot the cumulative distribution of the Matched Stars (Image and Catalog) and of the Stars in the Catalog inside the box image

        plt.figure(2)
        host = host_subplot(111)
        host.set_xlabel("Magnitude")
        host.set_ylabel("Number of Times")
        plt.plot(Mag_Box_Cat, np.arange(len(Mag_Box_Cat)), marker = '^',c='blue',markersize = 3, linestyle='None')
        plt.plot(sorted(Mag_Matched_Cat), np.arange(len(Mag_Matched_Cat)),marker = 'o', c='green',markersize = 3,linestyle='None')
        plt.plot(Mag_Matched_Image, np.arange(len(Mag_Matched_Image)),marker = 'o', c='red',markersize = 3,linestyle='None')
        savefig(prod_dir+"MagnitudesDistribution.png")
        #plt.show()



        # plot the Catalog and Image overlap in the same coordinate system
        plt.figure(3)
        #Image
        plt.plot(X_Trans_Image,Y_Trans_Image,linestyle = 'none',marker = 'o',c='red',markersize = 2.2)
        #Catalog
        plt.plot(X_Cat,Y_Cat,linestyle = 'none',marker = '^',c='blue', markersize = 2)
        #Missing Stars from the Catalog
        plt.plot(X_Missing_Image,Y_Missing_Image,linestyle = 'none',marker = 'x',c='green', markersize = 2.4)

        savefig(prod_dir+"CatImageOverlap.pdf")
        #plt.show()


        # Plot Size_Goodness, Stat1_Goodness and Stat2_Goodness for all the Missing objects
        plt.figure(4)
        host = host_subplot(111)
        host.set_xlabel("Magnitude")
        host.set_ylabel("Goodness")
        plt.gca().invert_xaxis()
        plt.scatter(Mag_Missing_Image, Size_Goodness,  marker = 'o', color = 'g')
        plt.scatter(Mag_Missing_Image, Stat1_Goodness, marker = 'o', color = 'r')
        plt.scatter(Mag_Missing_Image, Stat2_Goodness, marker = 'o', color = 'b')
        #plt.errorbar(Star_Size_Mean,Mag_Image_Sorted_Mean,xerr = Star_Size_StDev,yerr = Mag_Image_Sorted_StDev, marker = 'o', color = 'g', linestyle = "None")
        plt.title("Statistical Analysis")
        grid(False)
        savefig(prod_dir+"Goodness.png")
        plt.legend( loc='lower left')
        #plt.show()

        #  Plot Size_Goodness, Stat1_Goodness and Stat2_Goodness for the Star-like objects candidates belonging to the Missing Stars
        plt.figure(5)
        host = host_subplot(111)
        host.set_xlabel("Magnitude")
        host.set_ylabel("Goodness")
        plt.gca().invert_xaxis()
        plt.scatter(Mag_Miss_Stars_Image, Size_Good_Stars,  marker = 'o', color = 'g')
        plt.scatter(Mag_Miss_Stars_Image, Stat1_Good_Stars, marker = 'o', color = 'r')
        plt.scatter(Mag_Miss_Stars_Image, Stat2_Good_Stars, marker = 'o', color = 'b')
        plt.ylim(-2, 2)
        plt.title("Statistical Analysis of Possible Stars")
        grid(False)
        savefig(prod_dir+"PossibleStarsGoodness.png")
        plt.legend( loc='lower left')
        #plt.show()


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

        #RefpixX = (C*E - D*A)/(D*B - A*E) - 1/D
        
        RefpixX = (B*F-C*E)/(A*E-B*D)
        RefpixY = (A*F - D*C)/(D*B - A*E)
        #RefpixY = -((A*F-C*D)/(A*E-B*D))




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


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
