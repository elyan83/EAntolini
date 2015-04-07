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
from scipy.optimize import leastsq

sys.path.append('/ocs/modules')
import astro
import cfg
import image
import util
import pyfits


#------------------------------------------------------------------------------
# Implementation of Functions
#



#Create Catalog

def ourstars(ra,dec,filename,filedir):
    
    
    if os.path.exists(filedir+filename) == False :


        #nx = cfg.ccd_nx         # number of columns in simulated image
        #ny = cfg.ccd_ny         # number of rows in simulated image
        binf = 1                # binning factor
        fwhm = 1
        resolution = 5          # Number of images per FWHM for star trails
        band = 'V'              # Default band
        ra_rate = 0             # Default is untrailed image
        dec_rate = 0
        track_length = 0
        dx = 0
        dy = 0
        #ra = args[0]
        #dec = args[1]
        #filename = args[2]
        #filedir  = args[3]

        nadd = int(resolution*track_length/(binf*fwhm))+1      # no of images in track
        nadd2 = nadd/2
    
        nx = int(cfg.ccd_nx/binf)
        ny = int(cfg.ccd_ny/binf)
 
    
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
    
        print("Catalog Size arcsec = ra : " +str((ra_max-ra_min)*3600)+" Dec : "+str((dec_max-dec_min)*3600)+"\n")
    
        # Determine x and y increments for each image
        dx /= nadd
        dy /= nadd
    
        # Add stars. The origin is in the SE corner of the image. Dec is
        # proportional to y and ra is proportional to -x.
        # First get a list of stars from the USNO catalog.
        stars = astro.usno(cfg.catalog_dir+'usno185r',ra_min,ra_max,dec_min,dec_max)

    
        # Print out the Catalog
        fileout= open(filedir+filename, 'a+')
    
        if not stars:
            print('no stars in range',ra_min,ra_max,dec_min,dec_max)
            exit(0)
        else:
            print(len(stars),'stars')
            #fileout.write(str(len(stars))+"\n")
    
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
        
            fileout.write('{:4.1f} {:s} {:6.8f} {:s} {:6.8f} {:s}'.format(flux," ",(star.ra)*15," ",star.dec,"\n"))

        else :
            
            print(filedir+filename + " already exists"+"\n")


# Create doAstrometry Function and call it in ImageAnalyzer.py

def CheckFile(Filename):
    if os.path.exists(Filename) == False :
        
            Exec = True
            FileName = open(Filename, 'w')
            
            print(Filename+" Has been created "+"\n")
                
    else:
                         
            Exec = False
            FileName = None
            print(Filename+" already exists"+"\n")



    return Exec, FileName




def point_in_poly(vertices,x,y):
    from matplotlib.path import Path as mpPath
    
    path = mpPath(vertices)
    point = (x,y)
    return path.contains_point(point)


def GenBoxImge(FileTriangle,FileCatalog,nx,ny,prod_dir):

    box_image_x        = np.empty(4)
    box_image_y        = np.empty(4)

    # Catalog Stars inside the box Image
    X_Box_Cat    = array('f')
    Y_Box_Cat    = array('f')
    F_Box_Cat    = array('f')
    Mag_Box_Cat  = array('f')

    #Create the Box Image Size in the Catalog Coordinate Sytem
    
    
    # Obtain transformation parameters
    
    with open (FileTriangle,"r" ) as fileHandle:
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


    #Generate the file which will contain all the Stars in Catalog inside the size Image

    Execute, FCatStarsBoxImage = CheckFile(prod_dir+'FCatStarsBoxImage.txt')

    if Execute == True :
        
        #Load the Transformed and Sorted Catalog informations
        X_All_Cat,Y_All_Cat,F_All_Cat  = np.loadtxt(FileCatalog,usecols=[0,1,2],unpack=True)

        for i in range(len(X_All_Cat)) :
            inside = point_in_poly(vertices,X_All_Cat[i],Y_All_Cat[i])
        
            if inside == 1 :
            
                X_Box_Cat.append(X_All_Cat[i])
                Y_Box_Cat.append(Y_All_Cat[i])
                F_Box_Cat.append(F_All_Cat[i])
                Mag_Box_Cat.append(- 20 - ((log10(F_All_Cat[i]))/0.4))
                FCatStarsBoxImage.write(str(X_All_Cat[i])+" "+str(Y_All_Cat[i])+" "+str(F_All_Cat[i])+" "+str(- 20 - ((log10(F_All_Cat[i]))/0.4))+"\n")

    if FCatStarsBoxImage != None :
        FCatStarsBoxImage.close()

    return lastline,vertices

'''
def MatchMissObj(distances,prod_dir,vertices,X_Image_Obj,Y_Image_Obj,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat,X_Trans_Image,Y_Trans_Image,Flux_Image,Size_Image,Stat1_Image,Stat2_Image):


    # Stars in the Image which r > 5*mean
    X_Missing_Image    = array('f')
    Y_Missing_Image    = array('f')
    F_Missing_Image    = array('f')
    Mag_Missing_Image  = array('f')
    Size_Missing_Image     = array('f')
    Stat1_Missing_Image    = array('f')
    Stat2_Missing_Image    = array('f')

    # Stars in the Catalog which r > 5*mean
    X_Missing_Cat    = array('f')
    Y_Missing_Cat    = array('f')
    F_Missing_Cat    = array('f')


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


    # Stars in the Catalog which r <= 5*mean
    X_Matched_Cat    = array('f')
    Y_Matched_Cat    = array('f')
    F_Matched_Cat    = array('f')
    Mag_Matched_Cat  = array('f')



    #Compute Mean and Sigma of the distances

    mean  = np.median(np.array(distances))
    sigma = np.std(np.array(distances))
    
    print("The Mean of the Distances between the stars in the catalog and the object in the image is :"+str(mean)+" +/- "+str(sigma)+"\n")
    
    
    
    
    Execute1,FMissingObjects  = CheckFile(prod_dir+'FMissingObjects.txt')
    Execute2,FMatchedStars    = CheckFile(prod_dir+'FMatchedStars.txt')
    
    
    
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
                
                
                if Execute1 == True :
                    
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
                
                if Execute2 == True :
                    
                    FMatchedStars.write(str(X_Image_Obj[i])+" "+str(Y_Image_Obj[i])+" "+str(X_Trans_Image[i])+" "+str(Y_Trans_Image[i])+" "+str(Flux_Image[i])+" "+str(- 20 - ((log10(Flux_Image[i]))/0.4))+" "+str(Size_Image[i])+" "+str(Stat1_Image[i])+" "+str(Stat2_Image[i])+" "+str(X_Trans_Cat[i])+" "+str(Y_Trans_Cat[i])+" "+str(F_Trans_Cat[i])+" "+str(- 20 - ((log10(F_Trans_Cat[i]))/0.4))+"\n")
    
    
    if FMissingObjects != None :
        FMissingObjects.close()
    
    if FMatchedStars != None :
        FMatchedStars.close()

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


    for lenpoint in range(0,len(Mag_Matched_Image),lendivide):
  
  
            Mag_Mean_Val =np.mean(Mag_Matched_Image[lenPoint:lenPoint+lenDivide])
            Size_Mean_Val = np.mean(Size_Matched_Image[lenPoint:lenPoint+lenDivide])
            Stat1_Mean_Val = np.mean(Stat1_Matched_Image[lenPoint:lenPoint+lenDivide])
            Stat2_Mean_Val = np.mean(Stat2_Matched_Image[lenPoint:lenPoint+lenDivide])
    
    
    
        
            Mag_Dev.append(math.sqrt(Std_Mag_Mean_Val/lendivide))
            Mag_Mean_Val = 0
            Size_Mean.append(math.sqrt(Size_Mean_Val/lendivide))
            Size_Mean_Val = 0
            Stat1_Mean.append(math.sqrt(Stat1_Mean_Val/lendivide))
            Stat1_Mean_Val = 0
            Stat2_Mean.append(math.sqrt(Stat2_Mean_Val/lendivide))
            Stat2_Mean_Val = 0
        
            Std_Mag_Mean_Val  = np.sum(math.pow(Mag_Matched_Image[lenPoint:lenPoint+lenDivide]-Mag_Mean[Count2],2))
            Std_Size_Mean_Val = np.sum(math.pow(Size_Matched_Image[lenPoint:lenPoint+lenDivide]-Size_Mean[Count2],2))
            Std_Stat1_Mean_Val =np.sum(math.pow(Stat1_Matched_Image[lenPoint,lenPoint+lenDivide]-Stat1_Mean[Count2],2))
            Std_Stat2_Mean_Val =np.sum(math.pow(Stat2_Matched_Image[lenPoint,lenPoint+lenDivide]-Stat2_Mean[Count2],2))
            
        
    
                                                            
            Mag_Dev.append(math.sqrt(np.mean(math.pow(Mag_Matched_Image[lenPoint:lenPoint+lenDivide]-Mag_Mean[Count2],2))
                                                                    
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


    Execute1, FGoodnessAll = CheckFile(prod_dir+'FGoodnessAll.txt')
    Execute2, FGoodnessStars = CheckFile(prod_dir+'FGoodnessStars.txt')

    for i in range(len(Size_Goodness)):
        if Execute1 == True :
            FGoodnessAll.write(str(X_Missing_Image[i])+" "+str(Y_Missing_Image[i])+" "+str(F_Missing_Image[i])+" "+str(Mag_Missing_Image[i])+" "+str(Size_Goodness[i])+" "+str(Stat1_Goodness[i])+" "+str(Stat2_Goodness[i])+"\n")
    
    #Find potential Stars contained in th Missing Objects
    
        if (Size_Goodness[i] <= 1) and (Size_Goodness[i] >= -1):
            if (Stat1_Goodness[i] <= 1) and (Stat1_Goodness[i] >= -1):
                if (Stat2_Goodness[i] <= 1) and (Stat2_Goodness[i] >= -1):
                
                    Size_Good_Stars.append(Size_Goodness[i])
                    Stat1_Good_Stars.append(Stat1_Goodness[i])
                    Stat2_Good_Stars.append(Stat2_Goodness[i])
                    Mag_Miss_Stars_Image.append(Mag_Missing_Image[i])
                
                    if Execute1 == True :
                        FGoodnessStars.write(str(X_Missing_Image[i])+" "+str(Y_Missing_Image[i])+" "+str(F_Missing_Image[i])+" "+str(Mag_Missing_Image[i])+" "+str(Size_Goodness[i])+" "+str(Stat1_Goodness[i])+" "+str(Stat2_Goodness[i])+"\n")


    if FGoodnessStars != None :
        FGoodnessStars.close()

    if FGoodnessAll != None :
        FGoodnessAll.close()

    return Size_Goodness,Stat1_Goodness,Stat2_Goodness,Size_Good_Stars,Stat1_Good_Stars,Stat2_Good_Stars,Mag_Miss_Stars_Image,Mag_Missing_Image,X_Missing_Image,Y_Missing_Image

'''

def MakePlot(ncanvas,ndata,title,FigName,xlab,ylab,invert,x1,y1,x2,y2,x3,y3,show):
    
     if os.path.exists(FigName) == False :
    
        plt.figure(ncanvas)
        host = host_subplot(111)
        host.set_xlabel(xlab)
        host.set_ylabel(ylab)

        if invert == True :
            plt.gca().invert_xaxis()

        if ndata == 1 or ndata == 2 or ndata == 3 :
            plt.scatter(x1,y1,marker = 'o', color = 'g')

        if ndata == 2 or ndata == 3 :
            plt.scatter(x2,y2,marker = 'o', color = 'r')

        if ndata == 3 :
            plt.scatter(x3,y3,marker = 'o', color = 'b')

        plt.title(title)
        grid(False)
        savefig(FigName)
        plt.legend( loc='lower left')
        print(FigName+" has been created "+"\n")

        if show :
            plt.show()
                
                
     else :

        print(FigName + " already exists"+"\n")


def GetRefPoint(lastline,ccd_field_centre,CentRA,CentDEC) :
    
    # Transform Reference Point from ICS to CCS

    CRVAL1 =  (ccd_field_centre[0]*lastline[1] + ccd_field_centre[1]*lastline[2] + lastline[3])/3600
    CRVAL2 =  (ccd_field_centre[0]*lastline[4] + ccd_field_centre[1]*lastline[5] + lastline[6])/3600

    #CentRA = (CenterRA*15)

    # Compute the Reference Pixels in arcseconds

    A = lastline[1]
    B = lastline[2]
    C = lastline[3]
    D = lastline[4]
    E = lastline[5]
    F = lastline[6]



    RefpixX = (B*F-C*E)/(A*E-B*D)
    RefpixY = (A*F - D*C)/(D*B - A*E)

    return RefpixX,RefpixY,A,B,C,D,E,F

def GetDate():
    tm = time.gmtime()
    date = str(tm[0])+"-"+str(tm[1])+"-"+str(tm[2])+"T"+str(tm[3])+":"+str(tm[4])+":"+str(tm[5])

    return date


def LeastSquareFit(u_Image,v_Image,x_Cat,y_Cat,A,B,C,D,E,F):
    
    
    U_Cat    = array('f') #Distance North from the center of  Catalog object in the ICS
    V_Cat    = array('f') #Distance East  from the center of  Catalog object in the ICS
    DeltaU   = array('f') #Distance North from the center of  Catalog object in the ICS
    DeltaV   = array('f') #Distance East  from the center of  Catalog object in the ICS

    #Transform x,y (Catalog Objects in CCS) into U,V (Catalog Objects in ICS)
    
    for j in range(len(x_Cat)):
        U_Cat.append(x_Cat[j]*(1/A+(B*D)/(A*E - B*D))-(B/(A*E - B*D))*((y_Cat[j]-F)+D*C) -C/A)
        DeltaU.append(U_Cat[j]-u_Image[j])
        V_Cat.append((1/(A*E-D*B))*(A*(y_Cat[j] - F)-D*(x_Cat[j]-C)))
        DeltaV.append(V_Cat[j]-v_Image[j])

    #Compute inverse Matrix
    

    uv_Image = v_Image*u_Image
    u2_Image = u_Image*u_Image
    u3_Image = u_Image*u_Image*u_Image
    u4_Image = u2_Image * u2_Image
    v2_Image = v_Image*v_Image
    v3_Image = v_Image*v_Image*v_Image
    v4_Image = v2_Image*v2_Image

    x11 = np.sum(u4_Image)
    x12 = x21 = np.sum(u3_Image*v_Image)
    x13 = x22 = x31 = np.sum(u2_Image*v2_Image)
    x23 = x32 = np.sum(v3_Image *u_Image)
    x33 = np.sum(v4_Image)

    fMatrix = np.array([[x11,x12,x13],[x21,x22,x23],[x31,x32,x33]])

    fMatrixInv = inv(fMatrix)

    q1_U = np.sum(u2_Image*DeltaU)
    q2_U = np.sum(uv_Image*DeltaU)
    q3_U = np.sum(v2_Image*DeltaU)
    
    q1_V = np.sum(u2_Image*DeltaV)
    q2_V = np.sum(uv_Image*DeltaV)
    q3_V = np.sum(v2_Image*DeltaV)

    a_U = fMatrixInv[0][0]*q1_U+fMatrixInv[0][1]*q2_U+fMatrixInv[0][2]*q3_U
    b_U = fMatrixInv[1][0]*q1_U+fMatrixInv[1][1]*q2_U+fMatrixInv[1][2]*q3_U
    c_U = fMatrixInv[2][0]*q1_U+fMatrixInv[2][1]*q2_U+fMatrixInv[2][2]*q3_U
    
    a_V = fMatrixInv[0][0]*q1_V+fMatrixInv[0][1]*q2_V+fMatrixInv[0][2]*q3_V
    b_V = fMatrixInv[1][0]*q1_V+fMatrixInv[1][1]*q2_V+fMatrixInv[1][2]*q3_V
    c_V = fMatrixInv[2][0]*q1_V+fMatrixInv[2][1]*q2_V+fMatrixInv[2][2]*q3_V

    print("Fit Paremeters Found"+"\n")

    print("A_0_2 : "+ str(a_U)+"\n"+"A_1_1 : "+str(b_U)+"\n"+"A_2_0 : "+str(c_U)+"\n"+"B_0_2 : "+str(a_V)+"\n"+"B_1_1 : "+str(b_V)+"\n"+"B_2_0 : "+str(c_V)+"\n")
    print ("\n")

    return a_U,b_U,c_U,a_V,b_V,c_V,DeltaU,DeltaV


def InitPar(par):
    
    Remove = array('b')
    
    
    for i in range(0,len(par)):
        if par[i] == 'e':
            Remove.append(i)
    
    stopRemove = Remove[len(Remove)-1]+1

    if par[0] != '-':
        
        par = par.replace(par[0:stopRemove-1],'1')
    else :
        par = par.replace(par[0:stopRemove-1],'-1')

    return float(par)
    



def LeastSquareFitComparison(u_Image,v_Image,DeltaU,DeltaV,p):


    tplInitial1 = (p[0] ,p[1],p[2])
    tplInitial2 = (p[3] ,p[4],p[5])


    funcQuad=lambda tpl,u_Image,v_Image : tpl[0]*u_Image**2+tpl[1]*u_Image*v_Image+tpl[2]*v_Image**2
    func=funcQuad
    ErrorFunc1=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaU
    ErrorFunc2=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaV
    
    tplFinal1,success=leastsq(ErrorFunc1,tplInitial1[:],args=(u_Image,v_Image,DeltaU))
    tplFinal2,success=leastsq(ErrorFunc2,tplInitial2[:],args=(u_Image,v_Image,DeltaV))
    
    print("Fit Parameters Comparison Done :"+"\n")
    print "quadratic fit of U and u" ,tplFinal1
    print "quadratic fit of V and v" ,tplFinal2



def CreateWCS(Self,Dir,Filename,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Nx,Ny,Aorder,Border,Z,Date):
    
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
    
    Self.header['A_ORDER'] = Aorder,'Polynomial order, axis 1'
    Self.header['A_0_2']     = Z[0],' '
    Self.header['A_1_1']     = Z[1],' '
    Self.header['A_2_0']     = Z[2],' '
    
    Self.header['B_ORDER'] = Border,'Polynomial order, axis 2'
    Self.header['B_0_2']     = Z[3],' '
    Self.header['B_1_1']     = Z[4],' '
    Self.header['B_2_0']     = Z[5],' '
    
    # Date
    Self.header['DATE'] = Date,'Date this file was created.'
    
    if os.path.exists(Dir+Filename) == False :
        
        Self.writeto(Dir+Filename)
        
        print(Dir+Filename + " has been created"+"\n")
    
    else :
        
        print(Dir+Filename + " already exists"+"\n")
    
    return Self

def CreateNewImage(Self,hlist,Dir,Filename,cenVAL1,cenVAL2,cenPIX1,cenPIX2,a,b,c,d,e,f,Nx,Ny,Aorder,Border,Z,Date):
    
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
    
    Self.append(('A_ORDER',Aorder,'Polynomial order, axis 2'),end=True)
    Self.append(('A_0_2',Z[0],' '),end=True)
    Self.append(('A_1_1',Z[1],' '),end=True)
    Self.append(('A_2_0',Z[2],' '),end=True)
    
    Self.append(('B_ORDER',Border,'Polynomial order, axis 2'),end=True)
    Self.append(('B_0_2',Z[3],' '),end=True)
    Self.append(('B_1_1',Z[4],' '),end=True)
    Self.append(('B_2_0',Z[5],' '),end=True)
    
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



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
