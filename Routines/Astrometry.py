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

    if nargs < 7:
        parser.print_help()
    
    else:
        
    
        
    
        
        # File Product Directory
        prod_dir  = args[0]
        nx        = float(args[1])
        ny        = float(args[2])
        SortedCat = args[3]
        FilePeakImage = args[4]
        X_Cat = args[5]
        Y_Cat = args[6]
        
        #print Y_Cat
  
  
        ##########   IMAGE AND CATALOG MATCH INFROMATION ######################


        distances = array('f')
        x         = array('f')


        # Catalog Stars inside the box Image in the CCS
        X_All_Cat     = array('f')
        Y_All_Cat     = array('f')
        F_All_Cat     = array('f')

        #Load the Transformed and Sorted Catalog informations
        X_All_Cat,Y_All_Cat,F_All_Cat  = np.loadtxt(SortedCat,usecols=[0,1,2],unpack=True)


        # Image Transformed Coordinates CCS
        X_Trans_Image     = array('f')
        Y_Trans_Image     = array('f')

        X_Trans_Image,Y_Trans_Image= np.loadtxt(prod_dir+'Newfile',usecols=[0,1],dtype=[('f0',float),('f1',float)], unpack = True)
        
        
        # Image Objects in the ICS
        X_Imge       = array('f')
        Y_Image      = array('f')
        Flux_Image   = array('f')
        Size_Image   = array('f')
        Stat1_Image  = array('f')
        Stat2_Image  = array('f')

        X_Image,Y_Image,Flux_Image,Size_Image,Stat1_Image,Stat2_Image = np.loadtxt(FilePeakImage,usecols=[0,1,2,4,5,6],unpack = True)


        # Catalog Transformed Coordinates and Flux

        X_Trans_Cat        = array('f')
        Y_Trans_Cat        = array('f')
        F_Trans_Cat        = array('f')
        X_Image_Obj        = array('f')
        Y_Image_Obj        = array('f')

        #Load matched infromations
        X_Image_Obj,Y_Image_Obj, distances,X_Trans_Cat,Y_Trans_Cat,F_Trans_Cat  = np.loadtxt(prod_dir+'Match',usecols=[2,3,10,11,12,13],unpack=True)


        # Catalog Stars inside the box Image
        X_Box_Cat    = array('f')
        Y_Box_Cat    = array('f')
        F_Box_Cat    = array('f')
        Mag_Box_Cat  = array('f')

        box_image_x        = np.empty(4)
        box_image_y        = np.empty(4)


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




        # Create the Box Image Size in the Catalog Coordinate Sytem


        # Obtain transformation parameters

        with open (prod_dir+'TriangleOut.txt',"r" ) as fileHandle:
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




    
        # Sort objects from shorter to further distances (11th column of Match)
        x = sorted(distances)





        #Compute Mean and Sigma of the distances

        mean  = np.median(np.array(distances))
        sigma = np.std(np.array(distances))
    
        print("The Mean of the Distances between the stars in the catalog and the object in the image is :"+str(mean)+" +/- "+str(sigma)+"\n")

 
        #Generate the file which will contain all the Stars in Catalog inside the size Image
        
        Execute, FCatStarsBoxImage = CheckFile(prod_dir+'FCatStarsBoxImage.txt')
        
        if Execute == True :

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


        '''
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
        '''


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

        '''

#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()
