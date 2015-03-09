#!/opt/local/bin/python2.7

import os
import sys
from sys import argv
import pyfits
import time
from optparse import OptionParser
from numpy.linalg import inv

from scipy import optimize
import algopy
from array import *
import numpy as np

sys.path.append('/ocs/modules')
import cfg


#------------------------------------------------------------------------------
# Implementation of Functions
#


def residuals(a,b,c,y,u,v):

    err = y - (a*u*v + b*v*v + c*u*u)

    return err

def peval(u,v,p):
    
    return p[0]*u*v+p[1]*v*v+p[2]*u*u



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
    parser.add_option('-q','--quiet',dest='quiet',action='store_true',default=False,help='quiet output')
    parser.add_option('-v','--verbose',dest='verbose',action='store_true',default=False,help='verbose output')
    (options,args) = parser.parse_args()
    nargs = len(args)
    
    if nargs < 1:
        parser.print_help()
    
    else:

        # File Product Directory
        prod_dir = '/Users/Elisa/c/EAntolini/ProductFiles/'
        
        
        # Stars in the Catalog which r <= 5*mean
        u_Image    = array('f') #Distance North from the center of  Image object in the ICS
        v_Image    = array('f') #Distance East  from the center of  Image object in the ICS
        x_Cat      = array('f') #Distance North from the center of  Catalog object in the CCS
        y_Cat      = array('f') #Distance East  from the center of  Catalog object in the CCS
        U_Cat    = array('f') #Distance North from the center of  Catalog object in the ICS
        V_Cat    = array('f') #Distance East  from the center of  Catalog object in the ICS
        DeltaU    = array('f') #Distance North from the center of  Catalog object in the ICS
        DeltaV    = array('f') #Distance East  from the center of  Catalog object in the ICS
        
        #Load matched infromations
        u_Image,v_Image, x_Cat,y_Cat = np.loadtxt(prod_dir+'FMatchedStars.txt',usecols=[0,1,2,3],unpack=True)



        # Obtain transformation parameters

        with open (prod_dir+'TriangleOut.txt',"r" ) as fileHandle:
            lineList = fileHandle.readlines()
        lastline = lineList[len(lineList)-1]
        lastline = lastline.split(" ")

        for i in range(1,7):
            lastline[i] = float(lastline[i])


        A = lastline[1]
        B = lastline[2]
        C = lastline[3]
        D = lastline[4]
        E = lastline[5]
        F = lastline[6]

        #Transform x,y (Catalog Objects in CCS) into U,V (Catalog Objects in ICS)
        for j in range(len(x_Cat)):
            U_Cat.append(x_Cat[j]*(1/A+(B*D)/(A*E - B*D))-(B/(A*E - B*D))*((y_Cat[j]-F)+D*C) -C/A)
            DeltaU.append(U_Cat[j]-u_Image[j])
            V_Cat.append((1/(A*E-D*B))*(A*(y_Cat[j] - F)-D*(x_Cat[j]-C)))
            DeltaV.append(V_Cat[j]-v_Image[j])

        #Compute inverse Matrix

        for i in range(len(u_Image)):
            fMatrix = np.array([[u_Image[i]*v_Image[i],u_Image[i]*u_Image[i],v_Image[i]*v_Image[i]],[u_Image[i]*v_Image[i],u_Image[i]*u_Image[i],v_Image[i]*v_Image[i]],[u_Image[i]*v_Image[i],u_Image[i]*u_Image[i],v_Image[i]*v_Image[i]]])

            fMatrixInv = inv(fMatrix)
            print fMatrixInv


#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

