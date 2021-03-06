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
import math
from scipy.optimize import leastsq

sys.path.append('/ocs/modules')
import cfg


#------------------------------------------------------------------------------
# Implementation of Functions
#





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
        
        ytrue = array('f')
        
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

        x11 = 0.0
        x12 = 0.0
        x13 = 0.0
        x21 = 0.0
        x22 = 0.0
        x23 = 0.0
        x31 = 0.0
        x32 = 0.0
        x33 = 0.0

        q1_U = 0.0
        q2_U = 0.0
        q3_U = 0.0

        q1_V = 0.0
        q2_V = 0.0
        q3_V = 0.0

        p = [1e-07 ,-1e-9, 1e-8]


        for i in range(len(u_Image)):
            
         x11 += pow(u_Image[i],4)
         x12 += pow(u_Image[i],3)*v_Image[i]
         x13 += pow(v_Image[i],2)*pow(u_Image[i],2)

         x21 += pow(u_Image[i],3)*v_Image[i]
         x22 += pow(v_Image[i],2)*pow(u_Image[i],2)
         x23 += pow(v_Image[i],3)*u_Image[i]

         x31 += pow(v_Image[i],2)*pow(u_Image[i],2)
         x32 += pow(v_Image[i],3)*u_Image[i]
         x33 += pow(v_Image[i],4)

         q1_U += pow(u_Image[i],2)*DeltaU[i]
         q2_U += v_Image[i]*u_Image[i]*DeltaU[i]
         q3_U += pow(v_Image[i],2)*DeltaU[i]

         q1_V += pow(u_Image[i],2)*DeltaV[i]
         q2_V += v_Image[i]*u_Image[i]*DeltaV[i]
         q3_V += pow(v_Image[i],2)*DeltaV[i]
        
        


        fMatrix = np.array([[x11,x12,x13],[x21,x22,x23],[x31,x32,x33]])

        fMatrixInv = inv(fMatrix)

        a_U = fMatrixInv[0][0]*q1_U+fMatrixInv[0][1]*q2_U+fMatrixInv[0][2]*q3_U
        b_U = fMatrixInv[1][0]*q1_U+fMatrixInv[1][1]*q2_U+fMatrixInv[1][2]*q3_U
        c_U = fMatrixInv[2][0]*q1_U+fMatrixInv[2][1]*q2_U+fMatrixInv[2][2]*q3_U

        a_V = fMatrixInv[0][0]*q1_V+fMatrixInv[0][1]*q2_V+fMatrixInv[0][2]*q3_V
        b_V = fMatrixInv[1][0]*q1_V+fMatrixInv[1][1]*q2_V+fMatrixInv[1][2]*q3_V
        c_V = fMatrixInv[2][0]*q1_V+fMatrixInv[2][1]*q2_V+fMatrixInv[2][2]*q3_V

        print a_U,b_U,c_U
        print a_V,b_V,c_V


        #Comparison with Leastsquare Function

        tplInitial1 = (1e-07 ,-1e-9, 1e-8)
        tplInitial2 = (-1e-13 ,1e-13, -1e-13)

        funcQuad=lambda tpl,u_Image,v_Image : tpl[0]*u_Image**2+tpl[1]*u_Image*v_Image+tpl[2]*v_Image**2
        func=funcQuad
        ErrorFunc1=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaU
        ErrorFunc2=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaV

        tplFinal1,success=leastsq(ErrorFunc1,tplInitial1[:],args=(u_Image,v_Image,DeltaU))
        tplFinal2,success=leastsq(ErrorFunc2,tplInitial2[:],args=(u_Image,v_Image,DeltaV))
        print "quadratic fit of U and u" ,tplFinal1
        print "quadratic fit of V and v" ,tplFinal2



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

