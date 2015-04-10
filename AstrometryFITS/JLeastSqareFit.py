#!/opt/local/bin/python2.7

import os
import sys
from sys import argv
import pyfits
import time
from optparse import OptionParser
from numpy.linalg import inv

from scipy import optimize
# import algopy
from array import *
import numpy as np
import math
from scipy.optimize import leastsq

# sys.path.append('/ocs/modules')
# import cfg


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

        time1 = time.time()

        for tt in range(0,100):
            u2_Image=u_Image*u_Image
            v2_Image=v_Image*v_Image
            uv_Image=u_Image*v_Image
            
            x11 = np.sum(u2_Image*u2_Image)
            x12 = np.sum(u2_Image*uv_Image)
            x13 = np.sum(v2_Image*u2_Image)

            x21 = x12
            x22 = x13
            x23 = np.sum(v2_Image*uv_Image)
            
            x31 = x13
            x32 = x23
            x33 = np.sum(v2_Image*v2_Image)

            q1_U = np.sum(u2_Image*DeltaU)
            q2_U = np.sum(uv_Image*DeltaU)
            q3_U = np.sum(v2_Image*DeltaU)
            
            q1_V = np.sum(u2_Image*DeltaV)
            q2_V = np.sum(uv_Image*DeltaV)
            q3_V = np.sum(v2_Image*DeltaV)

            fMatrix = np.array([[x11,x12,x13],[x21,x22,x23],[x31,x32,x33]])

            #Compute inverse Matrix

            fMatrixInv = inv(fMatrix)

            a_U = fMatrixInv[0][0]*q1_U+fMatrixInv[0][1]*q2_U+fMatrixInv[0][2]*q3_U
            b_U = fMatrixInv[1][0]*q1_U+fMatrixInv[1][1]*q2_U+fMatrixInv[1][2]*q3_U
            c_U = fMatrixInv[2][0]*q1_U+fMatrixInv[2][1]*q2_U+fMatrixInv[2][2]*q3_U

            a_V = fMatrixInv[0][0]*q1_V+fMatrixInv[0][1]*q2_V+fMatrixInv[0][2]*q3_V
            b_V = fMatrixInv[1][0]*q1_V+fMatrixInv[1][1]*q2_V+fMatrixInv[1][2]*q3_V
            c_V = fMatrixInv[2][0]*q1_V+fMatrixInv[2][1]*q2_V+fMatrixInv[2][2]*q3_V

        print ("The first part runs in %g\n" % (time.time()-time1))

        print a_U,b_U,c_U
        print a_V,b_V,c_V

        time1 = time.time()

        for tt in range(0,100):
            tplInitial1 = (1e-07 ,-1e-9, 1e-8)
            tplInitial2 = (-1e-13 ,1e-13, -1e-13)

            funcQuad=lambda tpl,u_Image,v_Image : tpl[0]*u_Image**2+tpl[1]*u_Image*v_Image+tpl[2]*v_Image**2
            func=funcQuad
            ErrorFunc1=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaU
            ErrorFunc2=lambda tpl,u_Image,v_Image,DeltaU: func(tpl,u_Image,v_Image)-DeltaV

            tplFinal1,success=leastsq(ErrorFunc1,tplInitial1[:],args=(u_Image,v_Image,DeltaU))
            tplFinal2,success=leastsq(ErrorFunc2,tplInitial2[:],args=(u_Image,v_Image,DeltaV))

        print ("The second part runs in %g\n" % (time.time()-time1))

        print "quadratic fit of U and u" ,tplFinal1
        print "quadratic fit of V and v" ,tplFinal2



#------------------------------------------------------------------------------
# Start program execution.
#
if __name__ == '__main__':
    main()

