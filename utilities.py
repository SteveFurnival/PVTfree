
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Utilities, i.e. commonly used functions
#
#=======================================================================

#!/usr/bin/python3

import numpy as NP

#========================================================================
#  Normalise a 1D-array; just return array
#========================================================================

def Norm(x) :

    xSum = NP.sum(x)
    xNor = NP.divide(x,xSum)

#== return array ======================================================

    return xNor

#========================================================================
#  Normalise a 1D-array; return array and sum
#========================================================================

def NormSum(x) :

    xSum = NP.sum(x)
    xNor = NP.divide(x,xSum)

#== return array ======================================================

    return xNor,xSum

#========================================================================
#  Linear Fit of Y versus X
#========================================================================

def linearFit(X,Y) :

    Ndeg  = 1
    coefs = NP.polyfit(X,Y,Ndeg)

    slope = coefs[0] ; inter = coefs[1]

    #print("linearFit: slope,inter {:10.3e} {:10.3e}".format(slope,inter))

    return slope,inter

#========================================================================
#  Write a Vector of type vTyp (Int or Real) to Named File
#========================================================================

def writeVector(fOut,vNam,vArr,vTyp) :

    if   vTyp.upper() == "I" : nFor = "{:3d} "
    else                     : nFor = "{:10.3e} "

    sOut = vNam + "\n"
    fOut.write(sOut)

    sOut = " "
    for sElm in vArr :
        sOut = sOut + nFor.format(sElm)
    sOut = sOut + "\n"
    fOut.write(sOut)

#== Return ============================================================        
        
    return

#========================================================================
#  Write a Matrix
#========================================================================

def writeMatrix(fOut,mNam,mMat,mTyp) :

    if   mTyp == "I" : nFor = "{:3d} "
    else             : nFor = "{:10.3e} "

    sOut = mNam + "\n"
    fOut.write(sOut)

    for vArr in mMat :
        sOut = " "
        for sElm in vArr :
            sOut = sOut + nFor.format(sElm)
        sOut = sOut + "\n"
        fOut.write(sOut)

#== Return ============================================================        
        
    return

#========================================================================
#  Single Eigenvalue GDEM (for Psat & 1-Sided Stab-Check)
#========================================================================

def GDEM1(res0,res1,clsIO) :

#--------------------------------------------------------------------
#  Zick recommended max eigV = 1.0E+06; Soreide page 56
#  SRF has set      max eigV = 1.0E+03; Can get crashes with large eigV
#--------------------------------------------------------------------

#   eMax = 1.0E+06
    eMax = 1.0E+03

    b01 = NP.dot(res0,res1)
    b11 = NP.dot(res1,res1)

    dNom = b11 - b01

    if dNom != 0.0 :
        eigV = min(abs(b11/dNom),eMax)
    else :
        eigV = 1.0

    #if qDeb :
    #    sOut = "GDEM1: b01,b11,dNom,eigV {:10.3e} {:10.3e} {:10.3e} {:10.3e}\n".format(b01,b11,dNom,eigV)
    #    fDeb.write(sOut)
    
    return eigV

#========================================================================
#  End of Module
#========================================================================


