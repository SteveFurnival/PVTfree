
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Utilities, i.e. commonly used functions and constants
#
#=======================================================================

#!/usr/bin/python3

import os
import sys

import numpy as NP

from math import exp,sqrt

import allData as AD

#========================================================================
#  Constants
#========================================================================

macEPS = 2.22E-16   #-- 'Double Precision' Minimum Number
mGDEM1 = 4          #-- SS Iteration/GDEM Step for Psat & Stab
mGDEM2 = 5          #-- SS Iteration/GDEM Step for Flash

Third  =  1.0/3.0   #-- Used in cubicSolver
Sqr3   = sqrt(3.0)

btuCon =   1.98588  #-- Gas Constant/Energy   [Btu/(lbmol.degR)]
denAir =   0.0765   #-- Density of Air at STC [lb/ft3]
denWat =  62.428    #-- Density of Pure Water [lb/ft3]
dF2dR  = 459.67     #-- 0 degF                [degR]
gasCon =  10.7315   #-- Gas Constant/Volume   [psia.ft3/(lbmol.degR)]
molAir =  28.97     #-- Mole Weight of Air    [lb/lbmol]
pStand =  14.6959   #-- Standard Pressure     [psia]
tStand = 519.67     #-- Std Temp (60 degF)    [degR]
volMol = 379.483    #-- Vol 1 lbmol Ideal Gas [scf]

lbcMwX = -0.5       #-- LBC Viscosity Mw-Exp
lbcPcX = -2.0/3.0   #-- LBC Viscosity Pc-Exp
lbcTcX =  1.0/6.0   #-- LBC Viscosity Tc-Exp

#========================================================================
#  Create docstring to use in Help-About
#========================================================================

def aboutPVTfree() :

    '''
--
--  PVTfree Equation of State (EOS) PVT Program
--    
--  Copyright (C) 2020
--    
--  PVTfree is free software: you can redistribute it and/or modify
--  it under the terms of the GNU General Public License as published by
--  the Free Software Foundation, either version 3 of the License, or
--  (at your option) any later version.
--    
--  PVTfree is distributed in the hope that it will be useful,
--  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  GNU General Public License for more details.
--    
--  You should have received a copy of the GNU General Public License
--  along with PVTfree.  If not, see <https://www.gnu.org/licenses/>.
--    
--  Steve Furnival, HoBoil Limited, Banchory, UK, December 2020
--    
--  email: steve.furnival@btinternet.com
--
    '''

    return

#========================================================================
#  Normalise a 1D-array; just return array
#========================================================================

def Norm(x) :

    xSum = NP.sum(x)
    xNor = NP.divide(x,xSum)

#== return Nprmalised Array only ======================================

    return xNor

#========================================================================
#  Normalise a 1D-array; return array and sum
#========================================================================

def NormSum(x) :

    xSum = NP.sum(x)
    xNor = NP.divide(x,xSum)

#== return Normalised Array and Un-Normalised Sum =====================

    return xNor,xSum

#========================================================================
#  Linear Fit of Y versus X
#========================================================================

def linearFit(X,Y) :

    nDeg  = 1
    coefs = NP.polyfit(X,Y,nDeg)

    slope = coefs[0]
    inter = coefs[1]

    #print("linearFit: slope,inter {:10.3e} {:10.3e}".format(slope,inter))

#== Return slope (m) and intercept (C) ================================    

    return slope,inter

#========================================================================
#  Write a Vector of type vTyp (Int or Real) to Named File
#========================================================================

def writeVector(fOut,vNam,vArr,vTyp) :

    if   vTyp.upper() == "I" : nFor = "{:3d} "
    else                     : nFor = "{:10.3e} "

    #sOut = vNam + "\n"
    #fOut.write(sOut)

    #sOut = " "

    sOut = vNam + " "
    for sElm in vArr :
        sOut = sOut + nFor.format(sElm)
    sOut = sOut + "\n"
    fOut.write(sOut)

#== No Return =========================================================
        
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

#== No Return =========================================================
        
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

#== Return Eigenvalue =================================================    
    
    return eigV

#========================================================================
#  Wilson K-Values
#========================================================================

def wilsonK(pRes,tRes,clsEOS) :

    nCom = clsEOS.nComp
    wilK = NP.zeros(nCom)

    for iC in range(nCom) :

        rPred = clsEOS.gPP("PC",iC)/pRes  #-- Reciprocal Reduced Pressure
        rTred = clsEOS.gPP("TC",iC)/tRes  #-- Reciprocal Reduced Temperature

        wilK[iC] = rPred*exp(5.3727*(1.0 + clsEOS.gPP("AF",iC))*(1.0 - rTred))

#== Return array ======================================================        

    return wilK

#========================================================================
#  Mole Fraction of C7+
#========================================================================

def moleFracC7P(Z,clsEOS) :

    molC7P = 0.0

    for iC in range(clsEOS.nComp) :
        molW = clsEOS.gPP("MW",iC)
        if molW > 90.0 : molC7P = molC7P + Z[iC]

#== Return Mole Weight ================================================        

    return molC7P

#========================================================================
#  Can a string be converted to a number?
#========================================================================

def isNum(strnG) :

    try :
        float(strnG)
        return True
    except ValueError :
        return False

#========================================================================
#  Convert a String to a Float, Regardless ...
#========================================================================

def S2F(strnG) :

    try :
        numB = float(strnG)
        return numB
    except ValueError :
        return 0.0

#========================================================================
#  To avoid plotting zero values, replace with 'nan' (Not A Number)
#  matplotlib filters out nan's
#========================================================================

def zero_to_nan(values) :

    return [float('nan') if x == 0.0 else x for x in values]

#========================================================================
#  Machine Epsilon
#========================================================================

def calcMacEPS() :

    EPS = 1.0

#-- Keep halving EPS until 1+EPS is indistinguishable from 1 ----------    

    while EPS + 1.0 > 1.0 : EPS = 0.5*EPS

    return 2.0*EPS

#======================================================================
#  Return the sign of a number (or zero if number = 0)
#======================================================================

def signNumber(dNum) :

    if   dNum > 0.0 : dSgn =  1.0
    elif dNum < 0.0 : dSgn = -1.0
    else            : dSgn =  0.0

#== Return dSgn =======================================================

    return dSgn

#========================================================================
#  Process Input File Name/Directory
#========================================================================

def procFileName(filN) :

    sPlt = filN.split("\\")     #-- Split File/Dir-Name by '\'
    nSpl = len(sPlt)

    path = ""
    for iSpl in range(nSpl-1) :
        path = path + sPlt[iSpl] + "\\"     #-- Drive & Sub-Directories

    rooX = sPlt[nSpl-1]         #-- Rootname and Extension
    sPlt = rooX.split(".")
    root = sPlt[0]

#== Return path and root-name =========================================    

    return path,root

#========================================================================
#  Create I/O Files (Input; Output, Save & Log)
#========================================================================

def createIOFiles(path,root,qTem) :

    iERR = 0

#-- Concatenate Path and RootName -----------------------------------    

    pathRoot = path + '\\' + root

#-- Input File Required? --------------------------------------------

    if qTem : pInp = None
    else    : pInp = pathRoot + ".dat"

#-- Output, Save & Log Files ----------------------------------------    

    pOut = pathRoot + ".out"
    pSav = pathRoot + ".sav"
    pLog = pathRoot + ".log"

#== If Input Created, does it exist AND can it be opened? =============

    if not qTem :
        if os.path.isfile(pInp) and os.access(pInp,os.R_OK) :
            print("createIOFiles: Opened Input File ",pInp," and is readable")
        else :
            print("createIOFiles: Input File ",pInp," is missing or is not readable")
            clsIO = None
            iERR  = -1
            return iERR,clsIO

#-- Create the Output Files (and Input, if Appropriate) -------------        

    if pInp == None : fInp = None
    else            : fInp = open(pInp,'r')    #-- Open as Read -Only
    
    fOut = open(pOut,'w')     #-- Open as Write-Only
    fSav = open(pSav,'w+')    #-- Open as Read/Write
    fLog = open(pLog,'w')     #--         Write-Only

#== Create the I/O Class ==============================================

    clsIO = AD.classInpOut(path)

#-- Store file objects ----------------------------------------------

    clsIO.setFiles(root,fInp,fOut,fSav,fLog)

#-- Other files not set ---------------------------------------------    

    clsIO.setQDEB(False)
    clsIO.setQREG(False)
    clsIO.setQ100(False)
    clsIO.setQIMX(False)
    clsIO.setQMOR(False)
    clsIO.setQVIP(False)
    clsIO.setQ300(False)
    clsIO.setQGEM(False)

#== Return Values: Error Flag and I/O Class ===========================

    return iERR,clsIO

#========================================================================
#  Update I/O Files
#========================================================================

def updateIOFiles(path,root,clsIO) :

#-- First, create new input file ------------------------------------

    pathRoot = path + '\\' + root

    clsIO.patR = pathRoot

    pInp = pathRoot + '.dat'
    fInp = open(pInp,'r')

    clsIO.fInp = fInp

#-- Existing Out, Sav & Log File Names ------------------------------

    pOut = clsIO.fOut.name
    pSav = clsIO.fSav.name
    pLog = clsIO.fLog.name

#-- New Out, Sav & Log File Names -----------------------------------

    nOut = pathRoot + '.out'
    nSav = pathRoot + '.sav'
    nLog = pathRoot + '.log'

#== Remove Files if Null and Replace, Else Copy and Replace ===========

    iOut = os.path.getsize(pOut)
    clsIO.fOut = processOldFileNewFile(iOut,clsIO.fOut,nOut)
            
    iSav = os.path.getsize(pSav)
    clsIO.fSav = processOldFileNewFile(iSav,clsIO.fSav,nSav)
        
    iLog = os.path.getsize(pLog)
    clsIO.fLog = processOldFileNewFile(iLog,clsIO.fLog,nLog)

#== Return updated clsIO object =======================================

    return clsIO

#========================================================================
#  Close a File
#========================================================================

def closeFile(fStr) :

    fStr.close()            #-- Close the file
    
    fStrNam = fStr.name     #-- Get its name

#-- If the file is zero-length, delete it ---------------------------    
    
    if os.path.getsize(fStrNam) == 0 : os.remove(fStrNam)

#-- No return value -------------------------------------------------

    return False

#========================================================================
#  Process Old File and New File
#========================================================================

def processOldFileNewFile(iSiz,fOld,pNew) :

    fNew = open(pNew,'w')           #-- Open new file

    if iSiz > 0 :                   #-- If old file non-null, copy to new
        if fOld.readable() :
            with fOld : fNew.write(fOld.read())

    fOld.close()                    #-- Close old and delete
    os.remove(fOld.name)

#== Return New File ===================================================

    return fNew

#========================================================================
#  Replace a substring between two markers, i.e. [abc] -> [xyz]
#========================================================================

def replaceSS(oldStr,lMark,rMark,newSS) :

    beg = oldStr.find(lMark) + 1
    end = oldStr.find(rMark)

    remStr = oldStr[beg:end]
    newStr = oldStr[:beg] + newSS + oldStr[end:]

#== Return 'removed' string and new sting =============================    

    return remStr,newStr

#========================================================================
#  End of Module
#========================================================================


