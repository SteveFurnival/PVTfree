
#=======================================================================
#
#  PVTfree Equation of State (EOS) PVT Program
#
#  Copyright (C) 2019
#
#  This file is part of PVTfree.
#
#  PVTfree is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PVTfree is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with PVTfree.  If not, see <https://www.gnu.org/licenses/>.
#
#  Steve Furnival, HoBoil Limited, Banchory, UK, June 2019
#
#  email: dev.pvtfree@gmail.com
#
#========================================================================

#!/usr/bin/python3

import os
import sys
from pathlib import Path

import numpy      as NP

import constants  as CO
import calcPhsPlt as CP
import calcProps  as CR
import readBlk    as RB
import readComp   as RC
import readGen    as RE
import readGrp    as RG
import readReg    as RR
import readSamp   as RS
import readExps   as RX
import units      as UN
import writeOut   as WO

from math import sqrt

#========================================================================
#  Equation of State Class
#========================================================================

class classEoS :

    KIJ = []  #-- Binary Interaction Parameters

    PRP = []  #-- Component Names, Mw's, Tc's, etc.
    
    nPR = 17  #-- Number of Arrays in PRP
              #-- List of Arrays in pMap below

    #  pMap     Array   Array                   Unit
    #  Index    Code    Name                    Type
    #
    #   0       CN      Name
    #   1       MW      Mole Weight             lb/lbmol
    #   2       TC      Critical Temperature    degR
    #   3       PC      Critical Pressure       psia
    #   4       VC      Critical Volume         ft3/lbmol
    #   5       ZC      Critical Z-Factor       dim
    #   6       AF      Acentric Factor         dim
    #   7       TB      Normal Boiling Point    degR
    #   8       SG      Specific Gravity        wrtAir=1
    #   9       PA      Parachor                dim
    #  10       SS      Volume Shift            dim
    #  11       MA      Omega-A Multiplier      dim
    #  12       MB      Omega-B Multiplier      dim
    #  13       CA      Ideal Gas Cp A-Coeff    Btu/lbmol.degR
    #  14       CB      Ideal Gas Cp B-Coeff    Btu/lbmol.degR^2
    #  15       CC      Ideal Gas Cp C-Coeff    Btu/lbmol.degR^3
    #  16       CD      Ideal Gas Cp D-Coeff    Btu/lbmol.degR^4
    
    pMap = {"CN": 0, "MW": 1, "TC": 2, "PC": 3, "VC": 4, \
            "ZC": 5, "AF": 6, "TB": 7, "SG": 8, "PA": 9, \
            "SS":10, "MA":11, "MB":12, "CA":13, "CB":14, \
            "CC":15, "CD":16                             }
    
    def __init__(self,EOS) :

        self.Tref = 536.67  #-- Ref Temp (for Ideal Gas Enthalpy) [25 degC -> degR]

        self.EOS  = EOS

        if self.EOS == "SRK" :          #-- User set Soave-Redlich-Kwong (SRK) EoS
            
            self.n1 = 1.0               #-- Coeffs in the Martin EoS
            self.n2 = 2.0
            
            self.OA = 0.4274802335403410        #-- OmegaA
            self.OB = 0.0866403499649577        #-- OmegaB

            self.w0 =  0.47979          #-- Coeffs in the alpha(AcenF)       
            self.w1 =  1.576
            self.w2 = -0.1925
            self.w3 =  0.025

        else :                          #-- By default, Peng-Robinson (PR) EoS
            
            sqr2    = 1.4142135623730951

            self.n1 = 2.0 + sqr2
            self.n2 = 2.0 - sqr2
            
            self.OA = 0.45723552892138225
            self.OB = 0.07779607390388849

            self.w0 =  0.379642
            self.w1 =  1.48503
            self.w2 = -0.164423
            self.w3 =  0.016666


        self.n3 = self.n1 + self.n2         #-- Useful Work Coeffs
        self.n4 = self.n1 * self.n2

        self.n5 = 1.0/(self.n1 - self.n2)

#-- Set the Number of Components ------------------------------------

    def setNComp(self,nComp) :

        self.NC = nComp  #-- Now dimension the arrays

        #self.PRP = NP.zeros((self.nPR,self.NC))
        self.KIJ = NP.zeros((self.NC,self.NC))
        self.LIJ = NP.zeros((self.NC,self.NC))
        
        self.PRP = [[0.0 for i in range(self.NC)] for j in range(self.nPR)]

    def incNComp(self,newNC) :

        nEXT = newNC - self.NC              #-- Number of Elements to extend by
        
        aEXT = NP.zeros(nEXT)
        bEXT = NP.zeros(newNC)

        for iPRP in range(self.nPR) :       #-- Extend the Parameter Arrays
            self.PRP[iPRP].extend(aEXT)

        for iC in range(self.NC) :          #-- Extend existings rows of KIJ
            self.KIJ[iC].extend(aEXT)
            self.LIJ[iC].extend(aEXT)
            
        for iC in range(self.NC,newNC) :    #-- Add extra rows of KIJ
            self.KIJ.append(bEXT)
            self.LIJ.append(bEXT)

        self.NC = newNC                     #-- Update N-Components
            
#-- Set the Number of Samples ---------------------------------------

    def setNSamp(self,nSamp) : self.NS = nSamp
            
#-- Set the Number of User Components -------------------------------

    def setNUser(self,nUser) : self.NU = nUser
            
#-- Set the Number of Pseudo-Components (in Splitting) --------------

    def setNPseu(self,nPseu) : self.NP = nPseu

#-- Set/Get a component Mole Weight, etc ----------------------------

    def sPP(self,typ,iC,valU) :

        iMap = self.pMap[typ]
        self.PRP[iMap][iC] = valU

    def gPP(self,typ,iC) :

        iMap = self.pMap[typ]
        valU = self.PRP[iMap][iC]

        return valU
             
#-- Set/Get a BIC ---------------------------------------------------

    def sIJ(self,iC,jC,KIJ) :

        self.KIJ[iC][jC] =       KIJ
        self.LIJ[iC][jC] = 1.0 - KIJ

    def gIJ(self,iC,jC) :

        valU = self.KIJ[iC][jC]

        return valU

    def gLIJ(self,iC,jC) :

        valU = self.LIJ[iC][jC]

        return valU

#-- Construct the no-(P,T) EoS Coefficients -------------------------    

    def eosCoefsNoPT(self) :

        nCom = self.NC

        self.aP = NP.zeros(nCom)
        self.aQ = NP.zeros(nCom)
        self.bI = NP.zeros(nCom)
        self.vS = NP.zeros(nCom)
        self.Mw = NP.zeros(nCom)

        OA0 = self.OA
        OB0 = self.OB

        iTc = self.pMap["TC"]
        iPc = self.pMap["PC"]
        iAF = self.pMap["AF"]
        iMA = self.pMap["MA"]
        iMB = self.pMap["MB"]
        iMw = self.pMap["MW"]
        iSS = self.pMap["SS"]

        for iC in range(nCom) :
            
            Tc = self.PRP[iTc][iC]
            Pc = self.PRP[iPc][iC]
            AF = self.PRP[iAF][iC]
            MA = self.PRP[iMA][iC]
            MB = self.PRP[iMB][iC]

            OA = OA0*MA
            OB = OB0*MB

            mAcF = self.w0 + AF*(self.w1 + AF*(self.w2 + AF*self.w3))

            P = 1.0 + mAcF
            Q =     - mAcF/sqrt(Tc)

            self.bI[iC] = OB*Tc/Pc

            aSqr = Tc*sqrt(OA/Pc)

            self.aP[iC] = aSqr*P
            self.aQ[iC] = aSqr*Q

            self.Mw[iC] = self.PRP[iMw][iC]
            self.vS[iC] = self.PRP[iSS][iC]

#-- Update the 1-Kij Array ------------------------------------------

        self.LIJ = NP.subtract(1.0,self.KIJ)

#-- Dummy Scalar, Vector and Matrix for EoS Coef Calculations -------

        self.dumS = 0.0
        self.dumV = NP.zeros( nCom )
        self.dumM = NP.zeros((nCom,nCom))

        return

#========================================================================
#  The cInpOut class handles input and output
#  Including graphics and debug
#========================================================================

class classInpOut :

    def __init__(self,rNam,patR,fInP,fOut,fSav) :

        self.rNam = rNam  #-- Rootname for the Input File
        self.patR = patR  #-- Path & Rootname

        self.fInP = fInP  #-- Input       File
        self.fOut = fOut  #-- Output      File
        self.fSav = fSav  #-- Save        File
        
        self.Pall = True  #-- Plot All Experiments/Observations

        self.Deb  = {"PSAT" : 0, "FLASH": 0, "STAB2": 0, "STAB1": 0, \
                     "RR"   : 0, "REG"  : 0, "BLACK": 0, "GROUP": 0, \
                     "CCE"  : 0, "CVD"  : 0, "DLE"  : 0, "SEP"  : 0, \
                     "SAT"  : 0, "FLS"  : 0, "SWL"  : 0, "GRD"  : 0, \
                     "ENV"  : 0, "TSAT" : 0}

        self.Tit1 = None
        self.Tit2 = None
        self.Tit3 = None

    def setQXL( self,qXL ) : self.qXL  = qXL

    def setQREG(self,qReg) : self.qReg = qReg
    def setFREG(self,fReg) : self.fReg = fReg

    def setQDEB(self,qDeb) : self.qDeb = qDeb
    def setFDEB(self,fDeb) : self.fDeb = fDeb

    def setQ100(self,q100) : self.q100 = q100
    def setF100(self,f100) : self.f100 = f100

    def setQIMX(self,qIMX) : self.qIMX = qIMX
    def setFIMX(self,fIMX) : self.fIMX = fIMX

    def setQMOR(self,qMOR) : self.qMOR = qMOR
    def setFMOR(self,fMOR) : self.fMOR = fMOR

    def setQVIP(self,qVIP) : self.qVIP = qVIP
    def setFVIP(self,fVIP) : self.fVIP = fVIP

    def setQ300(self,q300) : self.q300 = q300
    def setF300(self,f300) : self.f300 = f300

    def setQGEM(self,qGEM) : self.qGEM = qGEM
    def setFGEM(self,fGEM) : self.fGEM = fGEM

def PVTfree():

#========================================================================
#  Main Routine: Parses the User Dataset
#========================================================================

    iERR,clsIO = checkInputFile()

    if iERR < 0 : return iERR

#--------------------------------------------------------------------
#  Dictionaries to Hold Classes of Samples and Experiments
#  Dictionary of Regression Classes is held locally
#--------------------------------------------------------------------

    dicSAM = {}
    dicEXP = {}

#-- Initialise the Units class --------------------------------------

    clsUNI = UN.classUnits()

#---------------------------------------------------------------------
#  Get Machine Epsilon
#---------------------------------------------------------------------

    macEPS = calcMacEPS()

    CO.macEPS = macEPS      #-- Over-write the value set in constants.py

#---------------------------------------------------------------------
#  Parse Input File Line by Line
#---------------------------------------------------------------------

    fInP = clsIO.fInP
    fOut = clsIO.fOut

    iCnt  = 0

    with fInP :

        for curL in fInP:
    
            iCnt += 1

#-- Current Line in curL; split into List tokS ----------------------
    
            tokS = curL.split( )
            nTok = len(tokS)
      
#-- Blank or Comment? -----------------------------------------------

            if   nTok == 0           : pass #-- Blank line
            elif tokS[0][:2] == "--" : pass #-- Comment!

#-- Initial Section -------------------------------------------------

            elif tokS[0][:4].upper() == "INIT" :

                iERR,clsEOS,dicSAM = readInit(clsIO,dicSAM,clsUNI)
                if iERR < 0 : break
                
#-- Is this a DEF/ENDDEF fluid definition? --------------------------

            elif tokS[0][:3].upper() == "DEF" :

                iERR,clsEOS,dicSAM = readDef(clsIO,dicSAM,clsUNI)
                if iERR < 0 : break

#-- Read Experiments Definitions ------------------------------------

            elif tokS[0].upper() == "EXP" :

                iERR,dicEXP = RX.readExps(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
                if iERR < 0 : break
                    
#-- Read Regression Definitions -------------------------------------

            elif tokS[0].upper() == "REG" :

                iERR,clsEOS,dicSAM,dicEXP = \
                    RR.readReg(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
                if iERR < 0 : break

#-- Read Grouping Definitions ---------------------------------------

            elif tokS[0][:2].upper() == "GR" :

                EOS = clsEOS.EOS
                    
                clsEOSG = classEoS(EOS)  #-- New Class
                dicSAMG = {}             #-- New Dictionary
                    
                iERR,clsEOS,dicSAM,dicEXP = \
                    RG.readGroup(clsIO,clsEOS,clsEOSG,dicSAM,dicSAMG,dicEXP,clsUNI)
                if iERR < 0 : break

#-- Read Blackoil Definitions ---------------------------------------

            elif tokS[0][:2].upper() == "BL" :

                #clsIO.fSim = clsIO.rNam + "_BO.sim"
                iERR = RB.readBlack(clsIO,clsEOS,dicSAM,clsUNI)
                if iERR < 0 : break

#-- Read Compositional Definitions ----------------------------------

            elif tokS[0][:2].upper() == "CO" :

                clsIO.fSim = clsIO.rNam + "_EOS.sim"
                iERR = RC.readComp(clsEOS,dicSAM,clsIO,clsUNI)
                if iERR < 0 : break

#-- Read DEBUG controls ---------------------------------------------

            elif tokS[0][:3].upper() == "DEB" :

                iERR = RE.readDebug(clsIO)
                if iERR < 0 : break

#-- Read TITLE1 -----------------------------------------------------

            elif tokS[0].upper() == "TITLE1" : clsIO.Tit1 = buildTitle(tokS)                  

#-- Read TITLE2 -----------------------------------------------------

            elif tokS[0].upper() == "TITLE2" : clsIO.Tit2 = buildTitle(tokS)                  

#-- Read TITLE3 -----------------------------------------------------

            elif tokS[0].upper() == "TITLE3" : clsIO.Tit3 = buildTitle(tokS)                  

#-- Temporary STOP Card ---------------------------------------------

            elif tokS[0].upper() == "STOP"   : break
                    
            else                             : pass

#== Close All Output Files =============================================

    closeFiles(clsIO)

#== Return the Error Flag ==============================================    
                
    return iERR

#========================================================================
#  Close All Output Files and Delete the Empty Ones
#========================================================================

def closeFiles(clsIO) :

    if clsIO.qXL :
        CR = input("Click Return to close PVTfree.py when ready")

    fOut = clsIO.fOut
    fSav = clsIO.fSav

    if clsIO.qDeb : closeFile(clsIO.fDeb)
    if clsIO.qReg : closeFile(clsIO.fReg)
    if clsIO.q100 : closeFile(clsIO.f100)
    if clsIO.qIMX : closeFile(clsIO.fIMX)
    if clsIO.q300 : closeFile(clsIO.f300)
    if clsIO.qGEM : closeFile(clsIO.fGEM)
    if clsIO.qVIP : closeFile(clsIO.fVIP)

#-- Close the Files -------------------------------------------------    
    
    fOut.close()
    fSav.close()
    
#== No return values ==================================================
    
    return

def closeFile(fStr) :

    fStr.close()            #-- Close the file
    
    fStrNam = fStr.name     #-- Get its name

#-- If the file is zero-length, delete it ---------------------------    
    
    if os.path.getsize(fStrNam) == 0 : os.remove(fStrNam)

#-- No return value -------------------------------------------------

    return

#========================================================================
#  Build a TITLE string
#========================================================================

def buildTitle(tokS) :

    nTok = len(tokS)

    sTit = "--  "     #-- Eclipse-Style Comment
    
    iTok = 1          #-- First token was the TITLEx identifier
    
    while iTok < nTok :
        sTit = sTit + tokS[iTok] + " "
        iTok += 1

    return sTit

#========================================================================
#  Input/Output File Checking
#========================================================================

def checkInputFile() :

    iERR = 0

#== Has the rootname been given on the command line? ==================    

    nArg = len(sys.argv)
    print(nArg)

    if nArg > 1 :

#-- Yes, process the string to get Path & Root-Name -----------------        
        
        myFile = str(sys.argv[1])

        sPlt = myFile.split("\\")
        nSpl = len(sPlt)

        sPath = ""
        for iSpl in range(nSpl-1):
            sPath = sPath+ sPlt[iSpl] + "\\"

        sRoo = sPlt[nSpl-1]
        sPlt = sRoo.split(".")
        nSpl = len(sPlt)

        rootName = sPlt[0]

        pathR = Path(sPath + rootName)

    else :

#-- No, will need to get the root-name via user input ---------------        
        
        print(" ")
        
        rootName = input("Enter the Rootname of the Input File: ")
        
        dPath = os.path.dirname(os.path.realpath(__file__))
        pathR = Path(dPath) /rootName

#== Build the File Names ==============================================
        
    pathInp = pathR.with_suffix(".dat")
    pathOut = pathR.with_suffix(".out")
    pathSav = pathR.with_suffix(".sav")

    if os.path.isfile(pathInp) and os.access(pathInp,os.R_OK) :
        print("checkInputFile: Opened Input File ",pathInp," and is readable")
    else :
        print("checkInputFile: Input File ",pathInp," is missing or is not readable")
        clsIO = None
        iERR = -1
        return iERR,clsIO

    fInP = open(pathInp,'r')
    fOut = open(pathOut,'w')
    fSav = open(pathSav,'w')

#== Create the I/O Class ==============================================    

    clsIO = classInpOut(rootName,pathR,fInP,fOut,fSav)

#-- Is this instance being run from Excel front-end? ----------------    

    if nArg > 2 : qXL = True
    else        : qXL = False

    clsIO.setQXL(qXL)

#-- Other Files not set ---------------------------------------------

    clsIO.setQDEB(False)
    clsIO.setQREG(False)
    clsIO.setQ100(False)
    clsIO.setQIMX(False)
    clsIO.setQMOR(False)
    clsIO.setQVIP(False)
    clsIO.setQ300(False)
    clsIO.setQGEM(False)

    return iERR,clsIO

#=======================================================================
#  Read Initial Section
#=======================================================================

def readInit(clsIO,dicSAM,clsUNI) :

    iERR  = 0
    iCnt  = 0
    nSplt = 0

    fInP = clsIO.fInP
    fOut = clsIO.fOut

#----------------------------------------------------------------------
#  Parse Input File Line by Line
#----------------------------------------------------------------------

    for curL in fInP :
    
        iCnt += 1

#-- Current Line in curL; split into List tokS ----------------------
    
        tokS = curL.split()
        nTok = len(tokS)

#== Process Options ===================================================

        if   nTok                == 0      : pass   #-- Blank Line
        elif tokS[0][:2]         == "-- "  : pass   #-- Comment
        elif tokS[0][:4].upper() == "ENDI" : break  #-- ENDINIT k/w => Exit

#-- EOS K/W Read => Create the clsEOS to hold data ------------------
        
        elif tokS[0].upper()     == "EOS"  :

            EOS = tokS[1].upper()

            clsEOS = classEoS(EOS)

            print("Equation of State Specified as ",EOS)

            sCom = "--"

            WO.outputHeader(fOut,sCom,clsIO)
            WO.outputEOS(fOut,sCom,clsIO,clsEOS)

#-- SPLIT k/w -------------------------------------------------------
                    
        elif tokS[0][:3].upper() == "SPL"  :

            nSplt = int(tokS[1])
            #print("SPLIT: nSplt ",nSplt)
            if nSplt < 2 or nSplt > 5 :
                print("SPLIT: nSplt ",nSplt," Out of Range, 2 =< Nsplt =< 5 - Error")
                iERR = -1
                break
            else:
                clsEOS.setNPseu(nSplt)
                print("Plus Fraction will be split into ",nSplt," Pseudo-Components")

#-- SAMPLES k/w -----------------------------------------------------
                    
        elif tokS[0][:4].upper() == "SAMP" :

            if nSplt == 0 :
                nSplt = 1
                clsEOS.setNPseu(nSplt)
                
            iERR,dicSAM = RS.readSamp(clsIO,tokS,clsEOS,dicSAM,clsUNI)
            if iERR < 0 : break

#-- DEBUG k/w -------------------------------------------------------

        elif tokS[0][:3].upper() == "DEB"  :

            iERR = RE.readDebug(clsIO)
            if iERR < 0 : break

#== Return values =====================================================

    return iERR,clsEOS,dicSAM

#=======================================================================
#  Read Define Section (Full Compositional Description)
#=======================================================================

def readDef(clsIO,dicSAM,clsUNI) :

    iERR = 0
    iCnt = 0
    fInP = clsIO.fInP
    fOut = clsIO.fOut

    macEPS = calcMacEPS()   #-- Compute Machine Epsilon

#---------------------------------------------------------------------
#  Parse Input File Line by Line
#---------------------------------------------------------------------

    for curL in fInP :
    
        iCnt += 1

#-- Current Line in curL; split into List tokS
    
        tokS = curL.split( )
        nTok = len(tokS)
      
#-- Blank or Comment? -----------------------------------------------

        if nTok == 0 :  #-- Blank line
            pass
        else :

            if   tokS[0][:2] == "--"           : pass
            elif tokS[0][:4].upper() == "ENDD" : break
            elif tokS[0][:3].upper() == "DEB"  :
                iERR = RE.readDebug(clsIO)
            elif tokS[0].upper()     == "EOS"  :

#-- Equation of State -----------------------------------------------                    

                EOS = tokS[1].upper()

                clsEOS = classEoS(EOS)
                    
                print("Equation of State Specified as ",EOS)

                sCom = "--"

                WO.outputHeader(fOut,sCom,clsIO)
                WO.outputEOS(fOut,sCom,clsIO,clsEOS)

            elif tokS[0].upper() == "NCOMP" :

                nComp = int(tokS[1])
                clsEOS.NC = nComp
                clsEOS.setNComp(nComp)  #-- Dimension the Arrays
                    
            elif tokS[0].upper() == "NSAMP" :

                nSamp = int(tokS[1])
                clsEOS.NS = nSamp
                    
            elif tokS[0].upper() == "PROPS" :

                curL = next(fInP)  #-- Property Names
                curL = next(fInP)  #-- Property Units

                for iC in range(nComp) :
                    
                    curL = next(fInP)
                    tokS = curL.split( )
                    nTok = len(tokS)

                    if nTok != 17 :
                        print("DEFINE PROPS: Expecting 17 Columns, Only ",nTok," Read")
                        iERR = -1
                        break

                    sN =       tokS[ 0]
                    Mw = float(tokS[ 1])   ; Tc = float(tokS[ 2])
                    Pc = float(tokS[ 3])   ; Vc = float(tokS[ 4])
                    Zc = float(tokS[ 5])   ; AF = float(tokS[ 6])
                    Tb = float(tokS[ 7])   ; SG = float(tokS[ 8])
                    PA = float(tokS[ 9])   ; SS = float(tokS[10])
                    MA = float(tokS[11])   ; MB = float(tokS[12])
                    CA = float(tokS[13])   ; CB = float(tokS[14])
                    CC = float(tokS[15])   ; CD = float(tokS[16])

                    clsEOS.sPP("CN",iC,sN)
                    clsEOS.sPP("MW",iC,Mw) ; clsEOS.sPP("TC",iC,Tc)
                    clsEOS.sPP("PC",iC,Pc) ; clsEOS.sPP("VC",iC,Vc)
                    clsEOS.sPP("ZC",iC,Zc) ; clsEOS.sPP("AF",iC,AF)
                    clsEOS.sPP("TB",iC,Tb) ; clsEOS.sPP("SG",iC,SG)
                    clsEOS.sPP("PA",iC,PA) ; clsEOS.sPP("SS",iC,SS)
                    clsEOS.sPP("MA",iC,MA) ; clsEOS.sPP("MB",iC,MB)
                    clsEOS.sPP("CA",iC,CA) ; clsEOS.sPP("CB",iC,CB)
                    clsEOS.sPP("CC",iC,CC) ; clsEOS.sPP("CD",iC,CD)
                
            elif tokS[0].upper() == "BIP"   :

                curL = next(fInP)  #-- Component Names
                
                for iC in range(nComp) :
                    
                    curL = next(fInP)
                    tokS = curL.split( )
                    nTok = len(tokS)
                
                    if nTok != nComp + 1 :
                        print("DEFINE BIP: Expecting ",nComp+1," Columns, Only ",nTok," Read")
                        iERR = -1
                        break

                    iTok = 1
                    while iTok < nTok :
                        KIJ = float(tokS[iTok])
                        clsEOS.sIJ(iC,iTok-1,KIJ)
                        iTok += 1

            elif tokS[0].upper() == "SAMPLES" :
                
                curL = next(fInP)  #-- Sample Names
                tokS = curL.split()
                nTok = len(tokS)

                for iSamp in range(nSamp) :
                    sName = tokS[iSamp+1]
                    csSAM = RS.classSample(sName)
                    dicSAM[iSamp] = csSAM
                    dicSAM[iSamp].setNComp(nComp)
                
                for iC in range(nComp) :
                    
                    curL = next(fInP)
                    tokS = curL.split( )

                    if nTok != nSamp + 1 :
                        print("DEFINE SAMPLES: Expecting ",nSamp+1," Columns, Only ",nTok," Read")
                        iERR = -1
                        break
                    
                    for iSamp in range(nSamp) :
                        ZI = float(tokS[iSamp+1])
                        if ZI < macEPS : ZI = 1.0E-20  #-- Protect against Z = 0
                        dicSAM[iSamp].sZI(iC,ZI)

            else :
                pass

#========================================================================
#  Has any slop crept into the sample definitions?
#========================================================================

    for iSamp in range(nSamp) :
        sumT = 0.0
        for iC in range(nComp) :
            sumT = sumT + dicSAM[iSamp].gZI(iC)
        sumT = 1.0/sumT
        for iC in range(nComp) :
            zI = sumT*dicSAM[iSamp].gZI(iC)
            dicSAM[iSamp].sZI(iC,zI)

#== Back-calculate the C7+ Properties =================================

    backCalcPlusFracProps(clsEOS,dicSAM)
            
#========================================================================
#  Do we need to sort the components?  Most to Least Volatility
#========================================================================

    clsEOS,dicSAM = CR.sortComponents(clsEOS,dicSAM)

#========================================================================
#  Output the data
#========================================================================

    sTit = "Initialisation from a saved DEFINE"

    WO.outputProps(clsIO,clsEOS,dicSAM,sTit)

#========================================================================
#  Write Fluid Description to the SAV file
#========================================================================

    WO.outputSave(sTit,clsEOS,dicSAM,clsIO)

#======================================================================
#  Generate (Approximate) Phase Plots
#======================================================================

    CP.allSamplesPhasePlot(clsEOS,dicSAM,clsIO)

#== Return values =====================================================    

    return iERR,clsEOS,dicSAM

#========================================================================
#  From Full DEF, back-calculate the Plus Fraction Properties
#========================================================================

def backCalcPlusFracProps(clsEOS,dicSAM) :

    nCom = clsEOS.NC
    nSam = len(dicSAM)

#== For each sample ===================================================    

    for iSam in range(nSam) :

        clsSAM = dicSAM[iSam]

        zPlus = 0.0 ; mPlus = 0.0 ; sPlus = 0.0

#-- For components with Mw > C(6.5) ---------------------------------        

        for iC in range(nCom) :

            Mw = clsEOS.gPP("MW",iC) ; SG = clsEOS.gPP("SG",iC)

            if Mw > 90.0 :

                ZI    = clsSAM.gZI(iC)
                zPlus = zPlus + ZI
                mPlus = mPlus + ZI*Mw
                sPlus = sPlus + ZI*Mw/SG

#-- SG by inverse weighting -----------------------------------------                

        sPlus = mPlus/sPlus ; mPlus = mPlus/zPlus

        clsSAM.mPlus = mPlus
        clsSAM.sPlus = sPlus
        clsSAM.aPlus = 1.0      #-- Default value for ALPHA (info lost)

#== No return values ==================================================

    return

#========================================================================
#  Machine Epsilon
#========================================================================

def calcMacEPS() :

    EPS = 1.0

#-- Keep halving EPS until 1+EPS is indistinguishable from 1 ----------    

    while EPS + 1.0 > 1.0 : EPS = 0.5*EPS

    return 2.0*EPS

#========================================================================
#  Calls the Main Routine
#========================================================================

iERR = PVTfree()

#========================================================================
#  End of Module
#========================================================================

