
#========================================================================
#  PVTFree Main File
#  See Usage Details in aboutPVTfree in utilities.py
#========================================================================

#!/usr/bin/python3

import os
import sys

import numpy       as NP

import allData     as AD
import calcPhsPlt  as CF
import calcProps   as CP
import readBlk     as RB
import readComp    as RC
import readDatFile as RD
import readGen     as RE
import readGrp     as RG
import readReg     as RR
import readExps    as RX
import units       as UN
import utilities   as UT
import writeOut    as WO

from math import sqrt

#========================================================================
#  Main Routine: Parses the User Dataset
#========================================================================

def PVTfree():

    iERR,clsIO = checkInputFile()

    if iERR < 0 : return iERR

    clsIO.setDocStr(UT.aboutPVTfree.__doc__)  #-- DocString => clsIO

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

    macEPS = UT.calcMacEPS()

    UT.macEPS = macEPS      #-- Over-write the value set in constants.py

#---------------------------------------------------------------------
#  Parse Input File Line by Line
#---------------------------------------------------------------------

    fInp = clsIO.fInp
    fOut = clsIO.fOut

    iCnt  = 0

    for curL in fInp :
    
        iCnt += 1

#-- Current Line in curL; split into List tokS ----------------------
    
        tokS = curL.split( )
        nTok = len(tokS)
      
#-- Blank or Comment? -----------------------------------------------

        if   nTok == 0           : pass #-- Blank line
        elif tokS[0][:2] == "--" : pass #-- Comment!

#-- Initial Section -------------------------------------------------

        elif tokS[0][:4].upper() == "INIT" :

            iERR,clsEOS = RD.readINIT(dicSAM,clsUNI,clsIO)
            if iERR < 0 : break
                
#-- Is this a DEF/ENDDEF fluid definition? --------------------------

        elif tokS[0][:3].upper() == "DEF" :

            fDef = clsIO.fInp
            iERR,clsEOS,dicSAM = RD.readDEF(fDef,dicSAM,clsIO)
            if iERR < 0 : break

#-- Read Experiments Definitions ------------------------------------

        elif tokS[0].upper() == "EXP" :

            iERR = RX.readExps(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
            if iERR < 0 : break
                    
#-- Read Regression Definitions -------------------------------------

        elif tokS[0].upper() == "REG" :

            iERR,clsEOS,dicSAM,dicEXP = \
                RR.readReg(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
            if iERR < 0 : break

#-- Read Grouping Definitions ---------------------------------------

        elif tokS[0][:2].upper() == "GR" :

            iERR,clsEOS,dicSAM,dicEXP = \
                RG.readGroup(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
            if iERR < 0 : break

#-- Read Blackoil Definitions ---------------------------------------

        elif tokS[0][:2].upper() == "BL" :

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

#-- Read OPTION controls --------------------------------------------

        elif tokS[0][:3].upper() == 'OPT' :

            iERR = RE.readOption(clsIO)
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

    fOut = clsIO.fOut
    fSav = clsIO.fSav

    if clsIO.qDeb : clsIO.qDeb = UT.closeFile(clsIO.fDeb)
    if clsIO.qReg : clsIO.qReg = UT.closeFile(clsIO.fReg)
    if clsIO.q100 : clsIO.q100 = UT.closeFile(clsIO.f100)
    if clsIO.qIMX : clsIO.qIMX = UT.closeFile(clsIO.fIMX)
    if clsIO.q300 : clsIO.q300 = UT.closeFile(clsIO.f300)
    if clsIO.qGEM : clsIO.qGEM = UT.closeFile(clsIO.fGEM)
    if clsIO.qVIP : clsIO.qVIP = UT.closeFile(clsIO.fVIP)
    if clsIO.qMOR : clsIO.qMOR = UT.closeFile(clsIO.fMOR)

#-- Close the Files -------------------------------------------------    
    
    fOut.close()
    fSav.close()
    
#== No return values ==================================================
    
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

    if nArg > 1 :
        
        myFile = str(sys.argv[1])
        dPath,rootN = UT.procFileName(myFile)
        if dPath == "" : dPath = os.getcwd() + "\\"     #-- Null Path?

    else :

#== No, will need to get the root-name via user input =================
        
        print(" ")
        rootN = input("Enter the Rootname of the Input File: ")
        dPath = os.path.dirname(os.path.realpath(__file__))

#== Path and Rootname - Extensions Next ===============================        

    qTem = False
    
    iERR,clsIO = UT.createIOFiles(dPath,rootN,qTem)

    return iERR,clsIO

#========================================================================
#  Calls the Main Routine
#========================================================================

iERR = PVTfree()

#========================================================================
#  End of Module
#========================================================================

