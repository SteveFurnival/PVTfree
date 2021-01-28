
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Data for Component Grouping Module
#
#=======================================================================

#!/usr/bin/python3

from   copy     import deepcopy

import allData   as AD
import calcExps  as CX
import calcGrp   as CG
import plotBatch as PB
import readGen   as RG
import writeOut  as WO

#========================================================================
#  Read the New/Old Component Names
#========================================================================

def readGroup(clsIO,clsEOS0,dicSAM0,dicEXP0,clsUNI) :

    print("Reading User Data for Grouping")

    iERR  = 0
    iLine = 0
    
    nOld  = clsEOS0.nComp   #-- "Old" Number of Components in "Old" Class
    nNew  = 0
    nTst  = 0
    iSam  = 0               #-- By default, we will weight by Sample[0]
    qExp  = False

    fInp = clsIO.fInp

    dicGRP = {}             #-- Dictionary to hold the new Groups

#== Loop over lines in the Input File =================================    

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of GROUP Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok                == 0      : pass           #-- Blank line!
        elif tokS[0][:2]         == "--"   : pass           #-- Comment!
        elif tokS[0][:4].upper() == "ENDG" : break          #-- ENDGROUP => Exit
        elif tokS[0][:3].upper() == "DEB"  :                #-- DEBUG
             iERR = RG.readDebug(clsIO)
             if iERR < 0 : break
        elif tokS[0][:3].upper() == "OPT"  :                #-- OPTION
             iERR = RG.readOption(clsIO)
             if iERR < 0 : break
        elif tokS[0][:4].upper() == "RUNE" : qExp = True    #-- Run Experiments
        elif tokS[0][:2].upper() == "SA"   :                #-- SAMPLE for weighting

            samNam = tokS[1]
            iSam   = checkSamp4Grp(samNam,dicSAM0)
            if iSam < 0 :
                print("Sample ",samNam," Proposed for Weighting Not Found - Error")
                iERR = -1
                return iERR

        else :

#== Must be reading New Name and Old-Names ============================

            newNam = tokS[0]
            clsGRP = AD.classGRP(newNam)   #-- Create the New Group

            dicGRP[nNew] = clsGRP          #-- Add to the Dictionary of Groups
            
            nNew += 1                      #-- Increment the Counter

            oldNam = []

            iTok = 1
            while iTok < nTok :
                oldNam.append(tokS[iTok])
                nTst += 1                   #-- Test the number of "old" comps
                iTok += 1

            clsGRP.setOldNames(oldNam)

#----------------------------------------------------------------------
#  Process Data
#----------------------------------------------------------------------

    #print("nOld,nNew,nTst ",nOld,nNew,nTst)

    cOK = [False for i in range(nOld)]
    
    if nTst != nOld :
        print("Currently have ",nOld," Components defined: Only ",nTst," Read - Error")
        iERR = -1
        return iERR,clsEOS0,dicSAM0,dicEXP0

    for iNew in range(nNew) :
        
        oldNam = dicGRP[iNew].oldNam
        oldNum = []
        nGrp   = len(oldNam)
        
        for iGrp in range(nGrp) :
            sCom = oldNam[iGrp]
            iCom = checkComp4Grp(sCom,clsEOS0)
            if iCom < 0 :
                print("Old Component Name ",sCom," Not Found - Error")
                iERR = -1
                return iERR,clsEOS0,dicSAM0,dicEXP0
            else:
                oldNum.append(iCom)
                cOK[iCom] = True
                
        dicGRP[iNew].setOldNums(oldNum)

#-- Do we have all 'old' components used? ---------------------------
        
    for iC in range(nOld) :
        if not cOK[iC] :
            sCom = clsEOS0.gNM(iC)
            print("Component " + sCom + " Not Found Amongst 'Old' Components - Error")
            iERR = -1
            return iERR,clsEOS0,dicSAM0,dicEXP0
        
#==  Do we want experiments run before and after grouping? ============

    if qExp :
        dicEXPG = deepcopy(dicEXP0)
        sExt   = " - Before Grouping"
        CX.calcExps(sExt,clsEOS0,dicSAM0,dicEXP0,clsUNI,clsIO)

#== Create 'Group' Copies of EOS/Samples ==============================        

    EOS = clsEOS0.EOS

    clsEOSG = AD.classEOS(EOS)
    dicSAMG = {}

    clsEOSG.setEOSparms(EOS)

#======================================================================
#  Perform the Grouping
#======================================================================

    iERR,clsEOSG,dicSAMG = \
        CG.calcGroup(iSam,dicGRP,clsEOS0,dicSAM0)

    if iERR == 0 :
        print("Grouping Performed")
    else :
        print("Grouping Failed: Program Stopping")
        return iERR,clsEOS0,dicSAM0,dicEXP0

#-- Write New System to the Output File -----------------------------

    sTit = "Grouping From " + str(nOld) + " Components To " + str(nNew)

    WO.outputProps(sTit,clsEOSG,dicSAMG,clsIO)

#-- Save the new EoS/Sample Definitions -----------------------------    

    WO.outputSave(sTit,clsEOSG,dicSAMG,clsIO)

#----------------------------------------------------------------------
#  Re-run existing experiments?
#----------------------------------------------------------------------

    if qExp :
        sExt   = " - After  Grouping"
        for iExp in range(len(dicEXPG)) : dicEXPG[iExp].IsAct = True
        CX.calcExps(sExt,clsEOSG,dicSAMG,dicEXPG,clsUNI,clsIO)
        qReg = False
        PB.regPlots(clsIO,dicEXP0,dicEXPG,dicSAMG,qReg,clsUNI)

#----------------------------------------------------------------------
#  Swap Deepcopy New Class/Dictionary into Old
#----------------------------------------------------------------------

    clsEOS0 = deepcopy(clsEOSG)
    dicSAM0 = deepcopy(dicSAMG)

#----------------------------------------------------------------------
#  Perform deepcopy to put Grouped EOS/Samples into Main Sets
#----------------------------------------------------------------------

    if qExp : dicEXP0 = deepcopy(dicEXPG)

#== End of Routine ====================================================

    return iERR,clsEOS0,dicSAM0,dicEXP0

#======================================================================
#  Check a Sample Name is Valid
#======================================================================

def checkSamp4Grp(sName,dicSAM) :

    iSamp = -1
    nSamp = len(dicSAM)

    for iS in range(nSamp) :
        sThis = dicSAM[iS].sNam
        if sThis == sName :
            iSamp = iS
            break

    return iSamp

#======================================================================
#  Check a Component Name is Valid
#======================================================================

def checkComp4Grp(sName,clsEOS) :

    iComp = -1
    nComp = clsEOS.nComp

    for iC in range(nComp) :
        sThis = clsEOS.gNM(iC)
        if sThis == sName :
            iComp = iC
            break

    return iComp

#======================================================================
#  End of Module
#======================================================================
