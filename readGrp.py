
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Data for Component Grouping Module
#
#=======================================================================

#!/usr/bin/python3

import numpy    as NP

from   copy     import deepcopy
from   math     import sqrt

import calcEOS  as CE
import calcExps as CX
import genPlots as GP
import readGen  as RG
import readSamp as RS
import writeOut as WO

class classGRP :

    def __init__(self,newNam) :
        self.newNam = newNam

    def setOldNames(self,oldNam) :
        self.oldNam = []
        self.oldNam = oldNam

    def setOldNums(self,oldNum) :
        self.oldNum = []
        self.oldNum = oldNum

    def setSampNum(self,samNum) :
        self.samNum = samNum

#========================================================================
#  Read the New/Old Component Names
#========================================================================

def readGroup(clsIO,clsEOS0,clsEOSG,dicSAM0,dicSAMG,dicEXP0,clsUNI) :

    print("Reading User Data for Grouping")

    iERR  = 0
    iLine = 0
    
    nOld  = clsEOS0.NC   #-- "Old" Number of Components in "Old" Class
    nNew  = 0
    nTst  = 0
    iSam  = 0           #-- By default, we will weight by Sample[0]
    qExp  = False

    fInP = clsIO.fInP

    dicGRP = {}          #-- Dictionary to hold the new Groups

#== Loop over lines in the Input File =================================    

    for curL in fInP :

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
            clsGRP = classGRP(newNam)   #-- Create the New Group

            dicGRP[nNew] = clsGRP       #-- Add to the Dictionary of Groups
            
            nNew += 1                   #-- Increment the Counter

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
            sCom = clsEOS0.gPP("CN",iC)
            print("Component " + sCom + " Not Found Amongst 'Old' Components - Error")
            iERR = -1
            return iERR,clsEOS0,dicSAM0,dicEXP0
        
#==  Do we want experiments run before and after grouping? ============

    if qExp :
        dicEXPG = deepcopy(dicEXP0)
        sExt   = " - Before Grouping"
        CX.calcExps(sExt,clsEOS0,dicSAM0,dicEXP0,clsIO,clsUNI)

#======================================================================
#  Perform the Grouping
#======================================================================

    iERR = calcGroup(iSam,dicGRP,clsEOS0,clsEOSG,dicSAM0,dicSAMG)

    if iERR == 0 :
        print("Grouping Performed")
    else :
        print("Grouping Failed: Program Stopping")
        return iERR,clsEOS0,dicSAM0,dicEXP0

#-- Write New System to the Output File -----------------------------

    sTit = "Grouping From " + str(nOld) + " Components To " + str(nNew)

    #WO.outputGroup(clsIO,nOld,nNew)
    WO.outputProps(clsIO,clsEOSG,dicSAMG,sTit)

#-- Save the new EoS/Sample Definitions -----------------------------    

    sTit = "Grouping From " + str(nOld) + " Components To " + str(nNew)
    WO.outputSave(sTit,clsEOSG,dicSAMG,clsIO)

#----------------------------------------------------------------------
#  Re-run existing experiments?
#----------------------------------------------------------------------

    if qExp :
        sExt   = " - After  Grouping"
        for iExp in range(len(dicEXPG)) : dicEXPG[iExp].IsAct = True
        CX.calcExps(sExt,clsEOSG,dicSAMG,dicEXPG,clsIO,clsUNI)
        qReg = False
        GP.regPlots(clsIO,dicEXP0,dicEXPG,dicSAMG,qReg,clsUNI)

#----------------------------------------------------------------------
#  Perform deepcopy to put Grouped EOS/Samples into Main Sets
#----------------------------------------------------------------------

    clsEOS0 = deepcopy(clsEOSG)
    dicSAM0 = deepcopy(dicSAMG)
    if qExp :
        dicEXP0 = deepcopy(dicEXPG)

#== End of Routine ====================================================

    return iERR,clsEOS0,dicSAM0,dicEXP0

#========================================================================
#  Calculate Grouping
#========================================================================

def calcGroup(sNumW,dicGRP,clsEOS0,clsEOSG,dicSAM0,dicSAMG) :

    iERR = 0

    nOld = clsEOS0.NC
    nNew = len(dicGRP)
    nSam = len(dicSAM0)
    
    tStd = 519.67      #-- Standard Temperature [degR = 60 degF]

    #print("calcGroup: nNew ",nNew)

#----------------------------------------------------------------------
#  Create New EOS-Class and Dictionary of Samples
#----------------------------------------------------------------------

    clsEOSG.NC = nNew
    clsEOSG.NS = nSam

    clsEOSG.setNComp(nNew)   #-- Dimensions the Arrays

#== Samples ===========================================================

    for iSam in range(nSam) :

        clsSAM0 = dicSAM0[iSam]
        
        sOld    = clsSAM0.sName
        clsSAMG = RS.classSample(sOld)
        
        dicSAMG[iSam] = clsSAMG
        dicSAMG[iSam].setNComp(nNew)

        clsSAMG.mPlus = clsSAM0.mPlus
        clsSAMG.sPlus = clsSAM0.sPlus
        clsSAMG.aPlus = clsSAM0.aPlus

        sTot = 0.0

        for iNew in range(nNew) :

            dZI = 0.0

            oldNum = dicGRP[iNew].oldNum
            nGrp   = len(oldNum)

            for iGrp in range(nGrp) :
                iC  = oldNum[iGrp]
                dZI = dZI + dicSAM0[iSam].gZI(iC)

            dicSAMG[iSam].sZI(iNew,dZI)

            sTot = sTot + dZI

        if   sTot - 1.0 >  1.0E-08 :
            print("Composition of Sample " + sOld + " > 1.0: One Old Component Counted Twice?")
            iERR = -1
            return iERR
        elif sTot - 1.0 < -1.0E-08 :
            print("Composition of Sample " + sOld + " < 1.0: One Old Component Missing?")
            iERR = -1
            return iERR

#== Work Array for Coats-A Calculation: Whitson Eqn.(5.93) ==========

    aSqr = NP.zeros(nOld)

    for iOld in range(nOld) :
        
        MA = clsEOS0.gPP("MA",iOld)
        Tc = clsEOS0.gPP("TC",iOld)
        Pc = clsEOS0.gPP("PC",iOld)
        AF = clsEOS0.gPP("AF",iOld)

        tR = tStd/Tc
        aS = CE.sqrtAlpha(tR,AF,clsEOSG)

        aSqr[iOld] = sqrt(MA/Pc)*Tc*aS  #-- sqrt(OmgA*Tc^2*alfa/Pc]

        #print("iOld,aSqr {:2d} {:10.3e}".format(iOld,aSqr[iOld]))

#== Components ======================================================

    for iNew in range(nNew) :

        clsGRP = dicGRP[iNew]

        newNam = clsGRP.newNam

        clsEOSG.sPP("CN",iNew,newNam)

        oldNum = clsGRP.oldNum
        nGrp   = len(oldNum)

        zNew = dicSAMG[sNumW].gZI(iNew)

        MwNew = 0.0 ; TcNew = 0.0 ; PcNew = 0.0 ; VcNew = 0.0 ; ZcNew = 0.0
        AFNew = 0.0 ; TbNew = 0.0 ; SGNew = 0.0 ; PANew = 0.0 ; SSNew = 0.0
        CANew = 0.0 ; CBNew = 0.0 ; CCNew = 0.0 ; CDNew = 0.0

        aIGrp = 0.0
        bIGrp = 0.0

        for iGrp in range(nGrp) :

            iOld = oldNum[iGrp]
            zOld = dicSAM0[sNumW].gZI(iOld)

            MwOld = clsEOS0.gPP("MW",iOld) ; TcOld = clsEOS0.gPP("TC",iOld)
            PcOld = clsEOS0.gPP("PC",iOld) ; VcOld = clsEOS0.gPP("VC",iOld)
            ZcOld = clsEOS0.gPP("ZC",iOld) ; AFOld = clsEOS0.gPP("AF",iOld)
            TbOld = clsEOS0.gPP("TB",iOld) ; SGOld = clsEOS0.gPP("SG",iOld)
            PAOld = clsEOS0.gPP("PA",iOld) ; SSOld = clsEOS0.gPP("SS",iOld)
            CAOld = clsEOS0.gPP("CA",iOld) ; CBOld = clsEOS0.gPP("CB",iOld)
            CCOld = clsEOS0.gPP("CC",iOld) ; CDOld = clsEOS0.gPP("CD",iOld)

            MwNew = MwNew + zOld*MwOld ; TcNew = TcNew + zOld*TcOld
            PcNew = PcNew + zOld*PcOld ; VcNew = VcNew + zOld*VcOld
            ZcNew = ZcNew + zOld*ZcOld ; AFNew = AFNew + zOld*AFOld
            TbNew = TbNew + zOld*TbOld
            PANew = PANew + zOld*PAOld ; SSNew = SSNew + zOld*SSOld
            CANew = CANew + zOld*CAOld ; CBNew = CBNew + zOld*CBOld
            CCNew = CCNew + zOld*CCOld ; CDNew = CDNew + zOld*CDOld

            SGNew = SGNew + zOld*MwOld/SGOld  #-- Specific Gravity special!

            tRed  = tStd/TcOld
            sqAlf = CE.sqrtAlpha(tRed,AFOld,clsEOSG)
            MAOld = clsEOS0.gPP("MA",iOld)
            aIOld = sqrt(MAOld/PcOld)*TcOld*sqAlf

            MBOld = clsEOS0.gPP("MB",iOld)
            bIOld = MBOld*TcOld/PcOld

            aIGrp = aIGrp + zOld*aIOld
            bIGrp = bIGrp + zOld*bIOld

        SGNew = MwNew/SGNew                   #-- Done as inverse

        MwNew = MwNew/zNew ; TcNew = TcNew/zNew ; PcNew = PcNew/zNew
        VcNew = VcNew/zNew ; ZcNew = ZcNew/zNew ; AFNew = AFNew/zNew
        TbNew = TbNew/zNew ; PANew = PANew/zNew
        SSNew = SSNew/zNew ; CANew = CANew/zNew ; CBNew = CBNew/zNew
        CCNew = CCNew/zNew ; CDNew = CDNew/zNew

        clsEOSG.sPP("MW",iNew,MwNew) ; clsEOSG.sPP("TC",iNew,TcNew)
        clsEOSG.sPP("PC",iNew,PcNew) ; clsEOSG.sPP("VC",iNew,VcNew)
        clsEOSG.sPP("ZC",iNew,ZcNew) ; clsEOSG.sPP("AF",iNew,AFNew)
        clsEOSG.sPP("TB",iNew,TbNew) ; clsEOSG.sPP("SG",iNew,SGNew)
        clsEOSG.sPP("PA",iNew,PANew) ; clsEOSG.sPP("SS",iNew,SSNew)
        clsEOSG.sPP("CA",iNew,CANew) ; clsEOSG.sPP("CB",iNew,CBNew)
        clsEOSG.sPP("CC",iNew,CCNew) ; clsEOSG.sPP("CD",iNew,CDNew)

#-- Compute Coats OmegaB --------------------------------------------

        tRed  = tStd/TcNew
        sqAlf = CE.sqrtAlpha(tRed,AFNew,clsEOSG)
        TcSqA = TcNew*sqAlf
        
        aI    = TcSqA*TcSqA/PcNew
        bI    =       TcNew/PcNew

        aIGrp = aIGrp/zNew
        bIGrp = bIGrp/zNew
        
        OmgAI = aIGrp*aIGrp/aI
        OmgBI =       bIGrp/bI  #-- Whitson Eqn.(5.93)

        clsEOSG.sPP("MA",iNew,OmgAI)
        clsEOSG.sPP("MB",iNew,OmgBI)

#== Binaries ==========================================================

    for iNew in range(nNew) :

        iOldN = dicGRP[iNew].oldNum
        nIOld = len(iOldN)

        zINew = dicSAMG[sNumW].gZI(iNew)

        for jNew in range(iNew+1,nNew) :

            jOldN = dicGRP[jNew].oldNum
            nJOld = len(jOldN)

            zJNew = dicSAMG[sNumW].gZI(jNew)

            KIJ = 0.0
            ZIJ = 0.0

            for iGrp in range(nIOld) :

                iOld   = iOldN[iGrp]
                zIOld  = dicSAM0[sNumW].gZI(iOld)
                
                for jGrp in range(nJOld) :

                    jOld   = jOldN[jGrp]
                    zJOld  = dicSAM0[sNumW].gZI(jOld)
                    KIJOld = clsEOS0.gIJ(iOld,jOld)

                    ZIJ = ZIJ + zIOld*zJOld
                    KIJ = KIJ + zIOld*zJOld*KIJOld

            KIJ = KIJ/ZIJ

            clsEOSG.sIJ(iNew,jNew,KIJ)
            clsEOSG.sIJ(jNew,iNew,KIJ)

#== End of Routine ====================================================

    return iERR

#======================================================================
#  Check a Component Name is Valid
#======================================================================

def checkComp4Grp(sName,clsEOS) :

    iComp = -1
    nComp = clsEOS.NC

    for iC in range(nComp) :
        sThis = clsEOS.gPP("CN",iC)
        if sThis == sName :
            iComp = iC
            break

    return iComp

#======================================================================
#  Check a Sample Name is Valid
#======================================================================

def checkSamp4Grp(sName,dicSAM) :

    iSamp = -1
    nSamp = len(dicSAM)

    for iS in range(nSamp) :
        sThis = dicSAM[iS].sName
        if sThis == sName :
            iSamp = iS
            break

    return iSamp

#======================================================================
#  End of Module
#======================================================================
