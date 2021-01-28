
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculate Component Grouping Module
#
#=======================================================================

#!/usr/bin/python3

import numpy    as NP

from   math     import sqrt

import allData   as AD
import calcEOS   as CE
import utilities as UT

#========================================================================
#  Calculate Grouping
#========================================================================

def calcGroup(sNumW,dicGRP,clsEOS0,dicSAM0) :

    iERR = 0

    nOld = clsEOS0.nComp
    nNew = len(dicGRP)
    nSam = len(dicSAM0)
    
    tStd = UT.tStand      #-- Standard Temperature [519.67 degR = 60 degF]

    #print("calcGroup: nOld,nNew,nSam ",nOld,nNew,nSam)

#== Create 'Group' Copies of EOS/Samples ==============================        

    EOS = clsEOS0.EOS

    clsEOSG = AD.classEOS(EOS)
    dicSAMG = {}

    clsEOSG.setEOSparms(EOS)

#----------------------------------------------------------------------
#  Create New EOS-Class and Dictionary of Samples
#----------------------------------------------------------------------

    clsEOSG.nComp = nNew
    clsEOSG.nSamp = nSam

    clsEOSG.setNComp(nNew)   #-- Dimensions the Arrays

#== Samples ===========================================================

    for iSam in range(nSam) :

        clsSAM0 = dicSAM0[iSam]
        
        sOld    = clsSAM0.sNam
        clsSAMG = AD.classSample(sOld)
        
        dicSAMG[iSam] = clsSAMG
        dicSAMG[iSam].setIntComp(nNew,nNew)    #-- All "User" Components

        clsSAMG.iPlsMW = clsSAM0.iPlsMW
        clsSAMG.iPlsSG = clsSAM0.iPlsSG
        clsSAMG.iPlsAL = clsSAM0.iPlsAL

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

        #print("iSam,sTot ",iSam,sTot)

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

        clsEOSG.sNM(iNew,newNam)

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

    return iERR,clsEOSG,dicSAMG

#======================================================================
#  End of Module
#======================================================================
