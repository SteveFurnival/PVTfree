
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Data for Compositional Data Output
#
#=======================================================================

#!/usr/bin/python3

import datetime

import numpy     as NP

import allData   as AD
import calcGrad  as CG
import writeComp as WC

#========================================================================
#  Read the Data
#========================================================================

def readComp(clsEOS,dicSAM,clsIO,clsUNI) :

    iERR  = 0
    iLine = 0

    Tres  = None
    sSam  = None

    pRef = None
    Salt = None

    pRef = None
    dRef = None

    qDuni = None
    dMin  = None
    dMax  = None
    dInc  = None

    fInp = clsIO.fInp

#== Create the clsCMP class to hold relevant data ======================    

    clsCMP = AD.classComp()

#----------------------------------------------------------------------
#  Loop over lines until ENDC[OMP] found
#----------------------------------------------------------------------

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of COMP Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0 :                        #-- Blank line!
            pass
        elif tokS[0][:2]         == "--"   :    #-- Comment
            pass
        elif tokS[0][:4].upper() == "ENDC" :    #-- ENDCOMP read
            break

        elif tokS[0][:2].upper() == "UN"   :    #-- UNITS keyword [FIELD/METRIC]

            uTyp = tokS[1].upper()

            if uTyp[:1] == "M" : clsCMP.setOutUnit("MET")
            else               : clsCMP.setOutUnit("FLD")

        elif tokS[0][:2].upper() == "SI"   :    #-- Simulator Type [ECL/MOR]

            uTyp = tokS[1].upper()

            if   uTyp[:1] == "T" or \
                 uTyp[:1] == "M" : clsCMP.setSimType("TEM")
            elif uTyp[:1] == "E" : clsCMP.setSimType("ECL")
            elif uTyp[:1] == "C" or \
                 uTyp[:1] == "G" : clsCMP.setSimType("CMG")
            elif uTyp[:1] == "V" or \
                 uTyp[:1] == "N" : clsCMP.setSimType("VIP")

        elif tokS[0][:1].upper() == "W" :       #-- Properties to define Water

            iTok = 1
            while iTok < nTok:
                if   tokS[iTok][:1].upper() == "P" :    #-- Reference Pressure
                    iTok += 1
                    pRef  = float(tokS[iTok])
                    iTok += 1
                    wUni  =       tokS[iTok]
                elif tokS[iTok][:1].upper() == "S" :    #-- Salinity
                    iTok += 1
                    Salt  = float(tokS[iTok])
                    iTok += 1
                    sUni  =       tokS[iTok]
                iTok += 1

        elif tokS[0][:6].upper() == "COMREF" :

            iTok = 1
            while iTok < nTok :
                if   tokS[iTok][:1].upper() == "P" :    #-- PREF
                    iTok += 1
                    pRef  = float(tokS[iTok])
                    iTok += 1
                    pDuni =       tokS[iTok]
                elif tokS[iTok][:1].upper() == "D" :    #-- DREF
                    iTok += 1
                    dRef  = float(tokS[iTok])
                    iTok += 1
                    dDuni  =       tokS[iTok]
                iTok += 1

        elif tokS[0][:6].upper() == "COMDEP" :

            dUni  = tokS[1].lower()             #-- 1st Argument must be units
            qDuni = dUni in clsUNI.sTy["leng"]

            if qDuni :
                clsCMP.setDepUni(dUni.upper())
            else :
                print("First argument after COMDEP must be valid Length Units [ft, m] - Error")
                iERR = -1
                break

            iTok = 2
            while iTok < nTok :
                if   tokS[iTok][:4].upper() == "DMAX" :    #-- DMAX
                    iTok += 1
                    dMax  = float(tokS[iTok])
                elif tokS[iTok][:4].upper() == "DMIN" :    #-- DMIN
                    iTok += 1
                    dMin  = float(tokS[iTok])
                elif tokS[iTok][:4].upper() == "DINC" :    #-- DINC
                    iTok += 1
                    dInc  = float(tokS[iTok])
                iTok += 1

#-- Line must contains non-Table information ------------------------                

        else :  

            iTok = 0
            while iTok < nTok :
                if   tokS[iTok].upper() == "SAMP" :
                    iTok += 1
                    sSam = tokS[iTok]
                elif tokS[iTok].upper() == "TRES" :
                    iTok += 1
                    Tres = float(tokS[iTok])
                    iTok += 1
                    Tuni =       tokS[iTok]
                iTok += 1
        
#----------------------------------------------------------------------
#  Check Data Has Been Defined
#----------------------------------------------------------------------

    if Tres == None :
        print("No Reservoir Temperature Has Been Read via TRES argument - Error")
        iERR = -1
        return iERR
    else            :
        clsCMP.setTres(clsUNI.X2I(Tres,Tuni),Tuni)

#-- Valid Reservoir Sample to work on? ------------------------------        

    if sSam == None :
        print("No Sample Name Has Been Read via SAMP argument - Error")
        iERR = -1
        return iERR
    else :
        nSam = len(dicSAM)
        qFnd = False
        for iSam in range(nSam) :
            sArr = dicSAM[iSam].sNam
            if sSam.upper() == sArr.upper() :
                clsCMP.setSamp(iSam,sSam)
                qFnd = True
                break
        if not qFnd :
            print("Did Not Find Sample ",sSam," In Previously Defined Samples - Error")
            iERR = -1
            return iERR

#-- Have we read the WATER card with pRef and Salt? -----------------

    if pRef == None and Salt == None :

        clsCMP.setBrine(False)
        clsCMP.setPrefWat(14.6959)   #-- Standard Pressure [psia]
        clsCMP.setPrefUni("psia")
        clsCMP.setSaltWat(0.0)       #-- Pure Water
        clsCMP.setSaltUni("mfrac")

    else :

        clsCMP.setBrine(True)

        if Salt != None :
            clsCMP.setSaltUni(sUni)
            clsCMP.setSaltWat(clsUNI.X2I(Salt,sUni))
        else :
            clsCMP.setSaltWat(0.0)       #-- Pure Water
            clsCMP.setSaltUni("mfrac")
        
        if pRef != None :
            clsCMP.setPrefUni(wUni)
            clsCMP.setPrefWat(clsUNI.X2I(pRef,wUni))
        else :
            clsCMP.setPrefWat(14.6959)   #-- Standard Pressure [psia]
            clsCMP.setPrefUni("psia")
            
#== Composition versus Depth Calculation wanted =========================

    qRef = False

    if pRef != None and dRef != None :

        clsCMP.setDepRef(clsUNI.X2I(dRef,dDuni),dDuni)
        clsCMP.setPrsRef(clsUNI.X2I(pRef,pDuni),pDuni)

        qRef = True

    if qDuni != None :

        clsCMP.setDepths(dMin,dMax,dInc,qDuni)

        if qRef : clsCMP.qDep = True

#-- OK to do calculation ----------------------------------------------

    if clsCMP.qDep :
        ZMFvD = compGrad(clsEOS,dicSAM,clsCMP,clsUNI,clsIO)
    else :
        ZMFvD = None

#----------------------------------------------------------------------
#  Data Read OK: Write Output
#----------------------------------------------------------------------

    if clsCMP.tSim == "ECL" :
        WC.writeE300(ZMFvD,clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "TEM" :
        WC.writeMORE(ZMFvD,clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "VIP" :
        WC.writeVIP( ZMFvD,clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "CMG" :
        WC.writeGEM( ZMFvD,clsEOS,dicSAM,clsCMP,clsIO,clsUNI)

#========================================================================
#  End of Routine
#========================================================================

    return iERR

#=======================================================================
#  Driver for Compositional Table Calculations
#=======================================================================

def compGrad(clsEOS,dicSAM,clsCMP,clsUNI,clsIO) :

    tRef = clsCMP.Tres

    iSam = clsCMP.sNum
    nCom = clsEOS.nComp

    clsSAM = dicSAM[iSam]

    zRef = NP.zeros(nCom)
    for iC in range(nCom) : zRef[iC] = clsSAM.gZI(iC)

    dRef = clsCMP.dRef
    dMin = clsCMP.dMin
    dMax = clsCMP.dMax
    dInc = clsCMP.dInc
    pRef = clsCMP.pRef

    ZMFvD = []   #-- Will hold Composition vs. Depth

#== Get Reference Depth (log) Fugacity Coeffs =========================

    iNeu = 0

    fRef,dumP = \
        CG.setupCalcGRDCoefsP(iNeu,pRef,tRef,zRef,clsEOS)  #-- Dont need Ref dlnPhi/dP

#== Minimum (Top) Depth ===============================================

    iTop = 0
    zTyp = "Top"
    dTop = dMin - dRef

    pTop,zTop,qTop,sTop = \
        CG.calcTopBottom(zTyp,dTop,pRef,tRef,zRef,fRef,clsEOS,clsIO)

    if not qTop : iTop = - 1    #-- Is DEW point

    ZMFvD.append([dMin,sTop,zTop])

#== Maximum (Bottom) Depth ============================================    
        
    iBot = 0
    zTyp = "Bot"
    dBot = dMax - dRef

    pBot,zBot,qBot,sBot = \
        CG.calcTopBottom(zTyp,dBot,pRef,tRef,zRef,fRef,clsEOS,clsIO)

    if qBot : iBot =   1    #-- Is BUB point

#== Gas-Oil-Contact (GOC) Present? ====================================

    if iTop*iBot < 0 :
       qGOC = True
       dGOC,pGOC,sGOC,rGOC,zC1,zC7,zDGOC,zBGOC = \
            CG.calcGOC(dMin,dMax,tRef,dRef,pRef,zRef,fRef,clsEOS,clsIO)
    else :
        qGOC = False
        dGOC = 0.0

    clsCMP.setDepGOC(dGOC)

#----------------------------------------------------------------------
#  Use Depth-Increment to step from Top to beyond Bottom
#----------------------------------------------------------------------

    dNxt = dMin + dInc
    zTyp = "Stp"
    qIns = True

    while dNxt < dMax :

        dDep = dNxt - dRef

        pDep,zDep,qDep,sDep = \
            CG.calcTopBottom(zTyp,dDep,pRef,tRef,zRef,fRef,clsEOS,clsIO)

#== Insert GOC? =======================================================

        if qIns and qGOC and dNxt > dGOC :
            ZMFvD.append([dGOC,sGOC,zDGOC])
            qIns = False

#-- Add this depths results -----------------------------------------

        ZMFvD.append([dNxt,sDep,zDep])

#-- Increment the Depth ---------------------------------------------

        dNxt = dNxt + dInc

#== Add the Bottom Depth Values =======================================        
        
    ZMFvD.append([dMax,sBot,zBot])

#== Return ============================================================

    return ZMFvD

#========================================================================
#  End of Module
#========================================================================

