#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Perform Expeiments Module
#  
#========================================================================

#!/usr/bin/python3

import numpy        as NP
import numpy.linalg as LA

import calcEOS   as CE
import constants as CO
import calcFlash as CF
import calcSat   as CS
import readSamp  as RS
import readExps  as RX
import utilities as UT
import writeOut  as WO

from math  import log, exp, sqrt

#========================================================================
#  Calculate Experiments
#========================================================================

def calcExps(sExt,clsEOS,dicSAM,dicEXP,clsIO,clsUNI) :

    nExp = len(dicEXP)

    sLabl = "calcExps: Running Experiments  " + sExt
    print(sLabl)

#== Setup the (P,T)-Independent EoS Coefficients ======================

    clsEOS.eosCoefsNoPT()

#== Loop over experiments =============================================

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            xTyp = clsEXP.xName

            sOut = "Experiment: iExp,Type {:2d} ".format(iExp+1)
            sOut = sOut + xTyp

            print(sOut)

            if   xTyp == "CCE" : calcCCE(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "CVD" : calcCVD(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "DLE" : calcDLE(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SEP" : calcSEP(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "FLS" : calcFLS(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SAT" : calcSAT(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SWL" : calcSWL(iExp,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "GRD" : calcGRD(iExp,clsEOS,dicSAM,clsEXP,clsIO)

#== Write the experiments to the Print Output File (fOut) =============

    WO.outputExps(clsIO,dicEXP,dicSAM,clsUNI)
    
    sLabl = "calcExps: Experiments Completed" + sExt
    print(sLabl)

    return

#========================================================================
#  Constant Composition Expansion (CCE) Experiment
#========================================================================

def calcCCE(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sName

    typDEP = RX.expDEPvars.get("CCE")
    typOBS = RX.expOBSvars.get("CCE")

    iPR = typDEP.index("PRES")

    iRV = typOBS.index("RELV")
    iSL = typOBS.index("SLIQ")
    iZF = typOBS.index("ZFAC")
    iDO = typOBS.index("DENO")
    iUG = typOBS.index("VISG")
    iUO = typOBS.index("VISO")

#-- Load Feed Composition -------------------------------------------

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Saturation Pressure ---------------------------------------------

    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsSAM,clsIO)

    clsEXP.setPsatCal(pSat)

    iNeu =  0
    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,Z,clsEOS)

    if qBub : vEst = 0.0
    else :    vEst = 1.0

#======================================================================
#  Loop over user defined pressures
#======================================================================

    for iPrs in range(nPrs) :

        pRes = clsEXP.dDep[iPrs][iPR]

#== 1-Phase or 2-Phase? ===============================================

        if pRes > pSat :
            
            M1P,V1P,D1P,Z1P,U1P,dmS,dmS = CE.calcProps(iNeu,pRes,tRes,Z,clsEOS)
            #print("calcCCE[Z]: iP,pR,Mw,Vm,Ro,ZF,vS {:2d} {:8.2f} {:7.3f} {:8.4f} {:8.3f} {:8.4f} {:8.4f}".format(iPrs,pRes,M1P,V1P,D1P,Z1P,V1P))

            relV = V1P/Vsat
            
            if qBub : sLiq = 1.0
            else :    sLiq = 0.0

            Zgas = Z1P
            Doil = D1P
            Ugas = U1P
            Uoil = U1P
            
        else :  #-- Two-Phase => Flash the Fluid

            V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)
            vEst    = V

            Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
            Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

            Voil = (1.0-V)*Vliq
            Vgas =      V *Vvap

            Vtot = Voil + Vgas

            relV = Vtot/Vsat

            if clsEXP.sLCCE == "TOT" : sLiq = Voil/Vtot
            else :                     sLiq = Voil/Vsat

#== Store the data for this stage =====================================

        clsEXP.dCal[iPrs][iRV] = relV
        clsEXP.dCal[iPrs][iSL] = sLiq
        clsEXP.dCal[iPrs][iZF] = Zgas
        clsEXP.dCal[iPrs][iDO] = Doil
        clsEXP.dCal[iPrs][iUG] = Ugas
        clsEXP.dCal[iPrs][iUO] = Uoil

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Constant Volume Depletion (CVD) Experiment
#========================================================================

def calcCVD(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sName

    #print("calcCVD: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("CVD")
    typOBS = RX.expOBSvars.get("CVD")

    iPR = typDEP.index("PRES")

    iMR = typOBS.index("MREM")
    iSL = typOBS.index("SLIQ")
    iZF = typOBS.index("ZFAC")

#-- Load Feed Composition -------------------------------------------

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Saturation Pressure ---------------------------------------------
    
    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsSAM,clsIO)

    clsEXP.setPsatCal(pSat)

    iNeu =  0
    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,Z,clsEOS)

    if qBub : vEst = 0.0
    else :    vEst = 1.0

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zMol = 1.0
    zTot = 0.0

    for iPrs in range(nPrs) :

        pRes = clsEXP.dDep[iPrs][iPR]

#== 1-Phase or 2-Phase? ===============================================

        if pRes > pSat :
            
            M1P,V1P,D1P,Z1P,U1P,dmS,dmS = CE.calcProps(iNeu,pRes,tRes,Z,clsEOS)

            if qBub : sLiq = 1.0
            else :    sLiq = 0.0

            zTot = 0.0
            Zgas = Z1P
            
        else :

            V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)
            vEst    = V

            Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
            Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

            Voil = zMol*(1.0-V)*Vliq
            Vgas = zMol*     V *Vvap

            sLiq = Voil/Vsat

            Vtot = Voil + Vgas
            Vrem = Vtot - Vsat

            zRem = pRes*Vrem/(Zgas*CO.gasCon*tRes)

            for iC in range(nCom) :
                Z[iC] = (zMol*Z[iC] - zRem*Y[iC])/(zMol - zRem)

            zMol = zMol - zRem
            zTot = zTot + zRem

#-- Load Data -------------------------------------------------------

        clsEXP.dCal[iPrs][iMR] = zTot
        clsEXP.dCal[iPrs][iSL] = sLiq
        clsEXP.dCal[iPrs][iZF] = Zgas

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Differential Liberation (DLE) Experiment
#========================================================================

def calcDLE(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sName

    #print("calcCVD: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("DLE")
    typOBS = RX.expOBSvars.get("DLE")

    iPR = typDEP.index("PRES")

    iBO = typOBS.index("BO")
    iRS = typOBS.index("GOR")
    iDO = typOBS.index("DENO")
    iBT = typOBS.index("BT")
    iBG = typOBS.index("BG")
    iZF = typOBS.index("ZFAC")
    iGG = typOBS.index("GGRV")
    iUO = typOBS.index("VISO")

#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Saturation Pressure ---------------------------------------------
    
    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsSAM,clsIO)

    clsEXP.setPsatCal(pSat)

    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

    if qBub : vEst = 0.0
    else :    vEst = 1.0

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zMol = 1.0
    rTot = 0.0
    psTs = CO.pStand/CO.tStand

    for iPrs in range(nPrs) :

        pRes = clsEXP.dDep[iPrs][iPR]

#== 1-Phase or 2-Phase? ===============================================

        if pRes > pSat :
            
            M1P,V1P,D1P,Z1P,U1P,dmS,dmS = CE.calcProps(iLiq,pRes,tRes,Z,clsEOS)

            if qBub : sLiq = 1.0
            else :    sLiq = 0.0

            Voil = V1P
            Vrem = 0.0
            Doil = D1P
            Vtot = V1P
            Bgas = 0.0
            Zgas = 0.0
            gGrv = 0.0
            Uoil = U1P
            
        else :

            V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)
            vEst    = V

            Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
            Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

            Voil = zMol*(1.0-V)*Vliq
            Vgas = zMol*     V *Vvap

            zRem = pRes*Vgas/(Zgas*CO.gasCon*tRes)
            zMol = zMol - zRem

            Vrem = zRem*CO.volMol

            rTot = rTot + Vrem

            Z = NP.copy(X)

            Vtot = Voil + Vgas
            Bgas = psTs*Zgas*tRes/pRes
            gGrv = Mgas/CO.molAir

#-- Temporary Storage Ahead of Stock Tank Corrections ---------------            

        clsEXP.dCal[iPrs][iBO] = Voil
        clsEXP.dCal[iPrs][iRS] = Vrem
        clsEXP.dCal[iPrs][iDO] = Doil
        clsEXP.dCal[iPrs][iBT] = Vtot
        clsEXP.dCal[iPrs][iBG] = Bgas
        clsEXP.dCal[iPrs][iZF] = Zgas
        clsEXP.dCal[iPrs][iGG] = gGrv
        clsEXP.dCal[iPrs][iUO] = Uoil

#== Stock Tank Volume =================================================

    Mst,Vst,Dst,Zst,Ust,dmS,dmS = CE.calcProps(iLiq,CO.pStand,CO.tStand,Z,clsEOS)

    Vsto = zMol*Vst

#== "Normalise" via Stock Tank Volume =================================

    for iPrs in range(nPrs) :

        rTot = rTot - clsEXP.dCal[iPrs][iRS]

        clsEXP.dCal[iPrs][iBO] = clsEXP.dCal[iPrs][iBO]/Vsto
        clsEXP.dCal[iPrs][iRS] =                   rTot/Vsto

        clsEXP.dCal[iPrs][iBT] = clsEXP.dCal[iPrs][iBO] + \
                                clsEXP.dCal[iPrs][iBG]*(clsEXP.dCal[0][iRS] - clsEXP.dCal[iPrs][iRS])

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Separator Test (SEP) Experiment
#========================================================================

def calcSEP(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    if clsIO.Deb["SEP"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sName

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("SEP")
    typOBS = RX.expOBSvars.get("SEP")

    iPR = typDEP.index("PRES")
    iTR = typDEP.index("TEMP")

    iBO = typOBS.index("BO")
    iRS = typOBS.index("GOR")
    iDO = typOBS.index("DENO")
    iGG = typOBS.index("GGRV")

    Lsep = clsEXP.Lsep
    Vsep = clsEXP.Vsep

    #print("calcSEP: Lsep ",Lsep)
    #print("calcSEP: Vsep ",Vsep)

#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

    molZ = NP.zeros( nPrs)
    comZ = NP.zeros((nCom,nPrs))

    xSTO = NP.zeros(nCom) ; ySTG = NP.zeros(nCom)
    molO = NP.zeros(nPrs) ; molG = NP.zeros(nPrs)

    mSTO = 0.0 ; mSTG = 0.0

#-- First Stage is (pSat,tRes) and all feed goes to Stage-2 ---------    

    molZ[0] = 1.0 ; molZ[1] = 1.0
    for iC in range(nCom) :
        comZ[iC][0] = Z[iC] ; comZ[iC][1] = Z[iC]

#-- 1st-Stage is Assumed to (Psat,tRes) -----------------------------

    pObs = clsEXP.dDep[0][iPR]
    tRes = clsEXP.dDep[0][iTR]

    #print("calcSEP: pObs,tRes ",pObs,tRes)

#-- Saturation Pressure ---------------------------------------------
    
    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsSAM,clsIO)

    #print("CalcSEP: Psat {:10.3f}".format(pSat))

    clsEXP.setPsatCal(pSat)

    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

    vEst = 0.5

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zGas = 0.0
    zOil = 0.0
    zTot = 1.0

    for iPrs in range(1,nPrs) :

        pRes = clsEXP.dDep[iPrs][iPR]
        tRes = clsEXP.dDep[iPrs][iTR]

        #print("iPrs,pRes,tRes {:2d} {:10.3f} {:8.3f}".format(iPrs,pRes,tRes))

#-- 2-Phase Flash ---------------------------------------------------

        for iC in range(nCom) : Z[iC] = comZ[iC][iPrs]

        #zSum = NP.sum(Z)
        #Z    = NP.divide(Z,zSum)

        Z = UT.Norm(Z)

        #print("iSep,Z ",iPrs,Z)

        V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)

        if   V < 0.0 : V = 0.0  #-- calcFlash sets X=Z & Y=Z if 1-Phase
        elif V > 1.0 : V = 1.0
        
        vEst    = V

        zGas = molZ[iPrs]*     V
        zOil = molZ[iPrs]*(1.0-V)

        molG[iPrs] = zGas
        molO[iPrs] = zOil

        #print("iSep,molZ,molG,molO {:2d} {:10.5f} {:10.5f} {:10.5f}".format(iPrs,molZ[iPrs],molG[iPrs],molO[iPrs]))

#-- Properties of the Liquid and Vapour Output of this stage --------        

        Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
        Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

#-- Where is Liquid Output Going Next? ------------------------------

        kSepL = Lsep[iPrs]

        if kSepL == 0 :
            if mSTO + zOil > 0.0 :
                for iC in range(nCom) :
                    xSTO[iC] = (xSTO[iC]*mSTO + zOil*X[iC])/(mSTO + zOil)
                mSTO = mSTO + zOil
            else :
                xSTO = BP.copy(X)
        else :
            for iC in range(nCom) :
                comZ[iC][kSepL] = (comZ[iC][kSepL]*molZ[kSepL] + X[iC]*zOil)/(molZ[kSepL] + zOil)
            molZ[kSepL] = molZ[kSepL] + zOil
        
#-- Where is Vapour Output Going Next? ------------------------------        

        kSepV = Vsep[iPrs]

        if kSepV == 0 :
            if mSTG + zGas > 0.0 :
                for iC in range(nCom) :
                    ySTG[iC] = (ySTG[iC]*mSTG + zGas*Y[iC])/(mSTG + zGas)
                mSTG = mSTG + zGas
            else :
                ySTG = NP.copy(Y)
        else :
            for iC in range(nCom) :
                comZ[iC][kSepV] = (comZ[iC][kSepV]*molZ[kSepV] + Y[iC]*zGas)/(molZ[kSepV] + zGas)
            molZ[kSepV] = molZ[kSepV] + zGas     
        
        Voil = zOil*Vliq
        Vgas = zGas*CO.volMol

#-- Temporary storage for Bo & Rs ahead of divide by Vsto below -----        

        clsEXP.dCal[iPrs][iBO] = Voil
        clsEXP.dCal[iPrs][iRS] = Vgas
        clsEXP.dCal[iPrs][iDO] = Doil
        clsEXP.dCal[iPrs][iGG] = Mgas/CO.molAir

#== Stock Tank Oil Volume =============================================

    nTot = 0.0
    rTot = 0.0
    gTot = 0.0
    Vsto = clsEXP.dCal[nPrs-1][iBO]

    for iPrs in range(1,nPrs) :

        clsEXP.dCal[iPrs][iBO] = clsEXP.dCal[iPrs][iBO]/Vsto  #-- Bo
        clsEXP.dCal[iPrs][iRS] = clsEXP.dCal[iPrs][iRS]/Vsto  #-- GOR

        nTot = nTot + molG[iPrs]
        rTot = rTot + molG[iPrs]*clsEXP.dCal[iPrs][iRS]
        gTot = gTot + molG[iPrs]*clsEXP.dCal[iPrs][iGG]

#== Saturation Pressure Stage =========================================

    if nTot > 0.0 : gTot = gTot/nTot
    else :          gTot = 0.0

    #print("mSTG,mSTO {:10.5f} {:10.5f}".format(mSTG,mSTO))

    clsEXP.dCal[0][iBO] = Vsat/Vsto
    clsEXP.dCal[0][iRS] = mSTG*CO.volMol/Vsto
    clsEXP.dCal[0][iDO] = Dsat
    clsEXP.dCal[0][iGG] = gTot
    
#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  2-Phase Flash (FLS) Experiment
#========================================================================

def calcFLS(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    iLiq =  1
    iVap = -1

    sNam = dicSAM[nSam].sName

    #print("calcFLS: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("FLS")
    typOBS = RX.expOBSvars.get("FLS")

    iPR = typDEP.index("PRES")
    iTR = typDEP.index("TEMP")

    iZF = typOBS.index("ZFAC")
    iDO = typOBS.index("DENO")
    iVF = typOBS.index("VFRC")

#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = dicSAM[nSam].gZI(iC)

#======================================================================
#  Loop over user defined 'stages'
#======================================================================

    for iPrs in range(nPrs) :

        pRes = clsEXP.dDep[iPrs][iPR]
        tRes = clsEXP.dDep[iPrs][iTR]

#-- 2-Phase Flash ---------------------------------------------------        

        vEst    = 0.5
        V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)

#-- Properties of the Liquid and Vapour Output of this stage --------        

        Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
        Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

        clsEXP.dCal[iPrs][iZF] = Zgas
        clsEXP.dCal[iPrs][iDO] = Doil
        clsEXP.dCal[iPrs][iVF] = V

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Saturation Pressure (SAT) Experiment
#========================================================================

def calcSAT(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nTem = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sName

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("SAT")
    typOBS = RX.expOBSvars.get("SAT")

    iTR = typDEP.index("TEMP")

    iPS = typOBS.index("PSAT")
    iZF = typOBS.index("ZFAC")
    iDO = typOBS.index("DENO")

#-- Load Feed Composition -------------------------------------------

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = dicSAM[nSam].gZI(iC)

#======================================================================
#  Loop over user defined 'stages'
#======================================================================

    for iTem in range(nTem) :

        tRes = clsEXP.dDep[iTem][iTR]

        pSatO = clsEXP.dObs[iTem][iPS]

        if pSatO > 0.0 : pEst =  0.95*pSatO
        else           : pEst = -1.0
        
        qBub,pSat,Ksat = CS.calcPsat(pEst,tRes,clsEOS,clsSAM,clsIO)

        iLiq = 0
        MPs,VPs,DPs,ZPs,UPs,dmS,dmS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

        clsEXP.dCal[iTem][iPS] = pSat
        clsEXP.dCal[iTem][iZF] = ZPs
        clsEXP.dCal[iTem][iDO] = DPs

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Swelling Test (SWL) Experiment
#========================================================================

def calcSWL(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nInj = clsEXP.nSinj
    nAdd = clsEXP.nRow

    sNam = dicSAM[nSam].sName

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typDEP = RX.expDEPvars.get("SWL")
    typOBS = RX.expOBSvars.get("SWL")

    iMO = typDEP.index("MOLE")

    iPS = typOBS.index("PSAT")
    iVS = typOBS.index("VSWL")

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Load Feed Composition -------------------------------------------    

    clsSAM = dicSAM[nSam]
    clsINJ = dicSAM[nInj]

    clsWRK = RS.classSample("SWLwrk")
    clsWRK.setNComp(nCom)

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Load Injection Composition --------------------------------------    

    Y = NP.zeros(nCom)
    for iC in range(nCom) : Y[iC] = clsINJ.gZI(iC)

#-- Composition Array to be used ------------------------------------

    W = NP.zeros(nCom)

    mTot = 1.0
    mInj = 0.0
    
#======================================================================
#  Loop over user defined 'stages'
#======================================================================

    for iAdd in range(nAdd) :

        mInj = clsEXP.dDep[iAdd][iMO]

        mTot = mInj + 1.0

        #print("calcSWL: iAdd,mAdd,mInj,mTot ",iAdd,mThs,mInj,mTot)

        wDiv = 1.0/mTot
        for iC in range(nCom) :
            W[iC] = (Z[iC] + mInj*Y[iC])*wDiv
            clsWRK.sZI(iC,W[iC])

        pSatO = clsEXP.dObs[iAdd][iPS]

        if pSatO > 0.0 : pEst =  0.95*pSatO
        else           : pEst = -1.0
        
        qBub,pSat,Ksat = CS.calcPsat(pEst,tRes,clsEOS,clsWRK,clsIO)

        iLiq = 0
        MPs,VPs,DPs,ZPs,UPs,dmS,dmS = CE.calcProps(iLiq,pSat,tRes,W,clsEOS)

        if iAdd == 0 : Vref = VPs

        Vswl = mTot*VPs/Vref

        clsEXP.dCal[iAdd][iPS] = pSat
        clsEXP.dCal[iAdd][iVS] = Vswl

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Composition versus Depth (GRD) Experiment
#========================================================================

def calcGRD(iExp,clsEOS,dicSAM,clsEXP,clsIO) :

    if clsIO.Deb["GRD"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else:
        qDeb = False

#-- Various dimensions ----------------------------------------------    

    nCom = clsEOS.NC
    nSam = clsEXP.nSamp
    nDep = clsEXP.nRow
    nRef = clsEXP.nDref

#-- Pointers into the potential observed data (for calculated data) -

    typDEP = RX.expDEPvars.get("GRD")
    typOBS = RX.expOBSvars.get("GRD")

    iHG = typDEP.index("HEIG")

    iPR = typOBS.index("PRES")
    iPS = typOBS.index("PSAT")
    iDN = typOBS.index("DENS")
    iC1 = typOBS.index("ZC1")
    iCP = typOBS.index("ZC7+")

#-- Sample Name -----------------------------------------------------    

    clsSAM = dicSAM[nSam]
    sNam   = clsSAM.sName

    clsWRK = RS.classSample("GRDwrk")
    clsWRK.setNComp(nCom)

#-- Load Reference Composition --------------------------------------

    zRef = NP.zeros(nCom)
    for iC in range(nCom) : zRef[iC] = clsSAM.gZI(iC)

#-- Reference Depth & Pressure --------------------------------------    

    pRef = clsEXP.Pref
    dRef = clsEXP.Dref

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#== Saturation Pressure at Reference Depth ============================        

    pSatO = clsEXP.PsatO

    if pSatO > 0.0 : pObs =  0.95*pSatO
    else           : pObs = -1.0
        
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsSAM,clsIO)

    iNeu =  0
    iLiq =  1
    iVap = -1

    mSat,vSat,dSat,zSat,uSat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,zRef,clsEOS)

    zC1,zC7 = calcMoleFracGRD(zRef,clsEOS)

    clsEXP.dCal[nRef][iPR] = pRef
    clsEXP.dCal[nRef][iPS] = pSat
    clsEXP.dCal[nRef][iDN] = dSat
    clsEXP.dCal[nRef][iC1] = zC1
    clsEXP.dCal[nRef][iCP] = zC7

#-- Get Fugacity Coeffs and their Pressure-Derivatives --------------

    fRef,refP = setupCalcGRDCoefsP(iNeu,pRef,tRes,zRef,clsEOS)

#-- Above and Below any Possible GOC --------------------------------

    hAbv = None
    hBel = None

#--------------------------------------------------------------------
#  Above Reference Depth
#--------------------------------------------------------------------

    hLst  = dRef
    pSatR = pSat

    for iDep in range(nRef-1,-1,-1) :

        hDep = clsEXP.dDep[iDep][iHG]
        dDep = hDep - dRef

        pDep,zDep = calcStepGRD(dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        mDep,vDep,yDep,gDen,uDep,dumS,dumS = CE.calcProps(iNeu,pDep,tRes,zDep,clsEOS)
        zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

        for iC in range(nCom) : clsWRK.sZI(iC,zDep[iC])

        qSat,pSat,Ksat = CS.calcPsat(pDep,tRes,clsEOS,clsWRK,clsIO)

        if qSat != qBub :  #-- Going Up => Must Have Gone from Bubble -> Dew
            hAbv = hDep
            hBel = hLst
        
        hLst = hDep
        qBub = qSat

        clsEXP.dCal[iDep][iPR] = pDep
        clsEXP.dCal[iDep][iPS] = pSat
        clsEXP.dCal[iDep][iDN] = yDep
        clsEXP.dCal[iDep][iC1] = zC1
        clsEXP.dCal[iDep][iCP] = zC7

        if qDeb :
            if qBub : isBub =  1
            else    : isBub = -1
            sOut = "Above: iDep,hDep,pDep,pSat,Bub,zC1,zC7+,Dens {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:6.4f} {:6.4f} {:7.3f}\n".format(iDep,hDep,pDep,pSat,isBub,zC1,zC7,yDep)
            fDeb.write(sOut)
            #WO.writeArrayDebug(fDeb,Ksat,"Ksat")

#--------------------------------------------------------------------
#  Below Reference Depth
#--------------------------------------------------------------------

    hLst = dRef
    pSat = pSatR

    for iDep in range(nRef+1,nDep) :
        
        hDep = clsEXP.dDep[iDep][iHG]
        dDep = hDep - dRef

        pDep,zDep = calcStepGRD(dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        mDep,vDep,yDep,gDen,uDep,dumS,dumS = CE.calcProps(iNeu,pDep,tRes,zDep,clsEOS)
        zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

        for iC in range(nCom) : clsWRK.sZI(iC,zDep[iC])

        qSat,pSat,Ksat = CS.calcPsat(pDep,tRes,clsEOS,clsWRK,clsIO)

        if qSat != qBub :  #-- Going Down => Must Have Gone from Dew -> Bubble
            hAbv = hLst
            hBel = hDep
        
        hLst = hDep
        qBub = qSat

        clsEXP.dCal[iDep][iPR] = pDep
        clsEXP.dCal[iDep][iPS] = pSat
        clsEXP.dCal[iDep][iDN] = yDep
        clsEXP.dCal[iDep][iC1] = zC1
        clsEXP.dCal[iDep][iCP] = zC7

        if qDeb :
            if qBub : isBub =  1
            else    : isBub = -1
            sOut = "Below: iDep,hDep,pDep,pSat,Bub,zC1,zC7+,Dens {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:6.4f} {:6.4f} {:7.3f}\n".format(iDep,hDep,pDep,pSat,isBub,zC1,zC7,yDep)
            fDeb.write(sOut)
            #WO.writeArrayDebug(fDeb,Ksat,"Ksat")

#--------------------------------------------------------------------
#  Have we bracketed a GOC?
#--------------------------------------------------------------------

    #print("hAbv,hBel ",hAbv,hBel)

    if hAbv != None and hBel != None :

        calcGOCGRD(hAbv,hBel,tRes,dRef,pRef,zRef,fRef,clsEOS,clsIO)

#======================================================================
#  End of Routine
#======================================================================

    return

#========================================================================
#  Take a Step Away from the Reference Depth (GRD Experiment)
#  Whitson & Brule, SPE Phase Behavior Monograph Section 4.6.2
#========================================================================

def calcStepGRD(dDep,pDep,tRes,zRef,fRef,clsEOS,clsIO) :

    mxSS = 50

    nCom = clsEOS.NC
    iLiq = 0

    zDep = NP.zeros(nCom)   #-- Composition at this Depth
    yMol = NP.zeros(nCom)
    res0 = NP.zeros(nCom)
    res1 = NP.zeros(nCom)
    fBar = NP.zeros(nCom)

#-- Gravity Correction: Whitson Eqn.(4.96) --------------------------

    grvC = dDep/(144*CO.gasCon*tRes)

    #print("dDep,grvC {:10.3f} {:10.3e}".format(dDep,grvC))

    for iC in range(nCom) :
        zDep[iC] = zRef[iC]
        fBar[iC] = fRef[iC]*exp(grvC*clsEOS.gPP("MW",iC))  #-- W&B Eqn.(4.96)

#======================================================================
#  Main Iterative Loop
#======================================================================

    iSS   = 0     #-- Successive-Substitution Counter
    qConv = False

    while not qConv :

        iSS += 1

        if iSS > mxSS : break   #-- Not converged!

#-- Fugacity Coefs and their P-Derivatives of New Fluid -------------

        fDep,dFdP = setupCalcGRDCoefsP(iLiq,pDep,tRes,zDep,clsEOS)

#-- Mole Numbers: Whitson Eqn.(4.95) --------------------------------

        sumY = 0.0 ; dQdP = 0.0
        for iC in range(nCom) :
            res1[iC] = res0[iC]
            res0[iC] = fBar[iC]/fDep[iC]                    #-- W&B Eqn.(4.97a)
            yMol[iC] = zDep[iC]*res0[iC]                    #-- W&B Eqn.(4.95)
            sumY     = sumY   + yMol[iC]
            dQdP     = dQdP   + yMol[iC]*dFdP[iC]/fDep[iC]  #-- W&B Eqn.(4.100) - WRONG!!

        sumR = 1.0/sumY
        for iC in range(nCom) : res0[iC] = sumR*res0[iC]    #-- W&B Eqn.(4.97b)

        qFun = 1.0 - sumY                                   #-- W&B Eqn.(4.94)
        dPrs = -qFun/dQdP                                   #-- W&B Eqn.(4.99a)

        pDep = pDep + dPrs                                  #-- W&B Eqn.(4.99b)

#== Converged? ========================================================

        cTst = 0.0
        for iC in range(nCom) :
            cWrk = log(res0[iC]*zDep[iC]/yMol[iC])
            cTst = cTst + cWrk
        cTst = cTst*cTst

        if cTst < 1.0E-13 and abs(dPrs) < 1.0E-06 : break

#== SS or GDEM Update? W&B Eqn.(4.83) =================================

        if iSS % CO.mGDEM1 > 0 : eigV = 1.0
        else :
            
            res0 = NP.log(res0)         #-- Take logs of Residuals
            res1 = NP.log(res1)
                
            eigV = UT.GDEM1(res0,res1,clsIO)  #-- GDEM works on log-Residuals

            res0 = NP.exp(res0)         #-- Restore Residuals
            res1 = NP.exp(res1)

        #print("calcStepGRD: iSS,qFun,dQdP,eigV,dPrs,cTst {:2d} {:10.3e} {:10.3e} {:8.3f} {:10.3e} {:10.3e}".format(iSS,qFun,dQdP,eigV,dPrs,cTst))

        for iC in range(nCom) : yMol[iC] = yMol[iC]*pow(res0[iC],eigV)

#-- Mole Numbers to Composition -------------------------------------

        sumY = 0.0
        for iC in range(nCom) : sumY = sumY + yMol[iC]
        
        sumR = 1.0/sumY
        for iC in range(nCom) : zDep[iC] = sumR*yMol[iC]

#== Return the updated pressure and composition =======================

    #print("calcStepGRD: iSS,pDep {:2d} {:10.3f}".format(iSS,pDep))

    return pDep,zDep

#========================================================================
#  Calculate a Gas-Oil-Contact Depth, Given Above/Below Bounds
#========================================================================

def calcGOCGRD(hAbv,hBel,tRes,hRef,pRef,zRef,fRef,clsEOS,clsIO) :

    if clsIO.Deb["GRD"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.NC
    iNeu = 0
    pObs = -1.0

    zDep = NP.zeros(nCom)   #-- Composition at this Depth

    clsWRK = RS.classSample("GRDwrk")
    clsWRK.setNComp(nCom)

#--------------------------------------------------------------------
#  Simple Interval Halving
#--------------------------------------------------------------------

    iMid = 1

    while hBel - hAbv > 1.0e-01 :

        hMid = 0.5*(hAbv + hBel)
        dDep =      hMid - hRef

        pDep,zDep = calcStepGRD(dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        mDep,vDep,yDep,gDen,uDep,dumS,dumS = CE.calcProps(iNeu,pDep,tRes,zDep,clsEOS)
        zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

        for iC in range(nCom) : clsWRK.sZI(iC,zDep[iC])

        qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,clsEOS,clsWRK,clsIO)

        if qDeb :
            if qBub : iBub =  1
            else    : iBub = -1
            sOut = "GOC search: iMid,hAbv,hMid,hBel,iBub,pSat {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:10.3f}\n" \
                   .format(iMid,hAbv,hMid,hBel,iBub,pSat)
            fDeb.write(sOut)
            
        if qBub : hBel = hMid
        else    : hAbv = hMid

        iMid += 1

    return

#========================================================================
#  Setup and Calculate GRD Fug & dFug/dP Coefficients (with Vol Shift)
#========================================================================

def setupCalcGRDCoefsP(iLiq,pRes,tRes,XY,clsEOS) :

    nCom = clsEOS.NC
    iPhs = 0
    pByT = pRes/tRes

    qP = True   #-- Onlt want Pres-derivatives, not Temp or Comp
    qT = False
    qX = False

#-- Calculate the log(Fugacity) Coefficients ------------------------

    FiC,FiP,dumV,dumX = CE.calcPhaseFugPTX(iPhs,qP,qT,qX,pRes,tRes,XY,clsEOS)

#-- EoS B[i] Coefficients -------------------------------------------    

    BiC = pByT*clsEOS.bI

#== Augment with the Volume Shift Parameters ==========================

    for iC in range(nCom) :

        cCor = BiC[iC]*clsEOS.gPP("SS",iC)              #-- bi*p/(RT)*si = Bi*si

        FiC[iC] = pRes*XY[iC]*exp(FiC[iC])              #-- log(phi) -> Fugacity
        FiP[iC] =     FiC[iC]*(   FiP[iC] + 1.0/pRes)   #-- dFug/dP

        FiC[iC] = FiC[iC]*exp(-cCor)                    #-- Volume Shift Corr
        FiP[iC] = FiP[iC]*exp(-cCor)*(1.0 - cCor/pRes)  #-- Volume Shift Corr

#== Return values =====================================================        

    return FiC,FiP

#========================================================================
#  Calculate Z[C1] and Z[C7+] for GRD Experiment
#========================================================================

def calcMoleFracGRD(Z,clsEOS) :

    zC1 = 0.0 ; zC7 = 0.0

    for iC in range(clsEOS.NC) :
        if clsEOS.gPP("MW",iC) < 25.0 : zC1 = zC1 + Z[iC]
        if clsEOS.gPP("MW",iC) > 90.0 : zC7 = zC7 + Z[iC]

#-- Return C1 & C7+ moles -------------------------------------------        

    return zC1,zC7

#========================================================================
#  End of Module
#========================================================================


