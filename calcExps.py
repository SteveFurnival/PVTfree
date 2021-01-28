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

import allData   as AD
import allProp   as AP
import calcEOS   as CE
import calcFlash as CF
import calcGrad  as CG
import calcSat   as CS
import utilities as UT
import writeOut  as WO

#========================================================================
#  Calculate Experiments
#========================================================================

def calcExps(sExt,clsEOS,dicSAM,dicEXP,clsUNI,clsIO) :

    nExp = len(dicEXP)
    nCom = clsEOS.nComp

    sLabl = "calcExps: Running Experiments  " + sExt
    print(sLabl)

#== Setup the (P,T)-Independent EoS Coefficients ======================

    clsEOS.eosCoefsNoPT()

#== Loop over experiments =============================================

    qDif = False    #-- We don't have a previous solution

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            xTyp = clsEXP.xName

            sOut = "Experiment: iExp,Type {:2d} ".format(iExp+1)
            sOut = sOut + xTyp

            print(sOut)

            clsEXP.setPsatInfo(None,None,NP.empty(nCom))

            if   xTyp == "CCE" :    calcCCE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "CVD" :    calcCVD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "DLE" :    calcDLE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SEP" :    calcSEP(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "FLS" :    calcFLS(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SAT" :    calcSAT(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SWL" :    calcSWL(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "GRD" : CG.calcGRD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)

#== Write the experiments to the Print Output File (fOut) =============

    WO.outputExps(clsEOS,dicSAM,dicEXP,clsUNI,clsIO)
    
    sLabl = "calcExps: Experiments Completed" + sExt
    print(sLabl)

    return

#========================================================================
#  Constant Composition Expansion (CCE) Experiment
#========================================================================

def calcCCE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sNam

#== Independent & Calculated Data =====================================

    typIND = AP.classLIB().INDshrt.get("CCE")
    typCAL = AP.classLIB().CALshrt.get("CCE")

    iPR = typIND.index("PRES")

    iRV = typCAL.index("RELV") ; iSL = typCAL.index("SLIQ")
    iZF = typCAL.index("ZFAC") ; iDO = typCAL.index("DENO")
    iUG = typCAL.index("VISG") ; iUO = typCAL.index("VISO")

    iMO = typCAL.index("MWO" ) ; iMG = typCAL.index("MWG" )
    iDG = typCAL.index("DENG") ; iFT = typCAL.index("IFT" )
    iHO = typCAL.index("HO"  ) ; iHG = typCAL.index("HG"  )
    iCO = typCAL.index("CPO" ) ; iCG = typCAL.index("CPG" )

#-- Load Feed Composition -------------------------------------------

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Saturation Pressure ---------------------------------------------

    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0

    qBub,pSat,logK = clsEXP.getPsatInfo(qDif,nCom)
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

    if not qDif : clsEXP.setPsatInfo(qBub,pSat,NP.log(Ksat))

    clsEXP.setPsatCal(pSat)

    iNeu =  0
    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,Z,clsEOS)

    #if qBub : vEst = 0.0
    #else :    vEst = 1.0

    vEst = None

#======================================================================
#  Loop over user defined pressures
#======================================================================

    for iPrs in range(nPrs) :

        pRes = clsEXP.dInd[iPR][iPrs]

#== 1-Phase or 2-Phase? ===============================================

        if pRes > pSat :
            
            Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iNeu,pRes,tRes,Z,clsEOS)
            #print("calcCCE[Z]: iP,pR,Mw,Vm,Ro,ZF,vS {:2d} {:8.2f} {:7.3f} {:8.4f} {:8.3f} {:8.4f} {:8.4f}".format(iPrs,pRes,M1P,V1P,D1P,Z1P,V1P))

            relV = Vliq/Vsat
            
            if qBub : sLiq = 1.0
            else :    sLiq = 0.0

            Mgas = Moil ; Vvap = Vliq
            Dgas = Doil ; Zgas = Zoil
            Ugas = Uoil

            IFT  = 0.0  #-- Surface Tension [dyne/cm]

            HO,CPO,uJTO = \
                CE.calcEnthSpecHeat(iNeu,pRes,tRes,Z,clsEOS)  #-- Enthalpy & Spec Heat

            HG = HO ; CPG = CPO ; uJTG = uJTO        #-- Gas = Oil = 1-Phase!

            #print("P,uJTO,uJTG,IFT {:8.2f} {:10.3e} {:10.3e} {:8.4f}".format(pRes,uJTO,uJTG,IFT))
            
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

            if sLiq < 0.0 :
                crash = 1.0/0.0

            IFT = CE.calcIFT(Vliq,Vvap,X,Y,clsEOS)

            HO,CPO,uJTO = \
                CE.calcEnthSpecHeat(iLiq,pRes,tRes,X,clsEOS)  #-- Enthalpy & Spec Heat
            HG,CPG,uJTG = \
                CE.calcEnthSpecHeat(iVap,pRes,tRes,Y,clsEOS)  #-- Enthalpy & Spec Heat

            #print("P,uJTO,uJTg {:8.2f} {:10.3e} {:10.3e} {:8.4f}".format(pRes,uJTO,uJTG,IFT))

#== Store the data for this stage =====================================

        clsEXP.dCal[iRV][iPrs] = relV ; clsEXP.dCal[iSL][iPrs] = sLiq
        clsEXP.dCal[iZF][iPrs] = Zgas ; clsEXP.dCal[iDO][iPrs] = Doil
        clsEXP.dCal[iUG][iPrs] = Ugas ; clsEXP.dCal[iUO][iPrs] = Uoil

        clsEXP.dCal[iMO][iPrs] = Moil ; clsEXP.dCal[iMG][iPrs] = Mgas
        clsEXP.dCal[iDG][iPrs] = Dgas ; clsEXP.dCal[iFT][iPrs] = IFT
        clsEXP.dCal[iHO][iPrs] = HO   ; clsEXP.dCal[iHG][iPrs] = HG
        clsEXP.dCal[iCO][iPrs] = CPO  ; clsEXP.dCal[iCG][iPrs] = CPG

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Constant Volume Depletion (CVD) Experiment
#========================================================================

def calcCVD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    if clsIO.Deb["CVD"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sNam

    if qDeb :
        sOut = "calcCVD: nCom,nSam,nRow {:2d} {:2d} {:2d}\n".format(nCom,nSam,nPrs)
        fDeb.write(sOut)

#-- Independent & Calculated Data -----------------------------------    

    typIND = AP.classLIB().INDshrt.get("CVD")
    typCAL = AP.classLIB().CALshrt.get("CVD")

    iPR = typIND.index("PRES")

    iMR = typCAL.index("MREM") ; iSL = typCAL.index("SLIQ")
    iZF = typCAL.index("ZFAC")
    iMO = typCAL.index("MWO" ) ; iMG = typCAL.index("MWG" )
    iDO = typCAL.index("DENO") ; iDG = typCAL.index("DENG")
    iUO = typCAL.index("VISO") ; iUG = typCAL.index("VISG")

#-- Load Feed Composition -------------------------------------------

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Saturation Pressure ---------------------------------------------
    
    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,logK = clsEXP.getPsatInfo(qDif,nCom)
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

    if qDeb :
        sOut = "calcCVD: tRes,qBub,pSat {:8.2f} {:} {:10.3f}\n".format(tRes,qBub,pSat)
        fDeb.write(sOut)

    if not qDif : clsEXP.setPsatInfo(qBub,pSat,NP.log(Ksat))

    clsEXP.setPsatCal(pSat)

    iNeu =  0
    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,Z,clsEOS)

    #if qBub : vEst = 0.0
    #else :    vEst = 1.0

    vEst = None

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zMol = 1.0
    zTot = 0.0

    for iPrs in range(nPrs) :

        pRes = clsEXP.dInd[iPR][iPrs]

#== 1-Phase or 2-Phase? ===============================================

        if pRes > pSat :
            
            Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iNeu,pRes,tRes,Z,clsEOS)

            if qBub : sLiq = 1.0
            else :    sLiq = 0.0

            zTot = 0.0
            Mgas = Moil ; Dgas = Doil
            Zgas = Zoil ; Ugas = Uoil
            
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

            zRem = pRes*Vrem/(Zgas*UT.gasCon*tRes)  #-- Moles Removed

            for iC in range(nCom) :
                Z[iC] = (zMol*Z[iC] - zRem*Y[iC])/(zMol - zRem)

            zMol = zMol - zRem
            zTot = zTot + zRem

            if qDeb :
                sOut = "calcCVD: zRem {:10.3e}\n".format(zRem)
                fDeb.write(sOut)

#-- Load Data -------------------------------------------------------

        clsEXP.dCal[iMR][iPrs] = zTot ; clsEXP.dCal[iSL][iPrs] = sLiq
        clsEXP.dCal[iZF][iPrs] = Zgas
        clsEXP.dCal[iMO][iPrs] = Moil ; clsEXP.dCal[iMG][iPrs] = Mgas
        clsEXP.dCal[iDO][iPrs] = Doil ; clsEXP.dCal[iDG][iPrs] = Dgas
        clsEXP.dCal[iUO][iPrs] = Uoil ; clsEXP.dCal[iUG][iPrs] = Ugas

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Differential Liberation (DLE) Experiment
#========================================================================

def calcDLE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sNam

    #print("calcCVD: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typIND = AP.classLIB().INDshrt.get("DLE")
    typOBS = AP.classLIB().OBSshrt.get("DLE")

    iPR = typIND.index("PRES")

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
    
    qBub,pSat,logK = clsEXP.getPsatInfo(qDif,nCom)
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

    if not qDif : clsEXP.setPsatInfo(qBub,pSat,NP.log(Ksat))

    clsEXP.setPsatCal(pSat)

    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

    #if qBub : vEst = 0.0
    #else :    vEst = 1.0

    vEst = None

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zMol = 1.0
    rTot = 0.0
    psTs = UT.pStand/UT.tStand

    for iPrs in range(nPrs) :

        pRes = clsEXP.dInd[iPR][iPrs]

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

            zRem = pRes*Vgas/(Zgas*UT.gasCon*tRes)
            zMol = zMol - zRem

            Vrem = zRem*UT.volMol

            rTot = rTot + Vrem

            Z = NP.copy(X)

            Vtot = Voil + Vgas
            Bgas = psTs*Zgas*tRes/pRes
            gGrv = Mgas/UT.molAir

#-- Temporary Storage Ahead of Stock Tank Corrections ---------------            

        clsEXP.dCal[iBO][iPrs] = Voil ; clsEXP.dCal[iRS][iPrs] = Vrem
        clsEXP.dCal[iDO][iPrs] = Doil ; clsEXP.dCal[iBT][iPrs] = Vtot
        clsEXP.dCal[iBG][iPrs] = Bgas ; clsEXP.dCal[iZF][iPrs] = Zgas
        clsEXP.dCal[iGG][iPrs] = gGrv ; clsEXP.dCal[iUO][iPrs] = Uoil

#== Stock Tank Volume =================================================

    Mst,Vst,Dst,Zst,Ust,dmS,dmS = CE.calcProps(iLiq,UT.pStand,UT.tStand,Z,clsEOS)

    Vsto = zMol*Vst

#== "Normalise" via Stock Tank Volume =================================

    for iPrs in range(nPrs) :

        rTot = rTot - clsEXP.dCal[iRS][iPrs]

        clsEXP.dCal[iBO][iPrs] = clsEXP.dCal[iBO][iPrs]/Vsto
        clsEXP.dCal[iRS][iPrs] =                   rTot/Vsto

        clsEXP.dCal[iBT][iPrs] = clsEXP.dCal[iBO][iPrs] + \
                                 clsEXP.dCal[iBG][iPrs]*(clsEXP.dCal[iRS][0] - clsEXP.dCal[iRS][iPrs])

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Separator Test (SEP) Experiment
#========================================================================

def calcSEP(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    if clsIO.Deb["SEP"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sNam

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typIND = AP.classLIB().INDshrt.get("SEP")
    typOBS = AP.classLIB().OBSshrt.get("SEP")

    iPR = typIND.index("PRES")
    iTR = typIND.index("TEMP")

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

    pObs = clsEXP.dInd[iPR][0]
    tRes = clsEXP.dInd[iTR][0]

    #print("calcSEP: pObs,tRes ",pObs,tRes)

#-- Saturation Pressure ---------------------------------------------
    
    if clsEXP.PsatO > 0.0 : pObs = clsEXP.PsatO
    else                  : pObs = -1.0
    
    qBub,pSat,logK = clsEXP.getPsatInfo(qDif,nCom)
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

    #print("CalcSEP: Psat {:10.3f}".format(pSat))

    if not qDif : clsEXP.setPsatInfo(qBub,pSat,NP.log(Ksat))

    clsEXP.setPsatCal(pSat)

    iLiq =  1
    iVap = -1

    Msat,Vsat,Dsat,Zsat,Usat,dumS,dumS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

    #vEst = 0.5

    vEst = None

#======================================================================
#  Loop over user defined pressures
#======================================================================

    zGas = 0.0
    zOil = 0.0
    zTot = 1.0

    for iPrs in range(1,nPrs) :

        pRes = clsEXP.dInd[iPR][iPrs]
        tRes = clsEXP.dInd[iTR][iPrs]

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
        Vgas = zGas*UT.volMol

#-- Temporary storage for Bo & Rs ahead of divide by Vsto below -----        

        clsEXP.dCal[iBO][iPrs] = Voil
        clsEXP.dCal[iRS][iPrs] = Vgas
        clsEXP.dCal[iDO][iPrs] = Doil
        clsEXP.dCal[iGG][iPrs] = Mgas/UT.molAir

#== Stock Tank Oil Volume =============================================

    nTot = 0.0
    rTot = 0.0
    gTot = 0.0
    Vsto = clsEXP.dCal[iBO][nPrs-1]

    for iPrs in range(1,nPrs) :

        clsEXP.dCal[iBO][iPrs] = clsEXP.dCal[iBO][iPrs]/Vsto  #-- Bo
        clsEXP.dCal[iRS][iPrs] = clsEXP.dCal[iRS][iPrs]/Vsto  #-- GOR

        nTot = nTot + molG[iPrs]
        rTot = rTot + molG[iPrs]*clsEXP.dCal[iRS][iPrs]
        gTot = gTot + molG[iPrs]*clsEXP.dCal[iGG][iPrs]

#== Saturation Pressure Stage =========================================

    if nTot > 0.0 : gTot = gTot/nTot
    else :          gTot = 0.0

    #print("mSTG,mSTO {:10.5f} {:10.5f}".format(mSTG,mSTO))

    clsEXP.dCal[iBO][0] = Vsat/Vsto
    clsEXP.dCal[iRS][0] = mSTG*UT.volMol/Vsto
    clsEXP.dCal[iDO][0] = Dsat
    clsEXP.dCal[iGG][0] = gTot
    
#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  2-Phase Flash (FLS) Experiment
#========================================================================

def calcFLS(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nPrs = clsEXP.nRow

    iLiq =  1
    iVap = -1

    sNam = dicSAM[nSam].sNam

    #print("calcFLS: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typIND = AP.classLIB().INDshrt.get("FLS")
    typOBS = AP.classLIB().OBSshrt.get("FLS")

    iPR = typIND.index("PRES")
    iTR = typIND.index("TEMP")

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

        pRes = clsEXP.dInd[iPR][iPrs]
        tRes = clsEXP.dInd[iTR][iPrs]

#-- 2-Phase Flash ---------------------------------------------------        

        vEst    = 0.5
        V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)

#-- Properties of the Liquid and Vapour Output of this stage --------        

        Moil,Vliq,Doil,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
        Mgas,Vvap,Dgas,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

        clsEXP.dCal[iZF][iPrs] = Zgas
        clsEXP.dCal[iDO][iPrs] = Doil
        clsEXP.dCal[iVF][iPrs] = V

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Saturation Pressure (SAT) Experiment
#========================================================================

def calcSAT(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nTem = clsEXP.nRow

    clsSAM = dicSAM[nSam]

    sNam = clsSAM.sNam

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typIND = AP.classLIB().INDshrt.get("SAT")
    typOBS = AP.classLIB().OBSshrt.get("SAT")

    iTR = typIND.index("TEMP")

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

        tRes  = clsEXP.dInd[iTR][iTem]
        pSatO = clsEXP.dObs[iPS][iTem]

        if pSatO > 0.0 : pEst =  0.95*pSatO
        else           : pEst = -1.0
        
        qBub = None
        pSat = None
        logK = NP.empty(nCom)
    
        qBub,pSat,Ksat = \
            CS.calcPsat(pEst,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

        iLiq = 0
        MPs,VPs,DPs,ZPs,UPs,dmS,dmS = CE.calcProps(iLiq,pSat,tRes,Z,clsEOS)

        clsEXP.dCal[iPS][iTem] = pSat
        clsEXP.dCal[iZF][iTem] = ZPs
        clsEXP.dCal[iDO][iTem] = DPs

        #print("calcSAT: tRes,qBub,pSat ",tRes,qBub,pSat)

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  Swelling Test (SWL) Experiment
#========================================================================

def calcSWL(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nInj = clsEXP.nSinj
    nAdd = clsEXP.nRow

    sNam = dicSAM[nSam].sNam

    #print("calcSEP: nCom,nSam,nRow ",nCom,nSam,nPrs)

    typIND = AP.classLIB().INDshrt.get("SWL")
    typOBS = AP.classLIB().OBSshrt.get("SWL")

    iMO = typIND.index("MOLE")

    iPS = typOBS.index("PSAT")
    iVS = typOBS.index("VSWL")

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#-- Load Feed Composition -------------------------------------------    

    clsSAM = dicSAM[nSam]
    clsINJ = dicSAM[nInj]

    clsWRK = AD.classSample("SWLwrk")
    clsWRK.setIntComp(nCom,nCom)

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

        mInj = clsEXP.dInd[iMO][iAdd]

        mTot = mInj + 1.0

        #print("calcSWL: iAdd,mAdd,mInj,mTot ",iAdd,mThs,mInj,mTot)

        wDiv = 1.0/mTot
        for iC in range(nCom) :
            W[iC] = (Z[iC] + mInj*Y[iC])*wDiv
            clsWRK.sZI(iC,W[iC])

        pSatO = clsEXP.dObs[iPS][iAdd]

        if pSatO > 0.0 : pEst =  0.95*pSatO
        else           : pEst = -1.0

        qBub = None
        pSat = None
        logK = NP.empty(nCom)
        
        qBub,pSat,Ksat = \
            CS.calcPsat(pEst,tRes,qBub,pSat,logK,clsEOS,clsWRK,clsIO)

        iLiq = 0
        MPs,VPs,DPs,ZPs,UPs,dmS,dmS = CE.calcProps(iLiq,pSat,tRes,W,clsEOS)

        if iAdd == 0 : Vref = VPs

        Vswl = mTot*VPs/Vref

        clsEXP.dCal[iPS][iAdd] = pSat
        clsEXP.dCal[iVS][iAdd] = Vswl

#======================================================================
#  End of Module
#======================================================================

    return

#========================================================================
#  End of Module
#========================================================================


