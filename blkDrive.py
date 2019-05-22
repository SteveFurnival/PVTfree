
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Driver for Blackoil Table Generation
#  1.  Find Psat
#  2.  For P =< Psat, Simulate Depletion Process as requested by User
#                     i.e. CCE, CVD or DLE
#                     Then flash Res Liq & Res Vap thru Multi-Stage SEPS
#  3.  For P > Psat, Extend Tables
#
#  Table Extension uses Singh, K., Fevang, O. & Whitson, C.H.,
#                       "Consistent Black-Oil PVT Table Modification"
#                       SPE 109,596, (2007).
#               https://www.onepetro.org/conference-paper/SPE-109596-MS
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import blkE100   as B1
import blkCoefs  as BC
import blkRegEOS as BE
import blkIMEX   as BI
import blkMORE   as BM
import blkVIP    as BN
import blkOther  as BO
import blkProps  as BP
import blkRegVis as BV

import calcEOS   as CE
import calcFlash as CF
import constants as CO
import calcReg   as CR
import calcSat   as CS
import calcWater as CW
import readExps  as RX
import utilities as UT
import writeOut  as WO

from math import exp,log

#=======================================================================
#  Driver for Blackoil Table Calculations
#=======================================================================

def calcBlack(clsIO,clsBLK,clsEOS,dicSAM,clsUNI) :

#== Set Output File Name [rootName.sim] ===============================

    if   clsBLK.tSim == "CMG" :

        if not clsIO.qIMX :
            pathIMX = clsIO.patR + ".imex"
            fIMX    = open(pathIMX,'w')
            clsIO.setQIMX(True)
            clsIO.setFIMX(fIMX)
        else :
            fIMX = clsIO.fIMX
            
        fSim = clsIO.fIMX
        sCom = "**"
        
    elif clsBLK.tSim == "MOR" :

        if not clsIO.qMOR :
            pathMOR = clsIO.patR + ".mor"
            fMOR    = open(pathMOR,'w')
            clsIO.setQMOR(True)
            clsIO.setFMOR(fMOR)
        else :
            fMOR = clsIO.fMOR
            
        fSim = clsIO.fMOR
        sCom = "--"
        
    elif clsBLK.tSim == "VIP" :

        if not clsIO.qVIP :
            pathVIP = clsIO.patR + ".vip"
            fVIP    = open(pathVIP,'w')
            clsIO.setQVIP(True)
            clsIO.setFVIP(fVIP)
        else :
            fVIP = clsIO.fVIP
            
        fSim = clsIO.fVIP
        sCom = "C "
        
    elif clsBLK.tSim == "ECL" :
        
        if not clsIO.q100 :
            path100 = clsIO.patR + ".e100"
            f100    = open(path100,'w')
            clsIO.setQ100(True)
            clsIO.setF100(f100)
        else :
            f100 = clsIO.fDep
            
        fSim = clsIO.f100
        sCom = "--"

#-- Write Header Information to Simulator File ----------------------        

    WO.outputHeader(fSim,sCom,clsIO)
    WO.outputEOS(fSim,sCom,clsIO,clsEOS)

#== Initialisation ====================================================    

    print("calcBlack: Entered")

    qLiq = True  ; qVap = False
    iLiq =  1    ; iVap = -1

    iSam = clsBLK.sNum
    nCom = clsEOS.NC

    clsSAM = dicSAM[iSam]

    Z = NP.zeros(nCom)
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

    sNam = dicSAM[iSam].sName

    xTyp = clsBLK.xTyp
    tRes = clsBLK.Tres
    pSep = clsBLK.pSep
    tSep = clsBLK.tSep

#-- Has User Defined Non-Standard Wiring for the Stage Order? -------

    if clsBLK.qLsep : Lsep = clsBLK.Lsep
    else            : Lsep = []

    if clsBLK.qVsep : Vsep = clsBLK.Vsep
    else            : Vsep = []

#-- Brine Density at Standard (Stock Tank) Conditions ---------------    

    pStc = CO.pStand    #-- Standard Pressure    [psia]
    tStc = CO.tStand    #-- Standard Temperature [degR = 60.0 degF]
    mFrc = clsBLK.bSalt #-- Mass Fraction of Salt in the Brine

    dSTW,comW = CW.calcRoweChouDen(mFrc,tStc,pStc,clsUNI)  #-- Stock Tank Brine Density

    clsBLK.setDenSTW(dSTW)

#-- Brine Properties at Reference Pressure & Tres (for PVTW keyword) 

    CW.calcPVTW(clsBLK,clsIO,clsUNI)

#======================================================================
#  Saturation Pressure Stage: Find Psat
#======================================================================

    clsEOS.eosCoefsNoPT()

#-- This experiment at varying temperature so set when needed -------

    pMea = -1.0
    qBub,pSat,Ksat = CS.calcPsat(pMea,tRes,clsEOS,clsSAM,clsIO)

    clsBLK.Psat = pSat

    Xsat = NP.zeros(nCom)
    Ysat = NP.zeros(nCom)
    
    if qBub :
        Xsat = NP.copy(Z)
        Ysat = NP.multiply(Z,Ksat)
    else :
        Xsat = NP.divide(Z,Ksat)
        Ysat = NP.copy(Z)

    Xsat = UT.Norm(Xsat) ; Ysat = UT.Norm(Ysat)

#-- Properties of the Feed and Incipient Phases at Psat -------------    

    Mx,Vx,Dx,Zx,Ux,Cx,Wx = CE.calcProps(iLiq,pSat,tRes,Xsat,clsEOS)
    My,Vy,Dy,Zy,Uy,dS,dS = CE.calcProps(iVap,pSat,tRes,Ysat,clsEOS)

    UgDew = Uy

#== Flash the Psat Liquid & Vapour through the Separator Train ========    

    if qBub : Vsat = Vx
    else    : Vsat = Vy

    Doo,Dgo,Rs,Bo,Bd,ySTG = \
        BO.sepFlash(iLiq,pSep,tSep,Lsep,Vsep,Xsat,Vx,clsEOS,clsIO)

    Dog,Dgg,Rv,Bg,Bd,ySTG = \
        BO.sepFlash(iVap,pSep,tSep,Lsep,Vsep,Ysat,Vy,clsEOS,clsIO)

    Md,Vd,Dd,Zd,Ud,dS,dS = CE.calcProps(iVap,pSat,tRes,ySTG,clsEOS)

    Bd = BO.gasFVF(pSat,tRes,Zd)

    RsSat = Rs      #-- Store the Saturated GOR & CGR 
    RvSat = Rv

#-- Stock Tank Densities from Saturated Stage depending on Type -----    

    if qBub :
        dSTO = Doo
        dSTG = Dgo
    else :
        dSTO = Dog
        dSTG = Dgg

    dTab = []  #-- Stores all Saturated data
    pTab = []  #-- Will store all values of Pres (including those > Psat)
    
    dRow = [pSat,Rs,Bo,Ux,Cx,Wx,Rv,Bg,Uy,Bd,Ud]
    dTab.append(dRow)

#== Setup Pressure Nodes in the Tables ================================    

    pMax = clsBLK.Pmax
    pMin = clsBLK.Pmin
    pInc = clsBLK.Pinc

    pRes = pMax
    qIns = True
    iRow = 0

#======================================================================
#  Build pTab-Array
#======================================================================

    while pRes >= pMin :

        if pRes > pSat :                #-- Above Psat
            pTab.append(pRes)
            iRow += 1
        else :
            if qIns :                   #-- Insert Psat
                pTab.append(pSat)
                qIns = False
                iSat = iRow
                iRow += 1
            pTab.append(pRes)           #-- Below Psat
            iRow += 1
            
        pRes = pRes - pInc              #-- Decrement Pressure

#-- Total number of stages (pressure nodes) -------------------------

    nPrs = len(pTab)

    aTab = NP.zeros((nPrs,11))

    aTab[iSat][:] = dRow                #-- Store Psat data in right aTab row

#======================================================================
#  Depletion Stages
#======================================================================

    if qBub : vEst = 0.0
    else :    vEst = 1.0

    zTot = 1.0
    pRes = pMax
    iRow = 0
    qIns = True

    while pRes >= pMin :

#-- Above Psat ------------------------------------------------------        

        if pRes > pSat :
            iRow += 1
        else :

#-- Psat-Row? -------------------------------------------------------

            if qIns :
                iRow += 1
                qIns = False

#-- This experiment at varying temperature so set when needed -------

            V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)
            vEst    = V

            Moil,Vliq,dLiq,Zoil,Uoil,Coil,Viso = CE.calcProps(iLiq,pRes,tRes,X,clsEOS)
            Mgas,Vvap,dVap,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,pRes,tRes,Y,clsEOS)

            zGas = zTot*     V
            zOil = zTot*(1.0-V)

            Xsep = NP.copy(X)  #-- Local copies of (X,Y) to flash thru Seps
            Ysep = NP.copy(Y)

#-- Flash Liquid thru the Separator Train ---------------------------

            Doo,Dgo,Rs,Bo,Bd,ySTG = \
                BO.sepFlash(iLiq,pSep,tSep,Lsep,Vsep,Xsep,Vliq,clsEOS,clsIO)

#-- Flash Vapour thru the Separator Train ---------------------------            
            
            Dog,Dgg,Rv,Bg,Bd,ySTG = \
                BO.sepFlash(iVap,pSep,tSep,Lsep,Vsep,Ysep,Vvap,clsEOS,clsIO)

            #Ud,dumm = CE.calcLBCderv(tRes,dVap,ySTG,clsEOS)  #-- Dry Gas Visc at pRes

            Md,Vd,Dd,Zd,Ud,dS,dS = CE.calcProps(iVap,pRes,tRes,ySTG,clsEOS)

            Bd = BO.gasFVF(pRes,tRes,Zd)

#== Modify the Moles to Next Stage, Depending on Depletion Type =======

            Voil = zOil*Vliq
            Vgas = zGas*Vvap
            
            Vtot = Voil + Vgas
            Vrem = Vtot - Vsat

            zRem = pRes*Vrem/(Zgas*CO.gasCon*tRes)

            if   xTyp == "CCE" :
                pass
            elif xTyp == "CVD" :
                for iC in range(nCom) : Z[iC] = (zTot*Z[iC] - zRem*Y[iC])/(zTot - zRem)
                zTot = zTot - zRem
            elif xTyp == "DLE" :
                Z    = NP.copy(X)
                zTot = zTot - zGas

            dRow = [pRes,Rs,Bo,Uoil,Coil,Viso,Rv,Bg,Ugas,Bd,Ud]
            dTab.append(dRow)

            aTab[iRow][:] = dRow    #-- Store data in aTab row

            iRow += 1

#== Decrement Pressure and Continue ===================================

        pRes = pRes - pInc

#== Work Arrays for "Back-Up" Fits ====================================        

    nSat = len(dTab)

    X = NP.zeros(nSat)
    Y = NP.zeros(nSat)

#======================================================================
#  Slope and Intercept of the (assumed linear) Bo versus Rs 
#======================================================================

    for i in range(nSat) :
        X[i] = dTab[i][clsBLK.iRs]
        Y[i] = dTab[i][clsBLK.iBo]

    slope,inter = UT.linearFit(X,Y)

    clsBLK.BoS = slope
    clsBLK.BoI = inter

#======================================================================
#  Slope and Intercept of the (assumed exp) Muo versus Rs 
#======================================================================

    for i in range(nSat) :
        X[i] =     dTab[i][clsBLK.iRs]
        Y[i] = log(dTab[i][clsBLK.iUo])

    slope,logIn = UT.linearFit(X,Y)
    inter = exp(logIn)

    clsBLK.UoS = slope
    clsBLK.UoI = inter

#======================================================================
#  Slope and Intercept of the (assumed linear) co.p versus Rs 
#======================================================================

    for i in range(nSat) :
        X[i] = dTab[i][clsBLK.iRs]
        Y[i] = dTab[i][clsBLK.iCo]*dTab[i][clsBLK.iPr]

    slope,inter = UT.linearFit(X,Y)

    clsBLK.CoS = slope
    clsBLK.CoI = inter

#======================================================================
#  Stock Tank Oil and Gas Properties
#======================================================================

    mSTO,oGrv = BP.initOilProps(dSTO,clsBLK,clsIO)
    mSTG,gGrv = BP.initGasProps(dSTG,clsBLK,clsIO)

    clsBLK.mSTO = mSTO
    clsBLK.mSTG = mSTG

    cCon = (dSTO/mSTO)*(CO.gasCon*CO.tStand/CO.pStand)  #-- Singh et al Eqn.(14)

    clsBLK.Co = cCon

#== Output Header =====================================================    

    WO.outputHeaderBO(fSim,iSam,sNam,clsBLK,clsIO,clsUNI)

#======================================================================
#  Generate the Saturated Oil and Gas STO Mole Fractions
#======================================================================

    cTab = []  #-- Table used to store 'calculated' data

    for iPrs in range(nSat) :

        Pr = dTab[iPrs][clsBLK.iPr]
        Rs = dTab[iPrs][clsBLK.iRs] ; Bo = dTab[iPrs][clsBLK.iBo]
        Rv = dTab[iPrs][clsBLK.iRv] ; Bg = dTab[iPrs][clsBLK.iBg]

        denO = BP.denOil(dSTO,dSTG,Rs,Bo,clsIO)
        denG = BP.denGas(dSTO,dSTG,Rv,Bg,clsIO)

        xLiq = cCon/(cCon +     Rs)
        yLiq = cCon/(cCon + 1.0/Rv)

        Mliq = BP.phaseMw(mSTG,mSTO,xLiq,clsIO)
        Mvap = BP.phaseMw(mSTG,mSTO,yLiq,clsIO)

        Vliq = Mliq/denO
        Vvap = Mvap/denG

        cRow = [xLiq,yLiq,denO,denG,Vliq,Vvap]
        
        cTab.append(cRow)

        #print("Pr,xO,yO,Vl,Vv {:10.3f}{:8.5f}{:8.5f}{:8.4f}{:8.4f}".format(Pr,xLiq,yLiq,Vliq,Vvap))

#======================================================================
#  Check the Oil Viscosity versus Pressure trend is physical
#  It can go wrong for condensates
#======================================================================

    qMonV = BP.oilViscMono(dTab,clsBLK)
    qMonC = BP.oilCompMono(dTab,clsBLK)

#== Initialise Extension Method =========================================        

    clsBLK.RT = CO.gasCon*tRes     #-- RT

    nEOS = clsBLK.nEOS

    clsBLK.EOS1["sOil"] = 0.5
    clsBLK.EOS1["sGas"] = 0.5

#========================================================================
#  Fit the 2-Component (STO & STG) EoS Parameters to Saturated Data
#========================================================================

    sOil,sGas = \
        BE.regEoSData(dTab,cTab,dicSAM,clsEOS,clsBLK,clsIO)

#========================================================================
#  Fit the 2-Component (STO & STG) Viscosity Parameters to Saturated Data
#========================================================================

    rOil,rGas = \
        BV.regViscData(dTab,cTab,dicSAM,clsEOS,clsBLK,clsIO)

#========================================================================
#  Calculate Convergence Pressure: Re-Set Saturated (sOil,sGas)
#========================================================================

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)
    
    pCon = BO.convPressure(dTab,clsBLK,clsIO)

    if pCon < pMax : pCon = pMax

    KoSat,KgSat,mO,mG = BO.extendTable(pSat,pCon,dTab,clsBLK,clsIO)

    RsMax,RvMax = BO.calcRsRv(cCon,pCon,pSat,KoSat,KgSat,mO,mG)

#== Extended Data =====================================================

    eTab = []
    uLst = dTab[0][clsBLK.iUo]
    iRow = 0

    for iExt in range(len(pTab)) :
        
        if pTab[iExt] > pSat :
            
            pRes  = pTab[iExt]
            RTp   = clsBLK.RT/pRes
            
            Rs,Rv = BO.calcRsRv(cCon,pRes,pSat,KoSat,KgSat,mO,mG)
            
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs ,clsBLK,clsIO)
            Co,Vo = BO.calcComp(qLiq,RTp,cCon,dSTO,dSTG,Rs,Bo,Uo,clsBLK,clsIO)

            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv ,clsBLK,clsIO)
            Bd,Ud = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,0.0,clsBLK,clsIO)

#-- For condensates, Uo can go astray -------------------------------

            if Uo > uLst : Uo = clsBLK.UoI*exp(Rs*clsBLK.UoS)
            uLst = Uo
            
            eRow = [pRes,Rs,Bo,Uo,Co,Vo,Rv,Bg,Ug,Bd,Ud]
            eTab.append(eRow)

            aTab[iRow][:] = eRow
            iRow         += 1

#========================================================================
#  Oil and Gas Output depends on Simulator Type
#========================================================================

    if   clsBLK.tSim == "MOR" :
        BM.outMORE(fSim,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)
    elif clsBLK.tSim == "CMG" :
        BI.outIMEX(fSim,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)
    elif clsBLK.tSim == "VIP" :
        BN.outVIP( fSim,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)
    else :
        B1.outE100(fSim,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)

#========================================================================
#  End of Routine
#========================================================================

    print("calcBlack: Finished")

    return

#========================================================================
#  End of Module
#========================================================================
