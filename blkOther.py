
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Performs Other Calculations for Blackoil Table Generation
#
#  Including Separator Train Flash
#            Convergence Pressure Estimate
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import blkCoefs  as BC
import blkProps  as BP
import calcFlash as CF
import calcEOS   as CE
import utilities as UT

from math import exp,log,log10

#========================================================================
#  Separator Flash
#========================================================================

def sepFlash(iRat,pSep,tSep,Lsep,Vsep,Z,vRes,clsEOS,clsIO) :

    nSep = len(pSep)
    nCom = clsEOS.nComp

    vEst = 0.5
    iLiq =   1 ; iVap = -1
    zGas = 0.0 ; zOil = 0.0 ; zTot = 1.0

#-- Accumulate the moles of liquid and vapour -----------------------

    molZ = NP.zeros( nSep)
    comZ = NP.zeros((nCom,nSep))

    xSTO = NP.zeros(nCom) ; ySTG = NP.zeros(nCom)

    mSTO = 0.0 ; mSTG = 0.0

    molZ[0] = zTot
    for iC in range(nCom) : comZ[iC][0] = Z[iC]

#== Loop over the stages of the Train =================================    

    for iSep in range(nSep) :

        pRes = pSep[iSep] ; tRes = tSep[iSep]
        
        for iC in range(nCom) : Z[iC] = comZ[iC][iSep]

        Z = UT.Norm(Z)

#-- This experiment at varying temperature so set when needed -------

        V,K,X,Y = CF.calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO)
        
        if   V < 0.0 : V = 0.0  #-- calcFlash sets X=Z & Y=Z if 1-Phase
        elif V > 1.0 : V = 1.0
                
        vEst    = V

        zGas = molZ[iSep]*     V
        zOil = molZ[iSep]*(1.0-V)

        #print("iSep,zInp,xOut,yOut {:2d} {:8.5f} {:8.5f} {:8.5f}".format(iSep,molZ[iSep],zOil,zGas))

#-- Where is Liquid Output Going Next? ------------------------------

        kSepL = Lsep[iSep]

        if kSepL == 0 :
            if mSTO + zOil > 0.0 :
                for iC in range(nCom) :
                    xSTO[iC] = (xSTO[iC]*mSTO + zOil*X[iC])/(mSTO + zOil)
                mSTO = mSTO + zOil
            else :
                xSTO = NP.copy(X)
        else :
            for iC in range(nCom) :
                comZ[iC][kSepL] = (comZ[iC][kSepL]*molZ[kSepL] + X[iC]*zOil)/(molZ[kSepL] + zOil) 
            molZ[kSepL] = molZ[kSepL] + zOil
            #print("toLiq[kSepL],molZ {:2d} {:8.5f}".format(kSepL,molZ[kSep]))
            
#-- Where is Vapour Output Going Next? ------------------------------

        kSepV = Vsep[iSep]

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
            #print("toVap[kSepV],molZ {:2d} {:8.5f}".format(kSepV,molZ[kSep]))

        #print("Lsep,Vsep ",Lsep[iSep],Vsep[iSep])

#========================================================================
#  Process the data: Stock Tank Moles, Volumes and Ratios
#========================================================================

    xSum = NP.sum(xSTO) ; ySum = NP.sum(ySTG)

    Moil,Vliq,dLiq,Zoil,Uoil,dumS,dumS = CE.calcProps(iLiq,UT.pStand,UT.tStand,xSTO,clsEOS)
    Mgas,Vvap,dVap,Zgas,Ugas,dumS,dumS = CE.calcProps(iVap,UT.pStand,UT.tStand,ySTG,clsEOS)

    vGas = mSTG*UT.volMol
    vOil = mSTO*Vliq

    if iRat == 1 :
        GOR = vGas/vOil
        FVF = vRes/vOil
    else :
        GOR = vOil/vGas
        FVF = vRes/vGas

    #print("Vvap,vGas {:10.5f} {:10.5f}".format(Vvap,vGas))

#-- Calculate Dry Gas FVF & Viscosity (if appropriate) --------------

    if iRat == 1 :
        Bdry = 0.0
        Udry = 0.0
    else :
        Bdry = zGas*Vvap/vGas
        Udry = Ugas

#== Return Values =====================================================
        
    return dLiq,dVap,GOR,FVF,Bdry,ySTG

#========================================================================
#  Estimate the Convergence Pressure: Singh Appendix B
#========================================================================

def convPressure(dTab,clsBLK) :

    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG
    mSTO = clsBLK.mSTO
    mSTG = clsBLK.mSTG
    slop = clsBLK.BoS
    intr = clsBLK.BoI
    Psat = clsBLK.Psat

    RsSat = dTab[0][clsBLK.iRs]
    RvSat = dTab[0][clsBLK.iRv]

    KoSat = (RsSat + cCon)/(1.0/RvSat + cCon)   #-- Singh Eqn.(3)
    
    KgSat =  KoSat/(RsSat*RvSat)                #-- Singh Eqn.(4)

    RsPK = - cCon*log(KoSat)/log(KgSat)         #-- Singh Eqn.(B1)
    BoPK = slop*RsPK + intr                     #-- Bo vs Rs at Conv Pres

    DoPK = BP.denOil(dSTO,dSTG,RsPK,BoPK)
    xOPK = cCon/(cCon + RsPK)
    MoPK = BP.phaseMw(mSTG,mSTO,xOPK)
    VoPK = MoPK/DoPK

    aCPK = BC.Acoef(xOPK,clsBLK)
    bCPK = BC.Bcoef(xOPK,clsBLK)
    cCPK = BC.Ccoef(xOPK,clsBLK)
    
    uCPK = BC.Ucoef(VoPK,bCPK,cCPK)
    wCPK = BC.Wcoef(VoPK,bCPK,cCPK)

    pCon = clsBLK.RT/uCPK - aCPK/wCPK

#== Return Convergence Pressure =======================================

    #print("convPressure: pCon {:10.3f}".format(pCon))

    return pCon

#========================================================================
#  Extended (P > Psat) Values
#========================================================================

def extendTable(pSat,pCon,dTab,clsBLK) :

    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG
    
    pRat = log10(pSat/(pCon+1.0))  #-- 1 psi "shift" to prevent K = 1 exactly

    #print("Ps ,Pk  {:10.3e} {:10.3e}".format(pSat,pCon))

#-- Saturated Rs and Rv ---------------------------------------------    

    RsSat = dTab[0][clsBLK.iRs]
    RvSat = dTab[0][clsBLK.iRv]

    #print("Rss,Rvs {:10.3e} {:10.3e}".format(RsSat,RvSat))

#-- Saturated oil and gas K-values: Singh Eqn.(3) and (4) -----------

    KoSat = (RsSat + cCon)/(1.0/RvSat + cCon)
    KgSat =  KoSat/(RsSat*RvSat)

    #print("Kos,Kgs {:10.3e} {:10.3e}".format(KoSat,KgSat))

#-- Singh Eqn.(6) and (8) -------------------------------------------

    mO = log10(KoSat)/pRat
    mG = log10(KgSat)/pRat

    #print("mO ,mG  {:10.3e} {:10.3e}".format(mO,mG))

    return KoSat,KgSat,mO,mG

#========================================================================
#  Calculate (Rs,Rv)
#========================================================================

def calcRsRv(cCon,pRes,pSat,KoSat,KgSat,mO,mG) :

    Ko = KoSat*pow(pRes/pSat,mO)    #-- Singh Eqn.(5)
    Kg = KgSat*pow(pRes/pSat,mG)    #-- Singh Eqn.(6)

    Rs = cCon*(1.0-Ko)/(Kg-1.0)     #-- Singh Eqn.(9)
    Rv = Ko/(Kg*Rs)                 #-- Singh Eqn.(10)

    #print("p,Ko,Kg,Rs,Rv {:10.3f} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(pRes,Ko,Kg,Rs,Rv))

#== Return values =====================================================    
    
    return Rs,Rv

#========================================================================
#  Calculate Saturated Properties
#========================================================================

def calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs,clsBLK) :

    if qLiq :
        xO  = cCon/(cCon + Rs)
        BoC = dSTO + Rs*dSTG
    else :
        xO = cCon*Rs/(cCon*Rs + 1.0)
        BoC = dSTG + Rs*dSTO

    Bo = BP.calcPhaseFVF(qLiq,RTp,BoC,xO,clsBLK)
    dO = BoC/Bo
    uO = BP.calcLBCvisc(xO,dO,clsBLK)

#== Return values =====================================================    
    
    return Bo,uO

#========================================================================
#  Calculate Compressibility and Viscosibility by Difference
#========================================================================

def calcComp(qLiq,RTp,cCon,dSTO,dSTG,Rs,Bo,Uo,clsBLK) :

    pCur = clsBLK.RT/RTp    #-- Back-out the Current Pressure
    
    delP = pCur*0.0001      #-- Perturbation
    pPer = pCur + delP
    
    RTpP = clsBLK.RT/pPer

    BoP,UoP = calcSatProp(qLiq,RTpP,cCon,dSTO,dSTG,Rs ,clsBLK)

    dBdp = (BoP - Bo)/delP
    dUdp = (UoP - Uo)/delP

    Co = - dBdp/Bo
    Vo =   dUdp/Uo

#== Return values =====================================================

    return Co,Vo
    
#========================================================================
#  Gas FVF [rcf/scf]
#========================================================================

def gasFVF(pRes,tRes,Z) :

    Bg = (UT.pStand/UT.tStand)*Z*tRes/pRes

#== Return value ======================================================

    return Bg

#========================================================================
#  Set EOS1/VIS1 Parameters from sOil/sGas and rOil/rGas
#========================================================================

def setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK) :

    clsBLK.EOS1["sOil"] = sOil[iSat]
    clsBLK.EOS1["sGas"] = sGas[iSat]
    
    clsBLK.VIS1["rOil"] = rOil[iSat]
    clsBLK.VIS1["rGas"] = rGas[iSat]

#== No return value ===================================================

    return

#========================================================================
#  Lee & Gonzalez Gas Viscosity Correlation
#  see Whitson & Brule, Section 3.3.5
#========================================================================

def gasViscLee(Rv,Bg,Tran,clsBLK) :

    denLB = (clsBLK.dSTG + Rv*clsBLK.dSTO)/Bg  #-- Gas Density [lb/ft3]
    denGM = denLB/62.428                       #-- gm/cm3

    yO   = clsBLK.Co*Rv/(clsBLK.Co*Rv + 1.0)

    molG = yO*clsBLK.mSTO + (1.0-yO)*clsBLK.mSTG

    A1 = (9.379+0.01607*molG)*pow(Tran,1.5)/(209.2+19.26*molG+Tran)

    A2 = 3.448 + 986.4/Tran + 0.01009*molG

    A3 = 2.447 - 0.2224*A2

    #print("A1,A2,A3 {:10.3e} {:10.3e} {:10.3e}".format(A1,A2,A3))

#== Return Value [cP] =================================================

    muLee = 0.0001*A1*exp(A2*pow(denGM,A3))

    #print("Rv,yO,denGM,molG,muLee {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(Rv,yO,denGM,molG,muLee))

    return muLee

#========================================================================
#  End of Module
#========================================================================
