#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Michelsen "Fast" (Approximate) Phase Plot Module (CP)
#  Ref: Michelsen, M.L.,
#  "A Simple Method for Calculation of Approximate Phase Boundaries"
#  Fluid Phase Equilibria, (1994), vol. 98, p. 1
#
#========================================================================

#!/usr/bin/python3

import matplotlib.pyplot as PL
import numpy             as NP
import numpy.linalg      as LA
import scipy.linalg      as SA

import os

import calcEOS   as CE
import calcFlash as CF
import calcReg   as CR
import calcSat   as CS
import utilities as UT
import writeOut  as WO

from math import copysign,exp,log,sqrt

#========================================================================
#  Main Routine (All Samples)
#========================================================================

def allSamplesPhasePlot(clsEOS,dicSAM,clsIO) :

    if clsIO.Deb["ENV"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else                    :
        qDeb = False
        fDeb = ""

    nSam = len(dicSAM)
    nCom = clsEOS.nComp

    #Z = NP.zeros(nCom)

    for iSam in range(nSam) :

        clsSAM = dicSAM[iSam]

        sNam = clsSAM.sNam

        print("Generating Phase Plot for Sample ",iSam+1,sNam)

        BubT,BubP,DewT,DewP,FlsT,FlsP,pCri,tCri = \
            approxPhasePlotT(clsEOS,clsSAM,clsIO,qDeb,fDeb)

        nBub = len(BubT) ; nDew = len(DewT) ; nFls = len(FlsT)

        #print("nBub,nDew ",nBub,nDew)

        clsSAM.setTempPsat(BubT,BubP,DewT,DewP,FlsT,FlsP,pCri,tCri)

#== Plot the Phase Plot ===============================================

        if nBub > 0 or nDew > 0 or nFls > 0 or pCri != None :
            plotPhasePlot(clsSAM)
        else :
            clsSAM.qPhas = False
            print("No Phase Plot Generated for Sample ",sNam)

#== Return from routine ===============================================        

    return

#========================================================================
#  Iterate in P or T and Alpha at a Given T or P
#========================================================================

def approxPhasePlotT(clsEOS,clsSAM,clsIO,qDeb,fDeb) :

    BubT = [] ; BubP = []
    DewT = [] ; DewP = []
    FlsT = [] ; FlsP = []

    pCri = None ; tCri = None

    pMin = 1.0/UT.macEPS ; pMax = 0.0
    tMin = 1.0/UT.macEPS ; tMax = 0.0
    
#== Calculate the Psat at low temperature =============================

    nCom = clsEOS.nComp
    
    Z = NP.zeros(nCom)
    
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)
  
    tRes = 519.67  #-- 60 degF => degR
  
    pSat =   -1.0
    tInc =  -20.0
    pInc =  100.0
    alfa =    1.0
    
    eD   =  -1.0
    nItr =     0

#== Setup EoS Coefficients ============================================

    clsEOS.eosCoefsNoPT()       #-- No P/T-Dependency => bi, api, aqi

#== Initial Saturation Temperature ====================================
    
    pRes = 300.0     #-- Michelsen Suggests 20 bar
    
    qBub,tSat,Ksat = CS.calcTsat(pRes,clsEOS,clsSAM,clsIO)

    if tSat < 0.0 :
        print("Could not find Initial Tsat for Approx Phase Plot")
        return BubT,BubP,DewT,DewP,FlsT,FlsP,pCri,tCri

    pMin = min(pMin,pRes) ; pMax = max(pMax,pRes)
    tMin = min(tMin,tSat) ; tMax = max(tMax,tSat)

    tLst = tSat ; pLst = pRes ; aLst = alfa

    pSatR = pRes ; tResI = tSat     #-- Calculated Psat at initial tRes

    if qDeb :
        sOut = "Initial Tsat: qBub " + str(qBub) + "\n"
        fDeb.write(sOut)

    if qBub :
        BubT.append(tSat) ; BubP.append(pRes)
    else    :
        DewT.append(tSat) ; DewP.append(pRes)

    #degF = tRes - 459.67
    if qDeb :
        sOut = "pRes,nItr,tSat,alfa,eD           {:8.3f} {:2d} {:10.3f} {:8.5f} {:10.3e}\n".format(pRes,nItr,tSat,alfa,eD)
        fDeb.write(sOut)
        WO.writeArrayDebug(fDeb,Ksat,"Init Ksat")

    lnKref = NP.log(Ksat)    #-- Log(Reference K-Values)

#----------------------------------------------------------------------
#  Increment Pressure
#----------------------------------------------------------------------

    pSat = pRes + pInc
    tRes = tSat
    qPrs = False             #-- At given Pres, Iterate in Temp

    nStp = 1

    while pSat > 0.0 :

#== Update the (P or T,Alpha) at New (T or P) =========================

        nItr,pSat,tRes,alfa,eD = \
            iterPorTAlfa(qPrs,pSat,tRes,alfa,Z,lnKref,clsEOS,clsIO)

        if pSat <     0.0 or tRes <  400.0 or \
           pSat > 15000.0 or tRes > 2000.0    :
            #print("pSat,tRes {:10.3f} {:10.3f}".format(pSat,tRes))
            break

        pMin = min(pMin,pSat) ; pMax = max(pMax,pSat)
        tMin = min(tMin,tRes) ; tMax = max(tMax,tRes)

#== Test the Refinement Procedure =====================================        

        nStp += 1

        if nStp > 100 : break

        #if nStp == 5 : refineInPresTemp(alfa,pSat,tRes,Z,lnKref,clsEOS,clsIO)

#-- Store (T,p) for use later ---------------------------------------
        
        if not qBub and alfa > 0.0 :
            DewT.append(tRes)
            DewP.append(pSat)
        else                   :
            BubT.append(tRes)
            BubP.append(pSat)

#-- Estimate dP/dT and da/dT ----------------------------------------

        dTem = tRes - tLst ; dPrs = pSat - pLst ; dalf = alfa - aLst

        dPdT = dPrs/dTem ; dadT = dalf/dTem
        dTdP = dTem/dPrs ; dadP = dalf/dPrs

#-- Switch to iterating in (pRes,alfa) if dP/dT > 10 ----------------            

        if qPrs :
            if abs(dPdT) > 10.0 :
                qPrs = False
                tInc = copysign( 20.0,dTem)
                pInc = copysign(100.0,dPrs)
        else :
            if abs(dTdP) >  0.1 :
                qPrs = True
                tInc = copysign( 20.0,dTem)
                pInc = copysign(100.0,dPrs)

        #degF = tRes - 459.67

        if qDeb :
            if qPrs :
                sOut = "tRes,nItr,pSat,alfa,eD,dPdT,dadT {:8.3f} {:2d} {:10.3f} {:8.5f} {:10.3e} {:10.3e} {:10.3e}\n".format(tRes,nItr,pSat,alfa,eD,dPdT,dadT)
                fDeb.write(sOut)
            else :
                sOut = "tRes,nItr,pSat,alfa,eD,dTdP,dadP {:8.3f} {:2d} {:10.3f} {:8.5f} {:10.3e} {:10.3e} {:10.3e}\n".format(tRes,nItr,pSat,alfa,eD,dTdP,dadP)
                fDeb.write(sOut)

#-- Have we passed the Critical Point? ------------------------------

        if aLst*alfa < 0.0 :
            if qPrs :
                dTCr = - aLst/dadT
                tCri =   tLst + dTCr
                pCri =   pLst + dPdT*dTCr
            else :
                dPcr = - aLst/dadP
                pCri =   pLst + dPcr
                tCri =   tLst + dTdP*dPcr
            #print("Critical Point (P,T) {:10.3f} {:10.5f}".format(pCri,tCri))

#-- Update Temperature for Next Step and Predict (pRes,alfa) --------

        tLst = tRes ; pLst = pSat ; aLst = alfa

        if qPrs :
            tRes = tRes + tInc
            pSat = pSat + tInc*dPdT
            alfa = alfa + tInc*dadT
        else :
            pSat = pSat + pInc
            tRes = tRes + pInc*dTdP
            alfa = alfa + pInc*dadP

#----------------------------------------------------------------------
#  Vapour Fraction Line, if Critical Point Found
#----------------------------------------------------------------------

    pMid = 0.5*(pMin+pMax)
    tMid = 0.5*(tMin+tMax)

    #print("pMin,pMax {:10.3f} {:10.3f}".format(pMin,pMax))
    #print("tMin,tMax {:10.3f} {:10.3f}".format(tMin,tMax))
    #print("pMid,tMid {:10.3f} {:10.3f}".format(pMid,tMid))
    
    #FlsT,FlsP = vapFracLine(pMid,tMid,Z,clsEOS,clsIO,qDeb,fDeb)

    if pCri != None :
        FlsT,FlsP = vapFracLine(pSatR,tResI,Z,clsEOS,clsIO,qDeb,fDeb)

        #michCP(pCri,tCri,Z,clsEOS,clsIO)

#== End of Routine ====================================================

    return BubT,BubP,DewT,DewP,FlsT,FlsP,pCri,tCri

#========================================================================
#  Iterate in P or T and Alpha at a Given T or P
#========================================================================

def iterPorTAlfa(qPrs,pRes,tRes,alfa,Z,lnKref,clsEOS,clsIO) :

    conV = 1.0E-10
    fSSQ = 1.0
    nItr = 0

    if clsIO.Deb["ENV"] > 1 :
        fDeb = clsIO.fDeb
        qDeb = True
    else :
        qDeb = False

#== Iterate ===========================================================    

    while fSSQ > conV :

        nItr += 1

        pNew,tNew,aNew,f1,f2,eD,df1dX,df2dX,df1da,df2da = \
            updatePorTAlfa(qPrs,pRes,tRes,alfa,Z,lnKref,clsEOS,clsIO)

        dPrs = pNew - pRes
        dTem = tNew - tRes
        dAlf = aNew - alfa

        fSSQ = 0.5*sqrt(f1*f1 + f2*f2)

        if qDeb :
            if qPrs : sOut = "nItr,pRes,alfa,f1,f2,fSSQ,df1dX,df2dX,df1da,df2da,dPoT,dAlf {:2d} {:10.3f} {:10.5f} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.5f} {:10.5f}\n". \
                           format(nItr,pRes,alfa,f1,f2,fSSQ,df1dX,df2dX,df1da,df2da,dPrs,dAlf)
            else    : sOut = "nItr,tRes,alfa,f1,f2,fSSQ,df1dX,df2dX,df1da,df2da,dPoT,dAlf {:2d} {:10.3f} {:10.5f} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.5f} {:10.5f}\n". \
                           format(nItr,tRes,alfa,f1,f2,fSSQ,df1dX,df2dX,df1da,df2da,dTem,dAlf)
            fDeb.write(sOut)

        pRes = pNew ; tRes = tNew ; alfa = aNew

        if pRes < 0.0 or tRes < 0.0 or nItr > 10 : break

#== Return Values =====================================================

    return nItr,pRes,tRes,alfa,eD

#========================================================================
#  Update (P or T,Alpha) at given (T or P)
#========================================================================

def updatePorTAlfa(qPrs,pRes,tRes,alfa,Z,lnKref,clsEOS,clsIO) :

    if clsIO.Deb["ENV"] > 0 :
        fDeb = clsIO.fDeb
        qDeb = True
    else :
        qDeb = False

    perT = 1.0E-04
    sqrP = sqrt(perT)

#== (f1,f2) and their Pressure-Derivatives ============================

    if qPrs : iDer =  1
    else    : iDer = -1

    f1,f2,eD,df1dX,df2dX = calcF1F2(pRes,tRes,Z,lnKref,alfa,iDer,clsEOS,clsIO)

#== Perturb alpha and calculate df1da and df2da by difference =========

    if abs(alfa) < sqrP : dPer = perT*UT.signNumber(alfa)
    else                : dPer = perT*alfa
    
    aPer = alfa + dPer

    iDer = 0  #-- Don't need any derivative information from this call

    a1,a2,aP,da1dX,da2dX = calcF1F2(pRes,tRes,Z,lnKref,aPer,iDer,clsEOS,clsIO)

    df1da = (a1 - f1)/dPer
    df2da = (a2 - f2)/dPer

#== Invert (2x2) Matrix and Calculate (P,alfa)-Updates ================

    dPoT,dAlf = solve2x2(f1,f2,df1dX,df1da,df2dX,df2da)

    alfa = alfa + dAlf

    if qPrs : pRes = pRes + dPoT
    else    : tRes = tRes + dPoT

    if dPoT == 0.0 and dAlf == 0.0 : pRes = -1.0

#== Return values =====================================================    

    return pRes,tRes,alfa,f1,f2,eD,df1dX,df2dX,df1da,df2da

#========================================================================
#  (F1,F2) Functions
#========================================================================

def calcF1F2(pRes,tRes,zRef,lnKref,alfa,iDer,clsEOS,clsIO) :

    nCom = clsEOS.nComp
    iLiq = 0
    qXD  = False

#== Calculate the Incipient Phase Composition =========================

    yMol,yNor = incipientY(alfa,zRef,lnKref,clsEOS,clsIO)

#== Fugacity and its Pressure-Derivatives (if required) ===============

    if   iDer ==  1 :
        qPD = True  ; qTD = False
    elif iDer == -1 :
        qPD = False ; qTD = True
    else :
        qPD = False ; qTD = False

    fugZ,dZdP,dZdT,dumM = CE.calcPhaseFugPTX(iLiq,qPD,qTD,qXD,pRes,tRes,zRef,clsEOS)
    fugY,dYdP,dYdT,dumM = CE.calcPhaseFugPTX(iLiq,qPD,qTD,qXD,pRes,tRes,yNor,clsEOS)

#== Calculate Functions ===============================================

    f1 = 0.0 ; f2 = 0.0 ; eD = 0.0

    fR = NP.log(yNor) + fugY - NP.log(zRef) - fugZ  #-- Vector

    f1 = NP.dot(yNor,fR)
    f2 = NP.dot(zRef,fR)

    fR = NP.multiply(fR,fR)
    eD = NP.dot(zRef,fR)
    
    df1dX = 0.0 ; df2dX = 0.0

    if   iDer ==  1 :
        pDiff = dYdP - dZdP         #-- Vector
        df1dX = NP.dot(yNor,pDiff)
        df2dX = NP.dot(zRef,pDiff)
    elif iDer == - 1 :
        tDiff = dYdT - dZdT         #-- Vector
        df1dX = NP.dot(yNor,tDiff)
        df2dX = NP.dot(zRef,tDiff)

#== Return values =====================================================            

    return f1,f2,eD,df1dX,df2dX

#========================================================================
#  Solve 2x2 Matrix
#  |  df1dX  df1dY  | |  dX  |  =  |  -f1  |
#  |  df2dX  df2dY  | |  dY  |  =  |  -f2  |
#========================================================================

def solve2x2(f1,f2,df1dX,df1dY,df2dX,df2dY) :

#-- Determinant -----------------------------------------------------

    detM = df1dX*df2dY - df1dY*df2dX

    if abs(detM) < UT.macEPS :
        dX = 0.0 ; dY = 0.0 ; return dX,dY

#-- Reduce numerical "slop" and get right signs ---------------------

    f1 = - f1/detM
    f2 = - f2/detM

#== (X,Y)-Updates =====================================================

    dX =    df2dY*f1 - df1dY*f2
    dY =  - df2dX*f1 + df1dX*f2

#== Return information ================================================

    return dX,dY

#========================================================================
#  Given an Alpha, Estimate the Incipient Phase Composition
#========================================================================

def incipientY(alfa,zRef,lnKref,clsEOS,clsIO) :

    nCom = clsEOS.nComp

    logZ = NP.log(zRef)
    logK = NP.multiply(lnKref,alfa)

    logY = NP.add(logZ,logK)

    yMol = NP.exp(logY)     #-- Y = exp(logY)

    ySum = NP.sum(yMol)

    yNor = NP.divide(yMol,ySum)

#== Return the Incipient Phase Composition ============================

    return yMol,yNor

#========================================================================
#  Refine (p,T,K) for Saturated Line
#========================================================================

def refineInPresTemp(alfa,pRes,tRes,Z,lnKref,clsEOS,clsIO) :

    nCom = clsEOS.nComp

#-- Current Incipient Composition -----------------------------------

    yApp,yNor = incipientY(alfa,Z,lnKref,clsEOS,clsIO)

#-- Current K-Values ------------------------------------------------
    
    K    = NP.divide(yApp,Z)
    logK = NP.log(K)
    Ze   = NP.multiply(Z,logK)

    fcP2e = NP.dot(Ze,logK)

#----------------------------------------------------------------------
#  Generate (p,T,K) updates
#----------------------------------------------------------------------

    fSSQ = 1.0
    nItr = 0

    while fSSQ > 1.0E-12 :

        fSSQ,pRes,tRes,logK = \
            updatePTK(pRes,tRes,Z,Ze,logK,fcP2e,clsEOS,clsIO)

        nItr += 1

        if nItr > 2 : break

#== Re-Calculate Ref-K's from Eqn.(18) ==============================

    #print("Old: lnKref ",lnKref)

    aRec   = 1.0/alfa
    lnKref = NP.multiply(aRec,logK)

    #print("New: lnKref ",lnKref)

    return

#========================================================================
#  Update (p,T,K) for Saturated Points
#========================================================================

def updatePTK(pRes,tRes,Z,Ze,logK,fcP2e,clsEOS,clsIO) :

    nCom = clsEOS.nComp

    iNeu = 0
    qPD  = True ; qTD  = True ; qXD = False

#-- Trial composition -----------------------------------------------

    K = NP.exp(logK)
    Y = Z*K

#-- fc+1: Equation (11) ---------------------------------------------

    sumY = NP.sum(Y)
    fcP1 = sumY - 1.0

    #print("sumY {:10.5f}".format(sumY))

#-- fc+2: Equation (13) ---------------------------------------------

    fcP2 = - fcP2e + NP.dot(Ze,logK)

#-- Liquid and Vapour Fugacities and their P-Derivatives ------------    

    yNor = Y/sumY   #-- Vector

    fugZ,dZdP,dZdT,dumM = CE.calcPhaseFugPTX(iNeu,qPD,qTD,qXD,pRes,tRes,Z   ,clsEOS)
    fugY,dYdP,dYdT,dumM = CE.calcPhaseFugPTX(iNeu,qPD,qTD,qXD,pRes,tRes,yNor,clsEOS)

#== Construct Jacobian and RHS terms ==================================

    fRes = fugY - fugZ  #-- All Vector Operations
    dfdP = dYdP - dZdP
    dfdT = dYdT - dZdT
    fRes = logK + fRes

    sumf1 = NP.dot(Y,fRes) ; sumf2 = NP.dot(Ze,fRes)
    df1dP = NP.dot(Y,dfdP) ; df2dP = NP.dot(Ze,dfdP)
    df1dT = NP.dot(Y,dfdT) ; df2dT = NP.dot(Ze,dfdT)

#== RHS-Terms =========================================================

    fcP1 = fcP1 - sumf1
    fcP2 = fcP2 - sumf2

#== Invert (2x2) Matrix ===============================================

    dPrs,dTem = solve2x2(fcP1,fcP2,df1dP,df1dT,df2dP,df2dT)

    fSSQ = 0.5*(dPrs*dPrs + dTem*dTem)

    fDeb = clsIO.fDeb

    sOut = "fcP1,fcP2,fSSQ,dPrs,dTem {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}\n".\
          format(fcP1,fcP2,fSSQ,dPrs,dTem)
    fDeb.write(sOut)

#== P,T & logK-Updates ================================================

    pRes = pRes + dPrs
    tRes = tRes + dTem

    for iC in range(nCom) :
        dlnK = - fRes[iC] - dPrs*dfdP[iC] - dTem*dfdT[iC]
        sOut = "iC,lnK,dlnK {:2d} {:10.3e} {:10.3e}\n".format(iC,logK[iC],dlnK)
        fDeb.write(sOut)
        logK[iC] = logK[iC] + dlnK

#== Return information ================================================

    return fSSQ,pRes,tRes,logK

#========================================================================
#  Calculate a Vapour Fraction (Internal) Line
#========================================================================

def vapFracLine(pRes,tRes,Z,clsEOS,clsIO,qDeb,fDeb) :

    nCom = len(Z)

    FlsT = [] ; FlsP = []

    aVal = [0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.55,0.50, \
            0.45,0.40,0.35,0.30,0.25,0.20,0.15,0.10,0.05]

    nVal = len(aVal)

#----------------------------------------------------------------------
#  Conventional Flash to Find p -> beta = 0.5
#----------------------------------------------------------------------

    tRes = 0.5*tRes
    betO = 0.5

    tRes,K = iterateInitFlash(pRes,tRes,Z,betO,clsEOS,clsIO,qDeb,fDeb)

    if tRes < 0.0 : return FlsT,FlsP

    FlsT.append(tRes) ; FlsP.append(pRes)

    lnKref = NP.log(K)

#----------------------------------------------------------------------
#  Decrement alfa (from unity) towards Critical Point
#----------------------------------------------------------------------

    for iVal in range(nVal) :

        alfa = aVal[iVal]

        pRes,tRes = findVapFracPT(alfa,betO,pRes,tRes,Z,lnKref,clsEOS)

        #print("alfa,pRes,tRes {:8.5f} {:10.3f} {:10.3f}".format(alfa,pRes,tRes))

        if pRes < 0.0 or tRes < 0.0 : break

        FlsT.append(tRes) ; FlsP.append(pRes)

#-- Decrement alfa --------------------------------------------------

        alfa = alfa - 0.05

#== Return Information ================================================

    return FlsT,FlsP

#========================================================================
#  For given (alfa,beta), find the (P,T) on constant vapour fraction line
#========================================================================

def findVapFracPT(alfa,betO,pRes,tRes,Z,lnKref,clsEOS) :

    nCom = clsEOS.nComp

    iLiq =  1   ; iVap = -1
    qPD  = True ; qTD  = True ; qXD = False

#----------------------------------------------------------------------
#  Iterate in (P,T)
#----------------------------------------------------------------------

    fSSQ = 1.0

    while fSSQ > 1.0E-12 :

#== Solve Equation (20) to find theta =================================        

        theta,K,c = findTheta(betO,alfa,Z,lnKref)

#-- Liquid and Vapour Mole Fractions --------------------------------

        X,Y = calcLiqVapMoleFrac(betO,Z,K,c)

#== Compute the Liquid & Vapour Fugacities and their (P,T)-Derivatives        
        
        fugX,dXdP,dXdT,dumM = CE.calcPhaseFugPTX(iLiq,qPD,qTD,qXD,pRes,tRes,X,clsEOS)
        fugY,dYdP,dYdT,dumM = CE.calcPhaseFugPTX(iVap,qPD,qTD,qXD,pRes,tRes,Y,clsEOS)

        logX = NP.log(X)
        logY = NP.log(Y)

#== Residuals and their (P,T)-Derivatives =============================        

        dRes = logY + fugY - logX - fugX  #-- Vector
        dPrs =        dYdP        - dXdP  #-- Vector  
        dTem =        dYdT        - dXdT  #-- Vector

        f1    = NP.dot(Y,dRes) ; f2    = NP.dot(X,dRes)
        df1dP = NP.dot(Y,dPrs) ; df2dP = NP.dot(X,dPrs)
        df1dT = NP.dot(Y,dTem) ; df2dT = NP.dot(X,dTem)

#== Solve the 2x2 Matrix ==============================================

        dPrs,dTem = solve2x2(f1,f2,df1dP,df1dT,df2dP,df2dT)

        fSSQ = 0.5*(dPrs*dPrs + dTem*dTem)

        pRes = pRes + dPrs
        tRes = tRes + dTem

        if pRes < 0.0 or tRes < 0.0 : break

#== Return Information ================================================

    return pRes,tRes

#========================================================================
#  Find Temperature at given Pressure that gives required beta
#========================================================================

def iterateInitFlash(pRes,tRes,Z,betO,clsEOS,clsIO,qDeb,fDeb) :

    nCom = len(Z)
    
    iLiq =  1 ; iVap = -1

#----------------------------------------------------------------------
#  Full Flash at Initial Pressure Estimate
#----------------------------------------------------------------------

    betC,K,X,Y = CF.calcFlash(pRes,tRes,Z,betO,clsEOS,clsIO)

#----------------------------------------------------------------------
#  Iterate in Pressure
#----------------------------------------------------------------------

    resB = betC - betO
    nItr = 0

    while abs(resB) > 1.0E-12 :

#== Calculate Liquid & Vapour Mole Fractions ==========================

        if qDeb :
            sOut = "pRes,tRes,betC {:10.3f} {:10.3f} {:10.7f}\n".format(pRes,tRes,betC)
            fDeb.write(sOut)

        if pRes < 0.0 : break

        tRes,K,c = refineInTempOnly(betO,pRes,tRes,Z,K,clsEOS)

        betC = solveRR(betC,Z,c)

        resB = betC - betO

        nItr += 1

        if nItr > 10 : break

#== Return Information ================================================

    return pRes,K

#========================================================================
#  Refine Temperature to Get Required Vapour Fraction (betO)
#========================================================================

def refineInTempOnly(betO,pRes,tRes,Z,K,clsEOS) :

    nCom = clsEOS.nComp

    iLiq =  1
    iVap = -1

#-- c-Vector: c = 1/(1-K[i]) ----------------------------------------

    c    = calcCvalues(K)
    logK = NP.log(K)

#== Rachford-Rice Function G(beta) using required beta ================

    G = 0.0
    for iC in range(nCom) : G = G + Z[iC]/(betO - c[iC])

#-- Liquid and Vapour Mole Numbers and Mole Fractions ---------------        

    X,Y = calcLiqVapMoleFrac(betO,Z,K,c)

    sumX = NP.sum(X)           ; sumY = NP.sum(Y)
    recX = 1.0/sumX            ; recY = 1.0/sumY
    xNor = NP.multiply(recX,X) ; yNor = NP.multiply(recY,Y)

#-- Liquid and Vapour Fugacities and their P-Derivatives ------------

    qPD = False ; qTD = True ; qXD = False

    fugX,dumV,dXdT,dumM = CE.calcPhaseFugPTX(iLiq,qPD,qTD,qXD,pRes,tRes,xNor,clsEOS)
    fugY,dumV,dYdT,dumM = CE.calcPhaseFugPTX(iVap,qPD,qTD,qXD,pRes,tRes,yNor,clsEOS)

#== Form residual and various sums ====================================

    GF = 0.0 ; GFP = 0.0

    fRes = NP.zeros(nCom)
    dfdT = NP.zeros(nCom)

    for iC in range(nCom) :

        fRes[iC] = logK[iC] + fugY[iC] - fugX[iC]
        dfdT[iC] =            dYdT[iC] - dXdT[iC]
        
        xyz  = xNor[iC]*yNor[iC]/Z[iC]

        GF  = GF  + xyz*fRes[iC]
        GFP = GFP + xyz*dfdT[iC]

#== Find delP =========================================================

    Grhs =  G - GF

    if abs(GFP) > UT.macEPS : delT = Grhs/GFP
    else                    : delT = 0.0

#-- If update takes us too far, cut-back ... ------------------------    

    if tRes + delT < 0.0 : delT = 0.5*delT

    #print("tRes,Grhs,GFP,delT {:10.3f} {:10.3e} {:10.3e} {:10.3e}".format(tRes,Grhs,GFP,delT))

    tRes = tRes + delT

#== Calculate d(lnK): Equation (14) ===================================

    for iC in range(nCom) :
        dlnK     = - fRes[iC] - delT*dfdT[iC]
        logK[iC] =   logK[iC] + dlnK

#== Return Information ================================================

    K = NP.exp(logK)
    c = calcCvalues(K)

    return tRes,K,c

#========================================================================
#  Find Theta to solve Equation (20)
#========================================================================

def findTheta(beta,alfa,Z,lnKref) :

    theta = 0.0
    dThe  = 1.0

    while abs(dThe) > 1.0E-12 :
    
#-- Calculate K-Values [Eqn.19] and c-Values ------------------------

        K = calcKvalues(alfa,lnKref,theta)
        c = calcCvalues(K)

#-- Update for theta ------------------------------------------------

        dThe = rachfordRiceTheta(beta,alfa,Z,K,c)

        theta = theta + dThe

        #print("theta,dThe {:10.3e} {:10.3e}".format(theta,dThe))

#== Return Information ================================================

    return theta,K,c

#========================================================================
#  Rachford-Rice Equation and its Theta-Derivative
#========================================================================

def rachfordRiceTheta(beta,alfa,Z,K,c) :

    Nc = len(Z)

    GB = 0.0 ; dGBdH = 0.0
    
    for iC in range(Nc) :
        
        wrkF = beta - c[iC]
        wrkD = c[iC]/wrkF
        
        GB    =  GB   + Z[iC]/wrkF
        dGBdH = dGBdH + Z[iC]*wrkD*wrkD*K[iC]

#== Return Information ================================================

    dTheta = -GB/dGBdH

    return dTheta

#========================================================================
#  (Conventional) Solve Rachford-Rice Equation
#========================================================================

def solveRR(beta,Z,c) :

    funC = 1.0

    while abs(funC) > 1.0E-12 :

        funC,derV = funcDervRR(beta,Z,c)

        dBet = - funC/derV
        beta = beta + dBet

#== Return information ================================================

    return beta

#========================================================================
#  Rachford Rice Function and Derivative
#========================================================================

def funcDervRR(beta,Z,c) :

    nCom = len(Z)

    funC = 0.0 ; derV = 0.0

    for iC in range(nCom) :

        dNom = 1.0/(beta - c[iC])
        ZNom = Z[iC]*dNom

        funC = funC + ZNom
        derV = derV - ZNom*dNom

#== Return information ================================================

    return funC,derV

#========================================================================
#  Calculate Liquid & Vapour Mole Fractions
#========================================================================

def calcLiqVapMoleFrac(beta,Z,K,c) :

    Nc = len(Z)

    X = NP.zeros(Nc)
    Y = NP.zeros(Nc)

    for iC in range(Nc) :
        X[iC] = Z[iC]*c[iC]/(c[iC] - beta)
        Y[iC] = X[iC]*K[iC]

    sumX = NP.sum(X)
    sumY = NP.sum(Y)

    #print("sumX,sumY {:10.7f} {:10.7f}".format(sumX,sumY))

#== Return Information ================================================    

    return X,Y

#========================================================================
#  Calculate c-Values
#========================================================================

def calcCvalues(K) :

    Nc = len(K)

    c = NP.zeros(Nc)

    for iC in range(Nc) :
        wrkI = 1.0 - K[iC]
        if wrkI < 0.0 : c[iC] = 1.0/ wrkI
        else          : c[iC] = 1.0/(wrkI + 1.0E-30)  #-- Protect / zero

#== Return Information ================================================    

    return c

#========================================================================
#  Update K-Values
#========================================================================

def calcKvalues(alfa,lnKref,theta) :

    logK = NP.multiply(alfa,lnKref) + theta

    K = NP.exp(logK)

#== Return Information ================================================    

    return K

#========================================================================
#  Plot Phase Plot (Batch)
#========================================================================

def plotPhasePlot(clsSAM) :

    dR2dF =  459.67

#-- Generate plot in sub-directory GraphOut; doesn't exist - create it    

    graphOut = "GraphOut"

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating Phase Plot in Sub-Directory ",graphOut)

#-- Data for plot ---------------------------------------------------

    sNam = clsSAM.sNam

    titL = "Approx Phase Plot for Sample " + sNam

    xLab = "Temperature/[degF]"
    yLab = "Pressure/[psia]"

    figN = graphOut + "/PhasePlot_" + sNam + ".png"

    nBub = len(clsSAM.BubT)
    nDew = len(clsSAM.DewT)
    nFls = len(clsSAM.FlsT)

    pCri = clsSAM.pCri
    tCri = clsSAM.tCri

    if pCri == None or tCri == None : qCrt = False
    else                            :
        qCrt = True
        tCri = tCri - dR2dF

#== Bubble Point Line =================================================    

    if nBub > 0 :

        bLab = "Bub"
        
        tBub = NP.copy(clsSAM.BubT) - UT.dF2dR      #-- degR -> degF
        pBub = NP.copy(clsSAM.BubP)

        PL.plot(tBub,pBub,label=bLab,color='g')

#== Dew    Point Line =================================================    

    if nDew > 0 :

        dLab = "Dew"
        
        tDew = NP.copy(clsSAM.DewT) - UT.dF2dR      #-- degR -> degF
        pDew = NP.copy(clsSAM.DewP)

        PL.plot(tDew,pDew,label=dLab,color='r')

#== V = 0.5 Point Line ================================================

    if nFls > 0 :

        fLab = "V=0.5"
        
        tFls = NP.copy(clsSAM.FlsT) - UT.dF2dR      #-- degR -> degF
        pFls = NP.copy(clsSAM.FlsP)

        PL.plot(tFls,pFls,label=fLab,color='b')

#== Critical Point ====================================================

    if qCrt :

        tLab = "CP"

        tCrt = [tCri]
        pCrt = [pCri]

        PL.plot(tCrt,pCrt,label=tLab,marker='s',linestyle='None',color='k')

#== Format the image and plot to file =================================

    PL.title(titL)
    PL.xlabel(xLab)
    PL.ylabel(yLab)
    PL.legend()

    PL.savefig(figN)

    PL.close()

#== Return ============================================================    

    return

#========================================================================
#  End of Module
#========================================================================


