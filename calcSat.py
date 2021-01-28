
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculates Saturation Pressure & Temperature
#  
#========================================================================

#========================================================================
#  The majority of the routines here follow these references:
#  1. Whitson & Brule (W&B), "Phase Behavior", SPE Monograph, (2000)
#  2. Soreide, PhD Thesis, "Improved Phase Behavior Predictions of
#                           Petroleum Reservoir Fluids from a Cubic
#                           Equation of State",
#     NTNU, Trondheim, Norway, April 1989
#     http://www.ipt.ntnu.no/~curtis/courses/Theses/Ingolf-Soreide-NTH-PhD-Thesis.pdf
#========================================================================

#!/usr/bin/python3

import numpy      as NP

import calcEOS    as CE
import calcFlash  as CF
import calcReg    as CR
import calcStab   as CS
import readExps   as RX
import utilities  as UT
import writeOut   as WO

from math  import acos,copysign,cos,exp,log,sqrt

#========================================================================
#  Saturation Pressure Calculation
#========================================================================

def calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO) :

    nCom = clsEOS.nComp
    mPSA = 101

    if clsIO.Deb["PSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False
    
#== Pre-Sweep; do we have approx phase env and/or observed Psat? ======

    if qBub == None :
        qBub,p2PH,p1PH,logK = sweepPsat(pObs,tRes,clsEOS,clsSAM,clsIO)
    else :
        p1PH = 1.01*pSat
        p2PH = 0.99*pSat
        pSat =      p2PH

    pSat = p2PH
    p2OK = p2PH
    
#========================================================================
#  Do we have a valid Pressure Interval to search?
#========================================================================

    if qDeb :
        sOut = "calcPsat: tRes,p2PH,p1PH {:8.3f} {:10.3f} {:10.3f}\n".format(tRes,p2PH,p1PH)
        fDeb.write(sOut)
        K = NP.exp(logK)
        WO.writeArrayDebug(fDeb,K,"Post-sweepPsat: K")

    pDIF = p1PH - p2PH

    if p2PH >= p1PH :
        print("calcSat: can't find valid search range - Error")
        Ksat = NP.zeros(nCom)
        pSat = -1.0
        return qBub,pSat,Ksat

#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)
    
#== Test the logK-values have the right sign! =========================
    
    logL = -1.0*logK

    tol21,yMol,res0,dFdP = calcFugDerv(qBub,pSat,tRes,Z,logK,clsEOS)
    tol22,yMol,res0,dFdP = calcFugDerv(qBub,pSat,tRes,Z,logL,clsEOS)

    if qDeb :
        sOut = "Initial tol21,tol22 {:10.3e} {:10.3e}\n".format(tol21,tol22)
        fDeb.write(sOut)

    if tol21 > tol22 :
        logK = logL
        K    = NP.exp(logK)
        if qDeb : WO.writeArrayDebug(fDeb,K,"Post-tol2: K")
    
#----------------------------------------------------------------------
#  Main Loop
#----------------------------------------------------------------------

    iPSA = 0

    while iPSA < mPSA :

        iPSA += 1

#-- Get Fugacity Coeffs and their Pressure-Derivatives --------------

        if iPSA > 1 : res1 = res0  #-- Copy 'last' residual for GDEM

        tol2,yMol,res0,dFdP = calcFugDerv(qBub,pSat,tRes,Z,logK,clsEOS)

        resX = NP.exp(res0)

#== Is this a GDEM Step?  If so, calculate Single Eigenvalue ==========            

        if iPSA % UT.mGDEM1 > 0 : eigV = 1.0
        else                    : eigV = UT.GDEM1(res0,res1,clsIO)

#== Perform Update; Calculate Various Sums ============================

        yMol = yMol*resX                #-- W&B Eqn.(4.83) with lambda=1
        dQdp = NP.dot(yMol,dFdP)        #-- W&B Eqn.(4.86)
        qFun = NP.sum(yMol)             #-- W&B Eqn.(4.87a) 
        wrk0 = eigV*res0                #-- W&B Eqn.(4.83) with any lambda
        logK = logK + wrk0

#-- Test for Trivial Solution [all K's -> 1] ------------------------            

        triV = NP.dot(logK,logK)        #-- W&B Eqn.(4.88)

#== SS or GDEM step?  If GDEM, re-compute Yi afresh ===================        
        
        if eigV == 1.0 :
            sumY = qFun
        else :
            K    = NP.exp(logK)
            yMol = Z*K
            sumY = NP.sum(yMol)
            if qDeb :
                WO.writeArrayDebug(fDeb,K,"Post-GDEM: K")

        qFun =     1.0 - qFun
        tol1 = abs(1.0 - sumY)          #-- W&B Eqn.(4.87a)

#== Pressure Update (Newton) ==========================================

        if abs(dQdp) < UT.macEPS :
            print("calcSat: Can't update Psat [dQdp = 0] - Error")
            Ksat = NP.zeros(nCom)
            pSat = -1.0
            return qBub,pSat,Ksat

        dPrs = - qFun/dQdp

        pOld = pSat
        pSat = pSat + dPrs              #-- W&B Eqn.(4.85)

        if qDeb :
            sOut = "iPSA,qFun,dQdp,to11,tol2,triV,eigV,dPrs,pSat {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.4f}\n"\
                .format(iPSA,qFun,dQdp,tol1,tol2,triV,eigV,dPrs,pSat)
            fDeb.write(sOut)

        if pSat > p1PH or pSat < p2PH :
            #print("pSat>p1PH: p2PH,p1PH,pSat {:10.3f} {:10.3f} {:10.3f}".format(p2PH,p1PH,pSat))
            pSat,logK = rescuePsat(qBub,p2PH,p1PH,tRes,Z,clsEOS,clsIO)
            break

#== Converged? ========================================================

        if tol1 < 1.0E-12 and tol2 < 1.0E-08 :      #-- Hard to acheive 1E-13
            qCon = True
            break

#== Trivial Solution? =================================================

        if triV < 1.0E-04 :
            qTrv = True
            #print("calcPsat: Trivial Solution?")
            break

#== Pressure-Step too big? ============================================

        """
        if   pSat > p1PH :
            #pSat = pOld + 0.001
            if qDeb :
                sOut = "p2PH,p1PH,pSat {:10.3f} {:10.3f} {:10.3f}\n".format(p2PH,p1PH,pSat)
                fDeb.write(sOut)
        elif pSat < p2PH :
            pSat = pOld + 0.5*pDIF
        else :
            p2PH = pSat
            pDIF = p1PH - p2PH
        """

#========================================================================
#  Converged: Calculate K-Values
#========================================================================

    Ksat = NP.exp(logK)

#== Ensure the K-Values are the 'Right-Way-Up' ========================

    if Ksat[nCom-1] > 1.0 :
        qBub = False
        #qBub = not qBub
        Ksat = NP.divide(1.0,Ksat)
        if qDeb :
            sOut = "End-calcPsat: K-Values Flipped\n"
            fDeb.write(sOut)
    else :
        qBub = True
        if qDeb :
            sOut = "End-calcPsat: K-Values OK\n"
            fDeb.write(sOut)

    return qBub,pSat,Ksat

#========================================================================
#  Calculate Feed & Incipient Fugacities and their P-Derivatives
#========================================================================

def calcFugDerv(qBub,pSat,tRes,Z,logK,clsEOS) :

    iPhs = 0  #-- Force cubic solver to pick root with lowest GFE

    qP = True       #-- Don't Need Temp or Comp Derivatives, Just Pres
    qT = False
    qX = False
    
#-- Incipient Phase -------------------------------------------------        

    K    = NP.exp(logK)
    yMol = Z*K              #-- Vector!

#-- Get Fugacity Coeffs and their Pressure-Derivatives --------------

    yNor = UT.Norm(yMol)
        
    fugZ,dZdP,dumV,dumM = CE.calcPhaseFugPTX(iPhs,qP,qT,qX,pSat,tRes,Z   ,clsEOS)
    fugY,dYdP,dumV,dumM = CE.calcPhaseFugPTX(iPhs,qP,qT,qX,pSat,tRes,yNor,clsEOS)

#== Residual and Derivative Difference ================================

    res0 = fugZ - fugY - logK   #-- Vector!
    dFdP = dYdP - dZdP          #-- Vector!

#== Calculate tol2 =====================================================

    tol2 = calcTOL2(logK,res0)

#== Return values =====================================================

    return tol2,yMol,res0,dFdP

#========================================================================
#  Calculate tol2: W&B Eqn.(4.87) but note square inside the summation!
#========================================================================

def calcTOL2(logK,res0) :

    work = NP.absolute(logK)        #-- Absolute value (as will square!)
    work = NP.maximum(work,1.0E-03) #-- Protect against |logK| ~ 0
    work = NP.divide(res0,work)
    tol2 = NP.dot(work,work)

#== Return values =====================================================

    return tol2

#========================================================================
#  Initialisation for Psat Calculation
#  If Approx Phase Enevlope has been Calculated, use info
#  Else Sweep using Stab-Check from Low-Pressure Up until bracket Psat 
#========================================================================

def sweepPsat(pObs,tRes,clsEOS,clsSAM,clsIO) :

    nCom = clsEOS.nComp
    facT = 0.05         #-- Might want to make this User-Adjustable?

    if clsIO.Deb["PSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False
    
#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)

#== Do we have an Observed Value? =====================================

    if pObs > 0.0 : qPobs = True
    else          : qPobs = False

    if qDeb :
        sOut = "sweepPsat: pObs {:10.3f}\n".format(pObs)
        fDeb.write(sOut)

#== Have we got Approx Phase Envelope? ================================    

    qPhas = clsSAM.qPhas

    if qPhas :
        
        tCri = clsSAM.tCri      #-- Crtical Temp
        pBar = clsSAM.pBar      #-- Cricondenbar   (max-Pres)
        tHrm = clsSAM.tHrm      #-- Cricondentherm (max-Temp)
        
        if tCri == None : tCri = - 999.0  #-- No Cricondentherm => Is Gas

        if tRes < tCri : qBub = True
        else           : qBub = False

#-- Interpolate to find approx pSat at tRes -------------------------

        nBub = len(clsSAM.BubT)
        nDew = len(clsSAM.DewT)

        if       qBub and nBub > 0 : pInt = NP.interp(tRes,clsSAM.BubT,clsSAM.BubP)
        elif not qBub and nDew > 0 : pInt = NP.interp(tRes,clsSAM.DewT,clsSAM.DewP)
        else :
            pInt  = -1.0
            qPhas = False

        if qDeb :
            sOut = "tCri,pBar,tHrm,pInt {:10.4f} {:10.3f} {:10.4f} {:10.3f}\n".format(tCri,pBar,tHrm,pInt)
            fDeb.write(sOut)

    else : qBub = False        

#== Set Limits Depending on What Data We Have Available ===============

    qDat = True

    if   qPhas :
        p1PH = (1.0+facT)*pInt ; p2PH = (1.0-facT)*pInt ; pInc = p1PH - p2PH
    elif qPobs :
        p1PH = (1.0+facT)*pObs ; p2PH = (1.0-facT)*pObs ; pInc = p1PH - p2PH
    else :
        qDat = False

#----------------------------------------------------------------------
#  Data Available, Check They Appropriate Bounds
#----------------------------------------------------------------------

    mC7P = UT.moleFracC7P(Z,clsEOS)    #-- z(C7+)

    if qDat :

        if qDeb :
            sOut = "sweepPsat: Data available at tRes,mC7P {:8.3f} {:8.5f}\n".format(tRes,mC7P)
            fDeb.write(sOut)

        #print("p1PH,p2PH,pInc {:10.3f} {:10.3f} {:10.3f}".format(p1PH,p2PH,pInc))

        p1PH           = searchOnePhase(qBub,pInc,p1PH,tRes,Z,clsEOS,clsIO)
        iTyp,p2PH,logK = searchTwoPhase(qBub,pInc,p2PH,tRes,Z,clsEOS,clsIO)
        
    else :

#----------------------------------------------------------------------
#  No information to guide us, just sweep through range of P's
#----------------------------------------------------------------------

        if qDeb :
            sOut = "sweepPsat: No data available at tRes,mC7P {:8.3f} {:8.5f}\n".format(tRes,mC7P)
            fDeb.write(sOut)

        if mC7P > 0.125 : qBub = True   #-- > 12.5 mol% - Likely Volatile Oil
        else            : qBub = False  #-- Else Condensate

        iTyp,p2PH,p1PH,logK = boundPsat(tRes,Z,clsEOS,clsIO)

#== Return Information =================================================

    if qDeb :
        sOut = "sweepPsat: qBub,iTyp " + str(qBub) + " " + str(iTyp) + "\n"
        fDeb.write(sOut)

    return qBub,p2PH,p1PH,logK

#========================================================================
#  Search (upward in pressure if need be) to find 1-Phase state
#========================================================================

def searchOnePhase(qBub,pInc,pRes,tRes,Z,clsEOS,clsIO) :

    q1PH = False
    qSat = True

#----------------------------------------------------------------------
#  Increment upward until we find 1-Phase state
#----------------------------------------------------------------------

    while not q1PH :

        iTyp,triV,K = CS.twoSidedStabTest(qSat,pRes,tRes,Z,clsEOS,clsIO)

        if iTyp == 0 : break
        else         : pRes = pRes + pInc

#== Return value ======================================================

    return pRes

#========================================================================
#  Search (downward in pressure if need be) to find 2-Phase state
#========================================================================

def searchTwoPhase(qBub,pInc,pRes,tRes,Z,clsEOS,clsIO) :

    iTyp = 0
    qSat = True

#----------------------------------------------------------------------
#  Increment downward (in pressure) until we find 2-Phase state
#----------------------------------------------------------------------

    while iTyp == 0 :

        iTyp,triV,K = CS.twoSidedStabTest(qSat,pRes,tRes,Z,clsEOS,clsIO)

        if iTyp > 0 : break
        else        : pRes = pRes - pInc

#== Return Information ================================================

    logK = NP.log(K)

    return iTyp,pRes,logK

#========================================================================
#  Main Psat Calculation Failed even though we have bounded result
#  Refine the search space
#========================================================================

def rescuePsat(qBub,p2PH,p1PH,tRes,Z,clsEOS,clsIO) :

    if clsIO.Deb["PSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#== Using interval halfing ============================================

    pSat = 0.5*(p2PH+p1PH)
    qSat = True

    while abs(p1PH-p2PH) > 0.001 :

        iTyp,triV,K = CS.twoSidedStabTest(qSat,pSat,tRes,Z,clsEOS,clsIO)

        if iTyp == 0 :
            p1PH = pSat
        else         :
            p2PH = pSat
            logK = NP.log(K)

        pSat = 0.5*(p2PH+p1PH)

#== Return arguments ==================================================    
        
    return p2PH,logK

#========================================================================
#  Bound Psat using Stability Test
#========================================================================

def boundPsat(tRes,Z,clsEOS,clsIO) :

    if clsIO.Deb["PSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#== Limiting Values ===================================================

    qSat = True

    pMax = 15010.0
    pMin =    10.0

    iMax,trMx,Kmax = CS.twoSidedStabTest(qSat,pMax,tRes,Z,clsEOS,clsIO)
    iMin,trMn,Kmin = CS.twoSidedStabTest(qSat,pMin,tRes,Z,clsEOS,clsIO)

    n2SC = 2

    #print("boundPsat: iMax,iMin,trMn {:2d} {:2d} {:10.3e}".format(iMax,iMin,trMn))

#-- Test limiting values --------------------------------------------

    if iMax > 0 :
        print("Fluid at (Pres,Tres) = ({:10.3f},{:8.3f}) is 2-Phase: No-Psat".format(pMax,tRes))
        iTyp = -1
        logK = NP.zeros(clsEOS.nComp)
        return iTyp,pMin,pMax,logK
      
#----------------------------------------------------------------------
#  Interval Halving
#----------------------------------------------------------------------

    pInt = 0.5*(pMin + pMax)
    qLog = False

    while abs(pMax - pMin) > 100.0 :

        iTyp,triV,K = CS.twoSidedStabTest(qSat,pInt,tRes,Z,clsEOS,clsIO)

        n2SC += 1

        #print("iTyp,pMin,pInt,triV,pMax {:2d} {:10.3f} {:10.3f} {:10.3e} {:10.3f}".format(iTyp,pMin,pInt,triV,pMax))

        if iTyp > 0 :
            pMin = pInt
            logK = NP.log(K)
            qLog = True
        else        :
            pMax = pInt

        pInt = 0.5*(pMin + pMax)

#-- Ensure logK array is defined ------------------------------------        

    if not qLog : logK = NP.zeros(clsEOS.nComp)

#== Return values =====================================================

    return iTyp,pMin,pMax,logK

#========================================================================
#  Saturation Pressure Calculation
#========================================================================

def calcTsat(pRes,clsEOS,clsSAM,clsIO) :

    nCom = clsEOS.nComp
    mTSA = 101

    if clsIO.Deb["TSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#-- Load Feed Composition -------------------------------------------    

    Z = NP.zeros(nCom)
    
    for iC in range(nCom) : Z[iC] = clsSAM.gZI(iC)
    
#== Pre-Sweep; 2-Phase/1-Phase Temperature Bounds =====================

    t2PH,t1PH,logK = boundTsat(pRes,Z,clsEOS,clsIO)

    tRes = t2PH     #-- This will be our iteration parameter

    if qDeb :
        sOut = "calcTsat: pRes,t2PH,t1PH {:8.3f} {:10.3f} {:10.3f}\n".format(pRes,t2PH,t1PH)
        fDeb.write(sOut)
        K = NP.exp(logK)
        WO.writeArrayDebug(fDeb,K,"Post-boundTsat: K")

#----------------------------------------------------------------------
#  Main Iterative Loop
#----------------------------------------------------------------------

    iTSA = 0
    iPhs = 0
    qPD  = False
    qTD  = True
    qXD  = False

    K    = NP.exp(logK)
    yMol = Z*K

    while iTSA < mTSA :

        iTSA += 1

        if iTSA > 1 : res1 = res0  #-- Copy 'last' residual for GDEM

#== Log Fug Coefs and their T-derivatives for Feed and Trial Phase ====

        fugZ,dZdP,dZdT,dumM = CE.calcPhaseFugPTX(iPhs,qPD,qTD,qXD,pRes,tRes,Z,clsEOS)

#-- Get Normalised Trial Composition --------------------------------

        yNor = UT.Norm(yMol)

        fugY,dYdP,dYdT,dumM = CE.calcPhaseFugPTX(iPhs,qPD,qTD,qXD,pRes,tRes,yNor,clsEOS)

#== Residual and Derivative Difference ================================

        res0 = fugZ - fugY - logK   #-- Vector
        dFdT = dYdT - dZdT          #-- Vector

#== Calculate tol2 =====================================================

        tol2 = calcTOL2(logK,res0)

        resX = NP.exp(res0)

#== Is this a GDEM Step?  If so, calculate Single Eigenvalue ==========            

        if iTSA % UT.mGDEM1 > 0 : eigV = 1.0
        else                    : eigV = UT.GDEM1(res0,res1,clsIO)

#== Perform Update; Calculate Various Sums ============================

        yMol = yMol*resX                #-- W&B Eqn.(4.83) with lambda=1
        dQdT = NP.dot(yMol,dFdT)        #-- W&B Eqn.(4.86)
        qFun = NP.sum(yMol)             #-- W&B Eqn.(4.87a) 
        wrk0 = eigV*res0                #-- W&B Eqn.(4.83) with any lambda
        logK = logK + wrk0

#-- Test for Trivial Solution [all K's -> 1] ------------------------            

        triV = NP.dot(logK,logK)        #-- W&B Eqn.(4.88)

#== SS or GDEM step?  If GDEM, re-compute Yi afresh ===================        
        
        if eigV == 1.0 :
            sumY = qFun
        else :
            K    = NP.exp(logK)
            yMol = Z*K
            sumY = NP.sum(yMol)
            if qDeb :
                WO.writeArrayDebug(fDeb,K,"Post-GDEM: K")

        qFun =     1.0 - qFun
        tol1 = abs(1.0 - sumY)          #-- W&B Eqn.(4.87a)

#== Pressure Update (Newton) ==========================================

        dTem = - qFun/dQdT

        tOld = tRes
        tRes = tRes + dTem              #-- W&B Eqn.(4.85)

        if qDeb :
            sOut = "iTSA,qFun,dQdT,to11,tol2,triV,eigV,dTem,tRes {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.4f}\n"\
                .format(iTSA,qFun,dQdT,tol1,tol2,triV,eigV,dTem,tRes)
            fDeb.write(sOut)

#== Converged? ========================================================

        if tol1 < 1.0E-12 and tol2 < 1.0E-08 :      #-- Hard to acheive 1E-13
            qCon = True
            break

#== Trivial Solution? =================================================

        if triV < 1.0E-04 :
            qTrv = True
            print("calcPsat: Trivial Solution?")
            break

#========================================================================
#  Converged: Calculate K-Values
#========================================================================

    Ksat = NP.exp(logK)

#== Ensure the K-Values are the 'Right-Way-Up' ========================

    if Ksat[nCom-1] > 1.0 :
        qBub = False
        #qBub = not qBub
        #Ksat = NP.divide(1.0,Ksat)
        if qDeb :
            sOut = "End-calcTsat: K-Values NOT Flipped\n"
            fDeb.write(sOut)
    else :
        qBub = True
        if qDeb :
            sOut = "End-calcTsat: K-Values OK\n"
            fDeb.write(sOut)

#== Return values =====================================================

    return qBub,tRes,Ksat

#========================================================================
#  Bound Tsat using Stability Test
#========================================================================

def boundTsat(pRes,Z,clsEOS,clsIO) :

    if clsIO.Deb["TSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#== Limiting Values ===================================================

    qSat = True

    tMax = 2500.0
    tMin =  200.0

    iMax,trMx,Kmax = CS.twoSidedStabTest(qSat,pRes,tMax,Z,clsEOS,clsIO)
    iMin,trMn,Kmin = CS.twoSidedStabTest(qSat,pRes,tMin,Z,clsEOS,clsIO)

    n2SC = 2

    if qDeb :
        sOut = "boundTsat: iMax,iMin,trMn {:2d} {:2d} {:10.3e}\n".format(iMax,iMin,trMn)
        fDeb.write(sOut)

#-- Test limiting values --------------------------------------------

    if iMax > 0 :
        print("Fluid at (Pres,Tres) = ({:10.3f},{:8.3f}) is 2-Phase: No-Psat".format(pRes,tMax))
        iTyp = -1
        logK = NP.zeros(clsEOS.nComp)
        return tMin,tMax,logK
      
#----------------------------------------------------------------------
#  Interval Halving
#----------------------------------------------------------------------

    tInt = 0.5*(tMin + tMax)

    while abs(tMax - tMin) > 100.0 :

        iTyp,triV,K = CS.twoSidedStabTest(qSat,pRes,tInt,Z,clsEOS,clsIO)

        n2SC += 1

        if qDeb :
            sOut = "iTyp,pRes,tMin,tInt,triV,tMax {:2d} {:10.3f} {:10.3f} {:10.3f} {:10.3e}{:10.3f}\n".format(iTyp,pRes,tMin,tInt,triV,tMax)
            fDeb.write(sOut)

        if iTyp > 0 :
            tMin = tInt
            logK = NP.log(K)
        else        :
            tMax = tInt

        tInt = 0.5*(tMin + tMax)

#== Return values =====================================================

    return tMin,tMax,logK

#========================================================================
#  Bound Psat using Stability Test with No Constraining Data
#  Just search low to high until we have p2PH < pSat < p1PH
#========================================================================

def sweepPsatNoData(qBub,tRes,Z,clsEOS,clsIO) :

    iTyp,p2PH,p1PH,logK = boundPsat(tRes,Z,clsEOS,clsIO)

    if clsIO.Deb["PSAT"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    logK = NP.zeros(clsEOS.nComp)

    if qBub : iFed =  1     #-- Feed is Likely a Liquid
    else    : iFed = -1     #             Else a Vapour
        
#== Pressure Limits and # Increments ==================================

    pMIN =    10.0
    pMAX = 15010.0
    nINC = 15

#== set pSat = pMIN and calculate Increment ===========================

    pSat =  pMIN
    pINC = (pMAX - pSat)/float(nINC)
    p1PH = -1.0

#======================================================================
#  Increment in pressure until we find 1-Phase State
#======================================================================

    qSat = True

    for iStep in range(nINC+1) :

#== Can we split of another composition at this pressure? =============
        
        iThs,tThs,Kths = CS.twoSidedStabTest(qSat,pSat,tRes,Z,clsEOS,clsIO)

        if iThs > 0 : i2PH = 2
        else        : i2PH = 1

        if qDeb :
            sOut = "pSat,tRes,i2PH {:10.3f} {:8.3f} {:1d}\n".format(pSat,tRes,i2PH)
            fDeb.write(sOut)
   
        if iThs > 0 :
            iTyp = iThs
            p2PH = pSat
            logK = NP.log(Kths)
        else        :
            p1PH = pSat
            break
            
#== Increment pSat estimator and repeat ===============================
                
        pSat = pSat + pINC

#== Return Information ================================================    

    return iTyp,p2PH,p1PH,logK

#========================================================================
#  End of Module
#========================================================================


