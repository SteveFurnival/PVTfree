
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculates 2-Sided & 1-Sided Stability Test (ST)
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

import numpy          as NP
import scipy.optimize as SO

import calcEOS   as CE
import calcReg   as CR
import utilities as UT
import writeOut  as WO

from math import exp

#========================================================================
#  Two-Sided Stability Test
#  SS/GDEM [Opt=0] or BFGS [Opt=1]
#========================================================================

def twoSidedStabTest(qSat,pRes,tRes,Z,clsEOS,clsIO) :

    Opt1 = clsIO.Opt["STAB"]

    if Opt1 == 1 :
        iTyp,bstT,bstK = stab2SidedBFGS(qSat,pRes,tRes,Z,clsEOS,clsIO)
    else         :
        iTyp,bstT,bstK = twoSidedStabCheck(qSat,pRes,tRes,Z,clsEOS,clsIO)

#== Return Values =====================================================

    return iTyp,bstT,bstK

#========================================================================
#  Two-Sided Stability Test [Default SS/GDEM Version]
#========================================================================

def twoSidedStabCheck(qSat,pRes,tRes,Z,clsEOS,clsIO) :

    nCom = clsEOS.nComp
    
    qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
    qT   = False
    qX   = False

    if clsIO.Deb["STAB2"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#-- Feed Log(Fugacity) Coefficients ---------------------------------    

    iLiq = 0
    fugZ,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,Z,clsEOS)

#-- Wilson K-Values -------------------------------------------------    

    wilK = UT.wilsonK(pRes,tRes,clsEOS)        

#-- Trial Compositions ----------------------------------------------

    lnKV = NP.log(wilK)
    lnKL = -lnKV

#-- Assume Feed is Liquid, Try to Split-Off Vapour ------------------

    iVap = -1
    iSTV,sumV,gVap,tolV,triV,lnKV = \
        oneSidedStabCheck(iVap,pRes,tRes,Z,fugZ,lnKV,clsEOS,clsIO)

    qVap = gVap < 0.0 and triV > 1.0E-04

    if qDeb :
        sOut = "Vap-Like Split: pRes,tRes,iSTV,gVap,tolV,triV {:10.3f} {:8.3f} {:3d} {:10.3e} {:10.3e} {:10.3e}\n".format(pRes,tRes,iSTV,gVap,tolV,triV)
        fDeb.write(sOut)
        #KV = NP.exp(lnKV)
        #WO.writeArrayDebug(fDeb,KV,"twoSidedStab: Wilson-KV")        

#-- Assume Feed is Liquid, Try to Split-Off Vapour ------------------

    iLiq =  1
    iSTL,sumL,gLiq,tolL,triL,lnKL = \
        oneSidedStabCheck(iLiq,pRes,tRes,Z,fugZ,lnKL,clsEOS,clsIO)

    qLiq = gLiq < 0.0 and triL > 1.0E-04

    if qDeb :
        sOut = "Liq-Like Split: pRes,tRes,iSTL,gLiq,tolL,triL {:10.3f} {:8.3f} {:3d} {:10.3e} {:10.3e} {:10.3e}\n".format(pRes,tRes,iSTL,gLiq,tolL,triL)
        fDeb.write(sOut)
        #KL = NP.exp(lnKL)
        #WO.writeArrayDebug(fDeb,KL,"twoSidedStab: Wilson-KL")        

#== Process the 1-Sided Tests =========================================

    #print("2Stab: pRes,tRes,gLiq,triL,gVap,triV {:10.3f} {:8.3f} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(pRes,tRes,gLiq,triL,gVap,triV))
    
    if   qLiq and qVap :

        difK = lnKL - lnKV          #-- Vector
        triV = NP.dot(difK,difK)
        
        if   triV > 1.0E-04 :       #-- Are the solutions distinct?
            iTyp = -3            
            if qSat :
                if gLiq < gVap :    #-- Called from Psat  - Liquid Lower-GFE
                    bstK = NP.exp(lnKL)
                    bstT = triL
                    iTyp = -1
                else           :    #-- Called from Psat  - Vapour Lower-GFE
                    bstK = NP.exp(lnKV)
                    bstT = triV
                    iTyp = -2
            else    :               #-- Called from Flash - Take Ratio
                bstK = NP.exp(difK)
                bstT = 0.0
        elif lnKV[nCom-1] < 0.0 :   #-- logK-Vap 'Right' Way-Up
            iTyp = -2
            bstK = NP.exp(lnKV)
            bstT = triV
        else :                      #-- LogK-Liq must be 'Right' Way-Up
            iTyp = -1
            bstK = NP.exp(lnKL)
            bstT = triL
            
    elif qVap :                     #-- Only Vapour-Like Unstable
        iTyp = 2
        bstK = NP.exp(lnKV)
        bstT = triV
    elif qLiq :                     #-- Only Liquid-Like Unstable
        iTyp = 1
        bstK = NP.exp(lnKL)
        bstT = triL
    else :                          #-- Is Stable 1-Phase at this (p,T)
        iTyp = 0
        bstK = NP.zeros(nCom)
        bstT = 0.0

#== Ensure the K-Values are the "Right" Way Up ========================

    if not qSat and bstK[nCom-1] > 1.0 :
        bstK = NP.divide(1.0,bstK)
        #print("Stab2S: Flipped K-values at iTyp,pRes,tRes {:2d} {:10.3f} {:8.3f}".format(iTyp,pRes,tRes))

    iTyp = abs(iTyp)
    
    return iTyp,bstT,bstK

#========================================================================
#  One-Sided Stability Test
#========================================================================

def oneSidedStabCheck(iPhs,pRes,tRes,Z,fugZ,logK,clsEOS,clsIO) :

    mSTB = 101
    nCom = clsEOS.nComp
    
    qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
    qT   = False
    qX   = False
    
    res0 = NP.zeros(nCom)

    iSTB = 0

    if clsIO.Deb["STAB1"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    if qDeb > 0 :
        sOut = "1-Sided SC at pRes,tRes {:10.3f} {:8.3f}\n".format(pRes,tRes)
        K = NP.exp(logK)
        WO.writeArrayDebug(fDeb,K,"1SidedStab: K")
        WO.writeArrayDebug(fDeb,Z,"1SidedStab: Z")
        fDeb.write(sOut)

#== Successive Substitution/GDEM Loop =================================

    while iSTB < mSTB :

        iSTB += 1

#-- Mole Numbers ----------------------------------------------------

        K    = NP.exp(logK)
        yMol = Z*K
        
#-- Mole composition (normalised) -----------------------------------
            
        yNor,sumY = UT.NormSum(yMol)

        fugY,dumV,dumV,dumM = \
            CE.calcPhaseFugPTX(iPhs,qP,qT,qX,pRes,tRes,yNor,clsEOS)

        gStr = 1.0 - sumY

#-- Residuals -------------------------------------------------------

        res1 = res0
        res0 = fugZ - fugY - logK

        tolR =        NP.dot(res0,res0)
        gStr = gStr - NP.dot(yMol,res0)

#-- GDEM Step? ------------------------------------------------------

        if iSTB % UT.mGDEM1 > 0 : eigV = 1.0
        else                    : eigV = UT.GDEM1(res0,res1,clsIO)

        if eigV < 1.0 : eigV = 1.0

#-- Update logK with/without GDEM Acceleration ----------------------

        eVr0 = eigV*res0
        logK = logK + eVr0
        
        triV = NP.dot(logK,logK)

        if qDeb :
            sOut = "Stab1S: iSTB,sumY,gStr,tolR,triV,eigV {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:8.4f}\n" \
               .format(iSTB,sumY,gStr,tolR,triV,eigV)
            fDeb.write(sOut)
            #K = NP.exp(logK)
            #WO.writeArrayDebug(fDeb,K,"1SidedStab: K")

#-- Trivial or Converged? -------------------------------------------

        if   triV < 1.0E-04 :
            break
        elif tolR < 1.0E-12 :
            break

#== Return the Gstar ==================================================

    if qDeb > 0 :
        WO.writeArrayDebug(fDeb,K,"End 1S-SC: K")
        
    return iSTB,sumY,gStr,tolR,triV,logK

#========================================================================
#  Two-Sided Stability Test: BFGS Version
#========================================================================

def stab2SidedBFGS(qSat,pRes,tRes,Z,clsEOS,clsIO) :

    if clsIO.Deb["STAB2"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.nComp

    qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
    qT   = False
    qX   = False

    iNeu = 0 ; iLiq = 1 ; iVap = -1

#-- Feed Log(Fugacity) Coefficients ---------------------------------

    fugZ,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iNeu,qP,qT,qX,pRes,tRes,Z,clsEOS)

    hVec = NP.log(Z) + fugZ  #-- Vector!
    
#-- Wilson K-Values -------------------------------------------------    

    wilK = UT.wilsonK(pRes,tRes,clsEOS)        

#-- Trial Compositions ----------------------------------------------

    yVap = Z*wilK   #-- Vector: Vapour-Like Trial
    yLiq = Z/wilK   #-- Vector: Liquid-Like Trial

#== Create (lower,upper) bounds for the alpha's =====================    

    bndA = []
    for iC in range(nCom) : bndA.append((1.0E-20,None))

#== Test SciPy Minimize Routine =====================================

    aVap = 2.0*NP.sqrt(yVap)

    qConV,yVap = stabDriv(iVap,aVap,bndA,pRes,tRes,Z,hVec,clsEOS)

    aLiq = 2.0*NP.sqrt(yLiq)

    qConL,yLiq = stabDriv(iLiq,aLiq,bndA,pRes,tRes,Z,hVec,clsEOS)

    #print("qConL,yLiq ",qConL,yLiq)
    #print("qConV,yVap ",qConL,yVap)
    
#== Process the 1-Sided Tests =========================================

    if   qConL and qConV :

#-- If both trials unstable, take ratio unless they are the same ----
        
        iTyp = 3
        bstK = yVap/yLiq

        qSam = True
        for iC in range(nCom) :
            if abs(bstK[iC] - 1.0) > 1.0E-03 :
                qSam = False
                break

        if qSam :
            #print("iTyp=3 and All K = 1")
            if gStrV < gStrL : bstK = yVap/Z        #-- Vector
            else             : bstK =    Z/yLiq     #-- Vector
        
    elif qConV :
        iTyp = 2
        bstK = yVap/Z
    elif qConL :
        iTyp = 1
        bstK =    Z/yLiq
    else :
        iTyp = 0
        bstK = NP.ones(nCom)

#== Ensure the K-Values are the "Right" Way Up ========================

    if bstK[nCom-1] > 1.0 : bstK = 1.0/bstK     #-- Vector

    if qDeb : fDeb.write("Stab2: iTyp " + str(iTyp) + "\n")

    bstT = 0.0

    return iTyp,bstT,bstK

#========================================================================
#  Driver for SciPi.Optimize.Minimize Routine
#  Working with alfa = 2.0*sqrt(Y)
#  function returns g* and dg*/dalfa [jac=True]
#  Using (limited memory) BFGS with limits [bounds]
#========================================================================

def stabDriv(iPhs,alfa,bndA,pRes,tRes,Z,h,clsEOS) :

    if iPhs > 0 : sPhs = "Liq"
    else        : sPhs = "Vap"

#== Initialise the Minimizer Class (and set Trivialty Flag) ===========

    minimizer = Minimizer()

#== Call SciPy routine ================================================

    alfa = minimizer.minimize(iPhs,alfa,bndA,tRes,pRes,Z,h,clsEOS)

#== Process results to return =========================================

    if alfa[0] == 0.0 :
        qCon = False
        yMol = alfa
    else              :
        qCon = True
        yMol = 0.25*alfa*alfa

#== Return values =====================================================
    
    return qCon,yMol

#========================================================================
#  Based on SciPy Minimize TookTooLong
#========================================================================

class Minimizer :

    def __init__(self) : self.isTriv = False

    def minimize(self,iPhs,alfa,bndA,tRes,pRes,Z,h,clsEOS) :

        try :

            res = SO.minimize(fun      = self.stabMini,
                              x0       = alfa,
                              jac      = True,
                              method   = 'L-BFGS-B',
                              bounds   = bndA,
                              callback = self.callback,
                              args     = (pRes,tRes,Z,h,clsEOS))

            nIter = res.nit
            alfa  = res.x

        except :

            #print("Raised Exception")
            nIter = -1
            res   = None
            alfa  = NP.zeros(clsEOS.nComp)

        #print("iPhas,qTriv,nIter ",iPhs,self.isTriv,nIter)
        
        return alfa

    def callback(self,x) :
        if self.isTriv :
            raise Exception("Is Trivial")
            return True

#========================================================================
#  Function for SciPi.Optimize.Minimize Routine
#========================================================================

    def stabMini(self,alfa,*args) :

        iLiq = 0

        qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
        qT   = False
        qX   = False

#== Unpack the arguments in *args =====================================

        pRes   = args[0]
        tRes   = args[1]
        Z      = args[2]
        h      = args[3]
        clsEOS = args[4]

#== Un-Normalised & Normalised Mole Fractions =========================    

        Y = 0.25*alfa*alfa

        y,SY = UT.NormSum(Y)

#-- log(fugacity) coefficients of Trial phase -----------------------    

        fugY,dumV,dumV,dumM = \
            CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,y,clsEOS)

#== Gstar and Gradient-Vector =========================================

        logY = NP.log(Y)                    #-- Vector
        dGdY = logY + fugY - h              #-- Vector

        gStr = 1.0 - SY + NP.dot(Y,dGdY)    #-- Michelsen Eqn.(14)

        dGda = 0.5*alfa*dGdY                #-- Michelsen Eqn.(25)

#== Michelsen Approach to Triviality: Eqn.(23) & (24) =================

        YmZ  = Y - Z
        beta = NP.dot(YmZ,dGdY)
        r    = 2.0*gStr/beta

        if gStr < 1.0E-3 and abs(r-1.0) < 0.2 : self.isTriv = True

        #print("gStr,r,Triv {:10.3e} {:10.3e} {:d}".format(gStr,r,self.isTriv))

#== Return function and its derivative (as tuple) =====================

        return gStr,dGda

#========================================================================
#  End of Module
#========================================================================
