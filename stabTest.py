
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

import numpy    as NP

import calcEOS   as CE
import constants as CO
import utilities as UT
import writeOut  as WO

from math import exp

#========================================================================
#  Two-Sided Stability Test
#========================================================================

def twoSidedStabCheck(qSat,pRes,tRes,Z,clsEOS,clsIO) :

    nCom = clsEOS.NC
    
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

    wilK = wilsonK(pRes,tRes,clsEOS)        

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
                else           :    #-- Called from Psat  - Vapour Lower-GFE
                    bstK = NP.exp(lnKV)
                    bstT = triV
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
    nCom = clsEOS.NC
    
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

        if iSTB % CO.mGDEM1 > 0 : eigV = 1.0
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
#  Wilson K-Values
#========================================================================

def wilsonK(pRes,tRes,clsEOS) :

    nCom = clsEOS.NC
    wilK = NP.zeros(nCom)

    for iC in range(nCom) :

        rPred = clsEOS.gPP("PC",iC)/pRes  #-- Reciprocal Reduced Pressure
        rTred = clsEOS.gPP("TC",iC)/tRes  #-- Reciprocal Reduced Temperature

        wilK[iC] = rPred*exp(5.3727*(1.0 + clsEOS.gPP("AF",iC))*(1.0 - rTred))

    return wilK

#========================================================================
#  End of Module
#========================================================================
