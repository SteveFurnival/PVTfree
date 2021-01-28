#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Water Properties Module (CM)
#  
#========================================================================

#!/usr/bin/python3

from   math  import exp

import utilities as UT

#========================================================================
#  E100 Water Properties (PVTW Keyword)
#========================================================================

def calcPVTW(clsSIM,clsUNI,clsIO) :

#-- Brine Density at Surface Conditions -----------------------------    

    pSur = UT.pStand
    tSur = UT.tStand
    mFrc = clsSIM.bSalt

    dSTW,dumW = calcRoweChouDen(mFrc,tSur,pSur,clsUNI)

    clsSIM.dSTW = dSTW

#-- Brine Density at Reference Conditions ---------------------------    

    pRef = clsSIM.pRefW
    tRes = clsSIM.Tres
    mFrc = clsSIM.bSalt

    dWat,comW = calcRoweChouDen(mFrc,tRes,pRef,clsUNI)

#-- Bw(Ref) = Den(STC)/Den(Ref) = V(Ref)/V(STC) ---------------------

    BwRf = clsSIM.dSTW/dWat

#-- Brine Viscosity and Viscosibility [= (1/mu)(dmu/dp)] ------------    

    uWRf,uVCB = calcKestinVisc(mFrc,tRes,pRef,clsUNI)

#== Store results for Output ==========================================

    clsSIM.bRefW = BwRf
    clsSIM.uRefW = uWRf
    clsSIM.cRefW = comW
    clsSIM.vRefW = uVCB

    return

#========================================================================
#  Rowe and Chou Density Correlation
#  Input units are internal: mFrc - Mass Frac, tRes - degR, pRes - psia
#  R&C   units             : mFrc            , tRes - Kelv, pRes - kgf/cm2
#========================================================================

def calcRoweChouDen(mFrc,tRes,pRes,clsUNI) :

    pCon = 0.070307         #-- Pressure in kgf/cm2 (from Google!!)

    mFr2 = mFrc*mFrc
    tKel = clsUNI.I2X(tRes,"kelv")

    pkgF = pRes*pCon    
    pkg2 = pkgF*pkgF

    T1 = tKel
    
    T2 = T1*T1
    T3 = T1*T2

    a1,da1dT = rc_A1(T1,T2,T3)
    a2,da2dT = rc_A2(T1,T2,T3)
    a3,da3dT = rc_A3(T1,T2,T3)
    a4,da4dT = rc_A4(T1,T2,T3)
    a5,da5dT = rc_A5(T1,T2,T3)
    a6,da6dT = rc_A6(T1,T2,T3)
    a7,da7dT = rc_A7(T1,T2,T3)
    a8,da8dT = rc_A8(T1,T2,T3)

#== Density [kg/m3] ===================================================

    dRec = a1 - a2*pkgF - a3*pkg2 \
              + a4*mFrc + a5*mFr2 \
              - a6*pkgF*mFrc - a7*pkgF*mFr2 - 0.5*a8*pkg2*mFrc

    denW = 1000.0/dRec  #-- Density [kg/m3]

#== Derivatives [wrt Temp,Pres & Salinity] ============================

    dadT = da1dT - da2dT*pkgF - da3dT*pkg2 \
                 + da4dT*mFrc + da5dT*mFr2 \
                 - da6dT*pkgF*mFrc - da7dT*pkgF*mFr2 - 0.5*da8dT*pkg2*mFrc

    dDdT = - denW*dadT/dRec

    dDdP = - denW*(- a2 - 2.0*a3*pkgF - a6*mFrc -     a7*mFr2      -     a8*pkgF*mFrc)
    dDdS = - denW*(  a4 + 2.0*a5*mFrc - a6*pkgF - 2.0*a7*pkgF*mFrc - 0.5*a8*pkg2)

#== Return Values =====================================================    

    denW = clsUNI.X2I(denW,"kg/m3")
    
    dDdP = clsUNI.X2I(dDdP,"kg/m3")     #-- Pressure derivative: first kg/m3 -> lb/ft3
    dDdP = dDdP*pCon                #-- Then 1/kgf/cm2 -> 1/psia (inverse of way-in!)

    comP = dDdP/denW                #-- Compressibility [1/psia]

    return denW,comP

#========================================================================
#  Rowe and Chou A1-Coefficient
#========================================================================

def rc_A1(T1,T2,T3) :

    a1 =  5.916365E+00
    a2 = -1.035794E-02
    a3 =  9.270048E-06
    a4 = -1.127522E+03
    a5 =  1.006741E+05

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A2-Coefficient
#========================================================================

def rc_A2(T1,T2,T3) :

    a1 =  5.2049100E-03
    a2 = -1.0482101E-05
    a3 =  8.3285320E-09
    a4 = -1.1702939E+00
    a5 =  1.0227830E+02

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A3-Coefficient
#========================================================================

def rc_A3(T1,T2,T3) :

    a1 =  1.185470E-08
    a2 = -6.599143E-11
    a3 =  0.0
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A4-Coefficient
#========================================================================

def rc_A4(T1,T2,T3) :

    a1 = -2.51660E+00
    a2 =  1.11766E-02
    a3 = -1.70522E-05
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A5-Coefficient
#========================================================================

def rc_A5(T1,T2,T3) :

    a1 =  2.84851E+00
    a2 = -1.54305E-02
    a3 =  2.23982E-05
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A6-Coefficient
#========================================================================

def rc_A6(T1,T2,T3) :

    a1 = -1.48140E-03
    a2 =  8.29639E-06
    a3 = -1.24690E-08
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A7-Coefficient
#========================================================================

def rc_A7(T1,T2,T3) :

    a1 =  2.7141E-03
    a2 = -1.5391E-05
    a3 =  2.2655E-08
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A8-Coefficient
#========================================================================

def rc_A8(T1,T2,T3) :

    a1 =  6.2158E-07
    a2 = -4.0075E-09
    a3 =  6.5972E-12
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Kestin Viscosity Correlation
#========================================================================

def calcKestinVisc(mFrc,tRes,pRes,clsUNI) :

    molal = calcMolal(mFrc)
    tdegC = clsUNI.I2X(tRes,"degc")
    pBar  = clsUNI.I2X(pRes,"bara")

    pMPa  = 0.1*pBar

    mu0w20 = 1.002  #-- Pure water at 20 degC

    a = kestinA(molal)
    b = kestinB(molal)
    c = kestinC(tdegC)

    d = a + b*c

    mu0w = mu0w20*pow(10.0,c)

    muws = mu0w  *pow(10.0,d)

    pCor = 0.001*(0.8 + 0.01*(tdegC - 90.0)*exp(-0.25*molal))

    visW = muws*(1.0 + pCor*pMPa)   #--   Viscosity    [cP]

    dvdM = muws*pCor                #-- d[Viscosity/dp [cP/MPa]
    dvdp = dvdM/145.038             #--                [cP/psia]

    visC = dvdp/visW                #-- Viscosibility [= (1/mu)(dmu/dp) 1/psia]

    return visW,visC

#========================================================================
#  Kestin A-Sum
#========================================================================

def kestinA(M) :

    a1 =  3.324E-02
    a2 =  3.624E-03
    a3 = -1.879E-04

    A = M*(a1 + M*(a2 + M*a3))

    return A

#========================================================================
#  Kestin B-Sum
#========================================================================

def kestinB(M) :

    b1 = -3.96E-02
    b2 =  1.02E-02
    b3 = -7.02E-04

    B = M*(b1 + M*(b2 + M*b3))

    return B

#========================================================================
#  Kestin C-Sum
#========================================================================

def kestinC(degC) :

    c1 =  1.2378E+00
    c2 = -1.3030E-03
    c3 =  3.0600E-06
    c4 =  2.5500E-08

    T = 20.0 - degC

    sumC = T*(c1 + T*(c2 + T*(c3 + T*c4)))

    C = sumC/(96.0 + degC)

    return C

#========================================================================
#  Calculate Molality from Mass Fraction for Brine: W&B Eqn.(9.2c)
#  Where 10^6/PPM = 1/mFrac
#========================================================================

def calcMolal(mFrc) :

    if mFrc > 0.0 :
        molal = 17.1/(1.0/mFrc - 1.0)
    else :
        molal = 0.0

    return molal

#========================================================================
#  End of Module
#========================================================================

