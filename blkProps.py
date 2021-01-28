
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculates Properties for Blackoil Table Generation
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import blkCoefs  as BC
import calcEOS   as CE
import utilities as UT

from math import exp,log,log10,sqrt

#========================================================================
#  Calculate Phase FVF to re-fit (sOil,sGas)
#========================================================================

def calcPhaseFVF(qLiq,RTp,BoCon,xOil,clsBLK) :

    aOil = BC.Acoef(xOil,clsBLK)
    bOil = BC.Bcoef(xOil,clsBLK)
    cOil = BC.Ccoef(xOil,clsBLK)
    
    Moil = phaseMw(clsBLK.mSTG,clsBLK.mSTO,xOil)

    Aoil = aOil/(RTp*clsBLK.RT)
    Boil = bOil/ RTp

    Boi2 = Boil*Boil
    Boi3 = Boil*Boi2

    E2 = Boil - 1.0
    E1 = Aoil - 3.0*Boi2 -  2.0*Boil
    E0 = Boi3 +     Boi2 - Aoil*Boil

    nRoot,Zsmal,Zlarg = solveCubic(E2,E1,E0)    #-- Local copy of cubic-solver!!

    if qLiq : Voil = Zsmal*RTp - cOil
    else    : Voil = Zlarg*RTp - cOil

    Doil = Moil/Voil

    BoCal = BoCon/Doil

    return BoCal

#========================================================================
#  Reservoir Oil Density, Singh et al, SPE 109596, Eqn.(11)
#========================================================================

def denOil(dSTO,dSTG,Rs,Bo) :

    denO = (dSTO + Rs*dSTG)/Bo

    return denO

#========================================================================
#  Reservoir Gas Density, Singh et al, SPE 109596, Eqn.(19)
#========================================================================

def denGas(dSTO,dSTG,Rv,Bg) :

    denG = (dSTG + Rv*dSTO)/Bg

    return denG

#========================================================================
#  Phase Mole Weight
#========================================================================

def phaseMw(mSTG,mSTO,xOil) :

    xGas = 1.0 - xOil

    Mw = xGas*mSTG + xOil*mSTO

    return Mw

#========================================================================
#  Calculate Under-Saturated Oil Viscosity via Standing Correlation
#  See W+B Eqn.(3.130)
#========================================================================

def calcUndViscStand(pB,uB,p) :

    StanCon = 0.001*(0.024*pow(uB,1.6) + 0.038*pow(uB,0.56))

    u = uB + StanCon*(p - pB)

    return u

#========================================================================
#  2-Component (STO & STG) LBC Viscosity
#========================================================================

def calcLBCvisc(xOil,dOil,clsBLK) :

    dCof = BC.Rcoef(xOil,clsBLK)
    eCof = BC.Ecoef(xOil,clsBLK)
    uCof = BC.Mcoef(xOil,clsBLK)

    dRed = dOil/dCof

    sumV = sumVisc(dRed)

    vLBC = uCof +(pow(sumV,4) - 0.0001)/eCof

    return vLBC

#========================================================================
#  Derivative of LBC Viscosity wrt Critical Density Multipliers
#========================================================================

def calcLBCderv(xOil,dOil,clsBLK) :

    dCof = BC.Rcoef(xOil,clsBLK)
    eCof = BC.Ecoef(xOil,clsBLK)

    dRed = dOil/dCof

    sumV = sumVisc(dRed)
    sumD = sumDerv(dRed)

    dudR = 4.0*pow(sumV,3)*sumD*(-dRed/dCof)/eCof

    return dudR

#========================================================================
#  Fourth Power in Reduced Density for LBC Viscosity
#========================================================================

def sumVisc(dRed) :

#-- LBC Coefficients ------------------------------------------------    

    a0 =  0.1023
    a1 =  0.023364
    a2 =  0.058533
    a3 = -0.040758
    a4 =  0.0093324

    sumV = a0 + dRed*(a1 + dRed*(a2 + dRed*(a3 + dRed*a4)))

    return sumV

#========================================================================
#  Derivative of Fourth Power in Reduced Density for LBC Viscosity
#========================================================================

def sumDerv(dRed) :

#-- LBC Coefficients ------------------------------------------------    

    d1 =  0.023364      #-- d1 =   a1
    d2 =  0.117066      #-- d2 = 2*a2
    d3 = -0.122273      #-- d3 = 3*a3
    d4 =  0.0373296     #-- d4 = 4*a4

    sumD = d1 + dRed*(d2 + dRed*(d3 + dRed*d4))

    return sumD

#========================================================================
#  EoS (a,b) constant for Stock Tank Oil (STO)
#========================================================================

def initOilProps(dSTO,clsBLK) :

    OmegaA =   0.457235529  #-- PR Omega-A parameter
    OmegaB =   0.077796074  #-- PR Omega-B parameter

    clsBLK.dSTO = dSTO

    mSTO,oGrv = stockTankOilMw(dSTO)

#-- Singh Eqn.(32) and (33) -----------------------------------------    

    Toil =  608.0 + 364.0*log10(mSTO-71.2) + (2450.0*log10(mSTO     )-3800.0)*log10(oGrv)
    Poil = 1188.0 - 431.0*log10(mSTO-61.1) + (-852.0*log10(mSTO-53.7)+2310.0)*     (oGrv-0.8)

    Aoil = -0.039 + mSTO*(0.004 + 0.000003*mSTO)    #-- Singh Eqn.(37)

#-- Singh Eqn.(35) and (36) -----------------------------------------    

    if Aoil <= 0.49 : mAcf = lowAFmCoef(Aoil)
    else            : mAcf = higAFmCoef(Aoil)

    Tred = UT.tStand/Toil                           #-- Reduced Temperature

    alfa = pow((1.0 + mAcf*(1.0 - sqrt(Tred))),2)   #-- Singh Eqn.(34)

    aCof = OmegaA*pow(UT.gasCon*Toil,2)*alfa/Poil   #-- Singh Eqn.(29)                         
    bCof = OmegaB*    UT.gasCon*Toil        /Poil   #-- Singh Eqn.(30)
    sCof = 0.0                                      #-- Singh Eqn.(31)

    Voil = 21.573 + 0.015122*mSTO - 27.656*oGrv + 0.070615*mSTO*oGrv    #--Singh Eqn.(66)
    Zoil = Poil*Voil/(UT.gasCon*Toil)

    rCof = mSTO/Voil                                            #-- Density
    eCof = 5.35*pow(Toil,UT.lbcTcX)*pow(mSTO,UT.lbcMwX)*pow(Poil,UT.lbcPcX) #-- Singh (70)

    if Tred <= 1.5 :
        uCof = 34.10E-05*pow(      Tred        ,0.94 )/eCof     #-- Singh (71)
    else :
        uCof = 17.78E-05*pow((4.58*Tred - 1.67),0.625)/eCof     #-- Singh (72)

    #print("InitOilProps: aCof,bCof,sCof {:10.3e}{:10.3e}{:10.3e}".format(aCof,bCof,sCof))
    #print("initOilProps: rCof,eCof,uCof {:10.3e}{:10.3e}{:10.3e}".format(rCof,eCof,uCof))

    clsBLK.setOilCrit(Toil,Poil,Voil,Zoil)        

#========================================================================
#  End of Routine
#========================================================================

    clsBLK.EOS0["aOil"] = aCof
    clsBLK.EOS0["bOil"] = bCof
    clsBLK.EOS0["sOil"] = sCof

    clsBLK.VIS0["rOil"] = rCof
    clsBLK.VIS0["eOil"] = eCof
    clsBLK.VIS0["uOil"] = uCof

    return mSTO,oGrv

#========================================================================
#  EoS (a,b) constant for Stock Tank Gas (STG)
#========================================================================

def initGasProps(dSTG,clsBLK) :

    OmegaA =   0.457235529  #-- PR Omega-A parameter
    OmegaB =   0.077796074  #-- PR Omega-B parameter

    clsBLK.dSTG = dSTG

    gGrv = dSTG/UT.denAir
    mSTG = gGrv*UT.molAir

    if gGrv < 0.75 :
        Tgas = 168.0 + gGrv*(325.0 - gGrv*12.5)     #-- Singh Eqn.(41)
        Pgas = 667.0 + gGrv*( 15.0 - gGrv*37.5)     #-- Singh Eqn.(42)
    else :
        Tgas = 187.0 + gGrv*(330.0 - gGrv*71.5)     #-- Singh Eqn.(43)
        Pgas = 706.0 + gGrv*( 51.7 - gGrv*11.1)     #-- Singh Eqn.(44)

    Agas = 0.1637*gGrv - 0.0792                     #-- Singh Eqn.(48)

    if Agas <= 0.49 : mAcf = lowAFmCoef(Agas) #-- Singh Eqn.(46)
    else            : mAcf = higAFmCoef(Agas) #-- Singh Eqn.(47)

    Tred = UT.tStand/Tgas                           #-- Reduced Temperature

    alfa = pow((1.0 + mAcf*(1.0 - sqrt(Tred))),2)   #-- Singh Eqn.(45)

    aCof = OmegaA*pow(UT.gasCon*Tgas,2)*alfa/Pgas   #-- Singh Eqn.(38)       
    bCof = OmegaB*    UT.gasCon*Tgas        /Pgas   #-- Singh Eqn.(39)
    sCof = 0.0                                      #-- Singh Eqn.(40)

    Zgas = 0.2905 - 0.085*Agas      #-- Reid Correlation, Whitson Eqn.(5.66)
    Vgas = Zgas*UT.gasCon*Tgas/Pgas #-- Singh Eqn.(68)

    rCof = mSTG/Vgas
    eCof = 5.35*pow(Tgas,UT.lbcTcX)*pow(mSTG,UT.lbcMwX)*pow(Pgas,UT.lbcPcX) #-- Singh (73)

    if Tred <= 1.5 :
        uCof = 34.10E-05*pow(      Tred,        0.94 )/eCof     #-- Singh (74)
    else :
        uCof = 17.78E-05*pow((4.58*Tred - 1.67),0.625)/eCof     #-- Singh (75)

    #print("InitGasProps: aCof,bCof,sCof {:10.3e}{:10.3e}{:10.3e}".format(aCof,bCof,sCof))
    #print("initGasProps: rCof,eCof,uCof {:10.3e}{:10.3e}{:10.3e}".format(rCof,eCof,uCof))

    clsBLK.setGasCrit(Tgas,Pgas,Vgas,Zgas)        

#========================================================================
#  End of Routine
#========================================================================

    clsBLK.EOS0["aGas"] = aCof
    clsBLK.EOS0["bGas"] = bCof
    clsBLK.EOS0["sGas"] = sCof

    clsBLK.VIS0["rGas"] = rCof
    clsBLK.VIS0["eGas"] = eCof
    clsBLK.VIS0["uGas"] = uCof

    return mSTG,gGrv

#========================================================================
#  Low (<  0.49) Acentric Factor m-Coefficient [Singh Eqn.35 and 46)
#========================================================================

def lowAFmCoef(AF) :

    mAcf = 0.37464 + 1.54226*AF - 0.26992*AF*AF

    return mAcf

#========================================================================
#  High (>= 0.49) Acentric Factor m-Coefficient [Singh Eqn.36 and 47)
#========================================================================

def higAFmCoef(AF) :

    AF2 = AF*AF

    mAcf = 0.3796 + 1.485*AF - 0.1644*AF2 + 0.01667*AF2*AF

    return mAcf

#========================================================================
#  Stock Tank Oil Mole Weight
#  Singh et al, SPE 109596, Eqn.(15)
#========================================================================

def stockTankOilMw(dSTO) :

    oGrv = dSTO/UT.denWat
    oAPI = 141.5/oGrv - 131.5

    mSTO = 240.0 - 2.22*oAPI

    return mSTO,oGrv

#========================================================================
#  Stock Tank Gas Mole Weight
#  Singh et al, SPE 109596, Eqn.(16)
#========================================================================

def stockTankGasMw(dSTG) :

    gGrv = dSTG/UT.denAir
    mSTG = gGrv*UT.molAir

    return mSTG,gGrv

#========================================================================
#  Ensure Oil Viscosity is Monotonically Decreasing
#  Condensates can go wrong!
#========================================================================

def oilViscMono(dTab,clsBLK) :

    nSat = len(dTab)

    qMon = True
    uLst = 0.0

#-- Test for Monotonicity -------------------------------------------    

    for i in range(nSat) :
        Uo = dTab[i][clsBLK.iUo]
        if Uo > uLst :
            uLst = Uo
        else :
            qMon = False    #-- NOT monotonic
            break

#== Use the fit of ln(Uo) = A + B*Rs instead ==========================        

    if not qMon :

        for i in range(nSat) :
            
            Pr = dTab[i][clsBLK.iPr]
            Rs = dTab[i][clsBLK.iRs]
            
            Uo = clsBLK.UoI*exp(clsBLK.UoS*Rs)
            
            #print("i,Pr,Rs,Uo {:2d}{:10.3f}{:10.3f}{:10.5f}".format(i,Pr,Rs,Uo))

            dTab[i][clsBLK.iUo] = Uo  #-- Over-write the saturated data

#== No return value ===================================================

    return qMon

#========================================================================
#  Ensure Oil Viscosity is Monotonically Decreasing
#  Condensates can go wrong!
#========================================================================

def oilCompMono(dTab,clsBLK) :

    nSat = len(dTab)

    qMon = True
    cLst = 0.0

#-- Test for Monotonicity -------------------------------------------    

    for i in range(nSat) :
        Co = dTab[i][clsBLK.iCo]
        if Co > cLst :
            cLst = Co
        else :
            qMon = False    #-- NOT monotonic
            break

#== Use the fit of ln(Uo) = A + B*Rs instead ==========================        

    if not qMon :

        for i in range(nSat) :
            
            Pr = dTab[i][clsBLK.iPr]
            Rs = dTab[i][clsBLK.iRs]
            
            Co = (clsBLK.CoI + clsBLK.CoS)/Pr
            
            #print("i,Pr,Rs,Co {:2d}{:10.3f}{:10.3f}{:10.5f}".format(i,Pr,Rs,Co))

            dTab[i][clsBLK.iCo] = Co  #-- Over-write the saturated data

#== No return value ===================================================

    return qMon

#========================================================================
#  Calculate Under-Saturated Oil Viscosity via Standing Correlation
#  See W+B Eqn.(3.130)
#========================================================================

def calcUndViscStand(pB,uB,p) :

    StanCon = 0.001*(0.024*pow(uB,1.6) + 0.038*pow(uB,0.56))

    u = uB + StanCon*(p - pB)

    return u

#========================================================================
#  A-Coeff of Abdul-Majeed Under-Saturated Viscosity Correlation
#  See W+B Eqn.(3.132b)
#========================================================================

def abdulMajeedA(dSTO,Rs) :

    gSTO =  dSTO/UT.denWat      #-- Oil Specific Gravity
    oAPI = 141.5/gSTO - 131.5   #-- Oil API Gravity

    lnRs = log(5.615*Rs)

    A = 1.9311 - 0.89941*lnRs - 0.001194*oAPI*oAPI + 0.0092545*oAPI*lnRs

#== Return value ======================================================

    return A

#========================================================================
#  Abdul-Majeed Under-Saturated Viscosity Correlation
#  See W+B Eqn.(3.132a)
#========================================================================

def abdulMajeedVisc(A,P,Pb,Ub) :

    expA = A - 5.2106 + 1.11*log10(P - Pb)

    Uo = Ub + pow(10.0,expA)

#== Return value ======================================================

    return Uo

#========================================================================
#  Solve the Cubic x^3 + E2 x^2 + E1 x + E0 = 0
#  See Numerical Recipes Section 5.6 [E2 = a, E1 = b, E0 = c]
#  http://www.aip.de/groups/soe/local/numres/bookfpdf/f5-6.pdf
#========================================================================

def solveCubic(E2,E1,E0) :

    F2 = E2*UT.Third                  #--  a/3
    F1 = E1*UT.Third                  #--  b/3

    F22 = F2*F2                       #-- (a/3)^2

    Q =    F22 - F1                   #-- Eqn.(5.6.10)
    R = F2*F22 - F1*F2*1.5 + 0.5*E0   #-- Eqn.(5.6.10)

#-- Discriminant ----------------------------------------------------    

    D = R*R - Q*Q*Q

#== One Real Root =====================================================
    
    if   D > 0.0 :
        
        argA = NP.absolute(R) + NP.sqrt(D)             #-- Arg of Eqn.(5.6.13)
        A13  = NP.power(argA,UT.Third)
        A    = -A13*NP.sign(R)                         #-- Eqn.(5.6.15)
        
        if A == 0.0 : B = 0.0                          #-- Eqn.(5.6.16)
        else        : B = Q/A
        
        x1 = A + B                                     #-- Eqn.(5.6.17)
        x3 = x1

        nR = 1

#== Three Real Roots: Ignore Central Root as Unphysical in EoS ========
        
    elif D < 0.0 :
        
        sqrtQ = NP.sqrt(Q)
        theta = NP.arccos(R/(Q*sqrtQ))                 #-- Eqn.(5.6.11)
        cTh13 = NP.cos(UT.Third*theta)
        
        x1 = -2.0*sqrtQ*cTh13                          #-- Eqn.(5.6.12-x1)
        
#-- Next expression have used cos(A-B) = cosA.cosB + sinA.sinB where
#-- cos(-2.PI/3) = -1/2, sin(-2.PI/3) = -sqrt(3)/2 and sinA = sqrt(1-cosA^2)

        x3 = sqrtQ*(cTh13+NP.sqrt(3.0*(1.0-cTh13*cTh13)))

        nR = 3

#== Error in NR Text: see http://numerical.recipes/forum/showthread.php?t=865
        
    else :
        
        sqrtQ = NP.sqrt(Q)
        
        if   R > 0.0 :
            x1 = -2.0*sqrtQ
            x3 =      sqrtQ
        elif R < 0.0 :
            x1 =     -sqrtQ
            x3 =  2.0*sqrtQ
        else :
            x1 = 0.0
            x3 = 0.0

        nR = 3

#== Return Values =====================================================

    eLow = x1 - F2
    eHig = x3 - F2

    if nR > 1 and eLow < 0.0 : eLow = eHig

    return nR,eLow,eHig

#========================================================================
#  End of Module
#========================================================================
