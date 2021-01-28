
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculates EOS Coefficients and Log-Phi-Coefficients
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

import numpy     as NP

import utilities as UT

from math  import log,sqrt

import calcTest as CT

#========================================================================
#  Calculate EoS Coefficients and its derivatives
#========================================================================

def calcEoSCoefsT(p,T,x,qT,clsEOS) :

    nCom = clsEOS.nComp
    sqrT = sqrt(T)

#== b => B and a => A conversion ======================================

    pByT  =    p/T
    pByT2 = pByT/T

#== Bi & B ============================================================

    Bi = pByT*NP.copy(clsEOS.bI)
    B  =      NP.dot( x,Bi)

#== Form a(T) = p + sqrt(T)*q =========================================

    aI = clsEOS.aP + sqrT*clsEOS.aQ    #-- Vector!

#-- Work arrays -----------------------------------------------------    

    cI = NP.multiply(x,aI)      #-- Vector!
    dI = NP.dot(clsEOS.LIJ,cI)  #-- Vector where Lij = 1-Kij

#== Lij = 1 - Kij =====================================================

    Ai = pByT2*NP.multiply(aI,dI)
    A  = pByT2*     NP.dot(cI,dI)

#----------------------------------------------------------------------
#  Temperature Derivative Terms Requested
#----------------------------------------------------------------------

    if qT :

        cT = 0.5*pByT2/sqrT
        tC = 2.0/T

        eI = NP.multiply(x,clsEOS.aQ)
        fI = NP.dot(clsEOS.LIJ,eI)

        #for iC in range(nCom) :
        #    fI[iC] = NP.dot(clsEOS.LIJ[iC,:],eI)

        dAidT = cT*NP.multiply(aI,fI) + cT*NP.multiply(clsEOS.aQ,dI) - tC*Ai
        dAdT  = cT*     NP.dot(cI,fI) + cT*     NP.dot(       eI,dI) - tC*A

    else :

        dAdT  = clsEOS.dumS
        dAidT = clsEOS.dumV

#== Return values =====================================================

    return A,B,Ai,Bi,dAdT,dAidT

#========================================================================
#  EoS Aij Matrix
#========================================================================

def calcEoSCoefsX(p,T,clsEOS) :

    sqrT  = sqrt(tRes)
    pByT2 = p/(T*T)

#== Calculate a = p + sqrt(T)*q =======================================    

    aI = clsEOS.aP + sqrT*clsEOS.aQ    #-- Vector!

#== Build the Aij Matrix ==============================================

    Aij = pByT2*NP.outer(aI,aI)
    Aij =       NP.multiply(clsEOS.LIJ,Aij)  #-- Lij = 1 - Kij

#== Return the Matrix =================================================

    return Aij

#========================================================================
#  Calculate Phase Fugacity and (Optionally) its (P,T,X)-Derivatives
#========================================================================

def calcPhaseFugPTX(iPhs,qP,qT,qX,p,T,X,clsEOS) :

    nCom = clsEOS.nComp
    
#== Calculate the EoS Coefficients ====================================

    A,B,Ai,Bi,dAdT,dAidT = calcEoSCoefsT(p,T,X,qT,clsEOS)

#-- alpha[i] = A[i]/A and beta[i] = B[i]/B --------------------------    

    alfa = Ai/A     #-- Vector!
    beta = Bi/B     #-- Vector

#== Setup and Solve the Cubic EoS =====================================

    Eta,dFdE  = setupSolveCubicD(iPhs,A,B,clsEOS)

#-- (P,T)-Derivaties of Eta -----------------------------------------

    if qP or qT or qX :
        bWrk = 2.0*clsEOS.n4*B
        EtaA = - Eta/dFdE
        EtaB = -(clsEOS.n3*Eta*Eta +(bWrk - clsEOS.n3)*Eta - bWrk)/dFdE

#== Fugacity Coefficients =============================================

    n1B = clsEOS.n1*B
    n2B = clsEOS.n2*B

    sig1 = Eta + n1B
    sig2 = Eta + n2B

    C1 = -log(Eta)
    C2 =      Eta + B - 1.0
    C3 =  clsEOS.n5*A/B
    C4 =  log(sig2/sig1)

    DA = C3*C4
    DB = C2 - DA
    DA = 2.0*DA

    fugC = C1 + DA*alfa + DB*beta   #-- Vector!

#== Pressure Derivatives ==============================================

    if qP :

        dAdP = A/p
        dBdP = B/p

        dEdP = EtaA*dAdP + EtaB*dBdP

        dC1dP = -dEdP/Eta
        dC2dP =  dEdP + dBdP
        dC3dP =  0.0
        dC4dP = (dEdP + clsEOS.n2*dBdP)/sig2 - (dEdP + clsEOS.n1*dBdP)/sig1

        dDAdP =  C3*dC4dP + C4*dC3dP
        dDBdP =     dC2dP -    dDAdP
        dDAdP = 2.0*dDAdP

        dFdP = dC1dP + dDAdP*alfa + dDBdP*beta  #-- Vector!

    else :
        
        dFdP = clsEOS.dumV

#== Temperature Derivatives ===========================================

    if qT :

        dBdT = -B/T      #-- NB. dAdT (and dAi/dT returned by calcEoSCoefsT)

        dEdT = EtaA*dAdT + EtaB*dBdT

        dC1dT = -dEdT/Eta
        dC2dT =  dEdT + dBdT
        dC3dT =  C3*(dAdT/A - dBdT/B)
        dC4dT = (dEdT + clsEOS.n2*dBdT)/sig2 - (dEdT + clsEOS.n1*dBdT)/sig1

        dDAdT =  C3*dC4dT + C4*dC3dT
        dDBdT =     dC2dT -    dDAdT
        dDAdT = 2.0*dDAdT

        con1 = DA/A
        con2 = dAdT*con1

        dAlf = con1*dAidT - con2*alfa                   #-- Vector!

        dFdT = dC1dT + dDAdT*alfa + dDBdT*beta + dAlf   #-- Vector!

    else :
        
        dFdT = clsEOS.dumV

#== Composition Derivatives ===========================================

    if qX :

        AiJ = calcEoSCoefsX(p,T,clsEOS)

        EtaAC = EtaA*A
        EtaBC = EtaB*B

        C1A = - EtaAC/Eta ; C1B = - EtaBC/Eta ; C1T = C1A + C1B
        C2A =   EtaAC     ; C2B =   EtaBC + B
        C3A =   C3        ; C3B = - C3
        
        C4A =   EtaAC*(1.0/sig2 - 1.0/sig1)
        C4B = (EtaAC + n2B*B)/sig2 - (EtaBC + n1B)/sig1

        DAA = C3*C4A + C4*C3A ; DAB = C3*C4B + C4*C3B
        DBA =    C2A -    DAA ; DBB =    C2B -    DAB

        EIJ = 2.0* DA/A
        EAA = 2.0*(DAA - DA)
        EAB = 2.0* DAB
        EBA =      DBA
        EBB =      DBB - DB

        FA  = EIJ + EAA + EAB
        FB  =       EBA + EBB

        sFX = C1T +  FA*alfa +  FB*beta     #-- Sum[k] x[k]*dlogPhi[i]/dx[k]
        C1X =       C1A*alfa + C1B*beta     #-- Work Array for C1

        TC1T = NP.array([C1X]).T            #-- Transpose vector to Column (Matrix)

#-- dlogPhi[i]/dX[j] (Mole Number Derivative!!) ---------------------

        mAA = EAA*NP.outer(alfa,alfa)       #-- Matrix!!
        mAB = EAB*NP.outer(alfa,beta)       #-- Matrix!!
        mBA = EBA*NP.outer(beta,alfa)       #-- Matrix!!
        mBB = EBB*NP.outer(beta,beta)       #-- Matrix!!

        dFdX = TC1T + EIJ*Aij + mAA + mAB + mBA + mBB - SFX

        #for iC in range(nCom) :
        #    for jC in range(nCom) :
        #        dFdX[iC][jC] = C1X[jC] + EIJ*Aij[iC][jC] + \
        #                       EAA*alfa[iC]*alfa[jC] + \
        #                       EAB*alfa[iC]*beta[jC] + \
        #                       EBA*beta[iC]*alfa[jC] + \
        #                       EBB*beta[iC]*beta[jC] - sFX[iC]

    else :

        dFdX = clsEOS.dumM

#== Return Information ================================================

    return fugC,dFdP,dFdT,dFdX

#========================================================================
#  Setup and Solve the Cubic EoS and its Eta-Derivative
#========================================================================

def setupSolveCubicD(iLiq,A,B,clsEOS) :

#----------------------------------------------
#  If 3 Real-Roots, then
#  iLiq =  1 => Look for Liquid (smallest) Root
#  iLiq = -1 => Look for Vapour (largest)  Root
#  iLiq =  0 => Use min-GFE to select      Root
#  Else if 1 Real-Root, that's it!!
#----------------------------------------------

#-- Cubic Coefficients of f(x) = E3*x^3 + E2*x^2 + E1*x + E0 = 0 ----   

    E0 = -clsEOS.n4*B*B
    E2 =  clsEOS.n3*B - 1.0
    E1 = -E0 - E2 - 1.0 + A

    #print("setupSolveCubic: E2,E1,E0 ",E2,E1,E0)

    Eta = sortCubicRoots(iLiq,E2,E1,E0,A,B,clsEOS)

    #CT.numCubicSolver(E2,E1,E0,Eta)

    dFdE = 3.0*Eta*Eta + 2.0*E2*Eta + E1

    return Eta,dFdE

#========================================================================
#  Sort the Roots of the Cubic EoS
#========================================================================

def sortCubicRoots(iLiq,E2,E1,E0,A,B,clsEOS) :

#== Analytic Soution ==================================================
        
    nReal,rSmal,rLarg = solveCubic(E2,E1,E0)

#-- Sort Roots Depending on Number/Preferred Type -------------------

    if nReal == 1 :        #-- Take the Only Root
        rRoot = rSmal
    else :
        if   iLiq ==  1 :  #-- Liquid Root Required
            rRoot = rSmal
        elif iLiq == -1 :  #-- Vapour Root Required
            rRoot = rLarg
        else :             #-- Select Root on min-GFE
            ABn5  = clsEOS.n5*(A/B)
            n1B   = clsEOS.n1* B
            n2B   = clsEOS.n2* B
            Bm1   =            B - 1.0
            gSmal = (rSmal + Bm1) - log(rSmal) + \
                    ABn5*log((rSmal + n2B)/(rSmal + n1B))
            gLarg = (rLarg + Bm1) - log(rLarg) + \
                    ABn5*log((rLarg + n2B)/(rLarg + n1B))
            if gSmal < gLarg : rRoot = rSmal
            else             : rRoot = rLarg

            #print("nReal>1: eSmal,gSaml,eLarg,gLarg {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(rSmal,gSmal,rLarg,gLarg))
        
#== Smallest root -ve => discard in favour of +ve root ================            

    if nReal > 1 and rSmal < 0.0 :
        nReal = 1
        rSmal = rLarg

    return rRoot

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
#  Square-Root of the EoS alpha-Coefficient
#========================================================================

def sqrtAlpha(tRed,AF,clsEOS) :

    mAcF = calcMAcF(AF,clsEOS)

    sqrtAlfa = 1.0 + mAcF*( 1.0 - sqrt(tRed) )

#== Return values =====================================================                

    return sqrtAlfa

#========================================================================
#  m(Acf) Expansion
#========================================================================

def calcMAcF(AF,clsEOS) :

    mAcF = clsEOS.w0 + AF*( clsEOS.w1 + AF*( clsEOS.w2 + AF*clsEOS.w3 ) )

#== Return values =====================================================                

    return mAcF

#========================================================================
#  d[m(Acf)]/dw
#========================================================================

def calcdMdw(AF,clsEOS) :

    dMdw =  clsEOS.w1 + AF*( 2.0*clsEOS.w2 + 3.0*AF*clsEOS.w3 )

#== Return values =====================================================                

    return dMdw

#========================================================================
#  Calculate Key Properties and their Derivatives
#========================================================================

def calcProps(iLiq,pRes,tRes,X,clsEOS) :

    nCom = clsEOS.nComp
    qT   = False

#-- Calculate EoS Coefficients ----------------------------------------

    A,B,AiC,BiC,dumS,dumV = calcEoSCoefsT(pRes,tRes,X,qT,clsEOS)

#-- Setup and Solve the Cubic EoS (and its p-Derivative) --------------

    Eta,dFdE  = setupSolveCubicD(iLiq,A,B,clsEOS)

    dAdP,dBdP,dEdP = calcdEdP(pRes,A,B,Eta,dFdE,clsEOS)

#-- 2-Parameter (EoS) Z-Factor (and its p-derivative) & Molar Volume 

    ZF2  = Eta  + B
    Vm2  = ZF2*UT.gasCon*tRes/pRes
    dZdp = dEdP + dBdP

#-- Mole Weight and Volume Shift --------------------------------------

    cSI = UT.gasCon*NP.multiply(clsEOS.bI,clsEOS.vS)    #-- Vector!

    cS = NP.dot(X,cSI)
    Mw = NP.dot(X,clsEOS.Mw)

#-- 3-Parameter (volume shifted) Molar Volume, Density & Viscosity --

    Vm3 = Vm2 - cS
    den = Mw/Vm3

#-- Viscosity and its derivative wrt Density ------------------------    

    visC,dvdR = calcLBCderv(tRes,den,X,clsEOS)

#-- 3-Parameter Z-Factor --------------------------------------------

    ZF3 = pRes*Vm3/(UT.gasCon*tRes)

#-- Compressibility -------------------------------------------------

    comp = (Vm2/Vm3)*(1.0/pRes - dZdp/ZF2)

#-- d(density)/dp ---------------------------------------------------

    dDdp = den*comp

#-- d(viscosity)/dp = d(visc)/d(den) * d(den)/dp --------------------

    dvdp = dvdR*dDdp

#== Return Values =====================================================    

    return Mw,Vm3,den,ZF3,visC,comp,dvdp

#========================================================================
#  Calculate d(Eta)/dp
#========================================================================

def calcdEdP(pRes,A,B,Eta,dFdE,clsEOS) :

    dAdP = A/pRes
    dBdP = B/pRes

    bWrk = 2.0*clsEOS.n4*B

    EtaA = - Eta/dFdE
    EtaB = - (clsEOS.n3*Eta*Eta + (bWrk - clsEOS.n3)*Eta - bWrk)/dFdE

    dEdP = EtaA*dAdP + EtaB*dBdP

    return dAdP,dBdP,dEdP

#========================================================================
#  Lohrenz Bray Clark Viscosity [Whitson & Brule, Page 38]
#  And its derivative (wrt Density)
#========================================================================

def calcLBCderv(tRes,denS,X,clsEOS) :

#-- LBC Correlation Coefficients ------------------------------------    

    con0 =  0.1023
    con1 =  0.023364
    con2 =  0.058533
    con3 = -0.040758
    con4 =  0.0093324

#== Vector Calculations from here ... =================================

    sqM  = NP.sqrt(       clsEOS.Mw)    #-- sqrt(Mwi)
    tRed = NP.divide(tRes,clsEOS.Tc)    #-- tRedi = tRes/Tci

    TcPow = NP.power(clsEOS.Tc,UT.lbcTcX)   #-- Tci^( 1/6)
    PcPow = NP.power(clsEOS.Pc,UT.lbcPcX)   #-- Pci^(-2/3)
    MwPow = NP.power(clsEOS.Mw,UT.lbcMwX)   #-- Mwi^(-1/2)

    Xi = 5.35*TcPow*MwPow*PcPow

    Tpc = NP.dot(X,clsEOS.Tc)   #-- sum[i] xi*Tci
    Ppc = NP.dot(X,clsEOS.Pc)   #-- sum[i] xi*Pci
    Vpc = NP.dot(X,clsEOS.Vc)   #-- sum[i] xi*Vci
    Mpc = NP.dot(X,clsEOS.Mw)   #-- sum[i] xi*Mwi
    muD = NP.dot(X,sqM)         #-- sum[i] xi*sqrt(Mwi)

    mu1 = 4.58*tRed - 1.67
    
    mu1 = NP.where(mu1>0.0,mu1,0.0)     #-- Protect against mu1 < 0.0
    
    mu1 = 17.78E-05*NP.power(mu1 ,0.58)
    mu2 = 34.00E-05*NP.power(tRed,0.94)

    muI = NP.where(tRed>1.5,mu1,mu2)

    muI = NP.divide(muI,Xi)

    smI = NP.multiply(sqM,muI)

    muN = NP.dot(X,smI)

#== Scalar Calculations from here on ... ============================    

    mu0 = muN/muD
    XiR = 5.35*pow(Tpc,UT.lbcTcX)*pow(Mpc,UT.lbcMwX)*pow(Ppc,UT.lbcPcX)
    Rpc = Mpc/Vpc

    #print("calcLBCvisc: Tpc,Ppc,Vpc,Mpc ",Tpc,Ppc,Vpc,Mpc)

    rhoR = denS/Rpc     #-- Reduced Density

    rhoR2 = rhoR*rhoR
    rhoR3 = rhoR*rhoR2
    rhoR4 = rhoR*rhoR3

    sumR = con0 + con1*rhoR +     con2*rhoR2 +     con3*rhoR3 +     con4*rhoR4
    sumD =        con1      + 2.0*con2*rhoR  + 3.0*con3*rhoR2 + 4.0*con4*rhoR3

    visC = mu0 + (pow(sumR,4) - 1.0E-04)/XiR

    dvdR = 4.0*pow(sumR,3)*sumD/(rhoR*XiR)

#== Return function and its derivative ================================    

    return visC,dvdR

#========================================================================
#  IFT: W&B, Sec.3.5, Eqn.(3.136)
#========================================================================

def calcIFT(mVL,mVV,x,y,clsEOS) :

    nCom = clsEOS.nComp

#-- Mole Weighted Parachors -----------------------------------------

    gasT = 0.0
    oilT = 0.0

    for iC in range(nCom) :

        paraC = clsEOS.gPP("PA",iC)

        gasT = gasT + y[iC]*paraC
        oilT = oilT + x[iC]*paraC

#-- Inverse Molar Volume of Gas & Oil [Vm in gmol/cm3!!] ------------

    gasT = gasT/(62.428*mVV)
    oilT = oilT/(62.428*mVL)

#-- 1/4-Power of IFT ------------------------------------------------

    IFT14 = oilT - gasT

#== IFT [dyne/cm = mN/m] ==============================================

    IFT = IFT14**4

#== Return IFT ========================================================

    return IFT

#========================================================================
#  Enthalpy and Specific Heat
#========================================================================

def calcEnthSpecHeat(iPhs,p,T,x,clsEOS) :

    pByT  =    p/T
    pByT2 = pByT/T

    RT = UT.btuCon*T

#== Ideal Gas Contributions ===========================================

    CpIdl,HIdl = calcIdealGasTherm(T,x,clsEOS)

#== b & a EoS coeffs (and a-derivatives [1st & 2nd] wrt temperature) ==

    b,a,dadT,d2adT2 = calcEoSCoefsA(p,T,x,clsEOS)

#-- (A,B) Coeffs ----------------------------------------------------

    A = a*pByT2
    B = b*pByT

#-- Temperature Derivatives -----------------------------------------

    dBdT = -    B/T
    dAdT = -2.0*A/T + pByT2*dadT

#== Setup and Solve the Cubic EoS =====================================

    Eta,dFdE  = setupSolveCubicD(iPhs,A,B,clsEOS)

#-- d[eta]/dT -------------------------------------------------------

    Eta2 = Eta*Eta
    bWrk = 2.0*clsEOS.n4*B
    
    EtaA = - Eta/dFdE
    EtaB = - (clsEOS.n3*Eta2 + (bWrk - clsEOS.n3)*Eta - bWrk)/dFdE

    dEdT = EtaA*dAdT + EtaB*dBdT

    dZdT = dEdT + dBdT

#== Assemble Residual Enthalpy ========================================

    sig1 = Eta + clsEOS.n1*B
    sig2 = Eta + clsEOS.n2*B

    C2 = Eta + B - 1.0
    C3 = clsEOS.n5*A/B
    C4 = log(sig2/sig1)

    aD = dadT/a

    C5 = T*aD - 1.0

    HRes  = RT*(C2 + C3*C4*C5)

#== Residual Specific Heat ============================================

    dC2dT =     dEdT   + dBdT
    dC3dT = C3*(dAdT/A - dBdT/B)  #-- dB/dT = -1/T but leave as-is
    dC4dT = (dEdT + clsEOS.n2*dBdT)/sig2 - (dEdT + clsEOS.n1*dBdT)/sig1

    dC5dT = aD - T*aD*aD + T*d2adT2/a

    CpRes = HRes/T + RT*(dC2dT + C3*C4*dC5dT + C3*dC4dT*C5 + dC3dT*C4*C5)

#== Total Enthalpy & Specific Heat ====================================

    HTot  = HIdl  + HRes
    CpTot = CpIdl + CpRes

#== Joule-Thomson Coefficient =========================================

    uJT = RT*T*dZdT/(p*CpTot)

#== Return H & Cp = dH/dT =============================================

    return HTot,CpTot,uJT

#========================================================================
#  Calculate Ideal Gas Specific Heat and Enthalpy
#========================================================================

def calcIdealGasTherm(T,x,clsEOS) :

    nCom = clsEOS.nComp
    Tref = clsEOS.Tref

    T1 = T     ; T2 = T1*T1 ; T3 = T1*T2 ; T4 = T1*T3
    U1 = Tref  ; U2 = U1*U1 ; U3 = U1*U2 ; U4 = U1*U3

    tA =  T1 - U1      ;  tB = (T2 - U2)/2.0
    tC = (T3 - U3)/3.0 ;  tD = (T4 - U4)/4.0

#== Assemble sums =====================================================

    CpI = 0.0 ; HId = 0.0

    for iC in range(nCom) :
        
        CpA = clsEOS.gPP("CA",iC)
        CpB = clsEOS.gPP("CB",iC)
        CpC = clsEOS.gPP("CC",iC)
        CpD = clsEOS.gPP("CD",iC)

        CpT = CpA    + CpB*T1 + CpC*T2 + CpD*T3
        HIT = CpA*tA + CpB*tB + CpC*tC + CpD*tD

        CpI = CpI + x[iC]*CpT
        HId = HId + x[iC]*HIT        

    return CpI,HId

#========================================================================
#  b, a, da/dT & d2a/dT2 (for Enthalpy/Specific Heat)
#========================================================================

def calcEoSCoefsA(p,T,x,clsEOS) :

    nCom = clsEOS.nComp
    sqrT = sqrt(T)

#== b-Coefficient =====================================================

    b = NP.dot(x,clsEOS.bI)

#== Sub-terms =========================================================

    pI = NP.multiply(x,clsEOS.aP)
    qI = NP.multiply(x,clsEOS.aQ)

    #sI = NP.zeros(nCom)
    #tI = NP.zeros(nCom)

    #for iC in range(nCom) :
    #    sI[iC] = NP.dot(clsEOS.LIJ[iC,:],pI)
    #    tI[iC] = NP.dot(clsEOS.LIJ[iC,:],qI)

    sI = NP.dot(clsEOS.LIJ,pI)
    tI = NP.dot(clsEOS.LIJ,qI)

    aT1 = NP.dot(pI,sI)
    aT2 = NP.dot(pI,tI)
    aT3 = NP.dot(qI,tI)

#== a and its first/second derivative wrt Temperature =================    

    a      = aT1 + aT2* sqrT + aT3*T
    dadT   =       aT2/ sqrT + aT3
    d2adT2 = - 0.5*aT2/(sqrT*T)

#== Return values =====================================================

    return b,a,dadT,d2adT2

#========================================================================
#  Cubic Newton
#========================================================================

def cubicNewton(eta,E2,E1,E0) :

    D2 = 2.0*E2         #-- Work Scalar for 1st/2nd Derivative

#-- Function, 1st & 2nd Derivatives at current estimate -------------    

    f0 = cub0(eta,E2,E1,E0) ; f1 = cub1(eta,D2,E1) ; f2 = cub2(eta,D2)
    
#-- Halley-Newton [3rd Order] Update --------------------------------        
        
    nCub = 1
    dNew =    f0/ f1
    dHal = dNew/(1.0 - 0.5*dNew*f2/f1)

#----------------------------------------------------------------------
#  Iterate
#----------------------------------------------------------------------

    while abs(dHal) > 1.0E-15 :

        #print("nCub,eta,f0,f1,f2,dHal {:2d} {:12.5e} {:12.5e} {:12.5e} {:12.5e} {:12.5e}".format(nCub,eta,f0,f1,f2,dHal))

        eta  = eta - dHal

        f0 = cub0(eta,E2,E1,E0) ; f1 = cub1(eta,D2,E1) ; f2 = cub2(eta,D2)

        nCub += 1
        dNew =   f0/ f1
        dHal = dNew/(1.0 - 0.5*dNew*f2/f1)

#== Return solution ===================================================

    return eta,nCub

#========================================================================
#  Cubic Function and Derivatives; note D2 = 2*E2 in 1st & 2nd Derv
#========================================================================

def cub0(eta,E2,E1,E0) : return E0 + eta*(E1 + eta*(E2 + eta))
def cub1(eta,D2,E1   ) : return           E1 + eta*(D2 + eta*3.0)
def cub2(eta,D2      ) : return                     D2 + eta*6.0

#========================================================================
#  End of Module
#========================================================================
