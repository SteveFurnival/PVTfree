
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

import numpy      as NP

import constants  as CO

from math  import acos,copysign,cos,log,sqrt

#========================================================================
#  Calculate EoS Coefficients and its derivatives
#========================================================================

def calcEoSCoefsT(p,T,x,qT,clsEOS) :

    nCom = clsEOS.NC
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
    dI = NP.zeros(nCom)

#== Lij = 1 - Kij =====================================================

    for iC in range(nCom) :
        dI[iC] = NP.dot(clsEOS.LIJ[iC,:],cI)

    Ai = pByT2*NP.multiply(aI,dI)
    A  = pByT2*     NP.dot(cI,dI)

#----------------------------------------------------------------------
#  Temperature Derivative Terms Requested
#----------------------------------------------------------------------

    if qT :

        cT = 0.5*pByT2/sqrT
        tC = 2.0/T

        eI = NP.multiply(x,clsEOS.aQ)
        fI = NP.zeros(nCom)

        for iC in range(nCom) :
            fI[iC] = NP.dot(clsEOS.LIJ[iC,:],eI)

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

    Aij = bByT2*NP.outer(aI,aI)
    Aij =       NP.multiply(clsEOS.LIJ,Aij)  #-- Lij = 1 - Kij

#== Return the Matrix =================================================

    return Aij

#========================================================================
#  Calculate Phase Fugacity and (Optionally) its (P,T,X)-Derivatives
#========================================================================

def calcPhaseFugPTX(iPhs,qP,qT,qX,p,T,X,clsEOS) :

    nCom = clsEOS.NC
    
#== Calculate the EoS Coefficients ====================================

    A,B,Ai,Bi,dAdT,dAidT = calcEoSCoefsT(p,T,X,qT,clsEOS)

    if qX :          AiJ = calcEoSCoefsX(p,T,clsEOS)

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

    dFdE = 3.0*Eta*Eta + 2.0*E2*Eta + E1

    return Eta,dFdE

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
#  Sort the Roots of the Cubic EoS
#========================================================================

def sortCubicRoots(iLiq,E2,E1,E0,A,B,clsEOS) :

    nReal,rSmal,rLarg = solveCubic(E2,E1,E0)

    if nReal > 1 and rSmal < 0.0 :
        
        nReal = 1
        rSmal = rLarg

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
            Bn1   = clsEOS.n1*B
            Bn2   = clsEOS.n2*B
            Bm1   =           B - 1.0
            gSmal = (rSmal + Bm1) - log(rSmal) + ABn5*log((rSmal + Bn2)/(rSmal + Bn1))
            gLarg = (rLarg + Bm1) - log(rLarg) + ABn5*log((rLarg + Bn2)/(rLarg + Bn2))
            if gSmal < gLarg : rRoot = rSmal
            else             : rRoot = rLarg

    return rRoot

#========================================================================
#  Solve the Cubic x^3 + E2 x^2 + E1 x + E0 = 0
#  See Numerical Recipes Section 5.6 [E2 = a, E1 = b, E0 = c]
#  http://www.aip.de/groups/soe/local/numres/bookfpdf/f5-6.pdf
#========================================================================

def solveCubic(E2,E1,E0) :

    F2 = E2*CO.Third                  #--  a/3
    F1 = E1*CO.Third                  #--  b/3

    F22 = F2*F2                       #-- (a/3)^2

    Q = F22 - F1                      #-- Eqn.(5.6.10)
    R = F22*F2 - 1.5*F2*F1 + 0.5*E0   #-- Eqn.(5.6.10)

#-- Discriminant ----------------------------------------------------    

    D = R*R - Q*Q*Q

#== One Real Root =====================================================
    
    if   D > 0.0 :
        
        argA = abs(R) + sqrt(D)                        #-- Arg of Eqn.(5.6.13)
        A13  = pow(argA,CO.Third)
        A    = -copysign(A13,R)                        #-- Eqn.(5.6.15)
        
        if A == 0.0 : B = 0.0                          #-- Eqn.(5.6.16)
        else        : B = Q/A
        
        x1 = A + B                                     #-- Eqn.(5.6.17)
        x3 = x1

        nR = 1

#== Three Real Roots: Ignore Central Root as Unphysical in EoS ========
        
    elif D < 0.0 :
        
        sqrtQ = sqrt(Q)
        theta = acos(R/(Q*sqrtQ))                      #-- Eqn.(5.6.11)
        cosTh = cos(CO.Third*theta)
        
        x1 = -2.0*sqrtQ*cosTh                          #-- Eqn.(5.6.12-x1)
        
#-- Next expression have used cos(A-B) = cosA.cosB + sinA.sinB where
#-- cos(-2.PI/3) = -1/2, sin(-2.PI/3) = -sqrt(3)/2 and sinA = sqrt(1-cosA^2)

        x3 = sqrtQ*(cosTh+sqrt(3.0*(1.0-cosTh*cosTh)))

        nR = 3

#== Error in NR Text: see http://numerical.recipes/forum/showthread.php?t=865
        
    else :
        
        sqrtQ = sqrt(Q)
        
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
#  Calculate Key Properties and their Deriavtives
#========================================================================

def calcProps(iLiq,pRes,tRes,X,clsEOS) :

    nCom = clsEOS.NC
    qT   = False

#-- Calculate EoS Coefficients ----------------------------------------

    A,B,AiC,BiC,dumS,dumV = calcEoSCoefsT(pRes,tRes,X,qT,clsEOS)

#-- Setup and Solve the Cubic EoS (and its p-Derivative) --------------

    Eta,dFdE  = setupSolveCubicD(iLiq,A,B,clsEOS)

    dAdP,dBdP,dEdP = calcdEdP(pRes,A,B,Eta,dFdE,clsEOS)

#-- 2-Parameter (EoS) Z-Factor (and its p-derivative) & Molar Volume 

    ZF2  = Eta  + B
    Vm2  = ZF2*CO.gasCon*tRes/pRes
    dZdp = dEdP + dBdP

#-- Mole Weight and Volume Shift --------------------------------------

    cSI = CO.gasCon*NP.multiply(clsEOS.bI,clsEOS.vS)    #-- Vector!

    cS = NP.dot(X,cSI)
    Mw = NP.dot(X,clsEOS.Mw)

#-- 3-Parameter (volume shifted) Molar Volume, Density & Viscosity --

    Vm3 = Vm2 - cS
    den = Mw/Vm3

#-- Viscosity and its derivative wrt Density ------------------------    

    visC,dvdR = calcLBCderv(tRes,den,X,clsEOS)

#-- 3-Parameter Z-Factor --------------------------------------------

    ZF3 = pRes*Vm3/(CO.gasCon*tRes)

#-- Compressibility -------------------------------------------------

    comp = (Vm2/Vm3)*(1.0/pRes - dZdp/ZF2)

#-- d(density)/dp ---------------------------------------------------

    dDdp = den*comp

#-- d(viscosity)/dp = d(visc)/d(den) * d(den)/dp --------------------

    dvdp = dvdR*dDdp

#== Return Values =====================================================    

    return Mw,Vm3,den,ZF3,visC,comp,dvdp

#========================================================================
#  Lohrenz Bray Clark Viscosity [Whitson & Brule, Page 38]
#  And its derivative (wrt Density)
#========================================================================

def calcLBCderv(tRes,denS,X,clsEOS) :

    nCom = clsEOS.NC

#-- LBC Correlation Coefficients ------------------------------------    

    con0 =  0.1023
    con1 =  0.023364
    con2 =  0.058533
    con3 = -0.040758
    con4 =  0.0093324

#-- Low Pressure Component Viscosities ------------------------------

    muN = 0.0 ; muD = 0.0 ; Tpc = 0.0
    Ppc = 0.0 ; Vpc = 0.0 ; Mpc = 0.0

    for iC in range(nCom) :

        Tcrit = clsEOS.gPP("TC",iC)
        Pcrit = clsEOS.gPP("PC",iC)
        Vcrit = clsEOS.gPP("VC",iC)
        molWt = clsEOS.gPP("MW",iC)

        tRed = tRes/Tcrit

        Xi = 5.35*pow(Tcrit,CO.lbcTcX)*pow(molWt,CO.lbcMwX)*pow(Pcrit,CO.lbcPcX)

        if tRed > 1.5 : muI = 17.78E-05*pow((4.58*tRed-1.67),0.58)
        else          : muI = 34.00E-05*pow(      tRed      ,0.94)

        muI = muI/Xi

        muN = muN + X[iC]*sqrt(molWt)*muI
        muD = muD + X[iC]*sqrt(molWt)

        Tpc = Tpc + X[iC]*Tcrit
        Ppc = Ppc + X[iC]*Pcrit
        Vpc = Vpc + X[iC]*Vcrit
        Mpc = Mpc + X[iC]*molWt

    mu0 = muN/muD
    XiR = 5.35*pow(Tpc,CO.lbcTcX)*pow(molWt,CO.lbcMwX)*pow(Ppc,CO.lbcPcX)
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
#  End of Module
#========================================================================


