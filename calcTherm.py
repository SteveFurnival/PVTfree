
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculate Thermal Properties (Enthalpy, Entropy, ...)
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

import numpy        as NP
import numpy.linalg as LA

import calcEOS   as CE
import constants as CO

from math  import log,exp,sqrt

#========================================================================
#  Phase Enthalpy (and Specific Heat)
#========================================================================

def calcEnthalpy(iLiq,pRes,tRes,X,clsEOS) :

#== Ideal Gas Contributions ===========================================

    CpIdl,HIdl = calcIdealGasTherm(tRes,X,clsEOS)

    #print("CpIdl,HIdl  {:10.3e} {:10.3e}".format(CpIdl,HIdl))

    A,B,dAdT,dBdT,HD,dHDdT = calcEoSCoefsT(pRes,tRes,X,clsEOS)

    eta,dfde = setupSolveCubicD(iLiq,A,B,clsEOS)

    #print("eta ,defde {:10.3e} {:10.3e}".format(eta,dfde))

#-- d[eta]/dT -------------------------------------------------------

    eta2 = eta*eta
    bWrk = 2.0*clsEOS.n4*B
    
    etaA = - eta/dfde
    etaB = - (clsEOS.n3*eta2 + (bWrk - clsEOS.n3)*eta - bWrk)/dfde

    dedT = etaA*dAdT + etaB*dBdT

#== Assemble Residual Enthalpy ========================================

    sig1 = eta + clsEOS.n1*B
    sig2 = eta + clsEOS.n2*B

    C2 = eta + B - 1.0
    C3 = clsEOS.n5*A/B
    C4 = log(sig2/sig1)

    HRes  = CO.btuCon*tRes*(C2 + C3*HD*C4)

#== Terms for the Residual Specific Heat ==============================

    dC2dT = dedT + dBdT
    dC3dT = C3*(dAdT/A + 1.0/tRes)
    dC4dT = (dedT + clsEOS.n2*dBdT)/sig2 - (dedT + clsEOS.n1*dBdT)/sig1

    CpRes =      HRes/tRes + \
            CO.btuCon*tRes*(dC2dT + C3*HD*dC4dT + C3*C4*dHDdT + HD*C4*dC3dT)

    #print("Hres ,CpRes {:10.3e} {:10.3e}".format(HRes,CpRes))

#== Ideal + Residual ================================================

    HTot  = HIdl  + HRes
    CpTot = CpIdl + CpRes

    #print("HTot ,CpTot {:10.3e} {:10.3e}".format(HTot,CpTot))

    return HTot,CpTot

#========================================================================
#  Calculate Ideal Gas Specific Heat and Enthalpy
#========================================================================

def calcIdealGasTherm(tRes,X,clsEOS) :

    nCom = clsEOS.NC
    Tref = clsEOS.Tref

    T1 = tRes  ; T2 = T1*T1 ; T3 = T1*T2 ; T4 = T1*T3
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

        CpI = CpI + X[iC]*CpT
        HId = HId + X[iC]*HIT        

    return CpI,HId

#========================================================================
#  Calculate EoS Coefficients for Temperature-Dependent Applications
#========================================================================

def calcEoSCoefsT(pRes,tRes,X,clsEOS) :

    nCom = clsEOS.NC

#-- Work Arrays -----------------------------------------------------    

    aS = NP.zeros(nCom)
    bE = NP.zeros(nCom)

    xP = NP.zeros(nCom)
    xQ = NP.zeros(nCom)

    OmgA0 = clsEOS.OA
    OmgB0 = clsEOS.OB
    
    for iC in range(nCom) :

        Tc = clsEOS.gPP("TC",iC)
        Pc = clsEOS.gPP("PC",iC)
        AF = clsEOS.gPP("AF",iC)
        MA = clsEOS.gPP("MA",iC)
        MB = clsEOS.gPP("MB",iC)

        mAcF = calcMAcF(AF,clsEOS)

        aS[iC] = sqrt(OmgA0*MA/Pc)*Tc
        bE[iC] =      OmgB0*MB*Tc/ Pc

        xP[iC] = X[iC]*(1.0 + mAcF)
        xQ[iC] = X[iC]*(    - mAcF/sqrt(Tc))

#-- Include BIPs ----------------------------------------------------

    aPP = 0.0 ; aMX = 0.0 ; aQQ = 0.0 ; b   = 0.0

    aP = NP.zeros(nCom)
    aQ = NP.zeros(nCom)

    for iC in range(nCom) :

        aP[iC] = 0.0 ; aQ[iC] = 0.0
        
        for jC in range(nCom) :

            aij = aS[iC]*aS[jC]*(1.0 - clsEOS.gIJ(iC,jC))

            aP[iC] = aP[iC] + xP[jC]*aij
            aQ[iC] = aQ[iC] + xQ[jC]*aij

        aPP = aPP + xP[iC]*aP[iC]
        aMX = aMX + xP[iC]*aQ[iC]
        aQQ = aQQ + xQ[iC]*aQ[iC]

        b   = b   +  X[iC]*bE[iC]

#== a, da/dT & d2a/dT2 ================================================

    sqT = sqrt(tRes)

    a      = aPP + 2.0*aMX* sqT + aQQ*tRes
    dadT   =           aMX/ sqT + aQQ
    d2adT2 =     - 0.5*aMX/(sqT*tRes)

    HD   = tRes*dadT/a - 1.0

    dHDdT = dadT*(1.0 - tRes*dadT/a)/a + tRes*d2adT2/a

#== A & B =============================================================

    A = a*pRes/(tRes*tRes)
    B = b*pRes/ tRes

    #print("thermal: A    ,B     {:10.3e} {:10.3e}".format(A,B))

#== Derivatives =======================================================

    dAdT   = A*(  dadT/a - 2.0/tRes)
    d2AdT2 = A*(d2adT2/a - 4.0*dadT/(a*tRes) + 6.0/(tRes*tRes))

    dBdT   = -B/tRes

    #print("thermal: dAdT ,dBdT  {:10.3e} {:10.3e}".format(dAdT,dBdT))

#========================================================================
#  End of Routine
#========================================================================

    return A,B,dAdT,dBdT,HD,dHDdT

#========================================================================
#  m(Acf) Expansion
#========================================================================

def calcMAcF(AF,clsEOS) :

    mAcF = clsEOS.w0 + AF*( clsEOS.w1 + AF*( clsEOS.w2 + AF*clsEOS.w3 ) )

#== Return values =====================================================                

    return mAcF

#========================================================================
#  End of Module
#========================================================================


