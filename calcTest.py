
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see readMain.py for further details
#
#  Test Algorithms (not currently used) Module (CT)
#
#  Most of the code here is developed from the excellent PhD theis of
#  Martin Petitfrere, ED211: Laboratorie des Fluides Complexes et leurs
#                            Reservoirs
#  "EOS based simulations of thermal and compositional flows in
#   porous media"
#  Universte de Pau et des Pays de l'Adour, 12th Sept 2014
#  www.theses.fr/2014PAUU3036.pdf
#
#========================================================================

#!/usr/bin/python3

import numpy     as NP

import calcEOS   as CE
import calcFlash as CF
import calcReg   as CR
import calcStab  as CS
import utilities as UT
import writeOut  as WO

from math  import log,exp,sqrt

#========================================================================
#  Calculate Saturation Pressure using BFGS Stability Test
#========================================================================

def calcPsatBFGS(pEst,tRes,Z,clsEOS,clsIO) :

    nCom = clsEOS.nComp
    iNeu = 0

#-- Minimum and Maximum Pressure(in psia) to consider for Psat ------    

    p2PH =    10.0
    p1PH = 15010.0

#== The Sat-Type probably depends on the mol% of C7+ ==================

    mC7P = UT.moleFracC7P(Z,clsEOS)

    print("calcPsatBFGS: mC7P {:8.5f}".format(mC7P))

    if mC7P > 0.125 :   #-- Volatile oils have z(C7+) > 12.5%
        qBub = True
        iLiq =  1
        iVap = -1
    else :              #-- Else is Gas Condensate
        qBub = False
        iLiq = -1
        iVap =  1

    print("iLiq,iVap ",iLiq,iVap)

#== Initialise pSat and Liq/Vap Estimates =============================        

#----------------------------------------------------------------------
#  Pre-Sweep using trial Liquid using BFGS Stability Test
#----------------------------------------------------------------------

    qLOK = False
    qVOK = False
    nBel = 0

    pBel = []
    gBel = []

    if pEst == None : pSat = (2.0*p2PH + p1PH)/3.0  #-- Bias towards low-P
    else            : pSat =      pEst

    while (p1PH - p2PH) > 100.0 :

        print("p2PH,pSat,p1PH {:10.3f} {:10.3f} {:10.3f}".format(p2PH,pSat,p1PH))

#-- Wilson-K Values at current (p,T) --------------------------------

        wilK = UT.wilsonK(pSat,tRes,clsEOS)

        if qLOK : yLiq = NP.copy(yLOK)
        else    : yLiq = NP.divide(Z,wilK)
            
        if qVOK : yVap = NP.copy(yVOK)
        else    : yVap = NP.multiply(Z,wilK)
            
#-- Feed Log Fugacity Coefficients ----------------------------------

        fugZ,dumV,dumV,dumM = \
            CE.calcPhaseFugPTX(iNeu,qP,qT,qX,pSat,tRes,Z,clsEOS)

#-- Try to split a Liquid off what is assumed to be Vapour ----------        

        qConL,gStrL,yLiq = \
            stabCheckBFGS(iLiq,pSat,tRes,Z,fugZ,yLiq,clsEOS)

        if qConL :
            qLOK = True
            yLOK = NP.copy(yLiq)

#-- Try to split a Vapour off what is assumed to be Liquid ----------        

        qConV,gStrV,yVap = \
            stabCheckBFGS(iVap,pSat,tRes,Z,fugZ,yVap,clsEOS)

        if qConV :
            qVOK = True
            yVOK = NP.copy(yVap)

#-- Did either trial converge? --------------------------------------

        if qConL or qConV :
            p2PH = pSat
            pBel.append(pSat)
            if   qConL and qConV : gStr = min(gStrL,gStrV)
            elif qConL           : gStr =     gStrL
            elif qConV           : gStr =           gStrV
            gBel.append(gStr)
            nBel += 1
        else              :
            p1PH = pSat

#-- Two or more solutions below gStr = 0 => pSat --------------------

        if nBel >= 2 :
            dPrs = pBel[nBel-1] - pBel[nBel-2]
            dGst = gBel[nBel-1] - gBel[nBel-2]
            grad = dGst/dPrs
            intr = gBel[nBel-1] - pBel[nBel-1]*grad
            pSat = -intr/grad
        else :
            pSat = 0.5*(p2PH + p1PH)
        
    return qBub,pSat

#========================================================================
#  Test 'Numerical' Cubic Solver
#========================================================================

def numCubicSolver(E2,E1,E0,eta) :

    D2  = 2.0*E2
    E23 =     E2*UT.Third

    E232 = E23*E23
    E233 = E23*E232

    NR  = 0

#-- Symmetry Point ['N' in Nickalls Fig.1] --------------------------    

    eN = - E23
    yN = 2.0*E233 - E23*E1 + E0

    dN  = cub1(eN,D2,E1)
    
#== Two Turning Points (+ve) or Point of Inflection (-ve) =============    
    
    del2 = E232 - UT.Third*E1

    if del2 < 0.0 :

#-- Can we use 2nd Derv to estimate root ----------------------------

        disC = E1*E1 - 4.0*E0*abs(E2)
        
        est = abs(0.5*(-E1+sqrt(disC))/E2)

        if est < eN : est = - E0*eN/(yN-E0)

        lst,nC = cubicNewton(est,E2,E1,E0)

        #print("del2<0: eN,yN,dN,E0,est,lst,nC {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:2d}".format(eN,yN,dN,E0,est,lst,nC))

        NR += 1

    else :

#== Two-Turning Points ================================================

        del1 = sqrt(del2)

        lamb =  UT.Sqr3*del1
        h    = 2.0*del2*del1

        eTP1 = -E23 - del1
        eTP2 = -E23 + del1

        if   eTP1 > 0.0 :
            y1 = cub0(eTP1,E2,E1,E0)
            y2 = y1 - 2.0*h
            if y1 > 0.0 :
                m = (y1-E0)/eTP1
                C =  y1 - m*eTP1
                est = -C/m
                dif = est - eta
                #print("T1a: y1,est,dif {:10.3e} {:10.3e} {:10.3e}".format(y1,est,dif))
                NR += 1
            if y2 < 0.0 :
                NR += 1
                #print("T1b: y2 {:10.3e}".format(y2))
        elif eTP2 < 0.0 :
            if lamb < -eN :
                disC = E1*E1 - 4.0*E2*E0
                est  = 0.5*(-E1+sqrt(disC))/E2
                dif = est - eta
                #print("3a: TP2,eta,est,dif {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(TP2,eta,est,dif))
            else :
                est = sqrt(-E0/E2)
                dif = est - eta
                #print("3b: TP2,eta,est,dif {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(TP2,eta,est,dif))
        else :
            y1 = cub0(eTP1,E2,E1,E0)
            eA = eN + lamb
            yA = yN
            eB = eTP1 + 3.0*del1
            yB = y1
            m  = (yB-yA)/(eB-eA)
            C  = yB - m*eB
            est = - C/m
            dif = est - eta
            #print("2 :eA,yA,eB,yB,E0,E1,dif {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(eA,yA,eB,yB,E0,E1,dif))

#== No Return (For Now) ===============================================

    #if NR > 1 : print("NR>1: E2,E1,E0,del2 {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(E2,E1,E0,del2))

#========================================================================
#  End of Module
#========================================================================

