
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

import numpy    as NP

import calcEOS  as CE
import calcReg  as CR
import stabTest as ST
import writeOut as WO

from math  import log,exp,sqrt

#========================================================================
#  Calculate Saturation Pressure using BFGS Stability Test
#========================================================================

def calcPsatBFGS(pEst,tRes,Z,clsEOS,clsIO) :

    nCom = clsEOS.NC
    iNeu = 0

#-- Minimum and Maximum Pressure(in psia) to consider for Psat ------    

    p2PH =    10.0
    p1PH = 15010.0

#== The Sat-Type probably depends on the mol% of C7+ ==================

    mC7P = moleFracC7P(Z,clsEOS)

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

        wilK = wilsonK(pSat,tRes,clsEOS)

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
#  BFGS-Based 2-Phase Flash Calculation
#  Ref: Martin Petitfrere, PhD Thesis
#========================================================================

def calcFlashBFGS(pRes,tRes,Z,vEst,clsEOS,clsIO) :

    nCom = clsEOS.NC
    mFLS = 101
    mEPS = clsEOS.macEPS

#== Check we have an unstable solution at (p,T) using Stab-Check ======

    #iTyp,K = CS.twoSidedStabCheck(pRes,tRes,Z,clsEOS)

    iTyp,K = stab2SidedBFGS(pRes,tRes,Z,clsEOS)

    print("iTyp,pRes ",iTyp,pRes)

#-- No, is Single Phase ---------------------------------------------    

    if iTyp == 0 :
        X = NP.copy(Z) ; Y = NP.copy(Z)
        K = NP.ones(nCom)
        V = -1.0
        return V,K,X,Y

#== Is Two-Phase, Solve Rachford-Rice =================================

    V,X,Y = CS.solveRR(nCom,Z,K,vEst)

#== Alpha or Beta Formalism? ==========================================

    L = 1.0 - V

    if L > V : qAlf = False
    else     : qAlf = True

    if qAlf : xK = calcAlfaBFGS(Z,L,X)
    else    : xK = calcBetaBFGS(Z,V,Y)

#-- Min and max values ----------------------------------------------

    #Xmin = 1.0E-08
    Xmin = 0.0
    
    xLo = [Xmin for iC in range(nCom)]
    xHi = NP.pi*NP.sqrt(Z) - Xmin

#== Initial Function (GFE) and its Gradient ===========================

    fK,gK = flashFuncDervBFGS(xK,qAlf,pRes,tRes,Z,clsEOS)

    gE = NP.dot(gK,gK)

#== Initialise the BFGS invBk matrix as Identity matrix ===============

    Bk = NP.identity(nCom)

#----------------------------------------------------------------------
#  BFGS Loop; Petitfrere Algorithm 3.3, Page 108
#----------------------------------------------------------------------

    kBFGS = 0

    while gE > 1.0E-12 :

#== First step is SS ==================================================

        if kBFGS == 0 :

            sK = NP.multiply(-1.0,gK)

#== BFGS-step =========================================================        

        else :

            pK = NP.subtract(xK,xKM1)
            yK = NP.subtract(gK,gKM1)

            sK,Bk = updateFlashBFGS(pK,yK,Bk,gK)

#== Store Previous Array and Gradient =================================

        xKM1 = NP.copy(xK)
        gKM1 = NP.copy(gK)
        fKM1 = fK

#== Test step against limits ==========================================

        sK = limitFlashBFGS(xKM1,sK,xLo,xHi)

#== Line-Search to Fine New Array and Gradient ========================

        iSS,fK,xK,gK = \
            lineSearchFlash(kBFGS,qAlf,pRes,tRes,Z,fKM1,xKM1,sK,gKM1,clsEOS)

#== Reject BFGS step and take SS instead? =============================

        """
        if fK > fKM1 :
            
            sK = NP.multiply(-1.0,gKM1)
            xK = NP.add(xKM1,sK)

            fKSS,gK = flashFuncDervBFGS(xK,pRes,tRes,Z,clsEOS)

            fK1 = fK   - fKM1
            fK2 = fKSS - fKM1

            print("fK-fKM1,fKSS-fKM1 {:10.3e} {:10.3e}".format(fK1,fK2))

            fK = fKSS
        """

#-- Convergence Criteria --------------------------------------------            

        gE = NP.dot(gK,gK)

        S2 = (fKM1 - fK)/max(abs(fK),abs(fKM1),1.0)

        print("kBFGS,pRes,gE,S2 {:2d} {:10.3f} {:10.3e} {:10.3e}".format(kBFGS,pRes,gE,S2))

#== Increment counter =================================================

        kBFGS += 1

#== Calculate Vapour Fraction, Liquid & Vapour Moles ==================

    if qAlf : V,niL,niV = calcLiqVapMoles(Z,xK)
    else    : V,niL,niV = calcVapLiqMoles(Z,xK)

#== Mole Fractions of Liquid and Vapour ===============================

    X,Y = molNumToMolFrac(V,niL,niV)

    K = NP.divide(Y,X)

#== Return values =====================================================

    return V,K,X,Y

#========================================================================
#  Check proposed step is within limits, if not, reduce step
#========================================================================

def limitFlashBFGS(xK,sK,xLo,xHi) :

    nCom = len(xK)

#== Test Updated Array Against Lo/Hi Values ===========================

    for iC in range(nCom) :

        if xK[iC] + sK[iC] < xLo[iC] : sK[iC] = xLo[iC] - xK[iC]

        if xK[iC] + sK[iC] > xHi[iC] : sK[iC] = xHi[iC] - xK[iC]

#== Scale the step ====================================================

    #print("limitFlash: stpL {:10.3e}".format(stpL))

#== Return values =====================================================

    return sK

#========================================================================
#  BFGS Flash Update
#========================================================================

def updateFlashBFGS(pK,yK,Bk,gK) :

#-- Work array (zK) and dot-products --------------------------------    

    sPY = NP.dot(pK,yK)
    zK  = NP.matmul(Bk,yK)
    sYZ = NP.dot(yK,zK)

    con1 = (sPY + sYZ)/(sPY*sPY)
    con2 =        -1.0/ sPY

#-- Equation (3.146); BFGS inverse Matrix update formula ------------

    Bk1 = con1*NP.outer(pK,pK)
    Bk2 = con2*NP.outer(zK,pK)
    Bk3 = con2*NP.outer(pK,zK)

    Bk = NP.add(Bk,Bk1)
    Bk = NP.add(Bk,Bk2)
    Bk = NP.add(Bk,Bk3)

#-- Proposed step ---------------------------------------------------

    sK = -1.0*NP.matmul(Bk,gK)

#== Return Values =====================================================    

    return sK,Bk

#========================================================================
#  Calculate the Alpha-Parameter from Liquid Mole Fraction and L = 1 - V
#  Petitfrere Equation (3.179)
#========================================================================

def calcAlfaBFGS(Z,L,X) :

#-- Liquid Moles ----------------------------------------------------

    niL = NP.multiply(L,X)

#-- Arguments of the arcsin function --------------------------------    

    sqrZ = NP.sqrt(Z)
    sqrL = NP.sqrt(niL)

    sinA = NP.divide(sqrL,sqrZ)

#-- acrsin ----------------------------------------------------------    

    arcS = NP.arcsin(sinA)

#== Required array ====================================================

    alfa = 2.0*NP.multiply(sqrZ,arcS)

#== Return Array ======================================================

    return alfa

#========================================================================
#  Calculate the Beta-Parameter from Vapour Mole Fraction and V
#  Petitfrere Equation (3.189)
#========================================================================

def calcBetaBFGS(Z,V,Y) :

#-- Liquid Fraction and Liquid Moles --------------------------------

    niV = NP.multiply(V,Y)

#-- Arguments of the arcsin function --------------------------------    

    sqrZ = NP.sqrt(Z)
    sqrV = NP.sqrt(niV)

    sinA = NP.divide(sqrV,sqrZ)

#-- acrsin ----------------------------------------------------------    

    arcS = NP.arcsin(sinA)

#== Required array ====================================================

    beta = 2.0*NP.multiply(sqrZ,arcS)

#== Return Array ======================================================

    return beta

#========================================================================
#  Function (GFE) and its Gradient for BFGS-Flash Calculation
#========================================================================

def flashFuncDervBFGS(AorB,qAlf,pRes,tRes,Z,clsEOS) :

    nCom = clsEOS.NC

    iLiq = 1 ; iVap = -1

#== Calculate Vapour Fraction, Liquid & Vapour Moles ==================

    if qAlf : sumV,niL,niV = calcLiqVapMoles(Z,AorB)
    else    : sumV,niL,niV = calcVapLiqMoles(Z,AorB)

#== Mole Fractions of Liquid and Vapour ===============================

    X,Y = molNumToMolFrac(sumV,niL,niV)

    Xmin = 1.0E-20
    
    X = NP.maximum(X,Xmin)
    Y = NP.maximum(Y,Xmin)

    logX = NP.log(X)
    logY = NP.log(Y)

#== Log Fugacity Coefficients =========================================    

    fugX,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,X,clsEOS)
    fugY,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iVap,qP,qT,qX,pRes,tRes,Y,clsEOS)

#== Liquid and Vapour GFE =============================================

    logX = NP.add(logX,fugX)
    logY = NP.add(logY,fugY)

    gfeL = NP.dot(niL,logX)
    gfeV = NP.dot(niV,logY)

#== Total (2-Phase) GFE ===============================================

    gfeT = gfeL + gfeV

#== Gradient Vector (in alfa-space) ===================================

    grad = NP.subtract(logX,logY)

    aMul = NP.multiply(niL,niV)
    aMul = NP.divide(aMul,Z)
    aMul = NP.sqrt(aMul)

    if qAlf : gMul =  1.0
    else    : gMul = -1.0

    Grad = gMul*NP.multiply(aMul,grad)  #-- Equation (3.197)

#== Return Values =====================================================

    return gfeT,Grad

#========================================================================
#  Given (Z,alfa), calculate (V,X,Y)
#========================================================================

def calcLiqVapMoles(Z,alfa) :

    nCM1 = len(Z) - 1

#== Compute the Vapour Moles from beta & Z ============================
    
    sqrZ = 2.0*NP.sqrt(Z)
    argS =     NP.divide(alfa,sqrZ)
    sinA =     NP.sin(argS)
    sin2 =     NP.multiply(sinA,sinA)

    niL  =     NP.multiply(Z,sin2)
    sumL =     NP.sum(niL)

#-- Vapour Moles ----------------------------------------------------

    sumV = 1.0 - sumL
    niV  = NP.subtract(Z,niL)

    #print("Z,sqrZ,argS,sinA,sin2,niL,niV ",Z[nCM1],sqrZ[nCM1],argS[nCM1],sinA[nCM1],sin2[nCM1],niL[nCM1],niV[nCM1])

#== Return values =====================================================
    
    return sumV,niL,niV

#========================================================================
#  Given (Z,beta), calculate (V,X,Y)
#========================================================================

def calcVapLiqMoles(Z,beta) :

#== Compute the Vapour Moles from beta & Z ============================
    
    sqrZ = 2.0*NP.sqrt(Z)
    argS =     NP.divide(beta,sqrZ)
    sinA =     NP.sin(argS)
    sin2 =     NP.multiply(sinA,sinA)

    niV  =     NP.multiply(Z,sin2)
    sumV =     NP.sum(niV)

#-- Liquid Moles ----------------------------------------------------

    sumL = 1.0 - sumV
    niL  = NP.subtract(Z,niV)

#== Return values =====================================================
    
    return sumV,niL,niV

#========================================================================
#  Moles Fractions from Mole Numbers
#========================================================================

def molNumToMolFrac(V,niL,niV) :

    L = 1.0 - V
    
    X = NP.divide(niL,L)
    Y = NP.divide(niV,V)

#== Return Values =====================================================    

    return X,Y

#========================================================================
#  BFGS Line Search
#========================================================================

def lineSearchFlash(kBFGS,qAlf,pRes,tRes,Z,fun0,aOld,dAlf,Grad,clsEOS) :

    nCom = clsEOS.NC

    aNew = [0.0 for i in range(nCom)]

#-- Woolfe Conditions -----------------------------------------------

    c1W = 1.0E-04
    c2W = 0.9

#== if grd0 > 0.0, solution will go uphil, replace by SS iteration ====

    grd0 = NP.dot(dAlf,Grad)

    if grd0 >= -1.0E-12 :
        iSS =  1
        dAlf = NP.multiply(-1.0,Grad)
        aNew = NP.add(aOld,dAlf)
        fun1,Grad = flashFuncDervBFGS(aNew,qAlf,pRes,tRes,Z,clsEOS)
        grd1 = NP.dot(Grad,dAlf)
        #print("SS, not BFGS: kBFGS,fun1,grd0,grd1 {:3d} {:10.3e} {:10.3e} {:10.3e}".format(kBFGS,fun1,grd0,grd1))
        return iSS,fun1,aNew,Grad
    else :
        iSS = -1

    #print("LS_Flash: fun0,grd0 {:10.3e} {:10.3e}".format(fun0,grd0))

    lamB = 1.0  ; lamM = 1.0E-04
    lam2 = lamB ; fun2 = fun0

#== Find Step-Length based on reduction in gStar ======================

    while lamB > lamM :

#-- Current step-length ---------------------------------------------

        for iC in range(nCom) : aNew[iC] = aOld[iC] + lamB*dAlf[iC]

#-- New GFE, etc ----------------------------------------------------
    
        fun1,Grad = flashFuncDervBFGS(aNew,qAlf,pRes,tRes,Z,clsEOS)

        grd1 = NP.dot(Grad,dAlf)
        fTst = fun1 - fun0 - c1W*lamB*grd0

        #print("lsBFGS: kBFGS,lamB,fTst,grd1 {:3d} {:10.3e} {:10.3e} {:10.3e}".\
        #      format(kBFGS,lamB,fTst,grd1))

#== What to do? =======================================================

        if fTst <= 0.0 and grd1 < c2W*abs(grd0) :

#-- Sufficient Improvement: set alpha and exit ----------------------

            break

        else :

#-- Insufficient improvement: back-tracking required ----------------

            lamT = CR.updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)

        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)

#== Return information ================================================        

    return iSS,fun1,aNew,Grad

#========================================================================
#  Two-Sided Stability Test: BFGS Version
#========================================================================

def stab2SidedBFGS(qSat,pRes,tRes,Z,clsEOS,clsIO) :

    nCom = clsEOS.NC

    qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
    qT   = False
    qX   = False

    iNeu = 0 ; iLiq = 1 ; iVap = -1

#-- Feed Log(Fugacity) Coefficients ---------------------------------    

    fugZ,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iNeu,qP,qT,qX,pRes,tRes,Z,clsEOS)

    print("Z ",Z)

    logZ = NP.log(Z)
    hVec = logZ + fugZ  #-- Vector!
    
#-- Wilson K-Values -------------------------------------------------    

    wilK = ST.wilsonK(pRes,tRes,clsEOS)        

#-- Trial Compositions ----------------------------------------------

    yVap = Z*wilK   #-- Vector
    yLiq = Z/wilK   #-- Vector
    
#-- Assume Feed is Liquid, Try to Split-Off Vapour ------------------

    print("stab2S: Vapour-Split at pRes {:10.3f}".format(pRes))

    qConV,gStrV,yVap = stab1SidedBFGS(iVap,pRes,tRes,Z,hVec,yVap,clsEOS,clsIO)

#-- Assume Feed is Liquid, Try to Split-Off Vapour ------------------

    print("stab2S: Liquid-Split at pRes {:10.3f}".format(pRes))

    qConL,gStrL,yLiq = stab1SidedBFGS(iLiq,pRes,tRes,Z,hVec,yLiq,clsEOS,clsIO)

#== Process the 1-Sided Tests =========================================

    if   qConL and qConV :

#-- If both trials unstable, take ratio unless they are the same ----
        
        iTyp = 3
        bstK = yVap/yLiq

        qSam = True
        for iC in range(nCom) :
            if abs(bstK[iC] - 1.0) > 1.0E-03 : qSam = False

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

    #print("Stab2: iTyp ",iTyp)

    bstT = 0.0

    return iTyp,bstT,bstK            

#========================================================================
#  One-Sided Stability Test BFGS
#========================================================================

def stab1SidedBFGS(iPhs,pRes,tRes,Z,hVec,yMol,clsEOS,clsIO) :

    mBFGS = 101
    nCom  = clsEOS.NC

    qConv = False ; qTriv = False

    if clsIO.Deb["STAB1"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

#-- Work Variable (alfa) --------------------------------------------

    xK = 2.0*NP.sqrt(yMol)

#== Initial Function and Gradient =====================================

    gStr,rRat,gK = stabFuncDervBFGS(pRes,tRes,Z,hVec,xK,clsEOS)

    #print("bfgsSC: iBFGS,gStr,rRat      {:2d} {:10.3e} {:10.3e}".format(iBFGS,gStr,rRat))

#== BFGS Iteration ====================================================

    iBFGS = 0
    isTru = True

    while isTru :

#--------------------------------------------------------------------
#  Calculate Step: 1st Iteration Steepest Descent
#                  Subsequently, BFGS update
#--------------------------------------------------------------------

        if iBFGS == 0 :
            sK = - gK           #-- Vector
        else :
            
            pK = xK - xKM1      #-- Vector
            yK = gK - gKM1      #-- Vector

            sK = stabStepBFGS(pK,yK,gK)

#-- L2-norm of the update dAlf --------------------------------------                

        gThs = NP.dot(sK,sK)

#-- Store 'Last' alpha & gradient -----------------------------------

        xKM1 = xK               #-- Vector
        gKM1 = gK               #-- Vector
        fKM1 = gStr             #-- Vector

#-- Check the step is physical (alfa > 0) ---------------------------

        #sK = limitStabBFGS(xKM1,sK)

#== Line Search =======================================================        

        iSS,gStr,rRat,xK,gK = \
            stabLS_BFGS(iBFGS,pRes,tRes,Z,hVec,xKM1,sK,gKM1,gStr,rRat,clsEOS)

        fK = gStr

        if qDeb > 0 :
            sOut = "bfgsSC: pRes,iPhs,iBFGS,gStr,rRat,gThs,iSS {:10.3f} {:2d} {:2d} {:10.3e} {:10.3e} {:10.3e} {:2d}\n". \
              format(pRes,iPhs,iBFGS,gStr,rRat,gThs,iSS)
            fDeb.write(sOut)

#== Convergence Test ==================================================

        if gStr < 0.0 and gThs < 1.0E-10 :   #-- Michelsen Eqn.(30)
            qConv = True
            break

#== Triviality Test ===================================================

        if gStr > 0.0 and gStr < 1.0E-03 and \
           abs(rRat-1.0) < 0.2 :                #-- Michelsen Eqn.(31)
            qTriv = True
            break

#== Nocedal Test: Petitfrere Eqn.(3.205) ==============================

        S2 = (fKM1 - fK)/max(abs(fK),abs(fKM1),1.0)

        if S2 < 1.0E-13 :
            if gStr > 0.0 : qTriv = True
            else          : qConv = True
            break

        #print("S2 {:10.3e}".format(S2))

#-- Update the Counter ----------------------------------------------        

        iBFGS += 1

        if iBFGS > mBFGS :
            print("yMol ",yMol)
            crash = 1.0/0.0

#----------------------------------------------------------------------
#  Converged, Trivial or Nowhere to go
#----------------------------------------------------------------------

    if qConv : iCon =  1
    else     : iCon = -1

    if qTriv : iTrv =  1
    else     : iTrv = -1

    #print("bfgsSC: iPhs,iBFGS,pRes,gStr,rRat,gThs,iCon,iTrv {:2d} {:2d} {:10.3f} {:10.3e} {:10.3e} {:10.3e} {:2d} {:2d}".\
    #    format(iPhs,iBFGS,pRes,gStr,rRat,gThs,iCon,iTrv))

    yMol = 0.25*NP.multiply(xK,xK)

#== Return the Gstar ==================================================

    #print("bfgsSC: qConv,qTriv,qFlat ",qConv,qTriv,qFlat)

    #print("stab1s: final-yMol ",yMol)
    
    return qConv,gStr,yMol

#========================================================================
#  Check the proposed step doesn't generate alfa[i] < 0.0
#========================================================================

def limitStabBFGS(aOld,dAlf) :

    nCom = len(aOld)
    aMin = 2.0E-10

#== Test Step =========================================================

    stpL = 1.0
    
    for iC in range(nCom) :

        if aOld[iC] + dAlf[iC] < 0.0 :
            stpT = (aMin - aOld[iC])/dAlf[iC]
            stpL = min(stpT,stpL)
            #print("iC,aOld,dAlf,stpT,stpL {:2d} {:10.3e} {:10.3e} {:10.3f} {:10.3f}".format(iC,aOld[iC],dAlf[iC],stpT,stpL))

    #print("stpL {:10.3e}".format(stpL))
    
    #dAlf = NP.multiply(stpL,dAlf)

#== Return values =====================================================

    return dAlf

#========================================================================
#  BFGS Step for Stability Test
#  Michelsen suggests inverse Hessian Bk is replaced by Identity each
#  iteration => not Bk-update and step = step(gK,pK,yK) using dot products
#========================================================================

def stabStepBFGS(pK,yK,gK) :

    nCom = len(gK)

    #sK = NP.zeros(nCom)

#-- Various Dot-Products --------------------------------------------

    sPY = NP.dot(pK,yK)
    sYY = NP.dot(yK,yK)
    sGP = NP.dot(gK,pK)
    sGY = NP.dot(gK,yK)

    con11 = (sPY+sYY)/sPY
    con12 =      sGP /sPY

    con2 = sGP/sPY
    con3 = sGY/sPY

    yCon = con2
    pCon = con3 - con11*con12

    sK = - gK + yCon*yK + pCon*pK   #-- Vector!

#== Return values =====================================================

    return sK

#========================================================================
#  BFGS Line Search for Stability Test
#========================================================================

def stabLS_BFGS(iBFGS,pRes,tRes,Z,hVec,aOld,dAlf,Grad,fun0,bet0,clsEOS) :

    nCom = clsEOS.NC

    #aNew = NP.zeros(nCom)

    grd0 = NP.dot(dAlf,Grad)

#-- Wolfe Conditions ------------------------------------------------

    c1W = 1.0E-04
    c2W = 0.9

    lamB = 1.0
    lamM = 1.0E-04
    lam2 = lamB
    fun2 = fun0

#== if grd0 > 0.0, solution will go uphil, replace by SS iteration ====    

    if grd0 >= -1.0E-12 :
        iSS =  1
        fun1,rRat,aNew,Grad = replaceBFGSbySS(pRes,tRes,Z,hVec,aOld,Grad,clsEOS)
        grd1 = NP.dot(Grad,dAlf)
        #print("lsBFGS: iBFGS,fun0,grd0,lamB,fun1,grd1,rRat {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".\
        #      format(iBFGS,fun0,grd0,lamB,fun1,grd1,rRat))
        return iSS,fun1,rRat,aNew,Grad
    else :
        iSS = -1

#== Find Step-Length based on reduction in gStar ======================

    while lamB > lamM :

#-- Current step-length ---------------------------------------------

        #for iC in range(nCom) : aNew[iC] = aOld[iC] + lamB*dAlf[iC]

        aNew = aOld + lamB*dAlf     #-- Vector

#-- New gStar, etc --------------------------------------------------        
    
        fun1,rRat,Grad = stabFuncDervBFGS(pRes,tRes,Z,hVec,aNew,clsEOS)

        grd1 = NP.dot(Grad,dAlf)

        #print("lsBFGS: iBFGS,fun0,grd0,lamB,fun1,grd1,rRat {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".\
        #      format(iBFGS,fun0,grd0,lamB,fun1,grd1,rRat))

#== What to do? =======================================================

        if fun1 <= fun0 + c1W*lamB*grd0 and grd1 < c2W*abs(grd0) :

#-- Sufficient Improvement: set alpha and exit ----------------------

            break

        else :

#-- Insufficient improvement: back-tracking required ----------------

            lamT = CR.updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)

        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)

#== Return information ================================================

    return iSS,fun1,rRat,aNew,Grad

#========================================================================
#  Replace BFGS step by SS Iteration
#========================================================================

def replaceBFGSbySS(pRes,tRes,Z,hVec,aOld,Grad,clsEOS) :

    dAlf = - Grad

    aNew = NP.add(aOld,dAlf)

#-- New gStar, etc --------------------------------------------------        
    
    fun1,rRat,Grad = stabFuncDervBFGS(pRes,tRes,Z,hVec,aNew,clsEOS)

#== Return information ================================================        

    return fun1,rRat,aNew,Grad

#========================================================================
#  BFGS Function and Gradient for Stability Check
#========================================================================

def stabFuncDervBFGS(pRes,tRes,Z,hVec,alfa,clsEOS) :

    nCom = clsEOS.NC
    iLiq = 0

    qP   = False    #-- No Pres, Temp or Comp Derivatives Needed
    qT   = False
    qX   = False
    
    aMin = 2.0E-10
    alfa = NP.maximum(alfa,aMin)

    yMol = 0.25*NP.multiply(alfa,alfa)

#-- Mole composition (normalised) -----------------------------------

    sumY = NP.sum(yMol)
    yNor = yMol/sumY        #-- Vector

    #print("stabFD: yMol[Nc] ",yMol[nCom-1])
            
    fugY,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,yNor,clsEOS)

#== Gstar and Gradient-Vector =========================================

    #gStr = 1.0 - sumY
    #beta = 0.0

    logY = NP.log(yMol)
    wrkY = yMol - Z
    dGdY = logY + fugY - hVec

    gStr = 1.0 - sumY + NP.dot(yMol,dGdY)       #-- Michelsen Eqn.(14)

    Grad = 0.5*NP.multiply(alfa,dGdY)           #-- Michelsen Eqn.(25)
    wrkY =                 wrkY*dGdY            #-- Michelsen Eqn.(23)

    beta = NP.dot(wrkY,dGdY)

    #for iC in range(nCom) :
    #    
    #    dGdY = log(yMol[iC]) + fugY[iC] - hVec[iC]
    #    
    #    gStr = gStr + yMol[iC]*dGdY             #-- Michelsen Eqn.(14)
    #    
    #    Grad[iC] = 0.5*alfa[iC]*dGdY            #-- Michelsen Eqn.(25)
    #
    #    beta = beta + (yMol[iC] - Z[iC])*dGdY   #-- Michelsen Eqn.(23)

#== Return gStar and the Gradient-Vector ==============================

    rRat =  2.0*gStr/beta                       #-- Michelsen Eqn.(24)

    return gStr,rRat,Grad

#========================================================================
#  Mole Fraction of C7+ (local copy!)
#========================================================================

def moleFracC7P(Z,clsEOS) :

    nCom = clsEOS.NC

    molF = 0.0

    for iC in range(nCom) :
        molW = clsEOS.gPP("MW",iC)
        if molW > 90.0 : molF = molF + Z[iC]

#== Return Mole Weight ================================================        

    return molF

#========================================================================
#  End of Module
#========================================================================

