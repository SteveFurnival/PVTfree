
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Calculates the Two-Phase Flash
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

import calcEOS   as CE
import calcReg   as CR
import calcStab  as CS
import utilities as UT
import writeOut  as WO

#========================================================================
#  2-Phase Flash
#========================================================================

def calcFlash(pRes,tRes,Z,vEst,clsEOS,clsIO) :

    if clsIO.Deb['FLASH'] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else                      :
        qDeb = False

    nCom = clsEOS.nComp
    #mFLS = 101
    mFLS = 13           #-- SS/GDEM for max-13 steps, then BFGS
    qSat = False

    iLiq =  1 ; iVap = -1

#-- Work arrays to hold residuals/GDEM ------------------------------

    res0 = NP.zeros(nCom)
    res1 = NP.zeros(nCom)
    res2 = NP.zeros(nCom)

#-- Check if we have an unstable solution at (P,T): get K-Values ----

    iTyp,triV,K = CS.twoSidedStabTest(qSat,pRes,tRes,Z,clsEOS,clsIO)

    if qDeb :
        sOut = "calcFlash: tRes,pRes,iTyp {:10.3f} {:10.3f} {:2d}\n".format(tRes,pRes,iTyp)
        fDeb.write(sOut)
        WO.writeArrayDebug(fDeb,K,"Post-2SidedStab: K")

    if iTyp == 0 :
        #print("Fluid at (P,T) = ({:8.1f},{:8.2f}) is Single-Phase".format(pRes,tRes))
        V = -1.0 ; K = NP.ones(nCom) ; X = Z ; Y = Z
        return V,K,X,Y

#-- Will work with log(K) -------------------------------------------    

    logK = NP.log(K)

#======================================================================
#  Main Iterative Loop
#======================================================================

    iFLS =  0 ; iCut =  0 ; icSS = 0
    
    qPro = False ; qCon = False ; qTrv = False

    V = vEst

    while iFLS < mFLS :

        iFLS += 1 ; icSS += 1
    
#== Call the RR-Solver ================================================

        V,X,Y = solveRR(Z,K,V)

#== Calculate the 2-Phase Gibbs Free Energy (GFE) and dG/d(Liq Moles) =        

        gfe2,dGdL = calcGFE2(pRes,tRes,V,X,Y,clsEOS)

#== If last step was GDEM, ensure GFE has been reduced and V bounded ==

        qUPD = True

        if qPro :   #-- Last step was promoted, check required!

            if V >= 0.0 and V <= 1.0 :          #-- 0 < V < 1: Flash Physical
                if gfe2 > lstG : qCut = True    #-- but GFE has increased
                else           : qCut = False   #-- OK
            else : qCut = True                  #-- Negative Flash Teritory

            qPro = False

#----------------------------------------------------------------------
#  Cut Step [-ve Flash or GFE increased) or Just Take SS
#  In Newton (or Newton-Like) Scheme, Would Do Line Seach Here
#  but would need dG/dlnKi and we only have dG/dnLi
#  As SS+GDEM, Just Half the Proposed Step, But Only 3 Times (1/8th step)
#  Why 4?  Seems to Work!!
#----------------------------------------------------------------------

            if qCut :
                iCut += 1
                if iCut < 4 :
                    #print("Promotion Step Halved: iCut ",iCut)
                    iFLS -= 1 ; icSS -= 1
                    dlnK = 0.5*dlnK
                    logK = logK - dlnK
                    K    = NP.exp(logK)
                else :
                    logK = logK - dlnK  #-- Remove the last of the Update
                    logK = logK + res0  #-- And make conventional SS
                    K    = NP.exp(logK)
                    
                qUPD  = False

#== Updates, Tolerance, Triviality Test, if allowed ===================

        if qUPD :

            res2 = res1        #-- Store previous two steps for GDEM
            res1 = res0
            res0 = dGdL        #-- dG/dnLi = dlnKi
            
            tolR = NP.dot(res0,res0)    #-- W&B Eqn.(4.30)
            triV = NP.dot(logK,logK)    #-- W&B Eqn.(4.51)

            if tolR < 1.0E-11 :   #-- 1E-13 hard to achieve!
                qCon = True
                break

            if triV < 1.0E-04 :
                qTrv = True
                break

#-- GDEM or Regular-SS step? ----------------------------------------

            lstG = gfe2  #-- Store 'last' GFE in case GDEM makes increase

            if icSS % UT.mGDEM2 > 0 :
                qPro = False
                logK = logK + res0    #-- W&B Eqn.(4.48)
                K    = NP.exp(logK)
            else :
                qPro,logK,K,dlnK = GDEM2(tolR,triV,res0,res1,res2,logK,clsIO)
                icSS = 0
                iCut = 0

#== End of Main Loop ==================================================                

        if qDeb :
            sOut = "calcFlash: I,gfe2,tolR,triV,V {:2d},{:10.3e},{:10.3e},{:10.3e},{:10.3e}\n".format(iFLS,gfe2,tolR,triV,V)
            fDeb.write(sOut)

#== K-Values, Vapour Fraction and Liquid/Vapour Mole Fractions ========        

    K = NP.exp(logK)            #-- K-Values

    V,X,Y = solveRR(Z,K,V)      #-- Solve RR for V,X & Y
    
#======================================================================
#  Not Converged and not trivial, try to finish off with BFGS
#======================================================================

    if not qCon and not qTrv :
        V,K,X,Y = calcFlashBFGS(pRes,tRes,Z,V,X,Y,clsEOS)

#======================================================================
#  End of Routine: Return Values
#======================================================================

    return V,K,X,Y

#========================================================================
#  Calculate the 2-Phase GFE (and its gradient)
#========================================================================

def calcGFE2(pRes,tRes,V,X,Y,clsEOS) :

    iLiq =  1
    iVap = -1

    qP = False  #-- Don't need Pres, Temp or Comp LogPhi Derivatives
    qT = False
    qX = False

    xyMin = 1.0E-20

#-- Liquid Moles ----------------------------------------------------

    L = 1.0 - V

#-- Protect against zero (X,Y) --------------------------------------

    X = NP.maximum(X,xyMin)
    Y = NP.maximum(Y,xyMin)

#-- Log Mole Fractions ----------------------------------------------

    logX = NP.log(X)
    logY = NP.log(Y)

#-- Liquid and Vapour Log(Fugacity) Coefficients --------------------

    fugX,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,X,clsEOS)
    fugY,dumV,dumV,dumM = \
        CE.calcPhaseFugPTX(iVap,qP,qT,qX,pRes,tRes,Y,clsEOS)

#-- Log Mole Fractions now become 'Chemical Potential' --------------

    logX = logX + fugX  #-- Vector!
    logY = logY + fugY  #-- Vector!

#== 2-Phase GFE =======================================================

    gfeL = NP.dot(X,logX)
    gfeV = NP.dot(Y,logY)

    gfe2 = L*gfeL + V*gfeV

#== Gradient (wrt Liquid Moles) =======================================

    dGdL = logX - logY  #-- Vector!

#== Return Values =====================================================

    return gfe2,dGdL

#========================================================================
#  2-Point GDEM Acceleration for Flash Calculation
#  See Soreide PhD Thesis Section 2.5, pp. 53-57
#========================================================================

def GDEM2(tolR,triV,res0,res1,res2,logK,clsIO) :

    nCom = len(res0)
    qPro = True

    dlnK = NP.zeros(nCom)

#== Coefficients bij from Eqn. 2.53 ===================================

    b01 = NP.dot(res0,res1)
    b02 = NP.dot(res0,res2)
    b11 = NP.dot(res1,res1)
    b12 = NP.dot(res1,res2)
    b22 = NP.dot(res2,res2)

    #print("GDEM2: b01,b02,b11,b12,b22 {:10.3e},{:10.3e},{:10.3e},{:10.3e},{:10.3e}".format(b01,b02,b11,b12,b22))        

#== Soreide Eqn. 2.57 =================================================

    e01 = (b01 - b12)/b12
    e02 = (b02 - b12)/b12
    e11 = (b11 - b12)/b12
    e22 = (b22 - b12)/b12

    #print("GDEM2: e01,e02,e11,e22 {:10.3e},{:10.3e},{:10.3e},{:10.3e}".format(e01,e02,e11,e22))        

#== Soreide Eqn. 2.58 =================================================

    eNum = e11 + e22 + e11*e22

    if abs(eNum) < 1.0E-15 :  #-- Have seen eNum = 0 => regular SS

        mul0 = 1.0
        mul1 = 0.0

    else :
    
        eDiv = 1.0/eNum

        mu1 = (e02 - e01 - e22 - e01*e22)*eDiv

#== Soreide Eqn. 2.59 =================================================

        mu2 = (e01 - e02 - e11 - e02*e11)*eDiv

#== Soreide Eqn. 2.60 =================================================

        mu3 = (e11*e22 - e01*e22 - e02*e11)*eDiv

#-- Any negative multipliers, take a regular SS step

        if mu3 < 0.0 or mu3 > 1.0E+06 :
            qPro = False
            mul0 = 1.0
            mul1 = 0.0
        else :
            mul0 = 1.0/mu3
            mul1 = mu2/mu3

        if mul0 < 0.0 or mul1 < 0.0 :  #-- Regular SS
            qPro = False
            mul0 = 1.0
            mul1 = 0.0

        #print("GDEM2: mul0,mul1 {:10.3e} {:10.3e}".format(mul0,mul1))        

#-- Additional restriction if mu3 < 0.1 as mul0 > 10 ----------------

        if mu3 < 0.1 :
            qPro = False
            mul0 = 1.0
            mul1 = 0.0

#== Update K-Values with GDEM step ====================================

#
#  Note: in Soreide's Thesis, Eqn.(2.54), the update is
#
#       lnKi[n+1] = lnKi[n] + (1/mu3)*dlnKi[n] + (mu2/mu3)*dlnKi[n-1]
#
#  whereas in Whitson & Brule, Eqn.(4.49), the sign of dlnKi[n-1]
#  term is negated.  This appears to work better so is used here!!
#

    dlnK = mul0*res0 - mul1*res1   #-- Vector!
    logK =      logK +      dlnK   #-- Vector!

    qPos = NP.all(logK > 0.0)      #-- All logK's positive? 
    qNeg = NP.all(logK < 0.0)      #-- All logK's negative?

#== All Ks >1 or All Ks < 1?  Reject GDEM and use regular SS ==========

    if qPos or qNeg or logK[nCom-1] >= 0.0 :
        qPro = False
        logK = logK - dlnK - res0   #-- Vector!

#== Return K and logK =================================================

    K = NP.exp(logK)

    return qPro,logK,K,dlnK

#========================================================================
#  Rachford-Rice Solver
#  Nichita D.V. and Leibovici, C.F.
#  Fluid Phase Equilibria, 353 (2013), pp. 38-49
#
#  NB. The components must have been ordered most->least volatile
#      before entering this routine!
#========================================================================

def solveRR(Z,K,vEst) :

    N = len(Z)

    mRR = 21

    z1 = Z[0] ; zN = Z[N-1]

    tolRR  = 1.0E-15        #-- eps
    tolRR2 = tolRR*tolRR    #-- eps^2

#--------------------------------------------------------------------
#  c-Vector: c = 1/(1-Ki) - must protect against Ki = 1
#  Described in text after Eqn.(2)
#--------------------------------------------------------------------    

    wI = 1.0 - K                                    #-- wi = 1 - Ki
    
    wInot0 = NP.where(abs(wI) > tolRR2,wI,tolRR2)   #-- Protect by eps^2
    
    c = NP.divide(1.0,wInot0)                       #-- ci = 1/(1-Ki)

#-- d-Vector [see text after Eqn.(3)] -------------------------------

    c1 = c[0] ; cN = c[N-1]
    
    cDiv = 1.0/(cN - c1)            #-- Scalar

    cWrk = c1 - c                   #-- Vector
    d    = cWrk*cDiv                #-- d[i] = (c1 - c[i])/(cN - c1)

    dP1 = d + 1.0                   #-- Vector [in RR-Func & Derv]

#== Bounds for a and its function-values ==============================        

    aL = z1/(1.0 - z1)      #-- Eqn.(35)
    aR = (1.0 - zN)/zN      #-- Eqn.(38)

#== Has an estimate for V been provided? ==============================

    if vEst == None or vEst < 0.0 or vEst > 1.0 :
        a0 = bracketRR(N,Z,d,dP1,aL,aR)
    else :
        a0 = (vEst - c1)/(cN - vEst)
        
    a  = a0 ; da = 1.0

    if a < aL or a > aR :
        #print("solveRR: vEst,aL,a,aR,c1,cN {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(vEst,aL,a,aR,c1,cN))
        a = bracketRR(N,Z,d,dP1,aL,aR)

#----------------------------------------------------------------------
#  Main Iterative Loop
#----------------------------------------------------------------------

    iRR = 0

    while iRR < mRR :

        iRR += 1

#== Step-4: F(a) and F'(a) ============================================

        Fa,FPa = funcDervRR(N,Z,d,dP1,a)
        inWin  = False        

        #print("solveRR[4]: (I,a,Fa,FPa) = ({:2d},{:10.3e},{:10.3e},{:10.3e})".format(iRR,a,Fa,FPa))

#== Step-5: F'(a) > 0 =================================================

        if FPa <= 0.0 :
            
#== Step-6: Newton Update =============================================

            da   = - Fa/FPa
            aNew = a + da

            #print("solveRR[6]: (da,aNew) = ({:10.3e},{:10.3e})".format(da,aNew))

#== Step-7: Test Convergence ==========================================

            if abs(Fa) < tolRR :
                V = (c1 + a*cN)/(1.0 + a)
                break

#== Step-8: Within Bounds? ============================================

            if aNew >= aL and aNew <= aR :
                #print("solveRR[8]: aL >= aNew >= aR")
                inWin = True

#== Step-9: F(a) > 0 ==================================================

        if not inWin :

#-- Potential Error Trap --------------------------------------------

            if abs(a) < tolRR or abs(1.0+a) < tolRR :
                Vnum = c1 + a*cN
                if abs(a) < tolRR : V = Vnum
                else              : V = Vnum/tolRR
                break

#== Select the appropriate Convex function ============================            

            Ga  = (a+1.0)*Fa/a                          #-- Eqn.(18)
            GPa = Ga*(1.0/(a+1.0) - 1.0/a + FPa/Fa)     #-- Derv of above
            
            if Fa > 0.0 :
                da =  - Ga/GPa
                #print("solveRR[9t]: a,Ga,GPa,da {:10.3e} {:10.3e} {:10.3e}".format(a,Ga,GPa,da))
            else:
                Ha  = - a*Ga                            #-- Eqn.(22)
                HPa = -   Ga - a*GPa                    #-- Derv of above
                da  = - Ha/HPa
                #print("solveRR[9f]: a,Ha,HPa,da {:10.3e} {:10.3e} {:10.3e}".format(a,Ha,HPa,da))

            aNew = a + da

#-- Have seen flip-flop between Ga/Ha solutions with small-da ------

            if abs(da) < tolRR or abs(Fa) < tolRR :
                V = (c1 + a*cN)/(1.0 + a)               #-- Eqn.(5)
                break

#== Value for next iterate ============================================
                
        a = aNew

#== Calculate (X,Y) ===================================================

    #print("solveRR: iRR,V,vEst = {:2d} {:8.5f} {:8.5f}".format(iRR,V,vEst))

    Zc = Z*c        #-- Vector
    cV = c - V      #-- Vector

    X = Zc/cV       #-- Vector
    Y =  X*K        #-- Vector

    #print("solveRR: iRR ",iRR)

    XYmin = 1.0E-20         #-- if any values < 1E-20, set to 1E-20

    X = NP.maximum(X,XYmin) #-- to protect log(x) & log(y) in the
    Y = NP.maximum(Y,XYmin) #-- fugacity expressions
    
    return V,X,Y

#========================================================================
#  Rachford-Rice Function
#========================================================================

def funcRR(N,Z,d,dP1,a) :

    wI   = d + a*dP1    #-- Denominator of Eqn.(3)
    zByW = Z/wI         #-- Numerator   of Eqn.(3), less a added below

    F    = a*NP.sum(zByW)

#== Return Function ===================================================        

    return F

#========================================================================
#  Rachford-Rice Function & Derivative
#  Numpy formalism is over twice the speed of the explicit loop - ???
#========================================================================

def funcDervRR(N,Z,d,dP1,a) :

    wI   = d + a*dP1    #-- Denominator of Eqn.(3)
    zByW = Z/wI         #-- Numerator   of Eqn.(3), less a added below
    dByW = d/wI
    derV = zByW*dByW    #-- z*d/[d + a(1+d)]^2

    F    = a*NP.sum(zByW)
    FP   =   NP.sum(derV)

#== Return Function and Derivative ====================================        

    return F,FP

#========================================================================
#  Bracket Rachford-Rice Limits and estimate initial a = a0
#  See text in 3rd paragraph of Section 6
#========================================================================

def bracketRR(N,Z,d,dP1,aL,aR) :

    FaL = funcRR(N,Z,d,dP1,aL)
    FaR = funcRR(N,Z,d,dP1,aR)
    
    mG  = (FaR-FaL)/(aR-aL)
    a0  = aL - FaL/mG

    #print("solveRR[2]: (FL,FR,a0) = ({:10.3e},{:10.3e},{:10.3e})".format(FaL,FaR,a0))

    return a0

#========================================================================
#  BFGS-Based 2-Phase Flash Calculation
#  Ref: Martin Petitfrere, PhD Thesis
#========================================================================

def calcFlashBFGS(pRes,tRes,Z,V,X,Y,clsEOS) :

    nCom = clsEOS.nComp
    mFLS = 21
    mEPS = UT.macEPS

#== Alpha or Beta Formalism? ==========================================

    L = 1.0 - V

    if L > V : qAlf = False
    else     : qAlf = True

    if qAlf : xK = calcAlfaBFGS(Z,L,X)
    else    : xK = calcBetaBFGS(Z,V,Y)

#-- Min and max values ----------------------------------------------

    Xmin = 0.0
    
    xLo = NP.zeros(nCom)
    xHi = NP.pi*NP.sqrt(Z)

#== Initial Function (GFE) and its Gradient ===========================

    fK,gK = flashFuncDervBFGS(xK,qAlf,pRes,tRes,Z,clsEOS)

    gE = NP.dot(gK,gK)

#== Initialise the BFGS invBk matrix as Identity matrix ===============

    Bk = NP.identity(nCom)

#----------------------------------------------------------------------
#  BFGS Loop; Petitfrere Algorithm 3.3, Page 97
#----------------------------------------------------------------------

    kBFGS = 0

    while gE > 1.0E-12 :

#== First step is SS ==================================================

        if kBFGS == 0 :

            sK = -gK

#== BFGS-step =========================================================        

        else :

            pK = xK - xKM1
            yK = gK - gKM1

            sK,Bk = updateFlashBFGS(pK,yK,Bk,gK)

#== Store Previous Array and Gradient =================================

        xKM1 = xK
        gKM1 = gK
        fKM1 = fK

#== Test step against limits ==========================================

        sK = limitFlashBFGS(xKM1,sK,xLo,xHi)

#== Line-Search to Fine New Array and Gradient ========================

        iSS,fK,xK,gK = \
            lineSearchFlash(kBFGS,qAlf,pRes,tRes,Z,fKM1,xKM1,sK,gKM1,clsEOS)

#-- Convergence Criteria --------------------------------------------            

        gE = NP.dot(gK,gK)

        S2 = (fKM1 - fK)/max(abs(fK),abs(fKM1),1.0)

        #print("kBFGS,pRes,gE,S2 {:2d} {:10.3f} {:10.3e} {:10.3e}".format(kBFGS,pRes,gE,S2))

#== Increment counter =================================================

        kBFGS += 1

#== Calculate Vapour Fraction, Liquid & Vapour Moles ==================

    if qAlf : V,niL,niV = calcLiqVapMoles(Z,xK)
    else    : V,niL,niV = calcVapLiqMoles(Z,xK)

#== Mole Fractions of Liquid and Vapour ===============================

    X,Y = molNumToMolFrac(V,niL,niV)

    K = Y/X

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

    Bk = Bk + Bk1 + Bk2 + Bk3

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

    niL = L*X           #-- Vector

#-- Arguments of the arcsin function --------------------------------    

    sqrZ = NP.sqrt(Z)
    sqrL = NP.sqrt(niL)

    sinA = sqrL/sqrZ    #-- Vector

#-- acrsin ----------------------------------------------------------

    sinA = NP.minimum(sinA,1.0)  #-- Protect against elem = 1 + EPS

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

    niV = V*Y           #-- Vector

#-- Arguments of the arcsin function --------------------------------    

    sqrZ = NP.sqrt(Z)
    sqrV = NP.sqrt(niV)

    sinA = sqrV/sqrZ    #-- Vector

#-- acrsin ----------------------------------------------------------    

    sinA = NP.minimum(sinA,1.0)  #-- Protect against elem = 1 + EPS

    arcS = NP.arcsin(sinA)

#== Required array ====================================================

    beta = 2.0*NP.multiply(sqrZ,arcS)

#== Return Array ======================================================

    return beta

#========================================================================
#  Function (GFE) and its Gradient for BFGS-Flash Calculation
#========================================================================

def flashFuncDervBFGS(AorB,qAlf,pRes,tRes,Z,clsEOS) :

    nCom = clsEOS.nComp

    qP = False      #-- Don't need P, T or X LogPhi Derivatives
    qT = False
    qX = False

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

    fugX,dumV,dumV,dumM = CE.calcPhaseFugPTX(iLiq,qP,qT,qX,pRes,tRes,X,clsEOS)
    fugY,dumV,dumV,dumM = CE.calcPhaseFugPTX(iVap,qP,qT,qX,pRes,tRes,Y,clsEOS)

#== Liquid and Vapour GFE =============================================

    logX = logX + fugX  #-- Vector!
    logY = logY + fugY  #-- Vector!

    gfeL = NP.dot(niL,logX)
    gfeV = NP.dot(niV,logY)

#== Total (2-Phase) GFE ===============================================

    gfeT = gfeL + gfeV

#== Gradient Vector (in alfa-space) ===================================

    grad = logX - logY  #-- Vector!

    aMul = NP.sqrt(niL*niV/Z)

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
    niV  =   Z - niL    #-- Vector

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
    niL  =   Z - niV    #-- Vector

#== Return values =====================================================
    
    return sumV,niL,niV

#========================================================================
#  Moles Fractions from Mole Numbers
#========================================================================

def molNumToMolFrac(V,niL,niV) :

    L = 1.0 - V
    
    X = niL/L   #-- Vector
    Y = niV/V   #-- Vector

#== Return Values =====================================================    

    return X,Y

#========================================================================
#  BFGS Line Search
#========================================================================

def lineSearchFlash(kBFGS,qAlf,pRes,tRes,Z,fun0,aOld,dAlf,Grad,clsEOS) :

    nCom = clsEOS.nComp

    aNew = NP.zeros(nCom)

#-- Woolfe Conditions -----------------------------------------------

    c1W = 1.0E-04
    c2W = 0.9

#== if grd0 > 0.0, solution will go uphil, replace by SS iteration ====

    grd0 = NP.dot(dAlf,Grad)

    if grd0 >= -1.0E-12 :
        iSS =  1
        dAlf = -Grad            #-- Vector
        aNew =  aOld + dAlf     #-- Vector
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

        aNew = aOld + lamB*dAlf     #-- Vector

#-- New GFE, etc ----------------------------------------------------
    
        fun1,Grad = flashFuncDervBFGS(aNew,qAlf,pRes,tRes,Z,clsEOS)

        grd1 = NP.dot(Grad,dAlf)
        fTst = fun1 - fun0 - c1W*lamB*grd0

        #print("lsBFGS: kBFGS,lamB,fTst,grd1 {:3d} {:10.3e} {:10.3e} {:10.3e}".\
        #     format(kBFGS,lamB,fTst,grd1))

#== What to do?  Converged or Back-Tracking? ==========================

        if fTst <= 0.0 and grd1 < c2W*abs(grd0) : break
        else : lamT = CR.updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)

        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)

#== Return information ================================================        

    return iSS,fun1,aNew,Grad

#========================================================================
#  End of Module
#========================================================================

