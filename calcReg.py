#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Regression Module
#  
#========================================================================

#!/usr/bin/python3

from copy  import deepcopy
from numpy import linalg   as LA
from math  import log,exp,sqrt

import numpy     as NP

import calcGrad  as CG
import calcProps as CP
import calcExps  as CX
import plotBatch as PB
import utilities as UT
import writeOut  as WO
import writeReg  as WR

#========================================================================
#  Regression Driver
#========================================================================

def runRegression(mIter,qExp,clsEOS0,dicSAM0,dicEXP0,dicREG,clsUNI,clsIO) :

#== Set Output File Name [rootName.reg] ===============================

    if not clsIO.qReg :
        pathReg = clsIO.patR + ".reg"
        fReg    = open(pathReg,'w')
        clsIO.setQREG(True)
        clsIO.setFREG(fReg)
    else :
        fReg = clsIO.fReg

    sCom = "--"

    WO.outputHeader(fReg,sCom,clsIO)
    WO.outputEOS(fReg,sCom,clsIO,clsEOS0)

#== Number of Variables ===============================================

    nVar = len(dicREG)

#-- Initial (normalised) variables ----------------------------------

    xVec = NP.ones(nVar)

#--------------------------------------------------------------------
#  Run the experiments with the initial EoS and set of sample info
#--------------------------------------------------------------------

    qDif = False

    runRegExps(qExp,qDif,clsEOS0,dicSAM0,dicEXP0,clsIO)

#== Initial residuals vector ==========================================

    qWrt = True

    WR.writeInitResHead(clsIO)
    
    ssq0,regI = calcResSSQ(qWrt,dicEXP0,qExp,clsIO)

    nReg = len(regI)  #-- Number of items in Residual Vector

    iTer = 0
    #print("iVar,ssqI {:2d}{:10.5f}".format(iTer,ssq0))

#== Take a (deep) copy of the EoS class and Samples dictionary ========
    
    clsEOS1 = deepcopy(clsEOS0)
    dicSAM1 = deepcopy(dicSAM0)
    dicEXP1 = deepcopy(dicEXP0)

#========================================================================
#  Main Loop
#========================================================================

    qConv = False
    nIter = 1

    while not qConv :

#======================================================================
#  Build Jacobian by perturbing each variable in turn
#======================================================================

        jacO = calcJaco(clsIO,dicREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,dicEXP1,qExp,xVec,regI)

#======================================================================
#  Calculate the Gradient and Hessian
#======================================================================

        Grad = NP.dot(jacO,regI)
        Hess = NP.dot(jacO,jacO.T)  #-- Note Transpose!

#======================================================================
#  Calculate the update vector delX
#======================================================================

        delX = rotDisc(nIter,Grad,Hess,dicREG,dicSAM1,clsEOS1,clsIO)

        #WR.writeVector("Grad   ",Grad)
        #WR.writeVector("xVec[B]",xVec)
        #WR.writeVector("delX   ",delX)

#======================================================================
#  Check if Variables Bounded?
#======================================================================

        delX,Grad = checkBounds(nIter,clsEOS0,dicREG,dicSAM0,xVec,delX,Grad,Hess,clsIO)

#======================================================================
#   Line Search
#======================================================================

        ssq1,xVec,regI = lineSearch(clsIO,dicREG,dicSAM0,dicSAM1,clsEOS0,clsEOS1,dicEXP1,qExp,ssq0,xVec,delX,Grad)

        #WR.writeVector("xVec[A]",xVec)

#======================================================================
#  Progress in SSQ?
#======================================================================

        if abs((ssq1 - ssq0)/ssq1) < 1.0E-04 :
            print("Regression converged!")
            qConv = True
        else :
            ssq0 = ssq1

        nIter += 1

        if nIter > mIter :
            print("Regression - Reached Max-Iterations = ",mIter)
            break

#== Plot Before/After Results =========================================        

    qReg = True
    PB.regPlots(clsIO,dicEXP0,dicEXP1,dicSAM1,qReg,clsUNI)

#======================================================================
#  (Deep) copy the working EoS/Samples back to the main EoS/Samples
#======================================================================
    
    clsEOS0 = deepcopy(clsEOS1)
    dicSAM0 = deepcopy(dicSAM1)
    dicEXP0 = deepcopy(dicEXP1)

#======================================================================
#  Write Final Set of Residuals
#======================================================================

    WR.writeFinalResHead(clsIO)
    
    ssq0,regI = calcResSSQ(qWrt,dicEXP0,qExp,clsIO)

#======================================================================
#  Write New EoS & Experiments to Print File
#======================================================================

    iERR,sTyp = regVarType(dicREG)

    if   sTyp == "REG" : sExt = "Regular EOS Variables"
    elif sTyp == "PLS" : sExt = "Plus Fraction Variables"
    elif sTyp == "LBC" : sExt = "LBC Viscosity Variables"

    sTit = "Regression Using " + sExt
    
    WO.outputProps(sTit,clsEOS0,dicSAM0,clsIO)

#-- Only Print the Regession "Active" Experiments [qExp] ------------    

    nExp = len(dicEXP0)
    qSav = [True for i in range(nExp)]

    for iExp in range(nExp) :
        qSav[iExp]    = qExp[iExp]
        clsEXP0       = dicEXP0[iExp]
        clsEXP0.IsAct = qExp[iExp]

#-- Print the Experiments -------------------------------------------        
    
    WO.outputExps(clsEOS0,dicSAM0,dicEXP0,clsUNI,clsIO)
    
#-- Restore the IsAct Flag ------------------------------------------    

    for iExp in range(nExp) :
        clsEXP0       = dicEXP0[iExp]
        clsEXP0.IsAct = qSav[iExp]
    
    sTit = "Regression Using " + sTyp + " Variable Types"
    WO.outputSave(sTit,clsEOS0,dicSAM0,clsIO)

#== Classes/Dictionaries Impicitly Update => No Return Arguments ======    

    return clsEOS0,dicSAM0,dicEXP0

    #return

#========================================================================
#  Are variables currently bounded, or will they become bounded?
#========================================================================

def checkBounds(nIter,clsEOS,dicREG,dicSAM,xVec,delX,Grad,Hess,clsIO) :

    nVar = len(dicREG)
    stpM = 0.0
    xTst = 1.0E-06

    WR.outputRegVars(clsIO,clsEOS,dicREG,dicSAM,xVec,delX,nIter)

    #WR.writeVector("delX[BoB]",delX)

#--------------------------------------------------------------------
#  Already bounded variables?
#--------------------------------------------------------------------

    nBnd,delX,Grad = existBounded(nIter,dicREG,dicSAM,clsEOS,xVec,delX,Grad,Hess,clsIO)

    if nBnd >= 0 :
        
#--------------------------------------------------------------------
#  Check Variables that will go bounded
#--------------------------------------------------------------------

        delX = limitBounds(dicREG,xVec,delX)

#====================================================================== 
#  End of Routine
#======================================================================

    #WR.writeVector("delX[BoA]",delX)

    return delX,Grad

#====================================================================== 
#  Check for variables that are already bounded
#======================================================================

def existBounded(nIter,dicREG,dicSAM,clsEOS,xVec,delX,Grad,Hess,clsIO) :

    nVar = len(dicREG)
    nBnd = 0
    xTst = 1.0E-06

    iBnd = [0 for i in range(nVar)]

    for iVar in range(nVar) :

        clsREG = dicREG[iVar]

        xMin = clsREG.vMin
        xMax = clsREG.vMax

        #print("iVar,xVec,delX,xMin,xMax {:2d} {:8.5f} {:8.5f} {:8.5f} {:8.5f}".format(iVar,xVec[iVar],delX[iVar],xMin,xMax))

#== Are any of the variables already bounded? =========================        

        if   delX[iVar] >= 0.0 and abs(xVec[iVar] - xMax) < xTst :
            iBnd[iVar] =  1
            nBnd      +=  1
        elif delX[iVar] <= 0.0 and abs(xVec[iVar] - xMin) < xTst :
            iBnd[iVar] = -1
            nBnd      +=  1
        else :
            iBnd[iVar] =  0

#--------------------------------------------------------------------
#  All Variables Bounded?
#--------------------------------------------------------------------

    if nBnd == nVar :
        print("All Variables Bounded - Stop Regression")
        nBnd = -1
        return nBnd,delX,Grad

#--------------------------------------------------------------------
#  Zero Rows & Columns (except Diag=1) for Bounded Variables
#--------------------------------------------------------------------

    if nBnd > 0 :

        #print("nBnd,iBnd ",nBnd,iBnd)

        for iVar in range(nVar) :
            if iBnd[iVar] != 0 :
                Grad[iVar] = 0.0
                for jVar in range(nVar) :
                    Hess[iVar][jVar] = 0.0
                    Hess[jVar][iVar] = 0.0
                Hess[iVar][iVar] = 1.0

        #WR.writeMatrix("Hess ",Hess)

        nIter = -1
        delX = rotDisc(nIter,Grad,Hess,dicREG,dicSAM,clsEOS,clsIO)

        #WR.writeVector("delX[C]",delX)

#====================================================================== 
#  End of Routine
#======================================================================

    return nBnd,delX,Grad

#====================================================================== 
#  Check for variables that will become bounded by proposed (full) step
#======================================================================

def limitBounds(dicREG,xVec,delX) :

    nVar = len(dicREG)

    xTst = 1.0E-06
    stpM = 1.0/xTst

    for iVar in range(nVar) :

        clsREG = dicREG[iVar]

        xMin = clsREG.vMin
        xMax = clsREG.vMax

        xPro = xVec[iVar] + delX[iVar]  #-- Proposed change in variable

        #print("iV,xV,dX,Mn,xP,Mx {:2d} {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(iVar,xVec[iVar],delX[iVar],xMin,xPro,xMax))

        if   xPro > xMax : stpT = (xMax - xVec[iVar])/delX[iVar]
        elif xPro < xMin : stpT = (xMin - xVec[iVar])/delX[iVar]
        else             : stpT = 1.0

        #print("iVar,stpT,stpM {:2d} {:10.5f} {:10.3e}".format(iVar,stpT,stpM))

        stpM = min(stpM,stpT)
        
#====================================================================== 
#  End of Routine
#  If step < 1, scale-back the delX so alfa[0] = 1 in lineSearch
#======================================================================

    if stpM < 1.0 : delX = stpM*delX    #-- Vector!

    #print("newlyBounded: stpM {:10.5f}".format(stpM))

    return delX

#========================================================================
#  Rotational Discrimination
#  Soreide PhD Thesis, Section 4.2, pp.108-110
#========================================================================

def rotDisc(nIter,Grad,Hess,dicREG,dicSAM,clsEOS,clsIO) :

#== Digaonalise the Matrix, i.e. H = ST*D*S ===========================

    detM,iOrd,eVal,eVec = diagHess(Hess)  #-- Eqn.(4.13)

#== Some diagnostics to help users ==================================

    if dicREG != None and nIter > 0 :
        WR.diagnoseHess(nIter,detM,iOrd,eVal,eVec,dicREG,dicSAM,clsEOS,clsIO)

#== Calculate the p-Vector ==========================================

    pVec = -NP.dot(eVec.T,Grad)            #-- Eqn.(4.17): Note Tranpose!

#== Calculate the y-Vector ==========================================

    yVec = calcYvec(pVec,eVal,clsIO)      #-- Eqn.(4.18) & Algorithm 4-6

#== Calculate the [re-ordered] step =================================

    delX = NP.dot(eVec,yVec)              #-- Eqn.(4.19)

#======================================================================
#  End of Routine
#======================================================================

    return delX

#======================================================================
#  Build the Jacobian Matrix
#======================================================================

def calcJaco(clsIO,dicREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,dicEXP,qExp,xVec,regI) :

    nVar = len(dicREG)
    nReg = len(regI)
    qWrt = False

    dVar = 1.0E-04

#== Create Perturbed-Residual Jacobian, Gradient and Hessian Arrays ===

    regP = NP.zeros( nReg)
    jacO = NP.zeros((nVar,nReg))

#======================================================================
#  Loop over variables
#======================================================================

    qDif = True

    for iVar in range(nVar) :

        clsREG = dicREG[iVar]  #-- Current regression variable class

#-- Save the current value and perturb variable ---------------------        

        xSav = xVec[iVar]
        dPer = dVar*xVec[iVar]
        xPer = xVec[iVar] + dPer

        perturbVariable(xPer,clsREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,clsIO)

#--------------------------------------------------------------------
#  Run the experiments with the perturbed EoS
#--------------------------------------------------------------------

        runRegExps(qExp,qDif,clsEOS1,dicSAM1,dicEXP,clsIO)

#-- Calculate the Perturbed residuals and SSQ -----------------------

        ssqP,regP = calcResSSQ(qWrt,dicEXP,qExp,clsIO)

        #print("iVar,ssqP {:2d}{:10.5f}".format(iVar+1,ssqP))

#--------------------------------------------------------------------
#  Build the current variable's column in the Jacobian
#--------------------------------------------------------------------

        for iReg in range(nReg) :
            jacO[iVar][iReg] = (regP[iReg] - regI[iReg])/dPer

#-- Restore the current value ---------------------------------------

        perturbVariable(xSav,clsREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,clsIO)

#======================================================================
#  Return the Jacobian Matrix
#======================================================================

    return jacO

#======================================================================
#  Calculate y-Vector [Soreide Eqn.(4.18)]
#======================================================================

def calcYvec(pVec,dVec,clsIO) :

    #sOut = "calcYvec: Enter\n"
    #fDeb.write(sOut)

    nVec = len(pVec)
    
    yH   =   0.2
    Acon = 100.0

    bTst = False
    cTst = False
    jVec = nVec
    lVec = nVec

    expA = -0.5/log(Acon)  #--Constant in Soreide Eqn.(4.22)

    yVec = NP.zeros(nVec)

    for iVec in range(nVec) :
        
        if abs(dVec[iVec]) < UT.macEPS :
            yVec[iVec] = 0.0
        else :
            yVec[iVec] = pVec[iVec]/dVec[iVec]

        #sOut = "iV,dV,pV,yV {:2d} {:10.3e} {:10.3e} {:10.3e}\n".format(iVec,dVec[iVec],pVec[iVec],yVec[iVec])
        #fDeb.write(sOut)

        if not bTst :
            if dVec[iVec]/dVec[0] < 1.0E-08 or dVec[iVec] < 0.0 :
                bTst = True
                jVec = iVec

        if not cTst :
            if abs(yVec[iVec]) > yH :
                cTst = True
                lVec = iVec
                dlV  = dVec[iVec]
                plV  = pVec[iVec]
                ylV  = yH*UT.signNumber(plV)

    #sOut = "calcYvec: bTst,cTst,jVec,lVec {:10.3e} {:10.3e} {:2d} {:2d}\n".format(bTst,cTst,jVec,lVec)
    #fDeb.write(sOut)

    if bTst :
        for iVec in range(jVec,nVec) :
            yVec[iVec] = 0.0

    if cTst :
        for iVec in range(lVec,jVec) :
            rI = 1.0 - exp(expA*log(dlV/dVec[iVec]))
            sI = pow(Acon,rI)
            yVec[iVec] = ylV*sI*pVec[iVec]/plV
            #sOut = "iV,rI,sI,yV {:2d} {:10.3e} {:10.3e} {:10.3e}\n".format(iVec,rI,sI,yVec[iVec])
            #fDeb.write(sOut)

    #WR.writeVector("yVec",yVec)

    #sOut = "calcYvec: Exit\n"
    #fDeb.write(sOut)

#== Return yVec =======================================================

    return yVec

#======================================================================
#  Digaonalise the Matrix, i.e. H = ST*D*S
#======================================================================

def diagHess(Hess) :

    nVec = len(Hess)

    eVal = NP.zeros( nVec)
    eVec = NP.zeros((nVec,nVec))
    iOrd = NP.zeros( nVec,dtype='i4')

    eVal,eVec = LA.eig(Hess)

#-- Ensure any complex numbers (because of round-off) are ignored ---

    eVal = NP.real(eVal)
    eVec = NP.real(eVec)

#-- Determinant of the Hessian --------------------------------------

    detM = LA.det(Hess)

#-- And order them [large->small] -----------------------------------

    iOrd = eVal.argsort()[::-1]

    eVal = eVal[  iOrd]
    eVec = eVec[:,iOrd]

    #print("iOrd ",iOrd)
    #WR.writeVector("eVal",eVal)
    #WR.writeMatrix("eVec",eVec)

#== Return Determinant, Eigenvalues (eVal) and Eigenvectors (eVec) ====

    return detM,iOrd,eVal,eVec    

#========================================================================
#  Perturb a Variable in the EoS
#========================================================================

def perturbVariable(vMul,clsREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,clsIO) :

    sReg = {'TCRIT':'TC','PCRIT':'PC','ACENF':'AF','SHIFT':'SS',
            'MULTA':'MA','MULTB':'MB','VCRIT':'VC'}
            
    vNam = clsREG.vNam
    vTyp = clsREG.vTyp

    if   vTyp == "REG" or vTyp == "LBC" :
        
        sCod = sReg[vNam]
        nGrp = len(clsREG.vGrp)
        
        for iGrp in range(nGrp) :
            iC   = clsREG.vGrp[iGrp]
            vBas = clsEOS0.gPP(sCod,iC)
            vPer = vBas*vMul
            clsEOS1.sPP(sCod,iC,vPer)
            
    elif vTyp == "PLS" :
        
        iSam = clsREG.vSam
        
        if   vNam == "MPLUS" :
            vBas = dicSAM0[iSam].iPlsMW
            vPer = vBas*vMul
            dicSAM1[iSam].iPlsMW = vPer
        elif vNam == "SPLUS" :
            vBas = dicSAM0[iSam].iPlsSG
            vPer = vBas*vMul
            dicSAM1[iSam].iPlsSG = vPer
        elif vNam == "APLUS" :
            vBas = dicSAM0[iSam].iPlsAL
            vPer = vBas*vMul
            dicSAM1[iSam].iPlsAL = vPer

        sumPlusFrac(dicSAM1,clsEOS1)
        iERR = CP.splitPlus(dicSAM1,clsEOS1,clsIO)

    elif vTyp == "KIJ" :
        
        iOne = clsREG.vOne
        nGrp = len(clsREG.vGrp)
        
        for iGrp in range(nGrp) :
            jC   = clsREG.vGrp[iGrp]
            vBas = clsEOS0.gIJ(iOne,jC)
            vPer = vBas*vMul
            clsEOS1.sIJ(iOne,jC,vPer)
            clsEOS1.sIJ(jC,iOne,vPer)

#======================================================================
#  End of Routine
#======================================================================

    return

#========================================================================
#  Sum Pseudo-Component Moles Back to Plus Fraction (Comp nUser)
#========================================================================

def sumPlusFrac(dicSAM,clsEOS) :

    nSamp = len(dicSAM)
    nUser = clsEOS.nUser
    nComp = clsEOS.nComp

    for iSamp in range(nSamp) :
        sPlus = 0.0
        for iC in range(nUser-1,nComp) : sPlus = sPlus + dicSAM[iSamp].gZI(iC)
        dicSAM[iSamp].sZI(nUser-1,sPlus)
        for iC in range(nUser,nComp)   : dicSAM[iSamp].sZI(iC,0.0)            

#========================================================================
#  End of Routine
#========================================================================

    return 

#========================================================================
#  Line Search [Numerical Recipes Chapter 9.7]
#========================================================================

def lineSearch(clsIO,dicREG,dicSAM0,dicSAM1,clsEOS0,clsEOS1,dicEXP,qExp,fun0,xVec,delX,Grad) :

    nVar = len(delX)
    qWrt = False

    grd0 = NP.dot(delX,Grad)

    if grd0 >= 0.0 :
        print("lineSearch: grd0 >= 0")
        fun1 = fun0
        regI = [0.0 for i in range(1)]
        return fun1,xVec,regI

    lamB = 1.0      #-- Default first   step-length [Full Newton!]
    lamM = 1.0E-04  #-- Minimum allowed step-length
    lam2 = lamB
    fun2 = fun0

    qDif = False

    while lamB > lamM :

#-- Update variables ------------------------------------------------

        for iVar in range(nVar) :
            clsREG   =   dicREG[iVar]
            xNew = xVec[iVar] + lamB*delX[iVar]
            perturbVariable(xNew,clsREG,clsEOS0,clsEOS1,dicSAM0,dicSAM1,clsIO)

#-- Update C6+ component Acentric Factor and Volume Shifts ----------

        #adjustAFandSS(clsEOS1)

#--------------------------------------------------------------------
#  Run the experiments with the Perturbed EoS (cF)
#--------------------------------------------------------------------

        runRegExps(qExp,qDif,clsEOS1,dicSAM1,dicEXP,clsIO)

#-- And calculate the new SSQ ---------------------------------------            

        fun1,regI = calcResSSQ(qWrt,dicEXP,qExp,clsIO)

        print("lineSearch: fun0,grd0,lamB,fun1 {:10.5f} {:10.3e} {:10.5f} {:10.5f}".format(fun0,grd0,lamB,fun1))

#== What to do? =======================================================

        if fun1 <= fun0 + 1.0E-04*lamB*grd0 :

#-- Sufficent improvement: update variables and exit ----------------

            xVec = xVec + lamB*delX     #-- Vector!

            break

        else :

            lamT = updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)

        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)  #-- lamB > 0.1*lamT

#========================================================================
#  End of Routine
#========================================================================

    return fun1,xVec,regI

#========================================================================
#  Adjust C6+ Acentric Factor and Volume Shift
#========================================================================

def adjustAFandSS(clsEOS) :

    nCom = clsEOS.nComp

    for iC in range(nCom) :

        Mw = clsEOS.gPP("MW",iC)

        if Mw > 90.0 :
            AFo = clsEOS.gPP("AF",iC)
            SSo = clsEOS.gPP("SS",iC)
            AFn = CP.adjustAcentFac(clsEOS,iC)
            SSn =   CP.calcVolShift(clsEOS,iC)
            print("iC,AFo,AFn,SSo,SSn {:2d} {:7.4f} {:7.4f} {:7.4f} {:7.4f}".format(iC,AFo,AFn,SSo,SSn))

#== Nothing to return =================================================    

    return

#========================================================================
#  Line Search Lambda Back-Tracking Calculation
#========================================================================

def updateLambda(lamB,lam2,fun0,fun1,fun2,grd0) :

#-- First time thru, use Quadratic Fit (f0,g0,f1) -------------------    

    if   lamB == 1.0 :
        
        lamT = -0.5*grd0/(fun1 - fun0 - grd0)
        #print("Quad: lamT ",lamT)

#-- Else perform Cubic Back-Track -----------------------------------
        
    else :
        
        rhs1 = fun1 - fun0 - lamB*grd0
        rhs2 = fun2 - fun0 - lam2*grd0
        
        lmSQ = lamB*lamB
        l2SQ = lam2*lam2
        
        aCof = (      rhs1/lmSQ -      rhs2/l2SQ)/(lamB - lam2)
        bCof = (-lam2*rhs1/lmSQ + lamB*rhs2/l2SQ)/(lamB - lam2)
        
        if aCof == 0.0 : lamT = -0.5*grd0/bCof
        else :
            
            disC = bCof*bCof - 3.0*aCof*grd0
            
            if   disC <  0.0 : lamT = 0.5*lamB
            elif bCof <= 0.0 : lamT = (-bCof + sqrt(disC))/(3.0*aCof)
            else             : lamT = -grd0/(bCof + sqrt(disC))
            
        if lamT > 0.5*lamB   : lamT = 0.5*lamB

        #print("Cubi: lamT ",lamT)

#======================================================================
#  Return Step-Length Update
#======================================================================

    return lamT

#========================================================================
#  Run the Experiments
#========================================================================

def runRegExps(qExp,qDif,clsEOS,dicSAM,dicEXP,clsIO) :

    nExp = len(dicEXP)
    fOut = None

    clsEOS.eosCoefsNoPT()

#----------------------------------------------------------------------
#  Loop over active experiments
#----------------------------------------------------------------------

    qDif = False

    for iExp in range(nExp) :

        if qExp[iExp] :  #-- Is this experiment active in regression?

            clsEXP = dicEXP[iExp]

            xTyp = clsEXP.xName
            
            if   xTyp == "CCE" : CX.calcCCE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "CVD" : CX.calcCVD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "DLE" : CX.calcDLE(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SEP" : CX.calcSEP(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "FLS" : CX.calcFLS(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SAT" : CX.calcSAT(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "SWL" : CX.calcSWL(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)
            elif xTyp == "GRD" : CG.calcGRD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO)

#======================================================================
#  End of Routine
#======================================================================

    return 

#========================================================================
#  Regression Residual and Sum of Squares Error (SSQ)
#========================================================================

def calcResSSQ(qWrt,dicEXP,qExp,clsIO) :

    nExp = len(dicEXP)
    nReg = 0

    numTyp = {'CCE':0,'CVD':0,'DLE':0,'SEP':0, \
              'FLS':0,'SAT':0,'SWL':0,'GRD':0}

    regI = []

    if qWrt : regO = []

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        nRow = clsEXP.nRow
        nObs = clsEXP.nObs

        xTyp = clsEXP.xName     #-- Exp Type [see numTyp]

        numTyp[xTyp] += 1       #-- Number of Times this Type

        xOut = xTyp + str(numTyp[xTyp])

        #print("calcResSSQ: iExp,nObs,nRow,qExp ",iExp,nObs,nRow,qExp[iExp])

#======================================================================
#  Is this experiment active in regression?
#======================================================================

        if qExp[iExp] :

#== Psat for CCE, CVD, DLE & SEP ======================================

            if clsEXP.nRsat >= 0 :
                if clsEXP.PsatO > 0.0 and clsEXP.PsatW > 0.0 :
                    sObs = "PSAT"
                    sStg = "   PSAT"
                    dVal  = clsEXP.PsatW*(clsEXP.PsatO - clsEXP.PsatC)/clsEXP.PsatO
                    regI.append(dVal)
                    nReg += 1
                    if qWrt :
                        rowR = [xOut,sObs,sStg,clsEXP.PsatW,clsEXP.PsatO,clsEXP.PsatC]
                        regO.append(rowR)

#== Array-based observations ==========================================        

            for iObs in range(nObs) :

                sObs = clsEXP.hObs[iObs]

                if clsEXP.qObs[iObs] :

                    for iRow in range(nRow) :

                        sStg = "{:7.1f}".format(clsEXP.dInd[0][iRow])

                        dWei = clsEXP.dWei[iObs][iRow]
                        dObs = clsEXP.dObs[iObs][iRow]
                        dCal = clsEXP.dCal[iObs][iRow]

                        if dWei > 0.0 and dObs > 0.0 :
                            
                            if sObs == "SLIQ" or sObs == "MREM" : cDiv = 1.0
                            else                                : cDiv = 1.0/dObs

                            dVal = dWei*(dObs - dCal)*cDiv
                            regI.append(dVal)
                            nReg += 1

                            if qWrt :
                                rowR = [xOut,sObs,sStg,dWei,dObs,dCal]
                                regO.append(rowR)

#----------------------------------------------------------------------
#  Assemble the SSQ - Remember factor of 1/2
#----------------------------------------------------------------------

    dSSQ = 0.0
    for iReg in range(nReg) : dSSQ = dSSQ + regI[iReg]*regI[iReg]

    dSSQ = 0.5*dSSQ

    #print("numObs: nReg,dSSQ ",nReg,dSSQ)

#----------------------------------------------------------------------
#  Output the Regression Information to reg-file?
#----------------------------------------------------------------------

    if qWrt : WR.outRegInfo(nReg,dSSQ,regO,regI,clsIO)

#======================================================================
#  Return SSQ and Residual Vector
#======================================================================

    return dSSQ,regI

#========================================================================
#  Regression Variable "Type"; REG, PLUS or LBC
#========================================================================

def regVarType(dicREG) :

    sReg = " "
    iERR = 0

    nVar = len(dicREG)

    for iVar in range(nVar) :
        clsREG = dicREG[iVar]
        sTyp   = clsREG.vTyp
        if   sTyp == "REG" or sTyp == "KIJ" : sThs = "REG"
        elif sTyp == "PLS"                  : sThs = "PLS"
        elif sTyp == "LBC"                  : sThs = "LBC"
        #print("iVar,sTyp,sThs ",iVar,sTyp,sThs)
        if iVar == 0 : sReg = sThs
        else         :
            if sThs == sReg :
                pass
            else :
                print("Can't mix ",sThs," and ",sReg," variable types")
                iERR = -1

#== Return value ======================================================

    return iERR,sReg

#======================================================================
#  End of Module
#======================================================================

