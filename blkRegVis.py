
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Fits Viscosity Data for Blackoil Table Generation
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import sys

import blkCoefs  as BC
import blkProps  as BP
import calcReg   as CR
import utilities as UT

#========================================================================
#  Fit the 2-Component (STO & STG) Viscosity Parameters to Saturated Data
#========================================================================

def regViscData(dTab,cTab,dicSAM,clsEOS,clsBLK,clsIO) :

    nEOS = clsBLK.nEOS
    nSat = len(dTab)

    uOil = NP.zeros(nSat) ; uGas = NP.zeros(nSat)
    xOil = NP.zeros(nSat) ; yOil = NP.zeros(nSat)
    dOil = NP.zeros(nSat) ; dGas = NP.zeros(nSat)

    for iSat in range(nSat) :
        uOil[iSat] = dTab[iSat][clsBLK.iUo]
        uGas[iSat] = dTab[iSat][clsBLK.iUg]
        xOil[iSat] = cTab[iSat][clsBLK.iXo]
        yOil[iSat] = cTab[iSat][clsBLK.iYo]
        dOil[iSat] = cTab[iSat][clsBLK.iDo]
        dGas[iSat] = cTab[iSat][clsBLK.iDg]

    print("Start : Regress Viscosity")

    ssq0,res0 = calcResVISBoth(uOil,uGas,xOil,yOil,dOil,dGas,clsBLK,clsIO)

    ssqL = ssq0

    qConv = False
    nIter = 1
    
    dicREG = None

#----------------------------------------------------------------------
#  Iteration Loop
#----------------------------------------------------------------------

    while not qConv :

#== Jacobian Matrix ===================================================

        jacB = calcJacVISBoth(uOil,uGas,xOil,yOil,dOil,dGas,clsBLK)

#== Assemble the Gradient Vector and Hessian Matrix ===================

        Grad = NP.dot(jacB,res0)
        Hess = NP.dot(jacB,jacB.T)  #-- Note Transpose

#== Calculate the update vector delX ==================================

        delX = CR.rotDisc(nIter,Grad,Hess,dicREG,dicSAM,clsEOS,clsIO)

#== Initial Residual Vector ===========================================    

        ssqT,res0 = lineSearchVIS(ssqL,delX,Grad,uOil,uGas,xOil,yOil,dOil,dGas,clsBLK,clsIO)

        #print("regViscData: nIter,ssqT ",nIter,ssqT)

#======================================================================
#  Progress in SSQ?
#======================================================================

        if abs((ssqT - ssqL)/ssqT) < 1.0E-04 :
            #print("Regression converged! nIter,ssq0,ssqT {:2d}{:10.3e}{:10.3e}".format(nIter,ssq0,ssqT))
            qConv = True
        else :
            ssqL = ssqT

        nIter += 1

        if nIter > 5 :
            #print("Regression - Too Many Iterations")
            break

    #xVec = storeVIScoefs(clsBLK)
    #UT.writeVector(sys.stdout,'xVec',xVec,'D')

#========================================================================
#  Fit Saturated Viscosities, one pair at a time
#========================================================================

    nSat = len(uOil)

    rOil = NP.zeros(nSat) ; rGas = NP.zeros(nSat)

    for iSat in range(nSat) :

        uO = uOil[iSat] ; uG = uGas[iSat]
        xO = xOil[iSat] ; yO = yOil[iSat]
        dO = dOil[iSat] ; dG = dGas[iSat]

        nI,rO,rG = fitSatVisc(uO,uG,xO,yO,dO,dG,clsBLK,clsIO)

        rOil[iSat] = rO
        rGas[iSat] = rG

        #print("iS,pR,rO,rG {:2d} {:2d} {:10.3f} {:10.3e} {:10.3e}".format(iSat,nI,pR,rO,rG))

    print("Finish: Regress Viscosity")

#== Return values =====================================================    

    return rOil,rGas

#========================================================================
#  Calculate Residuals for Oil AND Gas Viscosity
#========================================================================

def calcResVISBoth(uOil,uGas,xOil,yOil,dOil,dGas,clsBLK,clsIO) :

    nSat = len(uOil)
    nRes = 2*nSat

    resB = NP.zeros(nRes)

#-- Calculate Residuals for each phase individually -----------------    

    ssqL,resL = calcResVISPhase('Liq',uOil,xOil,dOil,clsBLK,clsIO)
    ssqV,resV = calcResVISPhase('Vap',uGas,yOil,dGas,clsBLK,clsIO)

#-- And then combine them into one array (Oil, then Gas) ------------    

    sTot = ssqL + ssqV

    for iRes in range(nSat) :
        resB[iRes]      = resL[iRes]
        resB[iRes+nSat] = resV[iRes]

    #print("sTot {:10.3e}".format(sTot))

    return sTot,resB

#========================================================================
#  Calculate Residuals for Oil or Gas Viscosity
#========================================================================

def calcResVISPhase(sPhs,uObs,xOil,dOil,clsBLK,clsIO) :

    nSat = len(uObs)

    resP = NP.zeros(nSat)
    ssQ  = 0.0

    for iSat in range(nSat) :

        uCal = BP.calcLBCvisc(xOil[iSat],dOil[iSat],clsBLK)

        resP[iSat] = (uCal - uObs[iSat])/uObs[iSat]

        #print("Phase,iS,uO,uC,rP {:3s} {:2d} {:10.3e} {:10.3e} {:10.3e}".format(sPhs,iSat,uObs[iSat],uCal,resP[iSat]))
        
        ssQ = ssQ + resP[iSat]*resP[iSat]

#== Return SSQ and Residuals ==========================================

    ssQ = 0.5*ssQ

    #print("ssQ {:10.3e}".format(ssQ))

    return ssQ,resP

#========================================================================
#  Calculate Jacobian for Oil AND Gas Viscosity
#========================================================================

def calcJacVISBoth(uObs,uGas,xOil,yOil,dOil,dGas,clsBLK) :

    nSat = len(uObs)
    nRes = 2*nSat
    nEOS = clsBLK.nEOS

    jacB = NP.zeros((nEOS,nRes))

#-- Calculate Jacobians for each phase individually -----------------

    jacL = calcJacVISPhase(uObs,xOil,dOil,clsBLK)
    jacV = calcJacVISPhase(uGas,yOil,dGas,clsBLK)

#== And then "add" Oil & Gas Jacobians together =======================

    for iEOS in range(nEOS) :
        for iSat in range(nSat) :
            jSat = iSat + nSat
            jacB[iEOS][iSat] = jacL[iEOS][iSat]
            jacB[iEOS][jSat] = jacL[iEOS][iSat]

    return jacB

#========================================================================
#  Calculate Jacobian for Oil or Gas Viscosity
#========================================================================

def calcJacVISPhase(uObs,xOil,dOil,clsBLK) :

    nSat = len(uObs)
    nEOS = clsBLK.nEOS

    jacO = NP.zeros((nEOS,nSat))

#== Loop over Saturated Stages ========================================

    for iSat in range(nSat) :

        rCof = BC.Rcoef(xOil[iSat],clsBLK)
        eCof = BC.Ecoef(xOil[iSat],clsBLK)

        dRed = dOil[iSat]/rCof

        sumV = BP.sumVisc(dRed)
        sumD = BP.sumDerv(dRed)

        drduC = 1.0/uObs[iSat]

        duds = 4.0*pow(sumV,3)/eCof
        dsdr = sumD
        drdD = -dRed/rCof

        dudr =           duds*dsdr*drdD/ uObs[iSat]
        dudu =                      1.0/ uObs[iSat]
        dude = - (pow(sumV,4) - 0.0001)/(uObs[iSat]*eCof*eCof)

        drdRo = dudr*BC.dRdRo(xOil[iSat],clsBLK)
        drdRg = dudr*BC.dRdRg(xOil[iSat],clsBLK)

        drdEo = dude*BC.dEdEo(xOil[iSat],clsBLK)
        drdEg = dude*BC.dEdEg(xOil[iSat],clsBLK)

        drdUo = dudu*BC.dMdUo(xOil[iSat],clsBLK)
        drdUg = dudu*BC.dMdUg(xOil[iSat],clsBLK)

        jacO[0][iSat] = drdRo ; jacO[3][iSat] = drdRg
        jacO[1][iSat] = drdEo ; jacO[4][iSat] = drdEg
        jacO[2][iSat] = drdUo ; jacO[5][iSat] = drdUg

        #print("iS,j1,j2,j3,j4,j5,j6 {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(iSat,drdRo,drdEo,drdUo,drdRg,drdEg,drdUg))

    return jacO

#========================================================================
#  Line Search for Viscosity Parameters [Numerical Recipes Chapter 9.7]
#========================================================================

def lineSearchVIS(fun0,delX,Grad,uOil,uGas,xOil,yOil,dOil,dGas,clsBLK,clsIO) :

    nVar = len(delX)

    grd0 = NP.dot(delX,Grad)

    if grd0 > 0.0 :
        print("lineSearch: grd0 > 0")
        fun1 = fun0
        res1 = [0.0 for i in range(1)]
        return fun1,res1

    lamB = 1.0      #-- Default first   step-length [Full Newton!]
    lamM = 1.0E-04  #-- Minimum allowed step-length
    lam2 = lamB
    fun2 = fun0

#== Create Work Array to Hold Current EOS Multipliers =================

    xVec = storeVIScoefs(clsBLK)

#----------------------------------------------------------------------
#  Find the Lambda Parameter
#----------------------------------------------------------------------

    qConv = False

    while not qConv :

#-- Update variables ------------------------------------------------

        updateVIScoefs(xVec,lamB,delX,clsBLK,clsIO)
        
#--------------------------------------------------------------------
#  Run the experiments with the Perturbed EoS (cF)
#--------------------------------------------------------------------

        fun1,res1 = calcResVISBoth(uOil,uGas,xOil,yOil,dOil,dGas,clsBLK,clsIO)

        #print("lnsrch: fun0,grd0,lamB,fun1 {:10.5f} {:10.5f} {:10.3e} {:10.5f}".format(fun0,grd0,lamB,fun1))

#== What to do? =======================================================

        if fun1 <= fun0 + 1.0E-04*lamB*grd0 :

#-- Sufficent improvement: update variables and exit ----------------

            for iVar in range(nVar) :
                xVec[iVar] = xVec[iVar] + lamB*delX[iVar]
                
            restoreVIScoefs(xVec,clsBLK)
            break

        elif lamB < lamM :

            restoreVIScoefs(xVec,clsBLK)
            break

        else :

            lamT = CR.updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)
            
        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)  #-- lamB > 0.1*lamT

#== Return values =====================================================        

    return fun1,res1

#========================================================================
#  Work Array for Current Viscosity Parameters
#========================================================================

def storeVIScoefs(clsBLK) :

    nEOS = clsBLK.nEOS

    xVec = [0.0 for i in range(nEOS)]

    xVec[0] = clsBLK.VIS1["rOil"]
    xVec[1] = clsBLK.VIS1["eOil"]
    xVec[2] = clsBLK.VIS1["uOil"]
    xVec[3] = clsBLK.VIS1["rGas"]
    xVec[4] = clsBLK.VIS1["eGas"]
    xVec[5] = clsBLK.VIS1["uGas"]

    return xVec

#========================================================================
#  Restore the Viscosity Parameters
#========================================================================

def restoreVIScoefs(xVec,clsBLK) :

    nEOS = clsBLK.nEOS

    clsBLK.VIS1["rOil"] = xVec[0]
    clsBLK.VIS1["eOil"] = xVec[1]
    clsBLK.VIS1["uOil"] = xVec[2]
    clsBLK.VIS1["rGas"] = xVec[3]
    clsBLK.VIS1["eGas"] = xVec[4]
    clsBLK.VIS1["uGas"] = xVec[5]

    return

#========================================================================
#  Update the set of Viscosity Parameters
#========================================================================

def updateVIScoefs(xVec,lamB,delX,clsBLK,clsIO) :

    clsBLK.VIS1["rOil"] = xVec[0] + lamB*delX[0]
    clsBLK.VIS1["eOil"] = xVec[1] + lamB*delX[1]
    clsBLK.VIS1["uOil"] = xVec[2] + lamB*delX[2]
    clsBLK.VIS1["rGas"] = xVec[3] + lamB*delX[3]
    clsBLK.VIS1["eGas"] = xVec[4] + lamB*delX[4]
    clsBLK.VIS1["uGas"] = xVec[5] + lamB*delX[5]

    return

#========================================================================
#  Fit Saturated Viscosities [2 variables in 2 unknowns]
#========================================================================

def fitSatVisc(uOobs,uGobs,xOil,yOil,dOil,dGas,clsBLK,clsIO) :

    qConv = False
    nIter = 0

    while not qConv :

        nIter += 1

        rO,rG = resSatVisc(uOobs,uGobs,xOil,yOil,dOil,dGas,clsBLK)

        jOO,jOG,jGO,jGG = jacSatVisc(uOobs,uGobs,xOil,yOil,dOil,dGas,clsBLK)

        drO,drG = solveSat(rO,rG,jOO,jOG,jGO,jGG,clsBLK,clsIO)

        ssQ = drO*drO + drG*drG

        #if abs(drO) > 0.5*abs(clsBLK.VIS1["rOil"]) : drO = 0.5*clsBLK.VIS1["rOil"]*NP.sign(drO)
        #if abs(drG) > 0.5*abs(clsBLK.VIS1["rGas"]) : drG = 0.5*clsBLK.VIS1["rGas"]*NP.sign(drG)

        clsBLK.VIS1["rOil"] = clsBLK.VIS1["rOil"] + drO
        clsBLK.VIS1["rGas"] = clsBLK.VIS1["rGas"] + drG

        #print("nIter,ssQ {:2d} {:10.3e}".format(nIter,ssQ))

        if ssQ < 1.0E-12 :
            break

        if nIter > 25 :
            break

    return nIter,clsBLK.VIS1["rOil"],clsBLK.VIS1["rGas"]

#========================================================================
#  Residuals for Fit Saturated Viscosities [2 variables in 2 unknowns]
#========================================================================

def resSatVisc(uOobs,uGobs,xOil,yOil,dOil,dGas,clsBLK) :

    uOcal = BP.calcLBCvisc(xOil,dOil,clsBLK)
    uGcal = BP.calcLBCvisc(yOil,dGas,clsBLK)

    resO = (uOcal - uOobs)/uOobs
    resG = (uGcal - uGobs)/uGobs

    return resO,resG

#========================================================================
#  Jacobian for Fit Saturated Viscosities [2 variables in 2 unknowns]
#========================================================================

def jacSatVisc(uOobs,uGobs,xOil,yOil,dOil,dGas,clsBLK) :

    duOdR = BP.calcLBCderv(xOil,dOil,clsBLK)
    duGdR = BP.calcLBCderv(yOil,dGas,clsBLK)

    jacOO = duOdR*BC.dRdRo(xOil,clsBLK)/uOobs
    jacOG = duOdR*BC.dRdRg(xOil,clsBLK)/uOobs
    jacGO = duGdR*BC.dRdRo(yOil,clsBLK)/uGobs
    jacGG = duGdR*BC.dRdRg(yOil,clsBLK)/uGobs

    return jacOO,jacOG,jacGO,jacGG

#========================================================================
#  Solve for Fit Saturated FVF/Viscosities [2 variables in 2 unknowns]
#========================================================================

def solveSat(rO,rG,jOO,jOG,jGO,jGG,clsBLK,clsIO) :

    disc = jOO*jGG - jOG*jGO

    if abs(disc) > 1.0E-15 :
        drO = -(jGG*rO - jOG*rG)/disc
        drG = -(jOO*rG - jGO*rO)/disc
    else :
        drO = 0.0
        drG = 0.0

    ssq = drO*drO + drG*drG

    #print("rO,rG,dRo,dRg {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(clsBLK.VIS1["rOil"],clsBLK.VIS1["rGas"],drO,drG,ssq))

    return drO,drG

#========================================================================
#  End of Module
#========================================================================
