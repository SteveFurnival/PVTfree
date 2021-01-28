
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Fits the EOS data for Blackoil Table Generation
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import blkCoefs  as BC
import blkProps  as BP
import calcReg   as CR

#========================================================================
#  Fit the 2-Component (STO & STG) EoS Parameters to Saturated Data
#========================================================================

def regEoSData(dTab,cTab,dicSAM,clsEOS,clsBLK,clsIO) :

    nEOS = clsBLK.nEOS
    nSat = len(dTab)

    pObs = NP.zeros(nSat)
    xOil = NP.zeros(nSat) ; yOil = NP.zeros(nSat)
    vOil = NP.zeros(nSat) ; vGas = NP.zeros(nSat)

    for iSat in range(nSat) :
        pObs[iSat] = dTab[iSat][clsBLK.iPr]
        xOil[iSat] = cTab[iSat][clsBLK.iXo]
        yOil[iSat] = cTab[iSat][clsBLK.iYo]
        vOil[iSat] = cTab[iSat][clsBLK.iMo]
        vGas[iSat] = cTab[iSat][clsBLK.iMg]

#== Initial Residual Vector ===========================================

    print("Start : Regress EOS")

    ssq0,res0 = calcResEOSBoth(pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO)

    ssqL = ssq0

    qConv = False
    nIter = 1

    dicREG = None

#----------------------------------------------------------------------
#  Iteration Loop
#----------------------------------------------------------------------

    while not qConv :

#== Jacobian Matrix ===================================================

        jacB = calcJacEOSBoth(pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO)

#== Assemble the Gradient Vector and Hessian Matrix ===================

        Grad = NP.dot(jacB,res0)
        Hess = NP.dot(jacB,jacB.T)  #-- Note Transpose

#== Calculate the update vector delX ==================================

        delX = CR.rotDisc(nIter,Grad,Hess,dicREG,dicSAM,clsEOS,clsIO)

#== Initial Residual Vector ===========================================    

        ssqT,res0 = lineSearchEOS(ssqL,delX,Grad,pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO)
        
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

    print("Finish: Regress EOS")
    
    #print("aOil,aGas {:10.3e} {:10.3e}".format(clsBLK.EOS1["aOil"],clsBLK.EOS1["aGas"]))
    #print("bOil,bGas {:10.3e} {:10.3e}".format(clsBLK.EOS1["bOil"],clsBLK.EOS1["bGas"]))

#========================================================================
#  Fit Saturated FVF, one pair at a time
#========================================================================

    nSat = len(dTab)
    
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

    sOil = NP.zeros(nSat) ; sGas = NP.zeros(nSat)

    for iSat in range(nSat) :

        Pr = dTab[iSat][clsBLK.iPr]

        RTp = clsBLK.RT/Pr

        Rs = dTab[iSat][clsBLK.iRs] ; Rv = dTab[iSat][clsBLK.iRv]
        Bo = dTab[iSat][clsBLK.iBo] ; Bg = dTab[iSat][clsBLK.iBg]

        BoCon = dSTO + Rs*dSTG
        BgCon = dSTG + Rv*dSTO

        xO = xOil[iSat] ; yO = yOil[iSat]

        sO,sG = fitSatFVF(RTp,Bo,Bg,BoCon,BgCon,xO,yO,clsBLK,clsIO)

        sOil[iSat] = sO ; sGas[iSat] = sG

        #sOut = "iS,Pr,sO,sG {:2d} {:10.3f} {:10.3e} {:10.3e}".format(iSat,Pr,sO,sG)
        #print(sOut)

#== Return values =====================================================

    return sOil,sGas

#========================================================================
#  Calculate Residuals for Oil AND Gas EOS
#========================================================================

def calcResEOSBoth(pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO) :

    nSat = len(pObs)
    nRes = 2*nSat

    resB = NP.zeros(nRes)

#-- Calculate Residuals for each phase individually -----------------    

    ssqL,resL = calcResEOSPhase(pObs,xOil,vOil,clsBLK,clsIO)
    ssqV,resV = calcResEOSPhase(pObs,yOil,vGas,clsBLK,clsIO)

#-- And then combine them into one array (Oil, then Gas) ------------    

    sTot = ssqL + ssqV

    for iRes in range(nSat) :
        resB[iRes]      = resL[iRes]
        resB[iRes+nSat] = resV[iRes]

    #print("sTot {:10.3e}".format(sTot))

    return sTot,resB

#========================================================================
#  Calculate Residuals for Oil or Gas EOS
#========================================================================

def calcResEOSPhase(pObs,xOil,vOil,clsBLK,clsIO) :

    if clsIO.Deb["BLACK"] > 0 : qDeb = True
    else                      : qDeb = False

    nSat = len(pObs)
    RT   = clsBLK.RT

    resP = NP.zeros(nSat)
    ssQ  = 0.0

    #print("calcResEOSPhase: aO,bO,sO,aG,bG,sG {:10.3f} {:10.5f} {:10.5f} {:10.3f} {:10.5f} {:10.5f}".format(clsBLK.EOS1["aOil"],clsBLK.EOS1["bOil"],clsBLK.EOS1["sOil"],clsBLK.EOS1["aGas"],clsBLK.EOS1["bGas"],clsBLK.EOS1["sGas"]))

    for iSat in range(nSat) :

        aCof = BC.Acoef(xOil[iSat],clsBLK)
        bCof = BC.Bcoef(xOil[iSat],clsBLK)
        cCof = BC.Ccoef(xOil[iSat],clsBLK)

        uCof = BC.Ucoef(vOil[iSat],bCof,cCof)
        wCof = BC.Wcoef(vOil[iSat],bCof,cCof)

        pCal = RT/uCof - aCof/wCof

        resP[iSat] = (pCal - pObs[iSat])/pObs[iSat]

        if qDeb > 0 :
            print("iS,pO,pC,rP {:2d} {:10.3e} {:10.3e} {:10.3e}".format(iSat,pObs[iSat],pCal,resP[iSat]))
        
        ssQ = ssQ + resP[iSat]*resP[iSat]

#== Return SSQ and Residuals ==========================================

    ssQ = 0.5*ssQ

    #print("ssQ {:10.3e}".format(ssQ))

    return ssQ,resP

#========================================================================
#  Calculate Jacobian for Oil AND Gas EOS
#========================================================================

def calcJacEOSBoth(pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO) :

    nSat = len(pObs)
    nRes = 2*nSat
    nEOS = clsBLK.nEOS

    jacB = NP.zeros((nEOS,nRes))

#-- Calculate Jacobians for each phase individually -----------------

    jacL = calcJacEOSPhase(pObs,xOil,vOil,clsBLK,clsIO)
    jacV = calcJacEOSPhase(pObs,yOil,vGas,clsBLK,clsIO)

#== And then "add" Oil & Gas Jacobians together =======================

    for iEOS in range(nEOS) :
        for iSat in range(nSat) :
            jSat = iSat + nSat
            jacB[iEOS][iSat] = jacL[iEOS][iSat]
            jacB[iEOS][jSat] = jacL[iEOS][iSat]

    return jacB

#========================================================================
#  Calculate Jacobian for Oil or Gas EOS
#========================================================================

def calcJacEOSPhase(pObs,xOil,vOil,clsBLK,clsIO) :

    nSat = len(pObs)
    nEOS = clsBLK.nEOS

    RT   = clsBLK.RT

    jacO = NP.zeros((nEOS,nSat))

#== Loop over Saturated Stages ========================================

    for iSat in range(nSat) :

        aCof = BC.Acoef(xOil[iSat],clsBLK)
        bCof = BC.Bcoef(xOil[iSat],clsBLK)
        cCof = BC.Ccoef(xOil[iSat],clsBLK)

        uCof = BC.Ucoef(vOil[iSat],bCof,cCof)
        wCof = BC.Wcoef(vOil[iSat],bCof,cCof)

        pCal = RT/uCof - aCof/wCof

        drda = -1.0/(     wCof*pObs[iSat])
        drdu = - RT/(uCof*uCof*pObs[iSat])
        drdw = aCof/(wCof*wCof*pObs[iSat])

        dudb = -1.0 ; dudc =  1.0

        dwdb = 2.0*(vOil[iSat] + cCof - bCof)
        dwdc = 2.0*(vOil[iSat] + cCof + bCof)

        drdb = drdu*dudb + drdw*dwdb
        drdc = drdu*dudc + drdw*dwdc

        dAdAo = BC.dAdAo(xOil[iSat],clsBLK)
        dAdAg = BC.dAdAg(xOil[iSat],clsBLK)
        
        dBdBo = BC.dBdBo(xOil[iSat],clsBLK)
        dBdBg = BC.dBdBg(xOil[iSat],clsBLK)

        dCdSo = BC.dCdSo(xOil[iSat],clsBLK)
        dCdSg = BC.dCdSg(xOil[iSat],clsBLK)

        drdAo = drda*dAdAo ; drdAg = drda*dAdAg
        drdBo = drdb*dBdBo ; drdBg = drdb*dBdBg
        drdSo = drdc*dCdSo ; drdSg = drdc*dCdSg

        jacO[0][iSat] = drdAo ; jacO[3][iSat] = drdAg
        jacO[1][iSat] = drdBo ; jacO[4][iSat] = drdBg
        jacO[2][iSat] = drdSo ; jacO[5][iSat] = drdSg

        #print("iS,j1,j2,j3,j4,j5,j6 {:3d} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(iSat,drdAo,drdBo,drdSo,drdAg,drdBg,drdSg))

    return jacO

#========================================================================
#  Line Search for EOS Parameters [Numerical Recipes Chapter 9.7]
#========================================================================

def lineSearchEOS(fun0,delX,Grad,pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO) :

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

    xVec = storeEOScoefs(clsBLK)

#----------------------------------------------------------------------
#  Find the Lambda Parameter
#----------------------------------------------------------------------

    qConv = False

    while not qConv :

#-- Update variables ------------------------------------------------

        updateEOScoefs(xVec,lamB,delX,clsBLK)
        
#--------------------------------------------------------------------
#  Run the experiments with the Perturbed EoS (cF)
#--------------------------------------------------------------------

        fun1,res1 = calcResEOSBoth(pObs,xOil,yOil,vOil,vGas,clsBLK,clsIO)

        #print("lnsrch: fun0,grd0,lamB,fun1 {:10.5f} {:10.5f} {:10.3e} {:10.5f}".format(fun0,grd0,lamB,fun1))

#== What to do? =======================================================

        if fun1 <= fun0 + 1.0E-04*lamB*grd0 :

#-- Sufficent improvement: update variables and exit ----------------

            for iVar in range(nVar) :
                xVec[iVar] = xVec[iVar] + lamB*delX[iVar]
                
            restoreEOScoefs(xVec,clsBLK)
            break

        elif lamB < lamM :

            restoreEOScoefs(xVec,clsBLK)
            break

        else :

            lamT = CR.updateLambda(lamB,lam2,fun0,fun1,fun2,grd0)
            
        lam2 = lamB
        fun2 = fun1
        lamB = max(lamT,0.1*lamB)  #-- lamB > 0.1*lamT

#========================================================================
#  End of module
#========================================================================

    return fun1,res1

#========================================================================
#  Work Array for Current EOS Parameters
#========================================================================

def storeEOScoefs(clsBLK) :

    nEOS = clsBLK.nEOS

    xVec = NP.zeros(nEOS)

    xVec[0] = clsBLK.EOS1["aOil"]
    xVec[1] = clsBLK.EOS1["bOil"]
    xVec[2] = clsBLK.EOS1["sOil"]
    xVec[3] = clsBLK.EOS1["aGas"]
    xVec[4] = clsBLK.EOS1["bGas"]
    xVec[5] = clsBLK.EOS1["sGas"]

    return xVec

#========================================================================
#  Restore the EOS Parameters
#========================================================================

def restoreEOScoefs(xVec,clsBLK) :

    nEOS = clsBLK.nEOS

    clsBLK.EOS1["aOil"] = xVec[0]
    clsBLK.EOS1["bOil"] = xVec[1]
    clsBLK.EOS1["sOil"] = xVec[2]
    clsBLK.EOS1["aGas"] = xVec[3]
    clsBLK.EOS1["bGas"] = xVec[4]
    clsBLK.EOS1["sGas"] = xVec[5]

    return

#========================================================================
#  Update the set of EOS Parameters
#========================================================================

def updateEOScoefs(xVec,lamB,delX,clsBLK) :

    clsBLK.EOS1["aOil"] = xVec[0] + lamB*delX[0]
    clsBLK.EOS1["bOil"] = xVec[1] + lamB*delX[1]
    clsBLK.EOS1["sOil"] = xVec[2] + lamB*delX[2]
    clsBLK.EOS1["aGas"] = xVec[3] + lamB*delX[3]
    clsBLK.EOS1["bGas"] = xVec[4] + lamB*delX[4]
    clsBLK.EOS1["sGas"] = xVec[5] + lamB*delX[5]

    aOM = clsBLK.EOS1["aOil"]
    bOM = clsBLK.EOS1["bOil"]
    sOM = clsBLK.EOS1["sOil"]

    aGM = clsBLK.EOS1["aGas"]
    bGM = clsBLK.EOS1["bGas"]
    sGM = clsBLK.EOS1["sGas"]

    #print("aO,bO,sO,aG,bG,sG {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(aOM,bOM,sOM,aGM,bGM,sGM))

    return

#========================================================================
#  Fit Saturated FVF [2 variables in 2 unknowns]
#========================================================================

def fitSatFVF(RTp,BoObs,BgObs,BoCon,BgCon,xOil,yOil,clsBLK,clsIO) :

    qConv = False
    nIter = 0

    while not qConv :

        nIter += 1

        rO,rG = resSatFVF(RTp,BoObs,BgObs,BoCon,BgCon,xOil,yOil,clsBLK)

        #print("sOil,sGas {:10.3e} {:10.3e}".format(clsBLK.EOS1["sOil"],clsBLK.EOS1["sGas"] ))
        #print("rO  ,rG   {:10.3e} {:10.3e}".format(rO ,rG ))

        jOO,jOG,jGO,jGG = jacSatFVF(BoObs,BgObs,BoCon,BgCon,xOil,yOil,clsBLK)

        #print("jOO,jOG {:10.3e} {:10.3e}".format(jOO,jOG))
        #print("jGO,jGG {:10.3e} {:10.3e}".format(jGO,jGG))

        drO,drG = solveSat(rO,rG,jOO,jOG,jGO,jGG,clsBLK)

        #print("drO,drG {:10.3e} {:10.3e}".format(drO,drG))

        ssQ = 0.0

        ssQ = drO*drO + drG*drG

        clsBLK.EOS1["sOil"] = clsBLK.EOS1["sOil"] + drO
        clsBLK.EOS1["sGas"] = clsBLK.EOS1["sGas"] + drG

        #fDeb = clsIO.fDeb
        #sOut= "nIter,ssQ {:2d} {:10.3e}\n".format(nIter,ssQ)
        #fDeb.write(sOut)

        if ssQ < 1.0E-12 :
            break

        if nIter > 5 :
            break

    return clsBLK.EOS1["sOil"],clsBLK.EOS1["sGas"]

#========================================================================
#  Residuals for Fit Saturated FVF [2 variables in 2 unknowns]
#========================================================================

def resSatFVF(RTp,BoObs,BgObs,BoCon,BgCon,xOil,yOil,clsBLK) :

    qOil = True
    qGas = False

    BoCal = BP.calcPhaseFVF(qOil,RTp,BoCon,xOil,clsBLK)
    BgCal = BP.calcPhaseFVF(qGas,RTp,BgCon,yOil,clsBLK)

    resO = (BoCal - BoObs)/BoObs
    resG = (BgCal - BgObs)/BgObs

    #fDeb = clsIO.fDeb
    #sOutO = "BoC,BoO,rO {:10.3e} {:10.3e} {:10.3e}\n".format(BoCal,BoObs,resO)
    #sOutG = "BgC,BgO,rG {:10.3e} {:10.3e} {:10.3e}\n".format(BgCal,BgObs,resG)
    #fDeb.write(sOutO)
    #fDeb.write(sOutG)

    return resO,resG

#========================================================================
#  Jacobian for Fit Saturated FVF [2 variables in 2 unknowns]
#========================================================================

def jacSatFVF(BoObs,BgObs,BoCon,BgCon,xOil,yOil,clsBLK) :

    Moil = BP.phaseMw(clsBLK.mSTG,clsBLK.mSTO,xOil)
    Mgas = BP.phaseMw(clsBLK.mSTG,clsBLK.mSTO,yOil)

    jOO = - BoCon*BC.dCdSo(xOil,clsBLK)/(Moil*BoObs)
    jOG = - BoCon*BC.dCdSg(xOil,clsBLK)/(Moil*BoObs)

    jGO = - BgCon*BC.dCdSo(yOil,clsBLK)/(Mgas*BgObs)
    jGG = - BgCon*BC.dCdSg(yOil,clsBLK)/(Mgas*BgObs)

    #print("jOO,jOG {:10.3e} {:10.3e}".format(jOO,jOG))
    #print("jGO,jGG {:10.3e} {:10.3e}".format(jGO,jGG))

    return jOO,jOG,jGO,jGG

#========================================================================
#  Solve for Fit Saturated FVF/Viscosities [2 variables in 2 unknowns]
#========================================================================

def solveSat(rO,rG,jOO,jOG,jGO,jGG,clsBLK) :

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
