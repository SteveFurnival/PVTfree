#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Perform Composition versus Depth (Gradient) 'Experiment' Module
#  
#========================================================================

#!/usr/bin/python3

import numpy        as NP

import allData   as AD
import allProp   as AP
import calcEOS   as CE
import calcSat   as CS
import calcStab  as CT
import utilities as UT
import writeOut  as WO

from math  import log,exp

#========================================================================
#  Composition versus Depth (GRD) Experiment
#========================================================================

def calcGRD(iExp,qDif,clsEOS,dicSAM,clsEXP,clsIO) :

    if clsIO.Deb["GRD"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else:
        qDeb = False

    qSat = True

#-- Various dimensions ----------------------------------------------    

    nCom = clsEOS.nComp
    nSam = clsEXP.nSamp
    nDep = clsEXP.nRow - 1      #-- Extra Row for GOC, if present
    nRef = clsEXP.nDref

    nGOC = nDep
    
#-- Pointers into the potential observed data (for calculated data) -

    typIND = AP.classLIB().INDshrt.get("GRD")
    typOBS = AP.classLIB().OBSshrt.get("GRD")

    iHG = typIND.index("HEIG")

    iPR = typOBS.index("PRES")
    iPS = typOBS.index("PSAT")
    iDN = typOBS.index("DENS")
    iC1 = typOBS.index("ZC1")
    iCP = typOBS.index("ZC7+")

#-- Sample Name -----------------------------------------------------    

    clsSAM = dicSAM[nSam]
    sNam   = clsSAM.sNam

    clsWRK = AD.classSample("GRDwrk")
    clsWRK.setIntComp(nCom,nCom)

#-- Load Reference Composition --------------------------------------

    zRef = NP.zeros(nCom)
    for iC in range(nCom) : zRef[iC] = clsSAM.gZI(iC)

#-- Reference Depth & Pressure --------------------------------------    

    pRef = clsEXP.Pref
    hRef = clsEXP.Dref

#-- Reservoir Temperature -------------------------------------------

    tRes = clsEXP.Tres

#== Saturation Pressure at Reference Depth ============================        

    pSatO = clsEXP.PsatO

    if pSatO > 0.0 : pObs =  0.95*pSatO
    else           : pObs = -1.0
        
    qBub = None
    pSat = None
    logK = NP.empty(nCom)
    
    qBub,pSat,Ksat = CS.calcPsat(pObs,tRes,qBub,pSat,logK,clsEOS,clsSAM,clsIO)

    if qBub : iRef =  1     #-- Ref-State is BUB
    else    : iRef = -1     #--              DEW

    iNeu =  0
    iLiq =  1
    iVap = -1

    mSat,vSat,dSat,zSat,uSat,dumS,dumS = CE.calcProps(iNeu,pSat,tRes,zRef,clsEOS)

    zC1,zC7 = calcMoleFracGRD(zRef,clsEOS)

    clsEXP.dCal[iPR][nRef] = pRef
    clsEXP.dCal[iPS][nRef] = pSat
    clsEXP.dCal[iDN][nRef] = dSat
    clsEXP.dCal[iC1][nRef] = zC1
    clsEXP.dCal[iCP][nRef] = zC7

    for iC in range(nCom) : clsEXP.zCal[nRef][iC] = zRef[iC]

#-- Get Reference (log) Fugacity Coeffs -----------------------------

    fRef,dumP = setupCalcGRDCoefsP(iNeu,pRef,tRes,zRef,clsEOS)  #-- Dont need Ref dlnPhi/dP

#======================================================================
#  GOC Search, if appropriate
#======================================================================

    iTop = 0
    iBot = 0

    if nRef > 0    :
        
        hTop = clsEXP.dInd[iHG][0]
        dTop = hTop - hRef
        zTyp = "Top"

#-- Calculate a Pressure & Composition at Top Depth -----------------

        pTop,zTop,qTop,sTop = \
            calcTopBottom(zTyp,dTop,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        if not qTop : iTop = - 1    #-- Is DEW point
        
    else : hTop = None

    if nRef < nDep :
        
        hBot = clsEXP.dInd[iHG][nDep-1]
        dBot = hBot - hRef
        zTyp = "Bot"

#-- Calculate a Pressure & Composition at Bottom Depth --------------

        pBot,zBot,qBot,sBot = \
            calcTopBottom(zTyp,dBot,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        if qBot : iBot = 1          #-- Is BUB point
        
    else : hBot = None

#--------------------------------------------------------------------
#  Search for Gas-Oil-Contact?
#--------------------------------------------------------------------

    if iTop*iBot < 0 :

        if iRef == iTop : hTop = hRef
        else            : hBot = hRef
        
        print("Gas-Oil-Contact Detected in Depth Interval")
        
        hGOC,pGOC,sGOC,dGOC,zC1,zC7,zDew,zBub = \
            calcGOC(hTop,hBot,tRes,hRef,pRef,zRef,fRef,clsEOS,clsIO)
    else :
        print("No Gas-Oil-Contact Detected in Depth Interval")
        hGOC = 0.0 ; pGOC = 0.0 ; sGOC = 0.0
        dGOC = 0.0 ; zC1  = 0.0 ; zC7  = 0.0
        zDew = NP.zeros(nCom) ; zBub = NP.zeros(nCom)

    clsEXP.dInd[iHG][nGOC] = hGOC
    
    clsEXP.dCal[iPR][nGOC] = pGOC
    clsEXP.dCal[iPS][nGOC] = sGOC
    clsEXP.dCal[iDN][nGOC] = dGOC
    clsEXP.dCal[iC1][nGOC] = zC1
    clsEXP.dCal[iCP][nGOC] = zC7

    for iC in range(nCom) : clsEXP.zCal[nGOC][iC] = zDew[iC]

#-- Above and Below any Possible GOC --------------------------------

    hAbv = None
    hBel = None

#--------------------------------------------------------------------
#  Above Reference Depth
#--------------------------------------------------------------------

    zTyp = "Abv"

    for iDep in range(nRef-1,-1,-1) :

        hDep = clsEXP.dInd[iHG][iDep]
        dDep = hDep - hRef

        #print("Above: iDep,hDep ",iDep,hDep)

        pDep,zDep,qDep,sDep = \
            calcTopBottom(zTyp,dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        mDep,vDep,yDep,gDen,uDep,dumS,dumS = CE.calcProps(iNeu,pDep,tRes,zDep,clsEOS)
        zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

        clsEXP.dCal[iPR][iDep] = pDep
        clsEXP.dCal[iPS][iDep] = sDep
        clsEXP.dCal[iDN][iDep] = yDep
        clsEXP.dCal[iC1][iDep] = zC1
        clsEXP.dCal[iCP][iDep] = zC7

        for iC in range(nCom) : clsEXP.zCal[iDep][iC] = zDep[iC]

        if qDeb :
            if qBub : isBub =  1
            else    : isBub = -1
            sOut = "Above: iDep,hDep,pDep,pSat,Bub,zC1,zC7+,Dens {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:6.4f} {:6.4f} {:7.3f}\n".format(iDep,hDep,pDep,pSat,isBub,zC1,zC7,yDep)
            fDeb.write(sOut)
            WO.writeArrayDebug(fDeb,zDep,"zDep")

#--------------------------------------------------------------------
#  Below Reference Depth
#--------------------------------------------------------------------

    zTyp = "Bel"

    for iDep in range(nRef+1,nDep) :
        
        hDep = clsEXP.dInd[iHG][iDep]
        dDep = hDep - hRef

        #print("Below: iDep,hDep ",iDep,hDep)

        pDep,zDep,qDep,sDep = \
            calcTopBottom(zTyp,dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        mDep,vDep,yDep,gDen,uDep,dumS,dumS = CE.calcProps(iNeu,pDep,tRes,zDep,clsEOS)
        zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

        clsEXP.dCal[iPR][iDep] = pDep
        clsEXP.dCal[iPS][iDep] = sDep
        clsEXP.dCal[iDN][iDep] = yDep
        clsEXP.dCal[iC1][iDep] = zC1
        clsEXP.dCal[iCP][iDep] = zC7

        for iC in range(nCom) : clsEXP.zCal[iDep][iC] = zDep[iC]

        if qDeb :
            if qBub : isBub =  1
            else    : isBub = -1
            sOut = "Below: iDep,hDep,pDep,pSat,Bub,zC1,zC7+,Dens {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:6.4f} {:6.4f} {:7.3f}\n".format(iDep,hDep,pDep,pSat,isBub,zC1,zC7,yDep)
            fDeb.write(sOut)
            WO.writeArrayDebug(fDeb,zDep,"zDep")

#======================================================================
#  End of Routine
#======================================================================

    return

#========================================================================
#  Calculates Pressure and Composition at Top or Bottom Depth
#  Prelude to potential GOC calculation
#========================================================================

def calcTopBottom(zTyp,dDep,pDep,tRes,zDep,fRef,clsEOS,clsIO) :

    nCom = clsEOS.nComp
    qSat = True
    qBub = None
    pSat = None

#-- Calculate a Pressure & Composition at Bottom Depth --------------
        
    pDep,zDep = calcStepGRD(dDep,pDep,tRes,zDep,fRef,clsEOS,clsIO)
        
#-- Is the new Composition Stable? ----------------------------------        
        
    iDep,depT,kDep = CT.twoSidedStabTest(qSat,pDep,tRes,zDep,clsEOS,clsIO)

#-- Unstable?  Re-Initialise and Calculate New Bottom Composition & Pressure

    if iDep > 0 :

        if   iDep == 1 :    #-- Liquid-Like Fluid Unstable
            zDep = kDep/zDep
        elif iDep == 2 :    #-- Vapour-Like Fluid Unstable
            zDep = kDep*zDep
        else :
            print("Search at " & zTyp & " Depth Returned Unstable Fluid with iStab = 3")
            crash = 1.0/0.0

        zDep = UT.Norm(zDep)
        
        pDep,zDep = calcStepGRD(dDep,pDep,tRes,zDep,fRef,clsEOS,clsIO)

        logK = NP.log(kDep)

    else :
        logK = NP.zeros(nCom)

#== Now Calculate the Saturation Pressure at this Depth ===============

    clsWRK = AD.classSample("GRDwrk")
    clsWRK.setIntComp(nCom,nCom)

    for iC in range(nCom) : clsWRK.sZI(iC,zDep[iC])

    qSat,pSat,Ksat = CS.calcPsat(pDep,tRes,qBub,pSat,logK,clsEOS,clsWRK,clsIO)

    #print("calcTopBottom: dDep,pDep,iDep,qSat,pSat ",dDep,pDep,iDep,qSat,pSat)
    
#== Return Values =====================================================

    return pDep,zDep,qSat,pSat

#========================================================================
#  Take a Step Away from the Reference Depth (GRD Experiment)
#  Whitson & Brule, SPE Phase Behavior Monograph Section 4.6.2
#========================================================================

def calcStepGRD(dDep,pDep,tRes,zRef,fRef,clsEOS,clsIO) :

    mxSS = 50

    nCom = clsEOS.nComp
    iLiq = 0

    zDep = NP.zeros(nCom)   #-- Composition at this Depth
    yMol = NP.zeros(nCom)
    res0 = NP.zeros(nCom)
    res1 = NP.zeros(nCom)
    fBar = NP.zeros(nCom)

#-- Gravity Correction: Whitson Eqn.(4.96) --------------------------

    grvC = dDep/(144.0*UT.gasCon*tRes)

    #print("dDep,grvC {:10.3f} {:10.3e}".format(dDep,grvC))

    zDep = zRef                                      #-- Vec: initial z
    fBar = NP.multiply(fRef,NP.exp(grvC*clsEOS.Mw))  #-- Vec: W&B Eqn.(4.96)

#======================================================================
#  Main Iterative Loop
#======================================================================

    iSS   = 0     #-- Successive-Substitution Counter
    qConv = False

    while not qConv :

        iSS += 1

        if iSS > mxSS : break   #-- Not converged!

#-- Fugacity Coefs and their P-Derivatives of New Fluid -------------

        fDep,dFdP = setupCalcGRDCoefsP(iLiq,pDep,tRes,zDep,clsEOS)

#-- Mole Numbers: Whitson Eqn.(4.95) --------------------------------

        res1 = res0                 #-- Vec: Might be used in GDEM
        res0 = fBar/fDep            #-- Vec
        yMol = zDep*res0            #-- Vec: W&B Eqn.(4.95)
        sumY = NP.sum(yMol)         #-- Scalar
        res0 = res0/sumY            #-- Vec: W&B Eqn.(4.97)
        wrkQ = yMol*res0*dFdP/fDep
        dQdP = NP.sum(wrkQ)         #-- Scalar: W&B Eqn.(4.100)

        qFun = 1.0 - sumY           #-- W&B Eqn.(4.94)
        dPrs = -qFun/dQdP           #-- W&B Eqn.(4.99a)
        pDep = pDep + dPrs          #-- W&B Eqn.(4.99b)

#== Converged? ========================================================

        cWrk = res0*zDep/yMol       #-- Vec
        cWrk = NP.log(cWrk)
        cTst = NP.sum(cWrk)
        cTst = cTst*cTst

        if cTst < 1.0E-12 : break

#== SS or GDEM Update? W&B Eqn.(4.83) =================================

        if iSS % UT.mGDEM1 > 0 : eigV = 1.0
        else :
            
            res0 = NP.log(res0)         #-- Take logs of Residuals
            res1 = NP.log(res1)
                
            eigV = UT.GDEM1(res0,res1,clsIO)  #-- GDEM works on log-Residuals

            res0 = NP.exp(res0)         #-- Restore Residuals
            res1 = NP.exp(res1)

        #print("calcStepGRD: iSS,qFun,dQdP,eigV,dPrs,cTst {:2d} {:10.3e} {:10.3e} {:8.3f} {:10.3e} {:10.3e}".format(iSS,qFun,dQdP,eigV,dPrs,cTst))

        yMol = yMol*NP.power(res0,eigV)

#-- Mole Numbers to Composition -------------------------------------

        zDep = UT.Norm(yMol)

#== Return the updated pressure and composition =======================

    #print("calcStepGRD: iSS,pDep {:2d} {:10.3f}".format(iSS,pDep))

    return pDep,zDep

#========================================================================
#  Calculate a Gas-Oil-Contact Depth, Given Above/Below Bounds
#========================================================================

def calcGOC(hAbv,hBel,tRes,hRef,pRef,zRef,fRef,clsEOS,clsIO) :

    if clsIO.Deb["GRD"] > 0 :
        qDeb = True
        fDeb = clsIO.fDeb
    else :
        qDeb = False

    nCom = clsEOS.nComp
    iNeu = 0
    pObs = -1.0

    zDep = NP.zeros(nCom)   #-- Composition at this Depth
    zDew = NP.zeros(nCom)
    zBub = NP.zeros(nCom)

    clsWRK = AD.classSample("GRDwrk")
    clsWRK.setIntComp(nCom,nCom)

    qBub = None
    pSat = None
    logK = NP.empty(nCom)
    
#--------------------------------------------------------------------
#  Simple Interval Halving
#--------------------------------------------------------------------

    iMid = 1
    zTyp = "GOC"
    zDep = zRef
    pDep = pRef

    while hBel - hAbv > 0.1 :

        hMid = 0.5*(hAbv + hBel)
        dDep =      hMid - hRef

        pDep,zDep,qDep,sDep = \
            calcTopBottom(zTyp,dDep,pRef,tRes,zRef,fRef,clsEOS,clsIO)

        if qDeb :
            if qDep : iBub =  1
            else    : iBub = -1
            sOut = "GOC search: iMid,hAbv,hMid,hBel,iBub,pSat {:2d} {:10.3f} {:10.3f} {:10.3f} {:2d} {:10.3f}\n" \
                   .format(iMid,hAbv,hMid,hBel,iBub,sDep)
            fDeb.write(sOut)
            
        if qDep :
            hBel = hMid
            zBub = zDep
        else    :
            hAbv = hMid
            zDew = zDep

        iMid += 1

#== Return values =====================================================

    #print("GOC: hMid,pMid,sMid {:10.3f} {:10.3f} {:10.3f}".format(hMid,pDep,sDep))

    iNeu = 0

    mSat,vSat,dSat,zSat,uSat,dumS,dumS = CE.calcProps(iNeu,sDep,tRes,zDep,clsEOS)

    zC1,zC7 = calcMoleFracGRD(zDep,clsEOS)

    return hMid,pDep,sDep,dSat,zC1,zC7,zDew,zBub

#========================================================================
#  Setup and Calculate GRD Fug & dFug/dP Coefficients (with Vol Shift)
#========================================================================

def setupCalcGRDCoefsP(iLiq,pRes,tRes,XY,clsEOS) :

    nCom = clsEOS.nComp
    iPhs = 0
    pByT = pRes/tRes
    recP =  1.0/pRes

    qP = True   #-- Onlt want Pres-derivatives, not Temp or Comp
    qT = False
    qX = False

#-- Calculate the log(Fugacity) Coefficients ------------------------

    f2,f2P,dumV,dumX = CE.calcPhaseFugPTX(iPhs,qP,qT,qX,pRes,tRes,XY,clsEOS)

    f2X = NP.exp(f2)        #-- Vec: exp(lnPhi)

#-- EoS B[i] Coefficients -------------------------------------------    

    BI = pByT*clsEOS.bI     #-- Vec: Bi*(p/RT)

    cI = BI*clsEOS.vS       #-- Vec: Ci = si*Bi
    XI = NP.exp(cI)         #-- Vec: exp(Ci)

    cIP = cI*recP           #-- Vec: Ci/p

#== 3-Parameter Fugacity and its P-Derivative =========================

    f3 = pRes*XY*f2X*XI     #-- Vec: p*xi*exp(lnPhi)*exp(Ci)

    f3P = recP + f2P + cIP

    f3P = f3*f3P

#== Return values =====================================================        

    return f3,f3P

#========================================================================
#  Calculate Z[C1] and Z[C7+] for GRD Experiment
#========================================================================

def calcMoleFracGRD(Z,clsEOS) :

    zC1 = 0.0 ; zC7 = 0.0

    for iC in range(clsEOS.nComp) :
        if clsEOS.gPP("MW",iC) < 25.0 : zC1 = zC1 + Z[iC]
        if clsEOS.gPP("MW",iC) > 90.0 : zC7 = zC7 + Z[iC]

#-- Return C1 & C7+ moles -------------------------------------------        

    return zC1,zC7

#========================================================================
#  End of Module
#========================================================================


