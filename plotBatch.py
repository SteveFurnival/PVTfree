#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Graphics Module (GP)
#
#========================================================================

#!/usr/bin/python3

import numpy             as NP
import os
import os.path           as OP
import matplotlib.pyplot as PL

import allData   as AD
import allProp   as AP
import readExps  as RX
import utilities as UT

#========================================================================
#  Generate Plots
#========================================================================

def genPlots(clsIO,dicEXP,dicSAM,clsUNI) :

    #return

    nExp = len(dicEXP)
    iInd = 0

#--------------------------------------------------------------------
#   By default, write graphs to GraphOut sub-directory
#   If it doesn't already exist, create it!
#--------------------------------------------------------------------

    graphOut = "GraphOut"
    rootN    = clsIO.rNam

    outRoot  = graphOut + "/" + rootN + "_"

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating Plots in Sub-Directory ",graphOut)

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            nObs = clsEXP.nObs
            nRow = clsEXP.nRow
            nSam = clsEXP.nSamp
            sNam = dicSAM[nSam].sNam

            xTyp = clsEXP.xName

            titLS = AP.classLIB().OBSlong.get(xTyp)

            for iObs in range(nObs) :

                if clsEXP.qObs[iObs] and clsEXP.qPlt[iObs] :

                    xVec = NP.zeros(nRow)
                    yObs = NP.zeros(nRow)
                    yCal = NP.zeros(nRow)

                    for iRow in range(nRow) :

                        xVec[iRow] = clsUNI.I2X(clsEXP.dInd[iInd][iRow],clsEXP.uInd[iInd])
                        yObs[iRow] = clsUNI.I2X(clsEXP.dObs[iObs][iRow],clsEXP.uObs[iObs])
                        yCal[iRow] = clsUNI.I2X(clsEXP.dCal[iObs][iRow],clsEXP.uObs[iObs])

                    yObs = UT.zero_to_nan(yObs)
                    yCal = UT.zero_to_nan(yCal)

                    titL = xTyp + " " + titLS[iObs] + " of Sample " + sNam

                    sLobs = clsEXP.hObs[iObs] + "-Obs"
                    sLcal = clsEXP.hObs[iObs] + "-Cal"

                    xLab = clsEXP.hInd[iInd] + "/[" + clsEXP.uInd[iInd] + "]"
                    yLab = clsEXP.hObs[iObs] + "/[" + clsEXP.uObs[iObs] + "]"

                    figN = outRoot + xTyp + "_" + str(iExp+1) + "_" + clsEXP.hObs[iObs] + ".png"
        
                    PL.plot(xVec,yObs,label=sLobs,marker="s",linestyle='None',color='r')
                    PL.plot(xVec,yCal,label=sLcal,color='r')
                
                    PL.title(titL)
                
                    PL.xlabel(xLab)
                    PL.ylabel(yLab)
                
                    PL.legend()
                
                    PL.savefig(figN)
                
                    PL.close()

    return

#========================================================================
#  Generate Plots: Before/After Regression or Grouping
#========================================================================

def regPlots(clsIO,dicEXP0,dicEXP1,dicSAM1,qReg,clsUNI) :

    #return

    nExp = len(dicEXP0)
    iInd = 0
    if qReg :
        sReg = "Reg]"
        sExt = "_REG"
    else    :
        sReg = "Grp]"
        sExt = "_GRP"

#--------------------------------------------------------------------
#   By default, write graphs to GraphOut sub-directory
#   If it doesn't already exist, create it!
#--------------------------------------------------------------------

    graphOut = "GraphOut"
    rootN    = clsIO.rNam

    outRoot  = graphOut + "/" + rootN + "_"

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating Plots in Sub-Directory ",graphOut)

    for iExp in range(nExp) :

        clsEXP0 = dicEXP0[iExp]
        clsEXP1 = dicEXP1[iExp]

        if clsEXP1.IsAct :

            nObs = clsEXP0.nObs
            nRow = clsEXP0.nRow
            nSam = clsEXP1.nSamp
            sNam = dicSAM1[nSam].sNam

            xTyp = clsEXP0.xName

            titLS = AP.classLIB().OBSlong.get(xTyp)

            for iObs in range(nObs) :

                if clsEXP0.qObs[iObs] :

                    xVec = NP.zeros(nRow)
                    yObs = NP.zeros(nRow)
                    yBef = NP.zeros(nRow)
                    yAft = NP.zeros(nRow)

                    for iRow in range(nRow) :

                        xVec[iRow] = clsUNI.I2X(clsEXP0.dInd[iInd][iRow],clsEXP0.uInd[iInd])
                        yObs[iRow] = clsUNI.I2X(clsEXP0.dObs[iObs][iRow],clsEXP0.uObs[iObs])
                        yBef[iRow] = clsUNI.I2X(clsEXP0.dCal[iObs][iRow],clsEXP0.uObs[iObs])
                        yAft[iRow] = clsUNI.I2X(clsEXP1.dCal[iObs][iRow],clsEXP0.uObs[iObs])

#-- Switch zeros to nan's to prevent plotting -----------------------

                    yObs = UT.zero_to_nan(yObs)
                    yBef = UT.zero_to_nan(yBef)
                    yAft = UT.zero_to_nan(yAft)

                    sObs = clsEXP0.hObs[iObs]

                    titL = xTyp + " " + titLS[iObs] + " of Sample " + sNam

                    sLobs  = clsEXP0.hObs[iObs] + "-Obs"
                    sLcal0 = clsEXP0.hObs[iObs] + "-Cal[Bef-" + sReg
                    sLcal1 = clsEXP0.hObs[iObs] + "-Cal[Aft-" + sReg

                    xLab = clsEXP0.hInd[iInd] + "/[" + clsEXP0.uInd[iInd] + "]"
                    yLab = clsEXP0.hObs[iObs] + "/[" + clsEXP0.uObs[iObs] + "]"

                    figN = graphFileName(outRoot,xTyp,iExp,sObs,sExt)

                    PL.plot(xVec,yObs,label=sLobs ,marker="s",linestyle='None',color='r')
                    PL.plot(xVec,yBef,label=sLcal0,color='r')
                    PL.plot(xVec,yAft,label=sLcal1,color='g')
                
                    PL.title(titL)
                
                    PL.xlabel(xLab)
                    PL.ylabel(yLab)
                
                    PL.legend()
                
                    PL.savefig(figN)
                
                    PL.close()

    return

#========================================================================
#  Blackoil Table Plots
#========================================================================

def blackPlots(dTab,eTab,fOil,fGas,clsBLK,clsUNI,clsIO) :

#-- Sample, Units and Sim-Type --------------------------------------    

    sNam = clsBLK.sNam  #-- Sample Name
    OutU = clsBLK.OutU  #-- FLD or MET
    sTyp = clsBLK.tSim  #-- Sim Type

#== Units =============================================================    

    if OutU[:3] == "MET" :
        pUni = "bara"
        RsUn = "sm3/sm3"
        BoUn = "rm3/sm3"
        UoUn = "cP"
        RvUn = "sm3/sm3"
        BgUn = "rm3/sm3"
        UgUn = "cP"
    else :
        pUni = "psia"
        RsUn = "Mscf/stb"
        BoUn = "rb/stb"
        UoUn = "cP"
        RvUn = "stb/Mscf"
        BgUn = "rb/Mscf"
        UgUn = "cP"

#-- All Plots use Pressure on X-Axis --------------------------------        

    iX   = clsBLK.iPr
    xUni = pUni
    xLab = "Pressure/[" + pUni + "]"

#== Which Keywords did we generate? ===================================

    iOil = -1
    iGas = -1

    if   sTyp[:3] == "ECL" :
        
        simT = "E100"
        oTyp = clsBLK.oTyp ; gTyp = clsBLK.gTyp
        
        if   oTyp == "PVTO" : iOil = 0
        elif oTyp == "PVDO" : iOil = 1
        elif oTyp == "PVCO" : iOil = 2
        
        if   gTyp == "PVTG" : iGas = 0
        elif gTyp == "PVDG" : iGas = 1
        
    elif sTyp[:3] == "TEM" :
        
        simT = "MORE"
        oTyp = clsBLK.oTyp ; gTyp = clsBLK.gTyp
        
        if   oTyp == "PVTO" : iOil = 0
        elif oTyp == "OPVD" : iOil = 1
        elif oTyp == "OPVT" : iOil = 2
        
        if   gTyp == "PVTG" : iGas = 0
        elif gTyp == "GPVT" : iGas = 1
        
    elif sTyp[:3] == "VIP" :
        
        simT = "VIP "
        oTyp = clsBLK.oTyp ; gTyp = clsBLK.gTyp
        
        if   oTyp == "BOOTAB" : iOil = 0
        elif oTyp == "BODTAB" : iOil = 1
        
        if   gTyp == "BOGTAB" : iGas = 0
        elif gTyp == "BDGTAB" : iGas = 1
        
    elif sTyp[:3] == "CMG" :
        
        simT = "IMEX"
        oTyp = clsBLK.oTyp ; gTyp = ""
        if   oTyp == "PVTG"    : iGas = 1
        elif oTyp == "PVT"     : iOil = 0 ; iGas = 1
        elif oTyp == "PVTCOND" : iOil = 1 ; iGas = 0
        elif oTyp == "PVTVO"   : iOil = 0 ; iGas = 1

    #print("blackPlots: simT,oTyp,gTyp,iOil,iGas ",simT,oTyp,gTyp,iOil,iGas)
        
#----------------------------------------------------------------------
#   By default, write graphs to GraphOut sub-directory
#   If it doesn't already exist, create it!
#----------------------------------------------------------------------

    graphOut = "GraphOut"
    rootN    = clsIO.rNam

    outRoot  = graphOut + "/" + rootN + "_"

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating BlackOil Plots in Sub-Directory ",graphOut)

#----------------------------------------------------------------------
#  GOR & CGR Plot
#----------------------------------------------------------------------

    if iOil == 0 and iGas == 0 :

        iY1 = clsBLK.iRs ; Y1Un = RsUn
        iY2 = clsBLK.iRv ; Y2Un = RvUn

        xX,Y1,Y2 = vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI)

#-- Create the GOR/CGR Plot -----------------------------------------
    
        titL  = simT + " GOR/CGR BlackOil Tables for Sample " + sNam
        y1Lab = "GOR/[" + RsUn + "]" ; y2Lab = "CGR/[" + RvUn + "]"
        s1Lab = "GOR"                ; s2Lab = "CGR"

        figN = outRoot + simT + "_GOR_CGR.png"    

        createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN)

    elif iOil == 0 :

        iY1 = clsBLK.iRs ; Y1Un = RsUn
        iY2 = -1         ; Y2Un = ""

        xX,Y1,Y2 = vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI)

        titL  = simT + " GOR BlackOil Tables for Sample " + sNam
        y1Lab = "GOR/[" + RsUn + "]" ; y2Lab = ""
        s1Lab = "GOR"                ; s2Lab = ""

        figN = outRoot + simT + "_GOR.png"    

        createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN)

    elif iGas == 0 :

        iY1 = clsBLK.iRv ; Y1Un = RvUn
        iY2 = -1         ; Y2Un = ""

        xX,Y1,Y2 = vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI)

        titL  = simT + " CGR BlackOil Tables for Sample " + sNam
        y1Lab = "CGR/[" + RvUn + "]" ; y2Lab = ""
        s1Lab = "CGR"                ; s2Lab = ""

        figN = outRoot + simT + "_CGR.png"    

        createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN)

#----------------------------------------------------------------------
#  Oil FVF & Viscosity
#----------------------------------------------------------------------

    if iOil >= 0 :

        iY1 = clsBLK.iBo ; Y1Un = BoUn
        iY2 = clsBLK.iUo ; Y2Un = UoUn

        xX,Y1,Y2 = vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI)

#-- Create the Oil FVF/Viscosity Plot -------------------------------
    
        titL  = simT + " Oil FVF & Visc for Sample " + sNam
        y1Lab = "Oil FVF/["  + BoUn + "]"
        y2Lab = "Oil Visc/[" + UoUn + "]"
        s1Lab = "Bo"
        s2Lab = "Uo"

        figN = outRoot + simT + "_Oil_FVF_VIsc.png"    

        createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN)

#----------------------------------------------------------------------
#  Gas FVF & Viscosity
#----------------------------------------------------------------------

    if iGas >= 0 :

        iY1 = clsBLK.iBg ; Y1Un = BgUn
        iY2 = clsBLK.iUg ; Y2Un = UgUn

        xX,Y1,Y2 = vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI)

#-- Create the Gas FVF/Viscosity Plot -------------------------------
    
        titL  = simT + " Gas FVF & Visc for Sample " + sNam
        y1Lab = "Gas FVF/["  + BgUn + "]"
        y2Lab = "Gas Visc/[" + UgUn + "]"
        s1Lab = "Bg"
        s2Lab = "Ug"

        figN = outRoot + simT + "_Gas_FVF_Visc.png"    

        createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN)

#----------------------------------------------------------------------
#  Under-Saturated Oil and Gas FVF & Viscosity
#----------------------------------------------------------------------

    if iOil == 0 : underSatOilPlot(simT,fOil,clsBLK,clsUNI,clsIO)

    if iGas == 0 : underSatGasPlot(simT,fGas,clsBLK,clsUNI,clsIO)

#== No return value ===================================================

    return

#========================================================================
#  Create the Blackoil Plot
#========================================================================

def createBlackPlot(xX,Y1,Y2,xLab,y1Lab,y2Lab,s1Lab,s2Lab,titL,figN) :

    fiG,aX1 = PL.subplots()  #-- Y1-Axis

    color = "tab:green"
    aX1.set_xlabel(xLab )
    aX1.set_ylabel(y1Lab,color=color)
    aX1.plot(xX,Y1,color=color,label=s1Lab)
    aX1.tick_params(axis="y",labelcolor=color)

    if y2Lab != "" :

        aX2 = aX1.twinx()        #-- Y2-Axis

        color = "tab:red"
        aX2.set_ylabel(y2Lab,color=color)
        aX2.plot(xX,Y2,color=color)
        aX1.plot(NP.nan,"-r",label=s2Lab)           #-- This adds Legend to AX1
        aX2.tick_params(axis="y",labelcolor=color)

    aX1.legend()             #-- Add legend to Y1-Axis
    PL.title(titL)
    fiG.tight_layout()
    PL.savefig(figN)
    PL.close()

#== No return value ===================================================

    return

#========================================================================
#  Create Vectors for Blackoil Plots
#========================================================================

def vecBlackPlot(dTab,eTab,iX,iY1,iY2,xUni,Y1Un,Y2Un,clsUNI) :

#== Dimension arrays ==================================================

    nSat = len(dTab)
    nExt = len(eTab)

    nVec = nSat + nExt

    iSat = nSat - 1
    iExt = nExt - 1
    iVec = 0

#== Saturated Data ====================================================

    xX = NP.zeros(nVec)
    Y1 = NP.zeros(nVec)
    Y2 = NP.zeros(nVec)

    while iSat >= 0 :

        xX[iVec] = clsUNI.I2X(dTab[iSat][iX ],xUni)
        Y1[iVec] = clsUNI.I2X(dTab[iSat][iY1],Y1Un)

        if iY2 > -1 : Y2[iVec] = clsUNI.I2X(dTab[iSat][iY2],Y2Un)

#-- Decrement iSat counter ------------------------------------------

        iSat -= 1
        iVec += 1

#== Extended Data =====================================================

    while iExt >= 0 :

        xX[iVec] = clsUNI.I2X(eTab[iExt][iX ],xUni)
        Y1[iVec] = clsUNI.I2X(eTab[iExt][iY1],Y1Un)
        if iY2 > -1 : Y2[iVec] = clsUNI.I2X(eTab[iExt][iY2],Y2Un)

#-- Decrement iExt counter ------------------------------------------

        iExt -= 1
        iVec += 1

#-- Switch zeros to nan's to prevent plotting -----------------------

    xX = UT.zero_to_nan(xX)
    Y1 = UT.zero_to_nan(Y1)
    Y2 = UT.zero_to_nan(Y2)

#== return vectors ====================================================

    return xX,Y1,Y2

#========================================================================
#  Under-Saturated Oil FVF & Viscosity Plot
#========================================================================

def underSatOilPlot(simT,fOil,clsBLK,clsUNI,clsIO) :

    sNam = clsBLK.sNam  #-- Sample Name
    OutU = clsBLK.OutU  #-- FLD or MET

    if OutU[:3] == "MET" :
        pUni = "bara"
        RsUn = "sm3/sm3"
        BoUn = "rm3/sm3"
        UoUn = "cP"
    else :
        pUni = "psia"
        RsUn = "Mscf/stb"
        BoUn = "rb/stb"
        UoUn = "cP"

#----------------------------------------------------------------------
#  Oil FVF
#----------------------------------------------------------------------

    titL = simT + " Undersaturated Oil FVF for Sample " + sNam
    xLab = "Pressure/[" + pUni + "]"
    yLab = "Oil FVF/["  + BoUn + "]"

    figN = "graphOut/" + clsIO.rNam + "_" + simT + "_UnderSatOilFVF.png"

#== Number of Under-Saturated Lines ===================================    

    nSet = len(fOil)

    iP  = 0
    iBo = 1
    iRs = 3

#== Loop over sets ====================================================

    for iSet in range(nSet) :

        fSet = fOil[iSet]
        nRow = len(fSet)

        xX = NP.zeros(nRow)
        Y1 = NP.zeros(nRow)

        for iRow in range(nRow) :

            dRow = fSet[iRow]

            xX[iRow] = clsUNI.I2X(dRow[iP ],pUni)
            Y1[iRow] = clsUNI.I2X(dRow[iBo],BoUn)

            if iRow == 0 :
                Rs  = clsUNI.I2X(dRow[iRs],RsUn)
                sRs = "Rs={:4.2f}".format(Rs)

#== Plot this line ====================================================            

        PL.plot(xX,Y1,label=sRs)

#== Finish the setup ==================================================        
                
    PL.title(titL)
    PL.xlabel(xLab)
    PL.ylabel(yLab)
    PL.legend()
    PL.savefig(figN)
    PL.close()

#----------------------------------------------------------------------
#  Oil Viscosity
#----------------------------------------------------------------------

    titL = simT + " Undersaturated Oil Viscosity for Sample " + sNam
    yLab = "Oil Visc/["  + UoUn + "]"

    figN = "graphOut/" + clsIO.rNam + "_" + simT + "_UnderSatOilVisc.png"

#== Number of Under-Saturated Lines ===================================    

    nSet = len(fOil)

    iUo = 2

#== Loop over sets ====================================================

    for iSet in range(nSet) :

        fSet = fOil[iSet]
        nRow = len(fSet)

        xX = NP.zeros(nRow)
        Y1 = NP.zeros(nRow)

        for iRow in range(nRow) :

            dRow = fSet[iRow]

            xX[iRow] = clsUNI.I2X(dRow[iP ],pUni)
            Y1[iRow] = clsUNI.I2X(dRow[iUo],UoUn)

            if iRow == 0 :
                Rs  = clsUNI.I2X(dRow[iRs],RsUn)
                sRs = "Rs={:4.2f}".format(Rs)

#== Plot this line ====================================================            

        PL.plot(xX,Y1,label=sRs)

#== Finish the setup ==================================================        
                
    PL.title(titL)
    PL.xlabel(xLab)
    PL.ylabel(yLab)
    PL.legend()
    PL.savefig(figN)
    PL.close()

#== No return value ===================================================

    return

#========================================================================
#  Under-Saturated Gas FVF & Viscosity Plot
#========================================================================

def underSatGasPlot(simT,fGas,clsBLK,clsUNI,clsIO) :

    sNam = clsBLK.sNam  #-- Sample Name
    OutU = clsBLK.OutU  #-- FLD or MET

    if OutU[:3] == "MET" :
        pUni = "bara"
        RvUn = "sm3/sm3"
        BgUn = "rm3/sm3"
        UgUn = "cP"
    else :
        pUni = "psia"
        if simT == "IMEX" : RvUn = "stb/scf"
        else              : RvUn = "stb/Mscf"
        BgUn = "rb/Mscf"
        UgUn = "cP"

    if simT == "IMEX" :
        xTxt = "Pressure/["
        xUni = pUni
        uLab = RvUn
        mulL = 1.0E+6
        sLab = "Rv={:5.1f}"
    else :
        xTxt = "CGR/["
        xUni = RvUn
        uLab = pUni
        mulL = 1.0
        sLab = "P={:5.0f}"

#----------------------------------------------------------------------
#  Gas FVF
#----------------------------------------------------------------------

    titL = simT + " Undersaturated Gas FVF for Sample " + sNam
    xLab = xTxt        + xUni + "]"
    yLab = "Gas FVF/[" + BgUn + "]"

    figN = "graphOut/" + clsIO.rNam + "_" + simT + "_UnderSatGasFVF.png"

#== Number of Under-Saturated Lines ===================================    

    nSet = len(fGas)

    iRv = 0
    iBg = 1
    iPr = 3

#== Loop over sets ====================================================

    for iSet in range(nSet) :

        fSet = fGas[iSet]
        nRow = len(fSet)

        xX = NP.zeros(nRow)
        Y1 = NP.zeros(nRow)

        for iRow in range(nRow) :

            dRow = fSet[iRow]

            xX[iRow] = clsUNI.I2X(dRow[iRv],xUni)
            Y1[iRow] = clsUNI.I2X(dRow[iBg],BgUn)

            if len(dRow) == 4 :
                Pr  = clsUNI.I2X(dRow[iPr],uLab)
                sPr = sLab.format(mulL*Pr)

#== Plot this line ====================================================            

        PL.plot(xX,Y1,label=sPr)

#== Finish the setup ==================================================        
                
    PL.title(titL)
    PL.xlabel(xLab)
    PL.ylabel(yLab)
    PL.legend()
    PL.savefig(figN)
    PL.close()

#----------------------------------------------------------------------
#  Gas Viscosity
#----------------------------------------------------------------------

    titL = simT + " Undersaturated Gas Viscosity for Sample " + sNam
    yLab = "Gas Visc/["  + UgUn + "]"

    figN = "graphOut/" + clsIO.rNam + "_" + simT + "_UnderSatGasVisc.png"

#== Number of Under-Saturated Lines ===================================    

    nSet = len(fGas)

    iUg = 2

#== Loop over sets ====================================================

    for iSet in range(nSet) :

        fSet = fGas[iSet]
        nRow = len(fSet)

        xX = NP.zeros(nRow)
        Y1 = NP.zeros(nRow)

        for iRow in range(nRow) :

            dRow = fSet[iRow]

            xX[iRow] = clsUNI.I2X(dRow[iRv],xUni)
            Y1[iRow] = clsUNI.I2X(dRow[iUg],UgUn)

            if len(dRow) == 4 :
                Pr  = clsUNI.I2X(dRow[iPr],uLab)
                sPr = sLab.format(mulL*Pr)

#== Plot this line ====================================================            

        PL.plot(xX,Y1,label=sPr)

#== Finish the setup ==================================================        
                
    PL.title(titL)
    PL.xlabel(xLab)
    PL.ylabel(yLab)
    PL.legend()
    PL.savefig(figN)
    PL.close()

#== No return value ===================================================

    return

#========================================================================
#  File Name for Graphics Output
#========================================================================

def graphFileName(sDir,sTyp,iExp,sObs,sExt) :

#-- Default File Name -----------------------------------------------    

    fNam = sDir + sTyp + "_" + str(iExp+1) + "_" + sObs + sExt + ".png"

#== If exists, modify until it doesn't ================================

    iExt = 2

    while OP.exists(fNam) :

        fNam  = sDir + sTyp + "_" + str(iExp+1) + "_" + sObs + sExt + "_" + str(iExt) + ".png"
        iExt += 1

#== Return value ======================================================

    return fNam
    
#========================================================================
#  End of Module
#========================================================================
