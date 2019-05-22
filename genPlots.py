#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Output Module (CM)
#
#========================================================================

#!/usr/bin/python3

import numpy             as NP
import os
import os.path           as OP
import matplotlib.pyplot as PL

import readExps  as RX

#========================================================================
#  Generate Plots
#========================================================================

def genPlots(clsIO,dicEXP,dicSAM,clsUNI) :

    nExp = len(dicEXP)
    iDep = 0

#--------------------------------------------------------------------
#   By default, write graphs to GraphOut sub-directory
#   If it doesn't already exist, create it!
#--------------------------------------------------------------------

    graphOut = "GraphOut"

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating Plots in Sub-Directory ",graphOut)

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            nObs = clsEXP.nObs
            nRow = clsEXP.nRow
            nSam = clsEXP.nSamp
            sNam = dicSAM[nSam].sName

            xTyp = clsEXP.xName

            titLS = RX.expOBSlong.get(xTyp)

            for iObs in range(nObs) :

                if clsEXP.qObs[iObs] and clsEXP.qPlt[iObs] :

                    xVec = NP.zeros(nRow)
                    yObs = NP.zeros(nRow)
                    yCal = NP.zeros(nRow)

                    for iRow in range(nRow) :

                        xVec[iRow] = clsUNI.I2X(clsEXP.dDep[iRow][iDep],clsEXP.uDep[iDep])
                        yObs[iRow] = clsUNI.I2X(clsEXP.dObs[iRow][iObs],clsEXP.uObs[iObs])
                        yCal[iRow] = clsUNI.I2X(clsEXP.dCal[iRow][iObs],clsEXP.uObs[iObs])

                    yObs = zero_to_nan(yObs)
                    yCal = zero_to_nan(yCal)

                    titL = xTyp + " " + titLS[iObs] + " of Sample " + sNam

                    sLobs = clsEXP.hObs[iObs] + "-Obs"
                    sLcal = clsEXP.hObs[iObs] + "-Cal"

                    xLab = clsEXP.hDep[iDep] + "/[" + clsEXP.uDep[iDep] + "]"
                    yLab = clsEXP.hObs[iObs] + "/[" + clsEXP.uObs[iObs] + "]"

                    figN = graphOut + "/" + xTyp + "_" + str(iExp+1) + "_" + clsEXP.hObs[iObs] + ".png"
        
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

    nExp = len(dicEXP0)
    iDep = 0
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

    if not os.path.exists(graphOut): os.makedirs(graphOut)

    print("Generating Plots in Sub-Directory ",graphOut)

    for iExp in range(nExp) :

        clsEXP0 = dicEXP0[iExp]
        clsEXP1 = dicEXP1[iExp]

        if clsEXP1.IsAct :

            nObs = clsEXP0.nObs
            nRow = clsEXP0.nRow
            nSam = clsEXP1.nSamp
            sNam = dicSAM1[nSam].sName

            xTyp = clsEXP0.xName

            titLS = RX.expOBSlong.get(xTyp)

            for iObs in range(nObs) :

                if clsEXP0.qObs[iObs] :

                    xVec = NP.zeros(nRow)
                    yObs = NP.zeros(nRow)
                    yBef = NP.zeros(nRow)
                    yAft = NP.zeros(nRow)

                    for iRow in range(nRow) :

                        xVec[iRow] = clsUNI.I2X(clsEXP0.dDep[iRow][iDep],clsEXP0.uDep[iDep])
                        yObs[iRow] = clsUNI.I2X(clsEXP0.dObs[iRow][iObs],clsEXP0.uObs[iObs])
                        yBef[iRow] = clsUNI.I2X(clsEXP0.dCal[iRow][iObs],clsEXP0.uObs[iObs])
                        yAft[iRow] = clsUNI.I2X(clsEXP1.dCal[iRow][iObs],clsEXP0.uObs[iObs])

#-- Switch zeros to nan's to prevent plotting -----------------------

                    yObs = zero_to_nan(yObs)
                    yBef = zero_to_nan(yBef)
                    yAft = zero_to_nan(yAft)

                    sObs = clsEXP0.hObs[iObs]

                    titL = xTyp + " " + titLS[iObs] + " of Sample " + sNam

                    sLobs  = clsEXP0.hObs[iObs] + "-Obs"
                    sLcal0 = clsEXP0.hObs[iObs] + "-Cal[Bef-" + sReg
                    sLcal1 = clsEXP0.hObs[iObs] + "-Cal[Aft-" + sReg

                    xLab = clsEXP0.hDep[iDep] + "/[" + clsEXP0.uDep[iDep] + "]"
                    yLab = clsEXP0.hObs[iObs] + "/[" + clsEXP0.uObs[iObs] + "]"

                    figN = graphFileName(graphOut,xTyp,iExp,sObs,sExt)

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
#  File Name for Graphics Output
#========================================================================

def graphFileName(sDir,sTyp,iExp,sObs,sExt) :

#-- Default File Name -----------------------------------------------    

    fNam = sDir + "/" + sTyp + "_" + str(iExp+1) + "_" + sObs + sExt + ".png"

#== If exists, modify until it doesn't ================================

    iExt = 2

    while OP.exists(fNam) :

        fNam  = sDir + "/" + sTyp + "_" + str(iExp+1) + "_" + sObs + sExt + "_" + str(iExt) + ".png"
        iExt += 1

#== Return value ======================================================

    return fNam
    
#========================================================================
#  To avoid plotting zero values, replace with 'nan' (Not A Number)
#========================================================================

def zero_to_nan(values) :

    return [float('nan') if x == 0.0 else x for x in values]

#========================================================================
#  End of Module
#========================================================================
