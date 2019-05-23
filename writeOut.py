#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Output Module (CM)
#
#========================================================================

#!/usr/bin/python3

import datetime
import numpy             as NP
import os
import os.path           as OP
import matplotlib.pyplot as PL

import calcReg   as CR
import readExps  as RX
import genPlots  as GP
import texttable as TT

#=======================================================================
#  Header, disclaimer, etc.
#=======================================================================

def outputHeader(fOut,sCom,clsIO) :

    sYear   = datetime.date.today().strftime("%Y")
    sMonth  = datetime.date.today().strftime("%B")
    sDayNam = datetime.date.today().strftime("%A")
    sDayNum = datetime.date.today().strftime("%d")
    
    sTime   = datetime.datetime.now().strftime("%H:%M:%S")

    fNam = fOut.name
    sTok = fNam.split(sep=".")
    nTok = len(sTok)

#          123456789012345678901234567890123456789012345678901234567890
    fOut.write("\n")
    fOut.write(sCom + "============================================================\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  PVTfree\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  Equation of State PVT Program\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  Copyright (C) 2019\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  Steve Furnival, HoBoil Limited, Banchory, UK, March 2019\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  PVTfree is free software: you can redistribute it and/or modify\n")
    fOut.write(sCom + "  it under the terms of the GNU General Public License as published by\n")
    fOut.write(sCom + "  the Free Software Foundation, either version 3 of the License, or\n")
    fOut.write(sCom + "  (at your option) any later version.\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  PVTfree is distributed in the hope that it will be useful,\n")
    fOut.write(sCom + "  but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
    fOut.write(sCom + "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
    fOut.write(sCom + "  GNU General Public License for more details.\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  You should have received a copy of the GNU General Public License\n")
    fOut.write(sCom + "  along with PVTfree.  If not, see <https://www.gnu.org/licenses/>.\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "============================================================\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "  Run on " + sDayNam + " "
                                  + sDayNum + " "
                                  + sMonth  + " "
                                  + sYear   + " at "
                                  + sTime   + "\n")
    fOut.write(sCom + "\n")
    fOut.write(sCom + "============================================================\n")
    fOut.write(sCom + "\n")

    outputTitle(fOut,sCom,clsIO)

    return

#=======================================================================
#  Output Titles
#=======================================================================

def outputTitle(fOut,sCom,clsIO) :

    if clsIO.Tit1 != None :
        fOut.write(sCom + "\n")
        fOut.write(clsIO.Tit1+"\n")
        if clsIO.Tit2 != None :
            fOut.write(clsIO.Tit2+"\n")
            if clsIO.Tit3 != None :
                fOut.write(clsIO.Tit3+"\n")
        fOut.write(sCom + "\n")
        fOut.write("\n")

    return

#=======================================================================
#  Output EOS Type
#=======================================================================

def outputEOS(fOut,sCom,clsIO,clsEOS) :

    if clsEOS.EOS == "SRK" :
        sEOS = "Soave-Redlich-Kwong (SRK)"
    else :
        sEOS = "Peng-Robinson (PR)"

    fOut.write("\n")
    fOut.write(sCom + "============================================================\n")
    fOut.write(sCom + "  Equation of State Selected = " + sEOS + "\n")
    fOut.write(sCom + "============================================================\n")
    fOut.write("\n")

    return

#=======================================================================
#  Output Header for Grouping Operation
#=======================================================================

def outputGroup(clsIO,nOld,nNew) :

    fOut = clsIO.fOut

    sLabl = "  Grouping From " + str(nOld) + " To " + str(nNew) + " Components\n"
    fOut.write("\n")
    fOut.write("============================================================\n")
    fOut.write(sLabl)
    fOut.write("============================================================\n")
    fOut.write("\n")

#=======================================================================
#  Output Components Properties, BIPs and Composition to File
#=======================================================================

def outputProps(clsIO,clsEOS,dicSAM,sTit) :

    fOut = clsIO.fOut

    nComp = clsEOS.NC
    nSamp = clsEOS.NS

    fOut.write("\n")
    fOut.write("#############################################\n")
    sOut = "## " + sTit + "\n"
    fOut.write(sOut)
    fOut.write("#############################################\n")
    fOut.write("\n")
    fOut.write("============================================================\n")
    fOut.write("  Component Properties\n")
    fOut.write("============================================================\n")

#-- Plus Fraction Properties ----------------------------------------

    fOut.write("\n")
    fOut.write("Plus Fraction (C7+) Properties\n")
    fOut.write("\n")

    sOut = "Sample    "
    for iSamp in range(nSamp) :
        sNam = dicSAM[iSamp].sName
        sNam = "{:>8}  ".format(sNam)
        sOut = sOut + sNam
    sOut = sOut + "\n"
    fOut.write(sOut)
    fOut.write("\n")

    sOut = "Mol Wt    "
    for iSamp in range(nSamp) :
        Mw  = dicSAM[iSamp].mPlus
        sMw = "{:8.3f}  ".format(Mw)
        sOut  = sOut + sMw + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

    sOut = "Spec Grv  "
    for iSamp in range(nSamp) :
        SG  = dicSAM[iSamp].sPlus
        sSG = "{:8.5f}  ".format(SG)
        sOut  = sOut + sSG + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

    sOut = "Alpha     "
    for iSamp in range(nSamp) :
        alf = dicSAM[iSamp].aPlus
        sAl = "{:8.5f}  ".format(alf)
        sOut  = sOut + sAl + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

#-- First Set of Properties -----------------------------------------    

    tab = TT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","f","f","f","f","f","f"]
    tabAlgn = ["l","r","r","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ,11 ,11 ]

    tab.set_cols_dtype(tabType)
    tab.set_cols_align(tabAlgn)
    tab.set_cols_width(tabColW)
    tab.set_precision(5)

    tOut = [[0.0 for jC in range(7)] for iC in range(nComp+1)]

#               12345678901--12345678901
    tabHead = ["Comp       \nName       ",
               "  Mol Wt   \n [lb/lbmol]",
               "  Tcrit    \n   [degR]  ",
               "  Pcrit    \n   [psia]  ",
               "  Vcrit    \n[ft3/lbmol]"    ,
               "  Zcrit    \n   [dim]   ",
               "  AcenF    \n   [dim]   "]

    tOut[0] = tabHead
    
    for iC in range(nComp) :

        cName = clsEOS.gPP("CN",iC)
        molWt = clsEOS.gPP("MW",iC) ; Tcrit = clsEOS.gPP("TC",iC)
        Pcrit = clsEOS.gPP("PC",iC) ; Vcrit = clsEOS.gPP("VC",iC)
        Zcrit = clsEOS.gPP("ZC",iC) ; AcenF = clsEOS.gPP("AF",iC)
        
        tRow = [cName,molWt,Tcrit,Pcrit,Vcrit,Zcrit,AcenF]

        tOut[iC+1] = tRow

    tab.add_rows(tOut)

    fOut.write("\n")
    fOut.write(tab.draw())
    fOut.write("\n")

#-- Second set of properties ------------------------------------------    
 
    tab = TT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","f","f","f","f","f","f"]
    tabAlgn = ["l","r","r","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ,11 ,11 ]

    tab.set_cols_dtype(tabType)
    tab.set_cols_align(tabAlgn)
    tab.set_cols_width(tabColW)
    tab.set_precision(5)

    tOut = [[0 for jC in range(7)] for iC in range(nComp+1)]

#               12345678901--12345678901
    tabHead = ["Comp       \nName       ",
               "  Tboil    \n   [degR]  ",
               " Spec Grav \n[wrt Wat=1]",
               " Parachor  \n   [dim]   ",
               " Vol Shift \n   [dim]   ",
               "Mult-OmegaA\n   [dim]   ",
               "Mult-OmegaB\n   [dim]   "]

    tOut[0] = tabHead
    
    for iC in range(nComp) :

        cName = clsEOS.gPP("CN",iC)
        Tboil = clsEOS.gPP("TB",iC) ; SpecG = clsEOS.gPP("SG",iC)
        ParaC = clsEOS.gPP("PA",iC) ; volSh = clsEOS.gPP("SS",iC)
        MultA = clsEOS.gPP("MA",iC) ; MultB = clsEOS.gPP("MB",iC)
        
        tRow = [cName,Tboil,SpecG,ParaC,volSh,MultA,MultB]

        tOut[iC+1] = tRow

    tab.add_rows(tOut)

    fOut.write("\n")
    fOut.write(tab.draw())
    fOut.write("\n")

#-- Third set of properties -------------------------------------------
 
    tab = TT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","e","e","e","e"]
    tabAlgn = ["l","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ]

    tab.set_cols_dtype(tabType)
    tab.set_cols_align(tabAlgn)
    tab.set_cols_width(tabColW)
    tab.set_precision(3)

    tOut = [[0 for jC in range(5)] for iC in range(nComp+1)]

#               12345678901--12345678901
    tabHead = ["Comp       \nName       ",
               "   IdCpA   \nBtu/lbmo.R ",
               "   IdCpB   \nBtu/lbmo.R2",
               "   IdCpC   \nBtu/lbmo.R3",
               "   IdCpD   \nBtu/lbmo.R4"]

    tOut[0] = tabHead
    
    for iC in range(nComp) :

        cName = clsEOS.gPP("CN",iC)
        IdCpA = clsEOS.gPP("CA",iC) ; IdCpB = clsEOS.gPP("CB",iC)
        IdCpC = clsEOS.gPP("CC",iC) ; IdCpD = clsEOS.gPP("CD",iC)
        
        tRow = [cName,IdCpA,IdCpB,IdCpC,IdCpD]

        tOut[iC+1] = tRow

    tab.add_rows(tOut)

    fOut.write("\n")
    fOut.write(tab.draw())
    fOut.write("\n")

#-- BICs ------------------------------------------------------------

    fOut.write("\n")
    fOut.write("============================================================\n")
    fOut.write("  Binary Iteraction Parameters\n")
    fOut.write("============================================================\n")

    nSETs = 0
    n1SET = 8
    nSumm = nComp

    nSETa = []

    while (nSumm-n1SET) > 0 :
        nSETs += 1
        nSETa.append(n1SET)
        nSumm = nSumm - n1SET

    nSETs += 1
    nSETa.append(nSumm)

    iHead = 0
    nOffS = 0

    for iS in range(nSETs) :
        
        tab = TT.Texttable()
        tab.set_deco(tab.HEADER)
    
        nBICs = nSETa[iS]

        #print("O/P BICs: iS,nBICs ",iS,nBICs)

        tabType = [] ; tabAlgn = [] ; tabColW = []
        tabHead = [] ; tOut    = []

        tabType.append("t")
        tabAlgn.append("l")
        tabColW.append( 8 )
#                       12345678
        tabHead.append("Comp    ")
        
        tOut = [[0 for jC in range(nBICs+1)] for iC in range(nComp+1)]

        for iC in range(nBICs) :

            cName = clsEOS.gPP("CN",iHead)
            iHead += 1
            
            tabType.append("f")
            tabAlgn.append("r")
            tabColW.append( 8 )
            tabHead.append(cName)

        tab.set_cols_dtype(tabType)
        tab.set_cols_align(tabAlgn)
        tab.set_cols_width(tabColW)
        tab.set_precision(5)

        tOut[0] = tabHead

        for iR in range(nComp) :

            cName = clsEOS.gPP("CN",iR)

            tRow = [cName]

            for jC in range(nBICs) :
                jOffS = nOffS + jC
                tRow.append(clsEOS.gIJ(iR,jOffS))

            tOut[iR+1] = tRow
    
        tab.add_rows(tOut)

        fOut.write("\n")
        fOut.write(tab.draw())
        fOut.write("\n")

        nOffS = nOffS + n1SET

#-- Compositions ------------------------------------------------------

    fOut.write("\n")
    fOut.write("============================================================\n")
    fOut.write("  Sample Compositions\n")
    fOut.write("============================================================\n")

    nSETs = 0
    n1SET = 7
    nSumm = nSamp

    nSETa = []

    while (nSumm-n1SET) > 0 :
        nSETs += 1
        nSETa.append(n1SET)
        nSumm = nSumm - n1SET

    nSETs += 1
    nSETa.append(nSumm)

    #print("nSETs,nSETa ",nSETs,nSETa)

    iHead = 0
    nOffS = 0

    for iS in range(nSETs) :
        
        tab = TT.Texttable()
        tab.set_deco(tab.HEADER)
    
        nSAMs = nSETa[iS]

        tabType = [] ; tabAlgn = [] ; tabColW = []
        tabHead = [] ; tOut    = []

        tabType.append("t")
        tabAlgn.append("l")
        tabColW.append( 8 )
#                       12345678
        tabHead.append("Com|Samp")
        
        tOut = [[0 for jC in range(nSAMs+1)] for iC in range(nComp+3)]

        for jS in range(nSAMs) :

            sName = dicSAM[iHead].sName
            iHead += 1
            
            tabType.append("f")
            tabAlgn.append("r")
            tabColW.append(10 )
            tabHead.append(sName)

        tab.set_cols_dtype(tabType)
        tab.set_cols_align(tabAlgn)
        tab.set_cols_width(tabColW)
        tab.set_precision(5)

        tOut[0] = tabHead

        for iR in range(nComp) :

            cName = clsEOS.gPP("CN",iR)

            tRow = [cName]

            for jS in range(nSAMs) :
                jOffS = nOffS + jS
                tRow.append(100.0*dicSAM[jOffS].gZI(iR))  #-- Output as %

            tOut[iR+1] = tRow

        tRow = [" "]
        for jS in range(nSAMs) : tRow.append(" ")

        tOut[iR+2] = tRow

        tRow = ["Total"]
        for jS in range(nSAMs) : tRow.append(100.0*dicSAM[nOffS + jS].gTot())

        tOut[iR+3] = tRow
    
        tab.add_rows(tOut)

        fOut.write("\n")
        fOut.write(tab.draw())
        fOut.write("\n")

        nOffS = nOffS + n1SET

    return

#========================================================================
#  Output the Experiments Results
#========================================================================

def outputExps(clsIO,dicEXP,dicSAM,clsUNI) :

    nExp = len(dicEXP)

    fOut = clsIO.fOut

#== Loop over experiments =============================================

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            xTyp = clsEXP.xName

            #print("Experiment: iExp,Type ",iExp+1,xTyp)

            nSam = clsEXP.nSamp
            sNam = dicSAM[nSam].sName

            if   xTyp == "CCE" :
                sLabl = "  Constant Composition Expansion (CCE) Experiment\n"
            elif xTyp == "CVD" :
                sLabl = "  Constant Volume Depletion (CVD) Experiment\n"
            elif xTyp == "DLE" :
                sLabl = "  Differential Liberation (DLE) Experiment\n"
            elif xTyp == "SEP" :
                sLabl = "  Separator Test (SEP) Experiment\n"
            elif xTyp == "FLS" :
                sLabl = "  2-Phase Flash (FLS) Experiment\n"
            elif xTyp == "SAT" :
                sLabl = "  Saturation Pressure (SAT) Experiment\n"
            elif xTyp == "SWL" :
                sLabl = "  Swelling Test (SWL) Experiment\n"
            elif xTyp == "GRD" :
                sLabl = "  Composition Versus Depth (GRD) Experiment\n"

            fOut.write("\n")
            fOut.write("============================================================\n")
            fOut.write(sLabl)
            fOut.write("============================================================\n")
            fOut.write("\n")

            sOut = "  Experiment Number "+str(iExp+1)+" run on Sample["+str(nSam+1)+"] = "+sNam+"\n"
            fOut.write(sOut)
            fOut.write("\n")

            if xTyp == "CCE" or xTyp == "CVD" or xTyp == "DLE" or \
               xTyp == "SWL" or xTyp == "GRD"                     :

                dOut = "{:8.3f}".format(clsUNI.I2X(clsEXP.Tres,clsEXP.Tuni))
                sOut = "  At Reservoir Temperature of               "+dOut+" "+clsEXP.Tuni+"\n"
                fOut.write(sOut)
                fOut.write("\n")

            if xTyp == "CCE" or xTyp == "CVD" or xTyp == "DLE" or \
               xTyp == "SEP"                                      :

                PsatU = clsEXP.PsatU
                dPobs = "{:10.3f}".format(clsUNI.I2X(clsEXP.PsatO,clsEXP.PsatU))
                dPcal = "{:10.3f}".format(clsUNI.I2X(clsEXP.PsatC,clsEXP.PsatU))

                sPobs = "  Observed   Saturation Pressure (Psat) = "+dPobs+" "+PsatU+"\n"
                sPcal = "  Calculated Saturation Pressure (Psat) = "+dPcal+" "+PsatU+"\n"

                fOut.write(sPobs)
                fOut.write(sPcal)
                fOut.write("\n")

            if xTyp == "SEP" : nSep = 1
            else             : nSep = 0

            nDep = clsEXP.nDep
            nObs = clsEXP.nObs
            nObU = clsEXP.nObU
            #nCol = clsEXP.nDep + 2*nObU
            nCol = nDep + nObU + nObs + nSep
            nRow = clsEXP.nRow

            #print("calcCCE: nCol,nRow ",nCol,nRow)

            tab = TT.Texttable()
            tab.set_deco(tab.HEADER)

            tabType = ["f" for i in range(nCol)]
            tabAlgn = ["r" for i in range(nCol)]
            tabColW = [ 11 for i in range(nCol)]

            tab.set_cols_dtype(tabType)
            tab.set_cols_align(tabAlgn)
            tab.set_cols_width(tabColW)
            tab.set_precision(5)
    
            tOut = [[0.0 for j in range(nCol)] for i in range(nRow+1)]

            tabHead = []

            for iDep in range(nDep) :
                sDep = clsEXP.hDep[iDep] + "\n" + clsEXP.uDep[iDep]
                tabHead.append(sDep)

            for iSep in range(nSep) :
                sSep = " LSep  VSep \n Stg#  Stg# "
                tabHead.append(sSep)
        
            for iObs in range(nObs) :
                if clsEXP.qObs[iObs] :
                    sObs = clsEXP.hObs[iObs] + "-OBS\n" + clsEXP.uObs[iObs]
                    tabHead.append(sObs)
                if clsEXP.qCal[iObs] :
                    sCal = clsEXP.hCal[iObs] + "-CAL\n" + clsEXP.uCal[iObs]
                    tabHead.append(sCal)

            tOut[0] = tabHead

#-- Loop over rows of Experiment ------------------------------------

            for iRow in range(nRow) :

                tRow = []

                for iDep in range(nDep) :
                    dDep = clsUNI.I2X(clsEXP.dDep[iRow][iDep],clsEXP.uDep[iDep])
                    tRow.append(dDep)

                for iSep in range(nSep) :
                    sSep = str(clsEXP.Lsep[iRow]) + "     " + \
                           str(clsEXP.Vsep[iRow])
                    tRow.append(sSep)

                for iObs in range(nObs) :
                    if clsEXP.qObs[iObs] :
                        dObs = clsUNI.I2X(clsEXP.dObs[iRow][iObs],clsEXP.uObs[iObs])
                        tRow.append(dObs)
                    if clsEXP.qCal[iObs] :
                        dCal = clsUNI.I2X(clsEXP.dCal[iRow][iObs],clsEXP.uCal[iObs])
                        tRow.append(dCal)

                tOut[iRow+1] = tRow
        
            tab.add_rows(tOut)

            fOut.write(tab.draw())
            fOut.write("\n")

            if nSep > 0 :
                qNonS = False
                for iSep in range(2,nRow) :
                    if clsEXP.Vsep[iSep] > 0 :
                        qNonS = True
                        exit
                if qNonS :
                    fOut.write("\n")
                    sOut = "NB. Separator Train Has Non-Standard Routing so GORs cannot be summed\n"
                    fOut.write(sOut)

#======================================================================
#  Write Plots
#======================================================================

    if clsIO.Pall : GP.genPlots(clsIO,dicEXP,dicSAM,clsUNI)

    return

#========================================================================
#  BlackOil Tables, Header
#========================================================================

def outputHeaderBO(fSim,iSam,sNam,clsBLK,clsIO,clsUNI) :

    tSim = clsBLK.tSim
    xTyp = clsBLK.xTyp

    tRes = clsBLK.Tres
    Tuni = clsBLK.tResU.upper()

    pSepU = clsBLK.pSepU.upper()
    tSepU = clsBLK.tSepU.upper()

    pSep = clsBLK.pSep
    tSep = clsBLK.tSep

#-- Simulator Type --------------------------------------------------

    sNull = "\n"
    
    if   tSim == "CMG" :
        sComm = "**"
        sHead = "**  (CMG) IMEX Blackoil Tables\n"
        sLine = "**============================================================\n"
        sPace = "**\n"
    elif tSim == "MOR" :
        sComm = "--"
        sHead = "--  (Tempest) MORE BlackOil Tables\n"
        sLine = "--============================================================\n"
        sPace = "--\n"
    elif tSim == "ECL" :
        sComm = "--"
        sHead = "--  ECLIPSE 100 (E100) BlackOil Tables\n"
        sLine = "--============================================================\n"
        sPace = "--\n"
    elif tSim == "VIP" :
        sComm = "C "
        sHead = "C  VIP/Nexus Blackoil Tables\n"
        sLine = "C ============================================================\n"
        sPace = "C  \n"
    
    fSim.write(sPace)
    fSim.write(sLine)
    fSim.write(sHead)
    fSim.write(sLine)
    fSim.write(sPace)
    fSim.write(sNull)

    if   xTyp == "CCE" : sLabl = sComm + " Generated using a Constant Composition Expansion (CCE) Experiment\n"
    elif xTyp == "CVD" : sLabl = sComm + " Generated using a Constant Volume Depletion (CVD) Experiment\n"
    elif xTyp == "DLE" : sLabl = sComm + " Generated using a Differential Liberation (DLE) Experiment\n"
    fSim.write(sPace)
    fSim.write(sLabl)

    tOut = clsUNI.I2X(tRes,Tuni)
    sOut = "{:7.3f} ".format(tOut)
    sLabl = sComm + " At a Reservoir Temperature of " + sOut + Tuni +"\n"
    fSim.write(sPace)
    fSim.write(sLabl)

    sLabl = sComm + " Using Sample Number[" + str(iSam+1) + "] called " + sNam + " \n"
    fSim.write(sPace)
    fSim.write(sLabl)

    nSep = len(pSep)

    sOutP = ""
    sOutT = ""

    for iSep in range(nSep) :
        dPrs = clsUNI.I2X(pSep[iSep],pSepU)
        dTem = clsUNI.I2X(tSep[iSep],tSepU)
        sPrs = "{:8.3f} ".format(dPrs)
        sTem = "{:8.3f} ".format(dTem)
        sOutP = sOutP + sPrs
        sOutT = sOutT + sTem

    sLabl = sComm + " With a " + str(nSep) + "-Stage Separator Train of:\n"
    fSim.write(sPace)
    fSim.write(sLabl)
    fSim.write(sPace)
    sLabl = sComm + " Pressures    [" + pSepU + "]: " + sOutP + "\n"
    fSim.write(sLabl)
    sLabl = sComm + " Temperatures [" + tSepU + "]: " + sOutT + "\n"
    fSim.write(sLabl)
    fSim.write(sPace)
    fSim.write(sNull)

#==  Units ============================================================

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    if   OutU == "MET" : uUnit = "METRIC"
    elif OutU == "SI"  : uUnit = "SI"
    else   :             uUnit = "FIELD"

    sLabl = sComm + " Units set for the Tables is " + uUnit +"\n"

    fSim.write(sPace)
    fSim.write(sLabl)
    fSim.write(sPace)
    fSim.write(sNull)

#==  Density Keyword ==================================================
    
    if tSim == "ECL" or tSim == "MOR" :
        outputECLDensity(OutU,fSim,clsBLK,clsUNI)

#== Water Properties (PVTW) ===========================================

    if tSim == "ECL" or tSim == "MOR" :
        outputECLPVTW(OutU,fSim,clsBLK,clsUNI)
    elif tSim == "VIP" :
        outputVIPPVTW(OutU,fSim,clsBLK,clsUNI)        

#== CMG DENSITY, Compressibility and Water Phase Properties ===========

    if tSim == "CMG" :
        outputCMGgeneral(OutU,fSim,clsBLK,clsUNI)

#== No return value ===================================================    

    return

#========================================================================
#  CMG IMEX General Output
#========================================================================

def outputCMGgeneral(OutU,fSim,clsSIM,clsUNI) :

    tRes = clsSIM.Tres

#-- Stock Tank Densities --------------------------------------------
    
    dOil = clsSIM.dSTO
    dWat = clsSIM.dSTW
    dGas = clsSIM.dSTG

#-- Unit Conversions ------------------------------------------------

    if OutU == "MET" :   #-- Metric-Equivalent
        sComp = "1/kPa"
        sDens = "kg/m3"
        sGfvf = "rm3/sm3"
        sGOR  = "sm3/sm3"
        sOfvf = "rm3/sm3"
        sPres = "kPa"
        sTemp = "degC"
        sVisc = "mPa.s"
    else :              #-- Field-Units
        sComp = "1/psi"
        sDens = "lb/ft3"
        sGfvf = "rb/scf"
        sGOR  = "scf/stb"
        sOfvf = "rb/stb"
        sPres = "psia"
        sTemp = "degF"
        sVisc = "cP"

#-- Headers ---------------------------------------------------------    

    fSim.write("**\n")
    fSim.write("**-----------------------------------------\n")
    fSim.write("**  Oil, Water and Gas Phase Properties\n")
    fSim.write("**\n")
    fSim.write("**  Units:  Compressibility " + sComp + "\n")
    fSim.write("**          Density         " + sDens + "\n")
    fSim.write("**          Gas FVF         " + sGfvf + "\n")
    fSim.write("**          GOR             " + sGOR  + "\n")
    fSim.write("**          Oil FVF         " + sOfvf + "\n")
    fSim.write("**          Pressure        " + sPres + "\n")
    fSim.write("**          Temperature     " + sTemp + "\n")
    fSim.write("**          Viscosity       " + sVisc + "\n")
    fSim.write("**-----------------------------------------\n")
    fSim.write("**\n")
    fSim.write("\n")

    tRes = clsUNI.I2X(tRes,sTemp)
    sOut = "{:7.3f}".format(tRes)

    if   clsSIM.oTyp == "PVTCOND" :
        fSim.write("*MODEL  *GASWATER_WITH_CONDENSATE\n")
    elif clsSIM.oTyp == "PVTG" :
        fSim.write("*MODEL  *GASWATER\n")
    else :
        fSim.write("*MODEL  *BLACKOIL\n")
    fSim.write("*TRES    " + sOut + "\n")
    fSim.write("\n")

    dOil = clsUNI.I2X(dOil,sDens)
    dGas = clsUNI.I2X(dGas,sDens)
    dWat = clsUNI.I2X(dWat,sDens)
    
    sOil = " {:10.3f}".format(dOil)
    sGas = " {:10.5f}".format(dGas)
    sWat = " {:10.3f}".format(dWat)

    fSim.write("*DENSITY  *OIL    " + sOil + "\n")
    fSim.write("*DENSITY  *GAS      " + sGas + "\n")
    fSim.write("*DENSITY  *WATER  " + sWat + "\n")
    fSim.write("\n")

#-- Compressibilities -----------------------------------------------

    pRefW = clsSIM.pRefW  #-- Reference Pressure        [psia   or bara]
    bRefW = clsSIM.bRefW  #--           FVF             [rb/stb or rm3/sm3]
    uRefW = clsSIM.uRefW  #--           Viscosity       [cP]
    cRefW = clsSIM.cRefW  #--           Compressibility [1/psi  or 1/bar]
    vRefW = clsSIM.vRefW  #--           Viscosibility   [1/psi  or 1/bar]

    pRefW = clsUNI.I2X(pRefW,sPres)
    cRefW = clsUNI.I2X(cRefW,sComp)
    vRefW = clsUNI.I2X(vRefW,sVisc)
    
    sPrf = " {:10.3f}".format(pRefW)
    sBrf = " {:10.5f}".format(bRefW)
    sUrf = " {:10.5f}".format(uRefW)
    sCrf = " {:10.3e}".format(cRefW)
    sVrf = " {:10.3e}".format(vRefW)

    fSim.write("*REFPW            " + sPrf + "\n")
    fSim.write("*BWI              " + sBrf + "\n")
    fSim.write("*VWI              " + sUrf + "\n")
    fSim.write("*CW               " + sCrf + "\n")
    fSim.write("*CVW              " + sVrf + "\n")
    fSim.write("\n")

#== No return values ==================================================

    return

#========================================================================
#  Write Eclipse DENSITY keyword
#========================================================================

def outputECLDensity(OutU,fSim,clsSIM,clsUNI) :

    if   OutU == "MET" : uUnit = "METRIC"
    elif OutU == "SI"  : uUnit = "SI"
    else   :             uUnit = "FIELD"

#== Did user provide the brine salinity? ==============================
    
    if clsSIM.setBr :

        if clsSIM.bSalt < 1.0E-06 :
            
            sLab1 = "--  User did not provide a value for Brine Salinity\n"
            sLab2 = "--  User will probably have to change 2nd argument of DENSITY k/w and most of PVTW k/w\n"

            fSim.write("--\n")
            fSim.write(sLab1)
            fSim.write(sLab2)
            fSim.write("--\n")
            fSim.write("\n")

        else :

            bSalt = clsSIM.bSalt
            uSalt = clsSIM.saltU

            bSalt = clsUNI.I2X(bSalt,uSalt)

            sLab1 = "--  Brine Salinity of " + str(bSalt) + " " + uSalt + "\n"

            fSim.write("--\n")
            fSim.write(sLab1)
            fSim.write("--\n")

    else :

        sLab1 = "--  User provided no Pref or Salt-Content to define Brine Properties\n"
        sLab2 = "--  Standard Pressure and Pure Water have been assumed\n"
        sLab3 = "--  User will probably have to change 2nd argument of DENSITY k/w and whole of PVTW k/w\n"
        fSim.write("--\n")
        fSim.write(sLab1)
        fSim.write(sLab2)
        fSim.write(sLab3)
        fSim.write("--\n")
        fSim.write("\n")

#==  Density Keyword ==================================================

    dOil = clsSIM.dSTO
    dWat = clsSIM.dSTW
    dGas = clsSIM.dSTG

    if uUnit == "METRIC" :
        sUni = "    KG/M3  "
        sCon = "kg/m3"
    else :
        sUni = "   LB/FT3  "
        sCon = "lb/ft3"

    dWat = clsUNI.I2X(dWat,sCon)

    if dOil < 0.0 : sOil = "    1*     "
    else :
        dOil = clsUNI.I2X(dOil,sCon)
        sOil = " {:10.3f}".format(dOil)
    
    if dGas < 0.0 : sGas = "    1*     "
    else :
        dGas = clsUNI.I2X(dGas,sCon)
        sGas = " {:10.5f}".format(dGas)
    
    sWat = " {:10.3f}".format(dWat)

    sLab1 = "--     ST Oil       ST Water     ST Gas\n"
    sLab2 = "--  " + sUni + "  " + sUni + "  " + sUni + "\n"
    sLab3 = "    " + sOil + "  " + sWat + "  " + sGas + "  /\n"

    fSim.write("--\n")
    fSim.write("--  Stock Tank Densities [Oil, Water, Gas]\n")
    fSim.write("--\n")
    fSim.write("\n")
    fSim.write("DENSITY\n")
    fSim.write(sLab1)
    fSim.write(sLab2)
    fSim.write(sLab3)
    fSim.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Output ECL PVTW Keyword
#========================================================================

def outputECLPVTW(OutU,fSim,clsSIM,clsUNI) :

    if   OutU == "MET" : uUnit = "METRIC"
    elif OutU == "SI"  : uUnit = "SI"
    else :               uUnit = "FIELD"

    pRefW = clsSIM.pRefW  #-- Reference Pressure        [psia   or bara]
    bRefW = clsSIM.bRefW  #--           FVF             [rb/stb or rm3/sm3]
    uRefW = clsSIM.uRefW  #--           Viscosity       [cP]
    cRefW = clsSIM.cRefW  #--           Compressibility [1/psi  or 1/bar]
    vRefW = clsSIM.vRefW  #--           Viscosibility   [1/psi  or 1/bar]

    uUni = "     cP    "
    if uUnit == "METRIC" :
        pUni  = "    BARA   "
        bUni  = "   RB/STB  "
        cUni  = "   1/BAR   "
        pRefW = clsUNI.I2X(pRefW,"bara")
        cRefW = clsUNI.I2X(cRefW,"1/bar")
        vRefW = clsUNI.I2X(vRefW,"1/bar")
    else :
        pUni = "    PSIA   "
        bUni = "  RM3/SM3  "
        cUni = "   1/PSI   "       

    sPrf = " {:10.3f}".format(pRefW)
    sBrf = " {:10.5f}".format(bRefW)
    sUrf = " {:10.5f}".format(uRefW)
    sCrf = " {:10.3e}".format(cRefW)
    sVrf = " {:10.3e}".format(vRefW)

    sLab1 = "--      Pref         BwRef        CwRef       UwRef        VwRef\n"
    sLab2 = "--  " + pUni + "  " + bUni + "  " + cUni + "  " + uUni + "  " + cUni + "\n"
    sLab3 = "    " + sPrf + "  " + sBrf + "  " + sCrf + "  " + sUrf + "  " + sVrf + "  /\n"

    sLab3 = "    %s  %s  %s  %s  %s  /\n" %(sPrf,sBrf,sCrf,sUrf,sVrf)

    fSim.write("--\n")
    fSim.write("--  Water Phase Properties\n")
    fSim.write("--\n")
    fSim.write("\n")
    fSim.write("PVTW\n")
    fSim.write(sLab1)
    fSim.write(sLab2)
    fSim.write(sLab3)
    fSim.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Output VIP PVTW Keyword
#========================================================================

def outputVIPPVTW(OutU,fSim,clsSIM,clsUNI) :

    if   OutU == "MET" : uUnit = "METRIC"
    elif OutU == "SI"  : uUnit = "SI"
    else :               uUnit = "FIELD"

    pRefW = clsSIM.pRefW  #-- Reference Pressure        [psia]
    bRefW = clsSIM.bRefW  #--           FVF             [rb/stb]
    uRefW = clsSIM.uRefW  #--           Viscosity       [cP]
    cRefW = clsSIM.cRefW  #--           Compressibility [1/psi]
    vRefW = clsSIM.vRefW  #--           Viscosibility   [1/psi]

    dSTW  = clsSIM.dSTW   #-- Water Density             [lb/ft3]          

    uUni = "     cP    "
    dUni = "  gm/cm3   "
    if uUnit == "METRIC" :
        pUni  = "    KPA   "
        bUni  = "   RB/STB  "
        cUni  = "   1/KPA   "
        pRefW = clsUNI.I2X(pRefW,"kpa")
        cRefW = clsUNI.I2X(cRefW,"1/kpa")
        vRefW = clsUNI.I2X(vRefW,"1/kpa")
    else :
        pUni = "    PSIA   "
        bUni = "  RB/STB   "
        cUni = "   1/PSI   "

    dSTW = clsUNI.I2X(dSTW,"gm/cm3")

    vRefW = vRefW*uRefW     #-- VIP/Nexus quotes dMu/dP = Mu*(1/Mu)*dMu/dP

    sPrf = " {:10.3f}".format(pRefW)
    sBrf = " {:10.5f}".format(bRefW)
    sUrf = " {:10.5f}".format(uRefW)
    sCrf = " {:10.3e}".format(cRefW)
    sVrf = " {:10.3e}".format(vRefW)
    sDen = " {:10.5f}".format(dSTW )

    sLab1 = "C              Pref        Dref       BwRef      CwRef     UwRef        VwRef\n"
    sLab2 = "C          " + pUni + "  " + dUni + "  " + bUni + "  " + cUni + "  " + uUni + "\n"
    sLab3 = "    " + sPrf + "  " + sBrf + "  " + sCrf + "  " + sUrf + "  " + sVrf + "  /\n"

    sLab3 = "      1   %s %s %s %s %s %s\n" %(sPrf,sDen,sBrf,sCrf,sUrf,sVrf)

#-- Unit Conversions ------------------------------------------------

    if OutU == "MET" :   #-- Metric-Equivalent
        sComp = "1/kPa"
        sDens = "gm/cm3"
        sGfvf = "rm3/sm3"
        sGOR  = "sm3/sm3"
        sOfvf = "rm3/sm3"
        sPres = "kPa"
        sTemp = "degC"
        sVisc = "cP"
    else :              #-- Field-Units
        sComp = "1/psi"
        sDens = "gm/cm3"
        sGfvf = "rb/scf"
        sGOR  = "scf/stb"
        sOfvf = "rb/stb"
        sPres = "psia"
        sTemp = "degF"
        sVisc = "cP"

#-- Headers ---------------------------------------------------------    

    fSim.write("C \n")
    fSim.write("C -----------------------------------------\n")
    fSim.write("C   Oil, Water and Gas Phase Properties\n")
    fSim.write("C \n")
    fSim.write("C   Units:  Compressibility " + sComp + "\n")
    fSim.write("C           Density         " + sDens + "\n")
    fSim.write("C           Gas FVF         " + sGfvf + "\n")
    fSim.write("C           GOR             " + sGOR  + "\n")
    fSim.write("C           Oil FVF         " + sOfvf + "\n")
    fSim.write("C           Pressure        " + sPres + "\n")
    fSim.write("C           Temperature     " + sTemp + "\n")
    fSim.write("C           Viscosity       " + sVisc + "\n")
    fSim.write("C -----------------------------------------\n")
    fSim.write("C \n")
    fSim.write("\n")

#== Water Phase Properties ============================================    

    fSim.write("C \n")
    fSim.write("C   Water Phase Properties\n")
    fSim.write("C \n")
    fSim.write("\n")
    fSim.write("PVTW  IPVTW   PBASEW        DWB        BWI       CW            VW        VWP\n")
    fSim.write(sLab3)
    fSim.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Write Regression Variables
#========================================================================

def outputRegVars(clsIO,clsEOS,dicREG,dicSAM,xVec,delX,nIter) :

    fReg = clsIO.fReg
    nVar = len(dicREG)

    sHead = "---------------------------------------\n"
    sLine = "  Regression Iteration (Variables): " + str(nIter) + "\n"

    fReg.write(sHead)
    fReg.write(sLine)
    fReg.write(sHead)
    fReg.write("\n")

#            12IV12vType12  xMin  12  xVec  12  xMax  12  delX  12  Bounded?  12Comps/Samps
    sLine = "  IV  vType    xMin      xVec      xMax      delX      Bounded?    Comps/Samp\n"
    sHead = "  --  -----  --------  --------  --------  --------  ------------  ----------\n"
    fReg.write(sHead)
    fReg.write(sLine)
    fReg.write(sHead)

    for iVar in range(nVar) :
        
        clsREG = dicREG[iVar]

        sNam = clsREG.vNam.ljust(5)
        sTyp = clsREG.vTyp

        xMin = clsREG.vMin
        xMax = clsREG.vMax

        sBou = "            "
                                            #         123456789012
        if   xVec[iVar] + delX[iVar] > xMax : sBou = "xMax Will Be"
        elif xVec[iVar] + delX[iVar] < xMin : sBou = "xMin Will Be"

        if   xVec[iVar] == xMin             : sBou = "xMin Already"
        elif xVec[iVar] == xMax             : sBou = "xMax Already"

        sVar =   "{:2d}".format(iVar+1)
        sMin = "{:8.5f}".format(xMin)
        sMax = "{:8.5f}".format(xMax)
        sVec = "{:8.5f}".format(xVec[iVar])
        sDlx = "{:8.5f}".format(delX[iVar])

        sExt = CR.regVarInfo(sTyp,clsREG,clsEOS,dicSAM)

        sLine = "  " + sVar + "  " + sNam + "  " + sMin + "  " + sVec + \
                "  " + sMax + "  " + sDlx + "  " + sBou + "  " + sExt + "\n"
        fReg.write(sLine)

    fReg.write(sHead)

    return

#========================================================================
#  Output SAV format
#========================================================================

def outputSave(sTit,clsEOS,dicSAM,clsIO) :

#-- Initialisation --------------------------------------------------    

    fSav = clsIO.fSav
    nCom = clsEOS.NC
    nSam = len(dicSAM)
    EOS  = clsEOS.EOS

#-- Define ----------------------------------------------------------    

    sLine = "DEFINE" + "\n"
#            123456789012345678901234567890123456789012345678901234567890 
    sHead = "------------------------------------------------------------\n"
    fSav.write("\n")
    fSav.write(sHead)
    fSav.write(sLine)
    fSav.write(sHead)

#-- Title -----------------------------------------------------------    

    sLine = "TITLE  " + sTit + "\n"
    fSav.write("\n")
    fSav.write(sLine)
    fSav.write("\n")

#-- EOS -------------------------------------------------------------    

    sLine = "EOS    " + EOS + "\n"
    fSav.write(sLine)
    fSav.write("\n")

#-- Number of Components and Samples --------------------------------    
    
    sLine = "NCOMP  " + str(nCom) + "\n"
    fSav.write(sLine)
    fSav.write("\n")

    sLine = "NSAMP  " + str(nSam) + "\n"
    fSav.write(sLine)
    fSav.write("\n")

#== Properties keyword ================================================    

    sLine = "PROPS\n"
    fSav.write(sLine)

#            12CompName12  MOLWT   12  TCRIT   12  PCRIT   12  VCRIT   
    sHed1 = "  COMP        MOLWT       TCRIT       PCRIT       VCRIT   " + \
                      "    ZCRIT       ACENF       TBOIL       SPECG   " + \
                      "    PARAC       VOLSH       MULTA       MULTB   " + \
                      "    IDCPA       IDCPB       IDCPC       IDCPD   " + "\n"
    sHed2 = "  NAME       lb/lbmol      degR        psia     ft3/lbmol " + \
                      "     dim         dim         degR      wrtWat=1 " + \
                      "     dim         dim         dim         dim    " + \
                      "  btu/lbmoR1  btu/lbmoR2  btu/lbmoR3  btu/lbmoR4" + "\n"
    fSav.write(sHed1)
    fSav.write(sHed2)

#----------------------------------------------------------------------
#  Component Properties
#----------------------------------------------------------------------

    for iC in range(nCom) :

        sNm = clsEOS.gPP("CN",iC).ljust(8)

        Mw = clsEOS.gPP("MW",iC) ; Tc = clsEOS.gPP("TC",iC)
        Pc = clsEOS.gPP("PC",iC) ; Vc = clsEOS.gPP("VC",iC)
        Zc = clsEOS.gPP("ZC",iC) ; AF = clsEOS.gPP("AF",iC)
        Tb = clsEOS.gPP("TB",iC) ; SG = clsEOS.gPP("SG",iC)
        Pa = clsEOS.gPP("PA",iC) ; SS = clsEOS.gPP("SS",iC)
        MA = clsEOS.gPP("MA",iC) ; MB = clsEOS.gPP("MB",iC)
        CA = clsEOS.gPP("CA",iC) ; CB = clsEOS.gPP("CB",iC)
        CC = clsEOS.gPP("CC",iC) ; CD = clsEOS.gPP("CD",iC)

        sMw = "  {:10.4f}".format(Mw) ; sTc = "  {:10.4f}".format(Tc)
        sPc = "  {:10.4f}".format(Pc) ; sVc = "  {:10.5f}".format(Vc)
        sZc = "  {:10.7f}".format(Zc) ; sAF = "  {:10.7f}".format(AF)
        sTb = "  {:10.4f}".format(Tb) ; sSG = "  {:10.7f}".format(SG)
        sPa = "  {:10.4f}".format(Pa) ; sSS = "  {:10.7f}".format(SS)
        sMA = "  {:10.7f}".format(MA) ; sMB = "  {:10.7f}".format(MB)
        sCA = "  {:10.3e}".format(CA) ; sCB = "  {:10.3e}".format(CB)
        sCC = "  {:10.3e}".format(CC) ; sCD = "  {:10.3e}".format(CD)

        sLin = "  " + sNm + sMw + sTc + sPc + sVc + sZc + \
                      sAF + sTb + sSG + sPa + sSS + \
                      sMA + sMB + sCA + sCB + sCC + sCD + "\n"

        fSav.write(sLin)

#----------------------------------------------------------------------
#  Binary Iteraction Parameters
#----------------------------------------------------------------------

    fSav.write("\n")
    sLine = "BIP\n"
    fSav.write(sLine)

    sLin = "  COMP    "
    for jC in range(nCom) :
        sNam = clsEOS.gPP("CN",jC)
        sLin = sLin + "  " + sNam.center(8)
    sLin = sLin + "\n"

    fSav.write(sLin)

    for iC in range(nCom) :
        sLin = "  "
        sNam = clsEOS.gPP("CN",iC)
        sLin = sLin + sNam.ljust(8)
        for jC in range(nCom) :
            KIJ  = clsEOS.gIJ(iC,jC)
            sIJ  = "  {:8.5f}".format(KIJ)
            sLin = sLin + sIJ
        sLin = sLin + "\n"
        fSav.write(sLin)

    fSav.write("\n")

#----------------------------------------------------------------------
#  Compositions
#----------------------------------------------------------------------

    fSav.write("\n")
    sLine = "SAMPLES\n"
    fSav.write(sLine)

    sLin = "  NAME  "
    for iSam in range(nSam) :
        sNam = dicSAM[iSam].sName
        sLin = sLin + sNam.rjust(10)
    sLin = sLin + "\n"

    fSav.write(sLin)

    for iC in range(nCom) :
        sLin = "  "
        sNam = clsEOS.gPP("CN",iC)
        sLin = sLin + sNam.ljust(8)
        for iSam in range(nSam) :
            dCom = dicSAM[iSam].gZI(iC)
            sCom = "  {:8.6f}".format(dCom)
            sLin = sLin + sCom
        sLin = sLin + "\n"
        fSav.write(sLin)

#-- EndDef ----------------------------------------------------------    

    sLine = "ENDDEF" + "\n"
    fSav.write("\n")
    fSav.write(sHead)
    fSav.write(sLine)
    fSav.write(sHead)
    fSav.write("\n")

    return

#========================================================================
#  Write an array to the Debug File
#========================================================================

def writeArrayDebug(fDeb,Array,ArrNam) :

    maxE = 6
    nArr = len(Array)
    nLin = int(nArr/maxE)

    if nArr % maxE > 0 : nLin += 1

    sOut = ArrNam + "\n"
    fDeb.write(sOut)

    nElm = 0

    for iLin in range(nLin) :
        sOut = ""
        for iElm in range(maxE) :
            if nElm > nArr-1 : break
            sThs = " {:10.3e}".format(Array[nElm])
            sOut = sOut + sThs
            nElm += 1
        sOut = sOut + "\n"
        fDeb.write(sOut)

#== No return values ==================================================

    return

#========================================================================
#  End of Module
#========================================================================
