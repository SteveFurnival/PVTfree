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

import plotBatch as PB
import readExps  as RX
import writeReg  as WR
import writeTT   as WT

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
    fOut.write(clsIO.sDocStr+"\n")
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

def outputProps(sTit,clsEOS,dicSAM,clsIO) :

    fOut = clsIO.fOut

    nComp = clsEOS.nComp
    nSamp = clsEOS.nSamp

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
        sNam = dicSAM[iSamp].sNam
        sNam = "{:>8}  ".format(sNam)
        sOut = sOut + sNam
    sOut = sOut + "\n"
    fOut.write(sOut)
    fOut.write("\n")

    sOut = "Mol Wt    "
    for iSamp in range(nSamp) :
        Mw  = dicSAM[iSamp].iPlsMW
        sMw = "{:8.3f}  ".format(Mw)
        sOut  = sOut + sMw + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

    sOut = "Spec Grv  "
    for iSamp in range(nSamp) :
        SG  = dicSAM[iSamp].iPlsSG
        sSG = "{:8.5f}  ".format(SG)
        sOut  = sOut + sSG + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

    sOut = "Alpha     "
    for iSamp in range(nSamp) :
        alf = dicSAM[iSamp].iPlsAL
        sAl = "{:8.5f}  ".format(alf)
        sOut  = sOut + sAl + "  "
    sOut = sOut + "\n"
    fOut.write(sOut)

#-- First Set of Properties -----------------------------------------    

    tab = WT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","f","f","f","f","f","f"]
    tabAlgn = ["l","r","r","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ,11 ,11 ]

    tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
    
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

        cName = clsEOS.gNM(iC)
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
 
    tab = WT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","f","f","f","f","f","f"]
    tabAlgn = ["l","r","r","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ,11 ,11 ]

    tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
    
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

        cName = clsEOS.gNM(iC)
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
 
    tab = WT.Texttable()
    tab.set_deco(tab.HEADER)
    
    tabType = ["t","e","e","e","e"]
    tabAlgn = ["l","r","r","r","r"]
    tabColW = [11 ,11 ,11 ,11 ,11 ]

    tab = setTabProps(tab,tabType,tabAlgn,tabColW,3)
    
    tOut = [[0 for jC in range(5)] for iC in range(nComp+1)]

#               12345678901--12345678901
    tabHead = ["Comp       \nName       ",
               "   IdCpA   \nBtu/lbmo.R ",
               "   IdCpB   \nBtu/lbmo.R2",
               "   IdCpC   \nBtu/lbmo.R3",
               "   IdCpD   \nBtu/lbmo.R4"]

    tOut[0] = tabHead
    
    for iC in range(nComp) :

        cName = clsEOS.gNM(iC)
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
        
        tab = WT.Texttable()
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

            cName = clsEOS.gNM(iHead)
            iHead += 1
            
            tabType.append("f")
            tabAlgn.append("r")
            tabColW.append( 8 )
            tabHead.append(cName)

        tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
    
        tOut[0] = tabHead

        for iR in range(nComp) :

            cName = clsEOS.gNM(iR)

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
        
        tab = WT.Texttable()
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

            sName = dicSAM[iHead].sNam
            iHead += 1
            
            tabType.append("f")
            tabAlgn.append("r")
            tabColW.append(10 )
            tabHead.append(sName)

        tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
    
        tOut[0] = tabHead

        for iR in range(nComp) :

            cName = clsEOS.gNM(iR)

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
#  Set 'tab' properties
#========================================================================

def setTabProps(tab,tabType,tabAlgn,tabColW,nPrec) :

    tab.set_cols_dtype(tabType)
    tab.set_cols_align(tabAlgn)
    tab.set_cols_width(tabColW)
    tab.set_precision(nPrec)

#== Return updated tab ================================================

    return tab

#========================================================================
#  Output the Experiments Results
#========================================================================

def outputExps(clsEOS,dicSAM,dicEXP,clsUNI,clsIO) :

    nExp = len(dicEXP)

    fOut = clsIO.fOut

#== Loop over experiments =============================================

    for iExp in range(nExp) :

        clsEXP = dicEXP[iExp]

        if clsEXP.IsAct :

            xTyp = clsEXP.xName

            #print("Experiment: iExp,Type ",iExp+1,xTyp)

            nSam = clsEXP.nSamp

            clsSAM = dicSAM[nSam]
            
            sNam = clsSAM.sNam
            nCom = clsSAM.nCom

            if   xTyp == "CCE" : sLabl = "  Constant Composition Expansion (CCE) Experiment\n"
            elif xTyp == "CVD" : sLabl = "  Constant Volume Depletion (CVD) Experiment\n"
            elif xTyp == "DLE" : sLabl = "  Differential Liberation (DLE) Experiment\n"
            elif xTyp == "SEP" : sLabl = "  Separator Test (SEP) Experiment\n"
            elif xTyp == "FLS" : sLabl = "  2-Phase Flash (FLS) Experiment\n"
            elif xTyp == "SAT" : sLabl = "  Saturation Pressure (SAT) Experiment\n"
            elif xTyp == "SWL" : sLabl = "  Swelling Test (SWL) Experiment\n"
            elif xTyp == "GRD" : sLabl = "  Composition Versus Depth (GRD) Experiment\n"

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

            nInd = clsEXP.nInd
            nObs = clsEXP.nObs
            nObU = clsEXP.nObU
            nCol = nInd + nObU + nObs + nSep
            nRow = clsEXP.nRow

            tab = WT.Texttable()
            tab.set_deco(tab.HEADER)

            tabType = ["f" for i in range(nCol)]
            tabAlgn = ["r" for i in range(nCol)]
            tabColW = [ 11 for i in range(nCol)]

            tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
        
            tOut = [[0.0 for j in range(nCol)] for i in range(nRow+1)]

            tabHead = []

            for iInd in range(nInd) :
                sInd = clsEXP.hInd[iInd] + "\n" + clsEXP.uInd[iInd]
                tabHead.append(sInd)

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

                for iInd in range(nInd) :
                    dInd = clsUNI.I2X(clsEXP.dInd[iInd][iRow],clsEXP.uInd[iInd])
                    tRow.append(dInd)

                for iSep in range(nSep) :
                    sSep = str(clsEXP.Lsep[iRow]) + "     " + \
                           str(clsEXP.Vsep[iRow])
                    tRow.append(sSep)

                for iObs in range(nObs) :
                    if clsEXP.qObs[iObs] :
                        dObs = clsUNI.I2X(clsEXP.dObs[iObs][iRow],clsEXP.uObs[iObs])
                        tRow.append(dObs)
                    if clsEXP.qCal[iObs] :
                        dCal = clsUNI.I2X(clsEXP.dCal[iObs][iRow],clsEXP.uCal[iObs])
                        tRow.append(dCal)

                tOut[iRow+1] = tRow
        
            tab.add_rows(tOut)

            fOut.write(tab.draw())
            fOut.write("\n")

#-- Extra comment for SEP if non-standard routing -------------------
            
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

#== Output Compositions for GRD experiment ============================

            if xTyp == "GRD" :

                fOut.write("\n")

                nCol = nRow + 1
                nRow = nCom

                tab = WT.Texttable()
                tab.set_deco(tab.HEADER)

                tabType = ["f" for i in range(nCol)]
                tabAlgn = ["r" for i in range(nCol)] ; tabAlgn[0] = "l"
                tabColW = [ 11 for i in range(nCol)]

                tab = setTabProps(tab,tabType,tabAlgn,tabColW,5)
        
                tOut = [[0.0 for j in range(nCol)] for i in range(nRow+1)]

                tabHead = []
                tabHead.append("Name")

                for iH in range(nCol-1) :
                    sHgt ="{:8.2f}".format(clsEXP.dInd[0][iH])
                    sCal = "H/[" + clsEXP.uInd[0] + "]\n" + sHgt
                    tabHead.append(sCal)

                tOut[0] = tabHead

#-- Loop over rows of Experiment ------------------------------------

                for iC in range(nRow) :

                    tRow = []
                    tRow.append(clsEOS.gNM(iC))

                    for iH in range(nCol-1) :
                        dCal = clsUNI.I2X(clsEXP.zCal[iH][iC],"perc")
                        tRow.append(dCal)

                    tOut[iC+1] = tRow

                tab.add_rows(tOut)

                fOut.write(tab.draw())
                fOut.write("\n")
        
#======================================================================
#  Write Plots
#======================================================================

    if clsIO.Pall : PB.genPlots(clsIO,dicEXP,dicSAM,clsUNI)

    return

#========================================================================
#  Output SAV format
#========================================================================

def outputSave(sTit,clsEOS,dicSAM,clsIO) :

    EOL = "\n"
#            123456789012345678901234567890123456789012345678901234567890 
    sHead = "------------------------------------------------------------\n"

#-- Initialisation --------------------------------------------------    

    fSav = clsIO.fSav
    nCom = clsEOS.nComp
    nSam = len(dicSAM)
    EOS  = clsEOS.EOS

    clsIO.nSav += 1    #-- Number of Times Sav has been written

#-- Define ----------------------------------------------------------    

    sLine  = "DEFINE  " + str(clsIO.nSav) + "  \n"
    setLin = [EOL,sHead,sLine,sHead]
    printLines(fSav,setLin)
    
#-- Title -----------------------------------------------------------    

    sLine  = "TITLE  " + sTit + EOL
    setLin = [EOL,sLine,EOL]
    printLines(fSav,setLin)
    
#-- EOS -------------------------------------------------------------    

    sLine  = "EOS    " + EOS + EOL
    setLin = [sLine,EOL]
    printLines(fSav,setLin)

#-- Number of Components and Samples --------------------------------    
    
    sLine  = "NCOMP  " + str(nCom) + EOL
    setLin = [sLine,EOL]
    printLines(fSav,setLin)

    sLine  = "NSAMP  " + str(nSam) + EOL
    setLin = [sLine,EOL]
    printLines(fSav,setLin)

#== Properties keyword ================================================    

    sLine = "PROPS\n"
#            12CompName12  MOLWT   12  TCRIT   12  PCRIT   12  VCRIT   
    sHed1 = "  COMP        MOLWT       TCRIT       PCRIT       VCRIT   " + \
                      "    ZCRIT       ACENF       TBOIL       SPECG   " + \
                      "    PARAC       VOLSH       MULTA       MULTB   " + \
                      "    IDCPA       IDCPB       IDCPC       IDCPD   " + EOL
    sHed2 = "  NAME       lb/lbmol      degR        psia     ft3/lbmol " + \
                      "     dim         dim         degR      wrtWat=1 " + \
                      "     dim         dim         dim         dim    " + \
                      "  btu/lbmoR1  btu/lbmoR2  btu/lbmoR3  btu/lbmoR4" + EOL

    setLin = [sLine,sHed1,sHed2]
    printLines(fSav,setLin)

#----------------------------------------------------------------------
#  Component Properties
#----------------------------------------------------------------------

    for iC in range(nCom) :

        sNm = clsEOS.gNM(iC).ljust(8)

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
                      sMA + sMB + sCA + sCB + sCC + sCD + EOL

        fSav.write(sLin)

#----------------------------------------------------------------------
#  Binary Iteraction Parameters
#----------------------------------------------------------------------

    fSav.write(EOL)
    sLine = "BIP\n"
    fSav.write(sLine)

    sLin = "  COMP    "
    for jC in range(nCom) :
        sNam = clsEOS.gNM(jC)
        sLin = sLin + "  " + sNam.center(8)
    sLin = sLin + EOL

    fSav.write(sLin)

    for iC in range(nCom) :
        sLin = "  "
        sNam = clsEOS.gNM(iC)
        sLin = sLin + sNam.ljust(8)
        for jC in range(nCom) :
            KIJ  = clsEOS.gIJ(iC,jC)
            sIJ  = "  {:8.5f}".format(KIJ)
            sLin = sLin + sIJ
        sLin = sLin + EOL
        fSav.write(sLin)

    fSav.write(EOL)

#----------------------------------------------------------------------
#  Compositions
#----------------------------------------------------------------------

    fSav.write(EOL)
    sLine = "SAMPLES\n"
    fSav.write(sLine)

    sLin = "  NAME  "
    for iSam in range(nSam) :
        sNam = dicSAM[iSam].sNam
        sLin = sLin + sNam.rjust(10)
    sLin = sLin + EOL

    fSav.write(sLin)

    for iC in range(nCom) :
        sLin = "  "
        sNam = clsEOS.gNM(iC)
        sLin = sLin + sNam.ljust(8)
        for iSam in range(nSam) :
            dCom = dicSAM[iSam].gZI(iC)
            sCom = "  {:8.6f}".format(dCom)
            sLin = sLin + sCom
        sLin = sLin + EOL
        fSav.write(sLin)

#-- EndDef ----------------------------------------------------------    

    sLine = "ENDDEF" + EOL
    setLin = [EOL,sHead,sLine,sHead,EOL]
    printLines(fSav,setLin)

    return

#========================================================================
#  Print a set of lines to a specified file
#========================================================================

def printLines(fOut,setLines) :

    for line in setLines : fOut.write(line)

    return

#========================================================================
#  Write an array to the Debug File
#========================================================================

def writeArrayDebug(fDeb,Array,ArrNam) :

    EOL = "\n"

    maxE = 6
    nArr = len(Array)
    nLin = int(nArr/maxE)

    if nArr % maxE > 0 : nLin += 1

    sOut = ArrNam + EOL
    fDeb.write(sOut)

    nElm = 0

    for iLin in range(nLin) :
        sOut = ""
        for iElm in range(maxE) :
            if nElm > nArr-1 : break
            sThs = " {:10.3e}".format(Array[nElm])
            sOut = sOut + sThs
            nElm += 1
        sOut = sOut + EOL
        fDeb.write(sOut)

#== No return values ==================================================

    return

#========================================================================
#  End of Module
#========================================================================
