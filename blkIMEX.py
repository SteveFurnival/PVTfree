
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Generates Blackoil Tables for IMEX
#
#=======================================================================

#!/usr/bin/python3

import blkOther  as BO
import constants as CO

def outIMEX(fIMX,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp

#== One Table plus Supplementary Keywords =============================

    outputHeader(fIMX,clsBLK,clsIO)

    outTable(fIMX,dTab,eTab,clsBLK,clsUNI,clsIO)
    
#======================================================================
#  PVT & PVTVO - Output BOT & VOT
#  PVTCOND     - Output BGUST & VGUST
#======================================================================

    if   oTyp == "PVTCOND" :
        genUnderSatGas(fIMX,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)
    elif oTyp == "PVT" or oTyp == "PVTVO" :
        genUnderSatOil(fIMX,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)
        
#== No return value ===================================================

    return

#========================================================================
#  Main Table
#========================================================================

def outTable(fIMX,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        if iSat == 0 : sExt = ""
        else         : sExt = ""

        outputRow(fIMX,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    sExt = ""

    while iExt >= 0 :

        outputRow(fIMX,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

#== Terminating Slash =================================================

    fIMX.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Write Saturated (up to 6-values) line of data
#========================================================================

def outputRow(fSim,iRow,fTab,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    Pr = fTab[iRow][clsBLK.iPr]
    Rs = fTab[iRow][clsBLK.iRs]
    Rv = fTab[iRow][clsBLK.iRv]
    Bo = fTab[iRow][clsBLK.iBo]
    Bg = fTab[iRow][clsBLK.iBg]
    Uo = fTab[iRow][clsBLK.iUo]
    Ug = fTab[iRow][clsBLK.iUg]
    Bd = fTab[iRow][clsBLK.iBd]
    Ud = fTab[iRow][clsBLK.iUd]

    if OutU == "MET" :
        Pr = clsUNI.I2X(Pr,"kpa")
        Uo = clsUNI.I2X(Uo,"mpa.s")
        Ug = clsUNI.I2X(Ug,"mpa.s")
        Ud = clsUNI.I2X(Ud,"mpa.s")
    else             :
        Rs = clsUNI.I2X(Rs,"scf/stb")
        Rv = clsUNI.I2X(Rv,"stb/scf")
        Bg = clsUNI.I2X(Bg,"rb/scf")
        Bd = clsUNI.I2X(Bd,"rb/scf")

    sPr = " {:10.3f}".format(Pr)
    sRs = " {:10.3f}".format(Rs)
    sRv = " {:10.3e}".format(Rv)
    sBo = " {:10.5f}".format(Bo)
    sBg = " {:10.3e}".format(Bg)
    sUo = " {:10.5f}".format(Uo)
    sUg = " {:10.5f}".format(Ug)
    sBd = " {:10.3e}".format(Bd)
    sUd = " {:10.3e}".format(Ud)

    if   oTyp == "PVTG" :
        sLabl = "    " + sPr + "  " + sBd + "  " + sUd + "\n"
    elif oTyp == "PVT" :
        sLabl = "    " + sPr + "  " + sRs + "  " + sBo + "  " + sBg + "  " + sUo + "  " + sUg + "\n"
    elif oTyp == "PVTVO" :
        sLabl = "    " + sPr + "  " + sRs + "  " + sRv + "  " + sBo + "  " + sBg + "  " + sUo + "  " + sUg + "  " + sBd + "  " + sUd + "\n"
    elif oTyp == "PVTCOND" :
        sLabl = "    " + sPr + "  " + sRs + "  " + sRv + "  " + sBo + "  " + sBg + "  " + sUo + "  " + sUg + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  Table Header
#========================================================================

def outputHeader(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    uUni = "     cP    "
    if OutU == "MET" :
        rUni = " KSM3/SM3  "
        pUni = "    BARA   "
        bUni = "  RM3/SM3  "
        cUni = "  1/BARA   "
    else :
        rUni = "  MSCF/STB "
        pUni = "    PSIA   "
        bUni = "   RB/STB  "
        cUni = "  1/PSIA   "

    if   oTyp == "PVTG" :
        sLab1 = "**      PGas       Gas FVF      Gas Visc\n"
    elif oTyp == "PVT"  :
        sLab1 = "**      PBub         GOR         Oil FVF      Gas FVF     Oil Visc     Gas Visc  \n"
    elif oTyp == "PVTVO" :
        sLab1 = "**      PBub         GOR          CGR         Oil FVF     Gas FVF      Oil Visc     Gas Visc      Dry FVF     Dry Visc\n"
    elif oTyp == "PVTCOND" :
        sLab1 = "**      PDew         GOR          CGR         Oil FVF     Gas FVF      Oil Visc     Gas Visc  \n"

    fSim.write("**\n")
    fSim.write("**  Fluid Properties\n")
    fSim.write("**\n")
    fSim.write("\n")
    fSim.write("*" + oTyp + "  *BG  1\n")
    fSim.write("\n")
    fSim.write(sLab1)
    fSim.write("\n")

    return

#========================================================================
#  Generate Undersatuarted Oil Data (BOT & VOT)
#========================================================================

def genUnderSatOil(fIMX,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    if OutU == "MET" :
        sPrs = "kpa"
        sFVF = "rm3/sm3"
    else :
        sPrs = "psia"
        sFVF = "rb/stb"
    sVis = "cp"

    fPrs = " {:10.3f}"
    fFVF = " {:10.5f}"
    fVis = " {:10.5f}"

    iP = 0 ; iB = 1 ; iU = 2

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    qLiq = True
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG
    pInc = clsBLK.Pinc

#== "Saturated" Data ==================================================    

    while iSat >= 0 :

        dOut = []

        Rs = dTab[iSat][clsBLK.iRs]

        jSat = iSat

        while jSat >= 0 :

            BO.setEoSViS(jSat,sOil,sGas,rOil,rGas,clsBLK)

            Pr  = dTab[jSat][clsBLK.iPr]
            RTp = clsBLK.RT/Pr

            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs,clsBLK,clsIO)

            dRow = [Pr,Bo,Uo]
            dOut.append(dRow)

#-- Decrement jSat counter ------------------------------------------

            jSat -= 1

#-- Extended Pressure Stages ----------------------------------------            

        BO.setEoSViS(0,sOil,sGas,rOil,rGas,clsBLK)

        jExt = nExt -1

        while jExt >= 0 :

            Pr  = eTab[jExt][clsBLK.iPr]
            RTp = clsBLK.RT/Pr

            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs,clsBLK,clsIO)

            dRow = [Pr,Bo,Uo]
            dOut.append(dRow)

#-- Decrement jExt counter ------------------------------------------

            jExt -= 1            

#== Write BOT/VOT for this sweep ======================================

        outUnderSat(fIMX,"BOT",dOut,iP,iB,sPrs,sFVF,fPrs,fFVF,clsUNI)

        outUnderSat(fIMX,"VOT",dOut,iP,iU,sPrs,sVis,fPrs,fVis,clsUNI)

#-- Decrement iSat counter ------------------------------------------

        iSat -= 1

#== "Extended" Data ===================================================

    BO.setEoSViS(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        dOut = []

        jExt = iExt

        while jExt >= 0 :

            Pr  = eTab[jExt][clsBLK.iPr]
            RTp = clsBLK.RT/Pr

            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs,clsBLK,clsIO)

            dRow = [Pr,Bo,Uo]
            dOut.append(dRow)

#-- Decrement jExt counter ------------------------------------------

            jExt -= 1

#-- Additional Row if last Row in the Extended Table ----------------

        if iExt == 0 :

            Pr  = Pr + pInc
            RTp = clsBLK.RT/Pr

            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rs,clsBLK,clsIO)

            dRow = [Pr,Bo,Uo]
            dOut.append(dRow)

#== Write BOT/VOT for this sweep ======================================

        outUnderSat(fIMX,"BOT",dOut,iP,iB,sPrs,sFVF,fPrs,fFVF,clsUNI)

        outUnderSat(fIMX,"VOT",dOut,iP,iU,sPrs,sVis,fPrs,fVis,clsUNI)

#-- Decrement iSat counter ------------------------------------------

        iExt -= 1

#== No return value ===================================================

    return

#========================================================================
#  Generate Undersatuarted Gas Data (BGUST & VGUST)
#========================================================================

def genUnderSatGas(fIMX,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    if OutU == "MET" :
        sPrs = "kpa"
        sFVF = "rm3/sm3"
    else :
        sPrs = "psia"
        sFVF = "rb/scf"
    sVis = "cp"

    fPrs = " {:10.3f}"
    fFVF = " {:10.3e}"
    fVis = " {:10.3e}"

    iP = 0 ; iB = 1 ; iU = 2

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    Rv0  = 0.0
    P0   = CO.pStand
    
    qVap = False
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

    while iSat >= 0 :

#-- Set EoS/Visc coefficients based on pressure ---------------------        

        BO.setEoSViS(iSat,sOil,sGas,rOil,rGas,clsBLK)

        Pr  = dTab[iSat][clsBLK.iPr]
        RTp = clsBLK.RT/Pr

        dOut = []

#-- Dry Gas Stage (Rv = 0.0) ----------------------------------------        

        Bg0,Ug0 = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv0,clsBLK,clsIO)

        dRow = [P0,Bg0,Ug0]
        dOut.append(dRow)

#== "Saturated Data" ==================================================

        jSat = nSat - 1

        while jSat >= iSat :

            Pr = dTab[jSat][clsBLK.iPr]
            Rv = dTab[jSat][clsBLK.iRv]

            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)

            dRow = [Pr,Bg,Ug]
            dOut.append(dRow)

#-- Decrement jSat counter ------------------------------------------

            jSat -= 1

#== Write BGUST/ZGUST for this sweep ==================================

        outUnderSat(fIMX,"BGUST",dOut,iP,iB,sPrs,sFVF,fPrs,fFVF,clsUNI)

        outUnderSat(fIMX,"VGUST",dOut,iP,iU,sPrs,sVis,fPrs,fVis,clsUNI)

#== Decrement iSat counter ============================================

        iSat -= 1        

#----------------------------------------------------------------------
#  "Extended" Data
#----------------------------------------------------------------------

    BO.setEoSViS(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        Pr  = eTab[iExt][clsBLK.iPr]
        RTp = clsBLK.RT/Pr

        dOut = []

#-- Dry Gas Stage (Rv = 0.0) ----------------------------------------        

        Bg0,Ug0 = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv0,clsBLK,clsIO)

        dRow = [P0,Bg0,Ug0]
        dOut.append(dRow)

#== "Saturated Data" ==================================================

        jSat = nSat - 1

        while jSat >= 0 :

            Pr = dTab[jSat][clsBLK.iPr]
            Rv = dTab[jSat][clsBLK.iRv]

            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)

            dRow = [Pr,Bg,Ug]
            dOut.append(dRow)

#-- Decrement jSat counter ------------------------------------------

            jSat -= 1

#== "Extended Data" ===================================================

        jExt = nExt - 1

        while jExt >= iExt :

            Pr = eTab[jExt][clsBLK.iPr]
            Rv = eTab[jExt][clsBLK.iRv]

            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)

            dRow = [Pr,Bg,Ug]
            dOut.append(dRow)

#-- Decrement jExt counter ------------------------------------------

            jExt -= 1

#== Write BGUST/ZGUST for this sweep ==================================

        outUnderSat(fIMX,"BGUST",dOut,iP,iB,sPrs,sFVF,fPrs,fFVF,clsUNI)

        outUnderSat(fIMX,"VGUST",dOut,iP,iU,sPrs,sVis,fPrs,fVis,clsUNI)

#-- Decrement iExt counter ------------------------------------------

        iExt -= 1

#== No return value ===================================================

    return

#========================================================================
#  Write a block of BOT/VOT or BGUST/VGUST data
#========================================================================

def outUnderSat(fIMX,sTit,dOut,i1,i2,s1,s2,f1,f2,clsUNI) :

    fIMX.write("*" + sTit + "  1  \n")
    fIMX.write("\n")

    nOut = len(dOut)

    for iOut in range(nOut) :
        
        X1  = dOut[iOut][i1]
        X2  = dOut[iOut][i2]
        
        X1  = clsUNI.I2X(X1,s1)
        X2  = clsUNI.I2X(X2,s2)
        
        sX1 = f1.format(X1)
        sX2 = f2.format(X2)

        sLab = "    " + sX1 + "  " + sX2 + "\n"

        fIMX.write(sLab)

    fIMX.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  End of Module
#========================================================================


