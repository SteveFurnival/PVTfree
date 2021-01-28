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

#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Blackoil Table Output Module (WB)
#
#========================================================================

#!/usr/bin/python3

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

    sim3 = tSim[:3]
    
    if   sim3 == "CMG" :
        sComm = "**"
        sHead = "**  (CMG) IMEX Blackoil Tables\n"
        sLine = "**============================================================\n"
        sPace = "**\n"
    elif sim3 == "TEM" :
        sComm = "--"
        sHead = "--  (Tempest) MORE BlackOil Tables\n"
        sLine = "--============================================================\n"
        sPace = "--\n"
    elif sim3 == "ECL" :
        sComm = "--"
        sHead = "--  ECLIPSE 100 (E100) BlackOil Tables\n"
        sLine = "--============================================================\n"
        sPace = "--\n"
    elif sim3 == "VIP" :
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

    if   OutU[:3] == "MET" : uUnit = "METRIC"
    elif OutU[:2] == "SI"  : uUnit = "SI"
    else                   : uUnit = "FIELD"

    sLabl = sComm + " Units set for the Tables is " + uUnit +"\n"

    fSim.write(sPace)
    fSim.write(sLabl)
    fSim.write(sPace)
    fSim.write(sNull)

#==  Density Keyword ==================================================
    
    if tSim[:3] == "ECL" or tSim[:3] == "TEM" :
        outputECLDensity(OutU,fSim,clsBLK,clsUNI)

#== Water Properties (PVTW) ===========================================

    if clsBLK.setBr :

        if sim3 == "ECL" or sim3 == "TEM" :
            outputECLPVTW(OutU,fSim,clsBLK,clsUNI)
        elif sim3 == "VIP" :
            outputVIPPVTW(OutU,fSim,clsBLK,clsUNI)        

#== CMG DENSITY, Compressibility and Water Phase Properties ===========

    if sim3 == "CMG" :
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

    if OutU[:3] == "MET" :   #-- Metric-Equivalent
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

    if clsSIM.setBr :

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

    if   OutU[:3] == "MET" : uUnit = "METRIC"
    elif OutU[:2] == "SI"  : uUnit = "SI"
    else                   : uUnit = "FIELD"

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

    if   OutU[:3] == "MET" : uUnit = "METRIC"
    elif OutU[:2] == "SI"  : uUnit = "SI"
    else                   :uUnit = "FIELD"

    pRefW = clsSIM.pRefW  #-- Reference Pressure        [psia   or bara]
    bRefW = clsSIM.bRefW  #--           FVF             [rb/stb or rm3/sm3]
    uRefW = clsSIM.uRefW  #--           Viscosity       [cP]
    cRefW = clsSIM.cRefW  #--           Compressibility [1/psi  or 1/bar]
    vRefW = clsSIM.vRefW  #--           Viscosibility   [1/psi  or 1/bar]

    uUni = "     cP    "
    if uUnit == "METRIC" :
        pUni  = "    BARA   "
        bUni  = "  RM3/SM3  "
        cUni  = "   1/BAR   "
        pRefW = clsUNI.I2X(pRefW,"bara")
        cRefW = clsUNI.I2X(cRefW,"1/bar")
        vRefW = clsUNI.I2X(vRefW,"1/bar")
    else :
        pUni = "    PSIA   "
        bUni = "   RB/STB  "
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

    if   OutU[:3] == "MET" : uUnit = "METRIC"
    elif OutU[:2] == "SI"  : uUnit = "SI"
    else                   : uUnit = "FIELD"

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

    if OutU[:3] == "MET" :   #-- Metric-Equivalent
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
#  End of Module
#========================================================================
