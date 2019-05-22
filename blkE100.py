#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Handles Output of E100 Tables
#
#========================================================================

#!/usr/bin/python3

import blkOther as BO
import blkProps as BP

def outE100(f100,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp
    gTyp = clsBLK.gTyp

#== Oil Keywords =====================================================

    outputHeaderOil(f100,clsBLK,clsIO)

    if   oTyp == "PVDO" :
        outPVDO(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
    elif oTyp == "PVCO" :
        outPVCO(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
    elif oTyp == "PVTO" :  #-- Must be PVTO
        outPVTO(f100,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)

#== Gas Keywords =====================================================    

    outputHeaderGas(f100,clsBLK,clsIO)

    if   gTyp == "PVDG" :
        outPVDG(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
    elif gTyp == "PVTG" :  #-- Must be PVTG
        outPVTG(f100,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)

#== No return values =================================================

    return

#========================================================================
#  PVDO Keyword
#========================================================================

def outPVDO(f100,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

#== No return value ===================================================

    return

#========================================================================
#  PVCO Keyword
#========================================================================

def outPVCO(f100,dTab,eTab,clsBLK,clsUNI,clsIO) :

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

        outputPVTOsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    sExt = ""

    while iExt >= 0 :

        outputPVTOsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

#== Terminating Slash =================================================

    f100.write("/\n")
    f100.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  PVTO Keyword
#========================================================================

def outPVTO(f100,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    qLiq = True
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG
    pInc = clsBLK.Pinc

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        if iSat == 0 : sExt = "  -- Psat"
        else         : sExt = "  -- Saturated"

        outputPVTOsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

        Pb = dTab[iSat][clsBLK.iPr]
        Rb = dTab[iSat][clsBLK.iRs]
        Ub = dTab[iSat][clsBLK.iUo]

#== "Saturated" Under-Saturated Data ==================================

        jSat = iSat - 1
        sExt = ""

        while jSat >= 0 :

            BO.setEoSVis(jSat,sOil,sGas,rOil,rGas,clsBLK)

            Pr    = dTab[jSat][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK,clsIO)
            if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

#-- Decrement the jSat counter --------------------------------------            

            jSat -= 1

#== "Extended" Under-Saturated Data ===================================

        jExt = nExt - 1

        BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

        while jExt >= 0 :

            if jExt == 0 : sExt = "  /"

            Pr    = eTab[jExt][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK,clsIO)
            if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

#-- Decrement the jExt counter --------------------------------------            

            jExt -= 1

#== Decrement the iSat counter ========================================

        f100.write("\n")

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        sExt = "  -- Extended"

        outputPVTOsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

        Pb = eTab[iExt][clsBLK.iPr]
        Rb = eTab[iExt][clsBLK.iRs]
        Ub = eTab[iExt][clsBLK.iUo]

#== "Extended" Under-Saturated Data ===================================

        jExt = iExt - 1

        while jExt >= 0 :

            if jExt == 0 : sExt = "  /"
            else         : sExt = ""

            Pr    = eTab[jExt][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK,clsIO)
            if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

#-- Decrement the jExt counter --------------------------------------            

            jExt -= 1

#== Extra Row =========================================================

        if iExt == 0 :
            Pr   = Pr + pInc
            RTp  = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK,clsIO)
            if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            sExt = "  /"
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

        f100.write("\n")

#== Decrement the iExt counter ========================================

        iExt -= 1

    f100.write("/\n")
    f100.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  PVDG Keyword
#========================================================================

def outPVDG(f100,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        if iSat == 0 : sExt = "  -- Psat"
        else         : sExt = "  -- Saturated"

        outputPVTGsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    sExt = "  -- Extended"

    while iExt >= 0 :

        outputPVTGsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#== Terminating Slash =================================================

    f100.write("/\n")
    f100.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  PVTG Keyword
#========================================================================

def outPVTG(f100,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    qVap = False
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        if iSat == 0 : sExt = "  -- Psat"
        else         : sExt = "  -- Saturated"

        outputPVTGsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

        Pr  = dTab[iSat][clsBLK.iPr]
        RTp = clsBLK.RT/Pr

#-- Re-fitted EoS/Visc data for under-saturated calculations --------        

        BO.setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK)

#== Undersaturated Data ===============================================

        jSat = iSat + 1
        sExt = ""

        while jSat < nSat :

            Rv = dTab[jSat][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#-- Decrement the jSat counter --------------------------------------

            jSat += 1

#-- Dry Gas (Rv = 0) stage ------------------------------------------            

        sExt = "  /\n"
        Rv   = 0.0
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)
        outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        sExt = "  -- Extended"

        outputPVTGsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

        Pr  = eTab[iExt][clsBLK.iPr]
        RTp = clsBLK.RT/Pr

#-- Re-fitted EoS/Visc data for under-saturated calculations --------        

        BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

#== Undersaturated Data ===============================================

        jExt = iExt + 1
        sExt = ""

        while jExt < nExt :

            Rv = eTab[jExt][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#-- Increment the jSat counter --------------------------------------

            jExt += 1

#== Undersaturated Data ===============================================

        jSat = 0
        sExt = ""

        while jSat < nSat :

            Rv = dTab[jSat][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#-- Decrement the jSat counter --------------------------------------

            jSat += 1

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#-- Dry Gas (Rv = 0) stage ------------------------------------------            

        sExt = "  /\n"
        Rv   = 0.0
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK,clsIO)
        outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#== Terminating Slash =================================================

    f100.write("/\n")
    f100.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  PVDO, PVCO or PVTO Header
#========================================================================

def outputHeaderOil(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    uUni = "     cP    "
    if OutU == "MET" :
        rUni = "  SM3/SM3  "
        pUni = "    BARA   "
        bUni = "  RM3/SM3  "
        cUni = "  1/BARA   "
    else :
        rUni = "  MSCF/STB "
        pUni = "    PSIA   "
        bUni = "   RB/STB  "
        cUni = "  1/PSIA   "

    if   oTyp == "PVDO" :
        sLab1 = "--      GOR        Oil FVF      Oil Visc\n"
        sLab2 = "--  " + rUni + "  " + bUni + "  " + uUni + "\n"
    elif oTyp == "PVCO" :
        sLab1 = "--      Pbub         GOR        Oil FVF      Oil Visc      Comp        Viscos\n"
        sLab2 = "--  " + pUni + "  " + rUni + "  " + bUni + "  " + uUni + "  " + cUni + "  " + cUni + "\n"
    else :
        sLab1 = "--      GOR          PBub       Oil FVF      Oil Visc\n"
        sLab2 = "--  " + rUni + "  " + pUni + "  " + bUni + "  " + uUni + "\n"

    fSim.write("--\n")
    fSim.write("--  Oil Phase Properties\n")
    fSim.write("--\n")
    fSim.write("\n")
    fSim.write(oTyp + "\n")
    fSim.write(sLab1)
    fSim.write(sLab2)
    fSim.write("\n")

    return

#========================================================================
#  PVDG or PVTG Header
#========================================================================

def outputHeaderGas(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    gTyp = clsBLK.gTyp

    uUni = "     cP    "
    if OutU == "MET" :
        pUni = "    BARA   "
        rUni = "  SM3/SM3  "
        bUni = "  RM3/SM3  "
    else :
        pUni = "    PSIA   "
        rUni = "  STB/MSCF "
        bUni = "   RB/MSCF "

    if gTyp == "PVDG" :
        sLab1 = "--      Pres       Gas FVF      Gas Visc\n"
        sLab2 = "--  " + pUni + "  " + bUni + "  " + uUni + "\n"
    else :
        sLab1 = "--      Pres         CGR        Gas FVF      Gas Visc\n"
        sLab2 = "--  " + pUni + "  " + rUni + "  " + bUni + "  " + uUni + "\n"

    fSim.write("--\n")
    fSim.write("--  Gas Phase Properties\n")
    fSim.write("--\n")
    fSim.write("\n")
    fSim.write(gTyp + "\n")
    fSim.write(sLab1)
    fSim.write(sLab2)
    fSim.write("\n")

    return

#========================================================================
#  Write Saturated (4-values) line of PVTO data, or
#                  (6-values) line of PVCO data
#========================================================================

def outputPVTOsat(fSim,iRow,fTab,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    Pr = fTab[iRow][clsBLK.iPr]
    Rs = fTab[iRow][clsBLK.iRs]
    Bo = fTab[iRow][clsBLK.iBo]
    Uo = fTab[iRow][clsBLK.iUo]

    if oTyp == "PVCO" :
        Co = fTab[iRow][clsBLK.iCo]
        Vo = fTab[iRow][clsBLK.iVo]
    else :
        Co = 0.0
        Vo = 0.0

    if OutU == "MET" :
        Pr = clsUNI.I2X(Pr,"bara")
        Co = clsUNI.I2X(Pr,"1/bar")
        Vo = clsUNI.I2X(Pr,"1/bar")
    else             :
        Rs = clsUNI.I2X(Rs,"mscf/stb")

    sRs = " {:10.5f}".format(Rs)
    sPr = " {:10.3f}".format(Pr)
    sBo = " {:10.5f}".format(Bo)
    sUo = " {:10.5f}".format(Uo)

    sLabl = "    " + sRs + "  " + sPr + "  " + sBo + "  " + sUo + sExt

    if oTyp == "PVCO" :
        sCo   = " {:10.3e}".format(Co)
        sVo   = " {:10.3e}".format(Vo)
        sLabl = sLabl + "  " + sCo + "  " + sVo + "\n"
    else :
        sLabl = sLabl + "\n"

    fSim.write(sLabl)
    
    return

#========================================================================
#  Write Saturated (4-values) line of PVTG data, or
#                  (3-values) line of PVDG data
#========================================================================

def outputPVTGsat(fSim,iRow,fTab,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    gTyp = clsBLK.gTyp

    Pr = fTab[iRow][clsBLK.iPr]

    if gTyp == "PVDG" :
        Rv = 0.0
        Bg = fTab[iRow][clsBLK.iBd]
        Ug = fTab[iRow][clsBLK.iUd]
    else :
        Rv = fTab[iRow][clsBLK.iRv]
        Bg = fTab[iRow][clsBLK.iBg]
        Ug = fTab[iRow][clsBLK.iUg]

    if OutU == "MET" :
        Pr   = clsUNI.I2X(Pr,"bara")
    else :
        Rv   = clsUNI.I2X(Rv,"stb/mscf")
        Bg   = clsUNI.I2X(Bg,"rb/mscf")

    sPr = " {:10.3f}".format(Pr)
    sRv = " {:10.5f}".format(Rv)
    sBg = " {:10.5f}".format(Bg)
    sUg = " {:10.5f}".format(Ug)

    if gTyp == "PVDG" :
        sLabl = "    " + sPr + "  " + sBg + "  " + sUg + sExt + "\n"
    else :
        sLabl = "    " + sPr + "  " + sRv + "  " + sBg + "  " + sUg + sExt + "\n"

    fSim.write(sLabl)
    
    return

#========================================================================
#  Write Under-Saturated (3-values) line of PVTO data
#========================================================================

def outputPVTOund(fSim,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    if OutU == "MET" :
        Pr   = clsUNI.I2X(Pr,"bara")

    sRs = "           "
    sPr = " {:10.3f}".format(Pr)
    sBo = " {:10.5f}".format(Bo)
    sUo = " {:10.5f}".format(Uo)

    sLabl = "    " + sRs + "  " + sPr + "  " + sBo + "  " + sUo + sExt + "\n"

    fSim.write(sLabl)
    
    return

#========================================================================
#  Write Under-Saturated (3-values) line of PVTG data
#========================================================================

def outputPVTGund(fSim,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    if OutU == "MET" :
        pass
    else :
        Rv = clsUNI.I2X(Rv,"stb/mscf")
        Bg = clsUNI.I2X(Bg,"rb/mscf")

    sPr = "           "
    sRv = " {:10.5f}".format(Rv)
    sBg = " {:10.5f}".format(Bg)
    sUg = " {:10.5f}".format(Ug)

    sLabl = "    " + sPr + "  " + sRv + "  " + sBg + "  " + sUg + sExt + "\n"

    fSim.write(sLabl)
    
    return

#========================================================================
#  End of Module
#========================================================================
