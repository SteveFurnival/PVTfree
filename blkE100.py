#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Handles Output of E100 Tables
#
#========================================================================

#!/usr/bin/python3

import blkGrad   as BG
import blkOther  as BO
import blkProps  as BP
import plotBatch as PB

def outE100(f100,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp
    gTyp = clsBLK.gTyp

    qDep = clsBLK.qDep

#== Oil Keywords =====================================================

    outputHeaderOil(f100,clsBLK,clsIO)

    if   oTyp == "PVDO" :
        outPVDO(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
        fOil = []
    elif oTyp == "PVCO" :
        outPVCO(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
        fOil = []
    elif oTyp == "PVTO" :
        fOil = outPVTO(f100,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)

#== Gas Keywords =====================================================    

    outputHeaderGas(f100,clsBLK,clsIO)

    if   gTyp == "PVDG" :
        outPVDG(f100,dTab,eTab,clsBLK,clsUNI,clsIO)
        fGas = []
    elif gTyp == "PVTG" :
        fGas = outPVTG(f100,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)

#== Composition versus Depth information? ============================

    if qDep :

        pass
    
#== Plot the Data ====================================================

    PB.blackPlots(dTab,eTab,fOil,fGas,clsBLK,clsUNI,clsIO)

#== No return values =================================================

    return

#========================================================================
#  PVDO Keyword
#========================================================================

def outPVDO(f100,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    sExt = ""

    while iSat >= 0 :

        outputPVTOsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputPVTOsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

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

    fOil = []   #-- Store Under-Saturated Oil Data for Plotting

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        fCol = []

        if iSat == 0 : sExt = "  -- Psat"
        else         : sExt = "  -- Saturated"

        outputPVTOsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

        Pb = dTab[iSat][clsBLK.iPr]
        Rb = dTab[iSat][clsBLK.iRs]
        Bb = dTab[iSat][clsBLK.iBo]
        Ub = dTab[iSat][clsBLK.iUo]

        dRow = [Pb,Bb,Ub,Rb]
        fCol.append(dRow)

        BO.setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK)

#== "Saturated" Under-Saturated Data ==================================

        jSat = iSat - 1
        sExt = ""

        while jSat >= 0 :

            Pr    = dTab[jSat][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)

            #if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Pr,Bo,Uo]
            fCol.append(dRow)

#-- Decrement the jSat counter --------------------------------------            

            jSat -= 1

#== "Extended" Under-Saturated Data ===================================

        jExt = nExt - 1

        while jExt >= 0 :

            if jExt == 0 : sExt = "  /"

            Pr    = eTab[jExt][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)
            
            #if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Pr,Bo,Uo]
            fCol.append(dRow)

#-- Decrement the jExt counter --------------------------------------            

            jExt -= 1

#== Decrement the iSat counter ========================================

        f100.write("\n")

        iSat -= 1

        fOil.append(fCol)
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        fCol = []

        sExt = "  -- Extended"

        outputPVTOsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

        Pb = eTab[iExt][clsBLK.iPr]
        Rb = eTab[iExt][clsBLK.iRs]
        Bb = eTab[iExt][clsBLK.iBo]
        Ub = eTab[iExt][clsBLK.iUo]

        dRow = [Pb,Bb,Ub,Rb]
        fCol.append(dRow)

#== "Extended" Under-Saturated Data ===================================

        jExt = iExt - 1

        while jExt >= 0 :

            if jExt == 0 : sExt = "  /"
            else         : sExt = ""

            Pr    = eTab[jExt][clsBLK.iPr]
            RTp   = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)
            
            #if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Pr,Bo,Uo]
            fCol.append(dRow)

#-- Decrement the jExt counter --------------------------------------            

            jExt -= 1

#== Extra Row =========================================================

        if iExt == 0 :
            
            Pr   = Pr + pInc
            RTp  = clsBLK.RT/Pr
            Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)

            #if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Pr)
            Uo = BP.calcUndViscStand(Pb,Ub,Pr)

            dRow = [Pr,Bo,Uo]
            fCol.append(dRow)

            sExt = "  /"
            outputPVTOund(f100,Pr,Bo,Uo,sExt,clsBLK,clsIO,clsUNI)

        f100.write("\n")

#== Decrement the iExt counter ========================================

        iExt -= 1

        fOil.append(fCol)
        
    f100.write("/\n")
    f100.write("\n")
        
#== Return value ======================================================

    return fOil

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
    Tres = clsBLK.Tres ; cCon = clsBLK.Co
    dSTO = clsBLK.dSTO ; dSTG = clsBLK.dSTG
    mSTO = clsBLK.mSTO ; mSTG = clsBLK.mSTG

    fGas = []   #-- Store Under-Saturated Gas Data for Plotting

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        fCol = []

        if iSat == 0 : sExt = "  -- Psat"
        else         : sExt = "  -- Saturated"

        outputPVTGsat(f100,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

        Pd  = dTab[iSat][clsBLK.iPr]
        Rd  = dTab[iSat][clsBLK.iRv]
        Bd  = dTab[iSat][clsBLK.iBg]
        Ud  = dTab[iSat][clsBLK.iUg]

        UG = BO.gasViscLee(Rd,Bd,Tres,clsBLK)

        dRow = [Rd,Bd,Ud,Pd]
        fCol.append(dRow)
        
        RTp = clsBLK.RT/Pd

#-- Re-fitted EoS/Visc data for under-saturated calculations --------        

        BO.setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK)

#== Undersaturated Data ===============================================

        jSat = iSat + 1
        sExt = ""

        while jSat < nSat :

            Rv = dTab[jSat][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Rv,Bg,Ug]
            fCol.append(dRow)

#-- Decrement the jSat counter --------------------------------------

            jSat += 1

#-- Dry Gas (Rv = 0) stage ------------------------------------------            

        sExt = "  /\n"
        Rv   = 0.0
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK)
        outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        fGas.append(fCol)
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        fCol = []
        sExt = "  -- Extended"

        outputPVTGsat(f100,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

        Pd  = eTab[iExt][clsBLK.iPr]
        Rd  = eTab[iExt][clsBLK.iRv]
        Bd  = eTab[iExt][clsBLK.iBg]
        Ud  = eTab[iExt][clsBLK.iUg]

        UG = BO.gasViscLee(Rd,Bd,Tres,clsBLK)

        dRow = [Rd,Bd,Ud,Pd]
        fCol.append(dRow)
        
        RTp = clsBLK.RT/Pd

#== Undersaturated Data ===============================================

        jExt = iExt + 1
        sExt = ""

        while jExt < nExt :

            Rv = eTab[jExt][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Rv,Bg,Ug]
            fCol.append(dRow)

#-- Increment the jSat counter --------------------------------------

            jExt += 1

#== Undersaturated Data ===============================================

        jSat = 0
        sExt = ""

        while jSat < nSat :

            Rv = dTab[jSat][clsBLK.iRv]
            Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK)
            outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

            dRow = [Rv,Bg,Ug]
            fCol.append(dRow)

#-- Decrement the jSat counter --------------------------------------

            jSat += 1

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#-- Dry Gas (Rv = 0) stage ------------------------------------------            

        sExt = "  /\n"
        Rv   = 0.0
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rv,clsBLK)
        outputPVTGund(f100,Rv,Bg,Ug,sExt,clsBLK,clsIO,clsUNI)

        dRow = [Rv,Bg,Ug]
        fCol.append(dRow)

        fGas.append(fCol)

#== Terminating Slash =================================================

    f100.write("/\n")
    f100.write("\n")
        
#== Return value ======================================================

    return fGas

#========================================================================
#  PVDO, PVCO or PVTO Header
#========================================================================

def outputHeaderOil(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    uUni = "     cP    "
    if OutU[:3] == "MET" :
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
    if OutU[:3] == "MET" :
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

    if OutU[:3] == "MET" :
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

    if   oTyp == "PVCO" :
        sCo   = " {:10.3e}".format(Co)
        sVo   = " {:10.3e}".format(Vo)
        sLabl = sLabl + "  " + sCo + "  " + sVo + "\n"
    elif oTyp == "PVDO" :
        sLabl = "    " + sPr + "  " + sBo + "  " + sUo + sExt
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

    if OutU[:3] == "MET" :
        Pr   = clsUNI.I2X(Pr,"bara")
        sRv  = " {:10.3e}".format(Rv)
        sBg  = " {:10.3e}".format(Bg)
    else :
        Rv   = clsUNI.I2X(Rv,"stb/mscf")
        Bg   = clsUNI.I2X(Bg,"rb/mscf")
        sRv  = " {:10.5f}".format(Rv)
        sBg  = " {:10.5f}".format(Bg)

    sPr = " {:10.3f}".format(Pr)
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

    if OutU[:3] == "MET" : Pr   = clsUNI.I2X(Pr,"bara")

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

    if OutU[:3] == "MET" :
        sRv = " {:10.3e}".format(Rv)
        sBg = " {:10.3e}".format(Bg)
    else :
        Rv = clsUNI.I2X(Rv,"stb/mscf")
        Bg = clsUNI.I2X(Bg,"rb/mscf")
        sRv = " {:10.5f}".format(Rv)
        sBg = " {:10.5f}".format(Bg)

    sPr = "           "
    sUg = " {:10.5f}".format(Ug)

    sLabl = "    " + sPr + "  " + sRv + "  " + sBg + "  " + sUg + sExt + "\n"

    fSim.write(sLabl)
    
    return

#========================================================================
#  End of Module
#========================================================================
