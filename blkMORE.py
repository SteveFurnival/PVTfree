
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Generates Black-Oil Output for Tempest-MORE
#
#========================================================================

#!/usr/bin/python3

import blkE100   as B1
import plotBatch as PB

def outMORE(fMOR,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp
    gTyp = clsBLK.gTyp

#== Oil Keywords ======================================================

    outputHeaderOil(fMOR,clsBLK,clsIO)

    if   oTyp == "OPVT" :
        outOPVT(fMOR,dTab,eTab,clsBLK,clsUNI,clsIO)
        fOil = []
    elif oTyp == "PVTO" :  #-- Must be PVTO
        fOil = B1.outPVTO(fMOR,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)

#== Gas Keywords ======================================================    

    outputHeaderGas(fMOR,clsBLK,clsIO)

    if   gTyp == "GPVT" :
        outGPVT(fMOR,dTab,eTab,clsBLK,clsUNI,clsIO)
        fGas = []
    elif gTyp == "PVTG" :  #-- Must be PVTG
        fGas = B1.outPVTG(fMOR,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)

#== Plot the Data ====================================================

    PB.blackPlots(dTab,eTab,fOil,fGas,clsBLK,clsUNI,clsIO)

#== No return value ===================================================

    return

#========================================================================
#  OPVT Table
#========================================================================

def outOPVT(fMOR,dTab,eTab,clsBLK,clsUNI,clsIO) :

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

        outputOPVTsat(fMOR,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    sExt = ""

    while iExt >= 0 :

        outputOPVTsat(fMOR,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

#== Terminating Slash =================================================

    fMOR.write("/\n")
    fMOR.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  GPVT Keyword
#========================================================================

def outGPVT(fMOR,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    sExt = ""

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        outputGPVTsat(fMOR,iSat,dTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputGPVTsat(fMOR,iExt,eTab,sExt,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#== Terminating Slash =================================================

    fMOR.write("/\n")
    fMOR.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  OPVT or PVTO Header
#========================================================================

def outputHeaderOil(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    oTyp = clsBLK.oTyp

    uUni = "     cP    "
    if OutU[:3] == "MET" :
        rUni = " KSM3/SM3  "
        pUni = "    BARA   "
        bUni = "  RM3/SM3  "
        cUni = "  1/BARA   "
    else :
        rUni = "  MSCF/STB "
        pUni = "    PSIA   "
        bUni = "   RB/STB  "
        cUni = "  1/PSIA   "

    if   oTyp == "OPVT" :
        sLab1 = "--      PBub       Oil FVF      Oil Visc      GOR         Comp         Viscos\n"
        sLab2 = "--  " + pUni + "  " + bUni + "  " + uUni + "  " + rUni + "  " + cUni + "  " + cUni + "\n"
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
#  GPVT or PVTG Header
#========================================================================

def outputHeaderGas(fSim,clsBLK,clsIO) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    gTyp = clsBLK.gTyp

    uUni = "     cP    "
    if OutU[:3] == "MET" :
        pUni = "    BARA   "
        rUni = " SM3/KSM3  "
        bUni = " RM3/KSM3  "
    else :
        pUni = "    PSIA   "
        rUni = "  STB/MSCF "
        bUni = "   RB/MSCF "

    if gTyp == "GPVT" :
        sLab1 = "--      Pres       Gas FVF      Gas Visc       CGR\n"
        sLab2 = "--  " + pUni + "  " + bUni + "  " + uUni + "  " + rUni + "\n"
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
#  Write Saturated (6-values) line of OPVT data
#========================================================================

def outputOPVTsat(fSim,iRow,fTab,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Pr = fTab[iRow][clsBLK.iPr]
    Rs = fTab[iRow][clsBLK.iRs]
    Bo = fTab[iRow][clsBLK.iBo]
    Uo = fTab[iRow][clsBLK.iUo]

    Co = fTab[iRow][clsBLK.iCo]
    Vo = fTab[iRow][clsBLK.iVo]

    if OutU[:3] == "MET" :
        Pr = clsUNI.I2X(Pr,"bara")
        Rs = clsUNI.I2X(Rs,"ksm3/sm3")
        Co = clsUNI.I2X(Co,"1/bar")
        Vo = clsUNI.I2X(Vo,"1/bar")
    else             :
        Rs = clsUNI.I2X(Rs,"mscf/stb")

    sRs = " {:10.5f}".format(Rs)
    sPr = " {:10.3f}".format(Pr)
    sBo = " {:10.5f}".format(Bo)
    sUo = " {:10.5f}".format(Uo)
    sCo = " {:10.3e}".format(Co)
    sVo = " {:10.3e}".format(Vo)

    sLabl = "    " + sPr + "  " + sBo + "  " + sUo + "  " + sRs + sExt + "  " + sCo + "  " + sVo + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  Write Saturated (4-values) line of GPVT data
#========================================================================

def outputGPVTsat(fSim,iRow,fTab,sExt,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Pr = fTab[iRow][clsBLK.iPr]
    Bg = fTab[iRow][clsBLK.iBd]
    Ug = fTab[iRow][clsBLK.iUd]
    Rv = fTab[iRow][clsBLK.iRv]

    if OutU[:3] == "MET" :
        Pr = clsUNI.I2X(Pr,"bara")
        Bg = clsUNI.I2X(Bg,"rm3/ksm3")
        Rv = clsUNI.I2X(Rv,"sm3/ksm3")
    else :
        Bg = clsUNI.I2X(Bg,"rb/mscf")
        Rv = clsUNI.I2X(Rv,"stb/mscf")

    sPr = " {:10.3f}".format(Pr)
    sBg = " {:10.5f}".format(Bg)
    sUg = " {:10.5f}".format(Ug)
    sRv = " {:10.5f}".format(Rv)

    sLabl = "    " + sPr + "  " + sBg + "  " + sUg + "  " + sRv + sExt + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  End of Module
#========================================================================

