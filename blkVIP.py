
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Generates Black-Oil Output for VIP or Nexus
#
#========================================================================

#!/usr/bin/python3

import blkOther  as BO
import blkProps  as BP
import plotBatch as PB
import utilities as UT

def outVIP(fVIP,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp
    gTyp = clsBLK.gTyp

#== Oil Keywords ======================================================

    outHeaderOil(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO)

    if   oTyp == "BODTAB" :
        outBODTAB(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO)
        fOil = []
    else :   #-- "BOOTAB" is default
        fOil = outBOOTAB(fVIP,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO)

#== Gas Keywords ======================================================    

    outHeaderGas(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO)

    if   gTyp == "BDGTAB" :
        outBDGTAB(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO)
        fGas = []
    else :   #-- "BOGTAB" is default
        fGas = outBOGTAB(fVIP,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO)

#== Plot the Data ====================================================

    PB.blackPlots(dTab,eTab,fOil,fGas,clsBLK,clsUNI,clsIO)

#== No return value ===================================================

    return

#========================================================================
#  Oil Table Header
#========================================================================

def outHeaderOil(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO) :

    oTyp = clsBLK.oTyp

    fVIP.write("C \n")
    fVIP.write("C Oil Phase Properties\n")
    fVIP.write("C \n")
    fVIP.write("\n")

#== Table Headers =====================================================

    dSTO = clsBLK.dSTO  #-- Stock Tank Oil Density     [lb/ft3]
    mSTO = clsBLK.mSTO  #-- Stock Tank Oil Mole Weight

    dSTO = clsUNI.I2X(dSTO,"gm/cm3")

    fVIP.write(oTyp + "  1\n")
    fVIP.write("     DOB         WRTO\n")

    sDens = " {:10.5f}".format(dSTO)
    sMolW = " {:10.3f}".format(mSTO)

    sLabl = sDens + "  " + sMolW + "\n"
    fVIP.write(sLabl)
    fVIP.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Gas Table Header
#========================================================================

def outHeaderGas(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO) :

    gTyp = clsBLK.gTyp

    fVIP.write("C \n")
    fVIP.write("C Gas Phase Properties\n")
    fVIP.write("C \n")
    fVIP.write("\n")

#== Table Headers =====================================================

    dSTG = clsBLK.dSTO  #-- Stock Tank Oil Density     [lb/ft3]

    dSTG = clsUNI.I2X(dSTG,"gm/cm3")

    fVIP.write(gTyp + "  1\n")
    fVIP.write("     DGB\n")

    sDens = " {:10.5f}".format(dSTG)

    sLabl = sDens + "\n"
    fVIP.write(sLabl)
    fVIP.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  BODTAB Table
#========================================================================

def outBODTAB(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

#== Column Headers ====================================================

    fVIP.write("    PSAT          BO            VO\n")

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        outputDeadOil(fVIP,iSat,dTab,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputDeadOil(fVIP,iExt,eTab,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

#== Terminating Slash =================================================

    fVIP.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  BOOTAB Table
#========================================================================

def outBOOTAB(fVIP,dTab,eTab,sOil,sGas,rOil,rGas,qMonV,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    qLiq = True
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

#== Column Headers ====================================================

    fVIP.write("C \n")
    fVIP.write("C  Saturated Data\n")
    fVIP.write("C \n")
    fVIP.write("\n")
    fVIP.write("      RS         PSAT         BO           VO\n")

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        outputLiveOil(fVIP,iSat,dTab,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputLiveOil(fVIP,iExt,eTab,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1

#== Terminating Slash =================================================

    fVIP.write("\n")

#----------------------------------------------------------------------
#  Undersaturated Data
#----------------------------------------------------------------------

    fOil = []
        
    fVIP.write("C \n")
    fVIP.write("C  Undersaturated Data\n")
    fVIP.write("C \n")
    fVIP.write("\n")

    lSat = 0
    lExt = 0

#== Saturated Data ====================================================    

    iSat = nSat - 1

    while iSat >= 0 :

        Ps = dTab[iSat][clsBLK.iPr]
        Rs = dTab[iSat][clsBLK.iRs]
        Us = dTab[iSat][clsBLK.iUo]
        
        BO.setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK)
        
        dOut = outputUnderSatOil(Ps,Rs,Us,dTab,eTab,lSat,lExt,qMonV,clsBLK,clsIO,clsUNI)
        fVIP.write("\n")
        lSat += 1

#-- Decrement the iSat counter --------------------------------------        

        iSat -= 1
        fOil.append(dOut)

#== Extended Data =====================================================

    iExt = nExt - 1

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :

        Ps = eTab[iExt][clsBLK.iPr]
        Rs = eTab[iExt][clsBLK.iRs]
        Us = eTab[iExt][clsBLK.iUo]

        dOut = outputUnderSatOil(Ps,Rs,Us,dTab,eTab,lSat,lExt,qMonV,clsBLK,clsIO,clsUNI)
        fVIP.write("\n")
        lExt += 1

#-- Decrement the iExt counter ========================================

        iExt -= 1
        fOil.append(dOut)

#== No return value ===================================================

    return fOil

#========================================================================
#  Blocks of Undersaturated Oil Data
#========================================================================

def outputUnderSatOil(Pb,Rb,Ub,dTab,eTab,lSat,lExt,qMonV,clsBLK,clsIO,clsUNI) :

#-- Initialisation --------------------------------------------------    

    OutU = clsBLK.OutU
    fVIP = clsIO.fVIP
    nSat = len(dTab)
    nExt = len(eTab)
    qLiq = True
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

    if OutU[:3] == "MET" :
        sPrs = "kpa"
        sGOR = "sm3/sm3"
        sFVF = "rm3/sm3"
    else :
        sPrs = "psia"
        sGOR = "scf/stb"
        sFVF = "rb/stb"
        
    nTot = nSat + nExt - lSat - lExt

    nCol = 4
    nMod =     nTot % nCol
    nRow = int(nTot / nCol)

    if nMod > 0 : nRow = nRow + 1

    iSat = nSat - 1 - lSat
    iExt = nExt - 1 - lExt

    dOut = []
    qFrs = True

#-- Saturated Data --------------------------------------------------    

    while iSat >= 0 :
        
        Ps = dTab[iSat][clsBLK.iPr]
        Rs = dTab[iSat][clsBLK.iRs]
        
        RTp = clsBLK.RT/Ps

        Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)

        if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Ps)
            
        if qFrs : dRow = [Ps,Bo,Uo,Rb] ; qFrs = False
        else    : dRow = [Ps,Bo,Uo]
        
        dOut.append(dRow)
            
#-- Decrement iSat --------------------------------------------------

        iSat -= 1

#-- Extended Data ---------------------------------------------------        

    while iExt >= 0 :
        
        Ps = eTab[iExt][clsBLK.iPr]
        Rs = eTab[iExt][clsBLK.iRs]
        
        RTp = clsBLK.RT/Ps

        Bo,Uo = BO.calcSatProp(qLiq,RTp,cCon,dSTO,dSTG,Rb,clsBLK)
        
        if not qMonV : Uo = BP.calcUndViscStand(Pb,Ub,Ps)
            
        if qFrs : dRow = [Ps,Bo,Uo,Rb] ; qFrs = False
        else    : dRow = [Ps,Bo,Uo]
        
        dOut.append(dRow)
            
#-- Decrement iSat --------------------------------------------------

        iExt -= 1

#== Output nRow's worth of Triplets ===================================

    iVal = 0
    qExt = False

    Pb  = clsUNI.I2X(Pb,sPrs)
    sPb = " {:10.3f}".format(Pb)

    for iRow in range(nRow) :
        sLin1 = "PSAT               "
        sLin2 = "P            "
        sLin3 = sPb + "  "
        for jCol in range(nCol) :
            Pr  = dOut[iVal][0]
            Bo  = dOut[iVal][1]
            Uo  = dOut[iVal][2]
            Pr  = clsUNI.I2X(Pr,sPrs)
            Bo  = clsUNI.I2X(Bo,sFVF)
            sPr = " {:10.3f}".format(Pr)
            sBo = " {:10.5f}".format(Bo)
            sUo = " {:10.5f}".format(Uo)
            sLin1 = sLin1 + sPr + "             "
            sLin2 = sLin2 + "    BO      " + "     VO     "
            sLin3 = sLin3 + sBo + "  " + sUo
            iVal += 1
            if iVal == nTot :
                qExt = True
                break
        sLin1 = sLin1 + "\n"
        sLin2 = sLin2 + "\n"
        sLin3 = sLin3 + "\n"
        fVIP.write(sLin1)
        fVIP.write(sLin2)
        fVIP.write(sLin3)
        if qExt : break

#== Return value ======================================================

    return dOut

#========================================================================
#  BDGTAB Keyword
#========================================================================

def outBDGTAB(fVIP,dTab,eTab,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

#== Column Headers ====================================================

    fVIP.write("    PSAT         BG           VG\n")

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        outputDeadGas(fVIP,iSat,dTab,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputDeadGas(fVIP,iExt,eTab,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#== Terminating Slash =================================================

    fVIP.write("\n")
        
#== No return value ===================================================

    return

#========================================================================
#  BOGTAB Keyword
#========================================================================

def outBOGTAB(fVIP,dTab,eTab,sOil,sGas,rOil,rGas,clsBLK,clsUNI,clsIO) :

    nSat = len(dTab)
    nExt = len(eTab)

    iSat = nSat - 1
    iExt = nExt - 1

    qVap = False
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

#== Column Headers ====================================================

    fVIP.write("C \n")
    fVIP.write("C  Saturated Data\n")
    fVIP.write("C \n")
    fVIP.write("\n")
    fVIP.write("    PSAT         RV           BG           VG\n")

#----------------------------------------------------------------------
#  Saturated Data
#----------------------------------------------------------------------

    while iSat >= 0 :

        outputLiveGas(fVIP,iSat,dTab,clsBLK,clsIO,clsUNI)

#== Decrement the iSat counter ========================================

        iSat -= 1
        
#----------------------------------------------------------------------
#  Extended Data
#----------------------------------------------------------------------

    while iExt >= 0 :

        outputLiveGas(fVIP,iExt,eTab,clsBLK,clsIO,clsUNI)

#== Decrement the iExt counter ========================================

        iExt -= 1
        
#== Terminating Slash =================================================

    fVIP.write("\n")
        
#----------------------------------------------------------------------
#  Undersaturated Data
#----------------------------------------------------------------------

    fGas = []

    fVIP.write("C \n")
    fVIP.write("C  Undersaturated Data\n")
    fVIP.write("C \n")
    fVIP.write("\n")

    Rv   = 0.0          #-- Dry Gas First
    Pr   = UT.pStand    #-- At Standard Pressure
    lSat = 0
    lExt = 0

    outputUnderSatGas(Rv,Pr,dTab,eTab,lSat,lExt,sOil,sGas,rOil,rGas,clsBLK,clsIO,clsUNI)
    fVIP.write("\n")

#== Saturated Data ====================================================    

    iSat = nSat - 1

    while iSat >= 0 :

        Rv = dTab[iSat][clsBLK.iRv]
        Pr = dTab[iSat][clsBLK.iPr]

        dOut = outputUnderSatGas(Rv,Pr,dTab,eTab,lSat,lExt,sOil,sGas,rOil,rGas,clsBLK,clsIO,clsUNI)
        fVIP.write("\n")
        lSat += 1

#-- Decrement the iSat counter --------------------------------------        

        iSat -= 1
        fGas.append(dOut)

#== Extended Data =====================================================

    iExt = nExt - 1

    while iExt >= 0 :

        Rv = eTab[iExt][clsBLK.iRv]
        Pr = eTab[iExt][clsBLK.iPr]
        
        dOut = outputUnderSatGas(Rv,Pr,dTab,eTab,lSat,lExt,sOil,sGas,rOil,rGas,clsBLK,clsIO,clsUNI)
        fVIP.write("\n")
        lExt += 1

#-- Decrement the iExt counter ========================================

        iExt -= 1
        fGas.append(dOut)

#== No return value ===================================================

    return fGas

#========================================================================
#  Blocks of Undersaturated Gas Data
#========================================================================

def outputUnderSatGas(Rd,Pd,dTab,eTab,lSat,lExt,sOil,sGas,rOil,rGas,clsBLK,clsIO,clsUNI) :

#-- Initialisation --------------------------------------------------    

    OutU = clsBLK.OutU
    fVIP = clsIO.fVIP
    nSat = len(dTab)
    nExt = len(eTab)
    qVap = False
    cCon = clsBLK.Co
    dSTO = clsBLK.dSTO
    dSTG = clsBLK.dSTG

    if OutU[:3] == "MET" :
        sPrs = "kpa"
        sCGR = "sm3/sm3"
        sFVF = "rm3/sm3"
    else :
        sPrs = "psia"
        sCGR = "stb/mscf"
        sFVF = "rb/mscf"
        
    nTot = nSat + nExt - lSat - lExt

    nCol = 4
    nMod =     nTot % nCol
    nRow = int(nTot / nCol)

    if nMod > 0 : nRow = nRow + 1

    iSat = nSat - 1 - lSat
    iExt = nExt - 1 - lExt

    dOut = []
    qFrs = True

#-- Saturated Data --------------------------------------------------    

    while iSat >= 0 :
        
        Pr  = dTab[iSat][clsBLK.iPr]
        Rv  = dTab[iSat][clsBLK.iRv]
        
        RTp = clsBLK.RT/Pr
        BO.setEoSVis(iSat,sOil,sGas,rOil,rGas,clsBLK)
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rd,clsBLK)

        if qFrs : dRow = [Rv,Bg,Ug,Rd] ; qFrs = False
        else    : dRow = [Rv,Bg,Ug]
        
        dOut.append(dRow)
            
#-- Decrement iSat --------------------------------------------------

        iSat -= 1

#-- Extended Data ---------------------------------------------------        

    BO.setEoSVis(0,sOil,sGas,rOil,rGas,clsBLK)

    while iExt >= 0 :
        
        Pr  = eTab[iExt][clsBLK.iPr]
        Rv  = eTab[iExt][clsBLK.iRv]
        
        RTp = clsBLK.RT/Pr
        Bg,Ug = BO.calcSatProp(qVap,RTp,cCon,dSTO,dSTG,Rd,clsBLK)
        
        if qFrs : dRow = [Rv,Bg,Ug,Rd] ; qFrs = False
        else    : dRow = [Rv,Bg,Ug]
        
        dOut.append(dRow)
            
#-- Decrement iSat --------------------------------------------------

        iExt -= 1

#== Output nRow's worth of Triplets ===================================

    iVal = 0
    qExt = False

    Rd  = clsUNI.I2X(Rd,sCGR)
    sRd = " {:10.5f}".format(Rd)

    for iRow in range(nRow) :
        sLin1 = "RVSAT              "
        sLin2 = "RV           "
        sLin3 = sRd + "  "
        for jCol in range(nCol) :
            Rv  = dOut[iVal][0]
            Bg  = dOut[iVal][1]
            Ug  = dOut[iVal][2]
            Rv  = clsUNI.I2X(Rv,sCGR)
            Bg  = clsUNI.I2X(Bg,sFVF)
            sRv = " {:10.5f}".format(Rv)
            sBg = " {:10.5f}".format(Bg)
            sUg = " {:10.5f}".format(Ug)
            sLin1 = sLin1 + sRv + "             "
            sLin2 = sLin2 + "    BG      " + "     VG     "
            sLin3 = sLin3 + sBg + "  " + sUg
            iVal += 1
            if iVal == nTot :
                qExt = True
                break
        sLin1 = sLin1 + "\n"
        sLin2 = sLin2 + "\n"
        sLin3 = sLin3 + "\n"
        fVIP.write(sLin1)
        fVIP.write(sLin2)
        fVIP.write(sLin3)
        if qExt : break

#== Return value ======================================================

    return dOut

#========================================================================
#  Write line of Dead Oil Data (3 values)
#========================================================================

def outputDeadOil(fSim,iRow,fTab,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Pr = fTab[iRow][clsBLK.iPr]
    Bo = fTab[iRow][clsBLK.iBo]
    Uo = fTab[iRow][clsBLK.iUo]

    if OutU[:3] == "MET" : Pr = clsUNI.I2X(Pr,"kpa")

    sPr = " {:10.3f}".format(Pr)
    sBo = " {:10.5f}".format(Bo)
    sUo = " {:10.5f}".format(Uo)

    sLabl = sPr + "  " + sBo + "  " + sUo + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  Write line of Dead Gas Data (3 values)
#========================================================================

def outputDeadGas(fSim,iRow,fTab,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Pr = fTab[iRow][clsBLK.iPr]
    Bd = fTab[iRow][clsBLK.iBd]
    Ud = fTab[iRow][clsBLK.iUd]

    if OutU[:3] == "MET" : Pr = clsUNI.I2X(Pr,"kpa")
    else             : Bd = clsUNI.I2X(Bd,"rb/mscf")

    sPr = " {:10.3f}".format(Pr)
    sBd = " {:10.5f}".format(Bd)
    sUd = " {:10.5f}".format(Ud)

    sLabl = sPr + "  " + sBd + "  " + sUd + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  Write line of Live Oil Data (4 values)
#========================================================================

def outputLiveOil(fSim,iRow,fTab,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Rs = fTab[iRow][clsBLK.iRs]
    Pr = fTab[iRow][clsBLK.iPr]
    Bo = fTab[iRow][clsBLK.iBo]
    Uo = fTab[iRow][clsBLK.iUo]

    if OutU[:3] == "MET" :
        Pr = clsUNI.I2X(Pr,"kpa")
    else :
        Rs = clsUNI.I2X(Rs,"scf/stb")

    sRs = " {:10.3f}".format(Rs)
    sPr = " {:10.3f}".format(Pr)
    sBo = " {:10.5f}".format(Bo)
    sUo = " {:10.5f}".format(Uo)

    sLabl = sRs + "  " + sPr + "  " + sBo + "  " + sUo + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  Write line of Live Gas Data (4 values)
#========================================================================

def outputLiveGas(fSim,iRow,fTab,clsBLK,clsIO,clsUNI) :

    OutU = clsBLK.OutU               #-- FLD (Field) or MET (Metric)

    Pr = fTab[iRow][clsBLK.iPr]
    Rv = fTab[iRow][clsBLK.iRv]
    Bg = fTab[iRow][clsBLK.iBg]
    Ug = fTab[iRow][clsBLK.iUg]

    if OutU[:3] == "MET" :
        Pr = clsUNI.I2X(Pr,"kpa")
    else :
        Rv = clsUNI.I2X(Rv,"stb/mscf")
        Bg = clsUNI.I2X(Bg,"rb/mscf")

    sPr = " {:10.3f}".format(Pr)
    sRv = " {:10.5f}".format(Rv)
    sBg = " {:10.5f}".format(Bg)
    sUg = " {:10.5f}".format(Ug)

    sLabl = sPr + "  " + sRv + "  " + sBg + "  " + sUg + "\n"

    fSim.write(sLabl)

#== No return value ===================================================    
    
    return

#========================================================================
#  End of Module
#========================================================================

