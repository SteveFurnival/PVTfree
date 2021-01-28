
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Calculates Composition Versus Depth Variation for Black Oil Tables
#
#========================================================================

#!/usr/bin/python3

import numpy    as NP

import calcGrad as CG

#=======================================================================
#  Driver for Blackoil Table Calculations
#=======================================================================

def blackGrad(clsEOS,dicSAM,clsBLK,clsUNI,clsIO) :

    tRef = clsBLK.Tres

    iSam = clsBLK.sNum
    nCom = clsEOS.nComp

    clsSAM = dicSAM[iSam]

    zRef = NP.zeros(nCom)
    for iC in range(nCom) : zRef[iC] = clsSAM.gZI(iC)

    dRef = clsBLK.dRef
    dMin = clsBLK.dMin
    dMax = clsBLK.dMax
    dInc = clsBLK.dInc
    pRef = clsBLK.pRef

    pBvD = []   #-- Will hold Calculated PB vs. Depth
    pDvD = []   #--                      PD vs. Depth

#== Get Reference Depth (log) Fugacity Coeffs =========================

    iNeu = 0

    fRef,dumP = \
        CG.setupCalcGRDCoefsP(iNeu,pRef,tRef,zRef,clsEOS)  #-- Dont need Ref dlnPhi/dP

#== Minimum (Top) Depth ===============================================

    iTop = 0
    zTyp = "Top"
    dTop = dMin - dRef

    pTop,zTop,qTop,sTop = \
        CG.calcTopBottom(zTyp,dTop,pRef,tRef,zRef,fRef,clsEOS,clsIO)

    if not qTop :
        iTop = - 1    #-- Is DEW point
        PaDT = [dMin,sTop]
        PaBT = [dMin,0.0 ]
    else :
        PaDT = [dMin,0.0 ]
        PaBT = [dMin,sTop]

#== Maximum (Bottom) Depth ============================================    
        
    iBot = 0
    zTyp = "Bot"
    dBot = dMax - dRef

    pBot,zBot,qBot,sBot = \
        CG.calcTopBottom(zTyp,dBot,pRef,tRef,zRef,fRef,clsEOS,clsIO)

    if qBot :
        iBot =   1    #-- Is BUB point
        PaDB = [dMax,0.0 ]
        PaBB = [dMax,sBot]
    else :
        PaDB = [dMax,sBot]
        PaBB = [dMax,0.0 ]

#== Gas-Oil-Contact (GOC) Present? ====================================

    if iTop*iBot < 0 :

        qGOC = True

        dGOC,pGOC,sGOC,rGOC,zC1,zC7,zDew,Zbub = \
            CG.calcGOC(dMin,dMax,tRef,dRef,pRef,zRef,fRef,clsEOS,clsIO)

        PaDCa = [dGOC     ,sGOC]  #-- Just above Contact => DEW
        PaBCa = [dGOC-0.01,sGOC]

        PaDCb = [dGOC+0.01,sGOC]  #-- Just below Contact => BUB
        PaBCb = [dGOC     ,sGOC]

    else :

        qGOC = False
        dGOC = 0.0
        sGOC = 0.0

#== Add the Top Depth Values before incrementing depths ===============        

    PaBT = [dMin,sGOC]

    pDvD.append(PaDT)
    pBvD.append(PaBT)

#----------------------------------------------------------------------
#  Use Depth-Increment to step from Top to beyond Bottom
#----------------------------------------------------------------------

    dNxt = dMin + dInc
    zTyp = "Stp"
    qIns = True

    while dNxt < dMax :

        dDep = dNxt - dRef

        pDep,zDep,qDep,sDep = \
            CG.calcTopBottom(zTyp,dDep,pRef,tRef,zRef,fRef,clsEOS,clsIO)

#== Insert GOC? =======================================================

        if qIns and qGOC and dNxt > dGOC :

            pDvD.append(PaDCa)
            pDvD.append(PaDCb)

            pBvD.append(PaBCa)
            pBvD.append(PaBCb)

            qIns = False

#-- Add this depths results -----------------------------------------

        if qDep :    #-- Is BUB point
            PaDD = [dNxt,sGOC]
            PaBD = [dNxt,sDep]
        else :
            PaDD = [dNxt,sDep]
            PaBD = [dNxt,sGOC]

        pDvD.append(PaDD)
        pBvD.append(PaBD)

#-- Increment the Depth ---------------------------------------------

        dNxt = dNxt + dInc

#== Add the Bottom Depth Values =======================================

    PaDB = [dMax,sGOC]
        
    pDvD.append(PaDB)
    pBvD.append(PaBB)

#----------------------------------------------------------------------
#  Write to the E100 File
#----------------------------------------------------------------------

    outE100dep(dGOC,pDvD,pBvD,clsBLK,clsUNI,clsIO)

#== Return ============================================================

    return

#========================================================================
#  Write Eclipse 100 (E100) Output
#========================================================================

def outE100dep(dGOC,pDvD,pBvD,clsBLK,clsUNI,clsIO) :

    OutU = clsBLK.OutU  #-- Output Units, FLD (Field) or MET (Metric)

    fSim = clsIO.f100   #-- E100 Output File

    dRef = clsBLK.dRef
    pRef = clsBLK.pRef

    if OutU[:3] == "MET" :
        pUni = "bara"
        dUni = "m"
        sPun = "bar "
        sDun = " m  "
        dRef = clsUNI.I2X(dRef,"m")
        pRef = clsUNI.I2X(pRef,"bara")
        dGOC = clsUNI.I2X(dGOC,"m")
        sFor = "{:8.3f}"
        sDrf = sFor.format(dRef)
        sPrf = sFor.format(pRef)
        sDgo = sFor.format(dGOC)
    else :
        pUni = "psia"
        dUni = "ft"
        sPun = "psia"
        sDun = " ft "
        sFor = "{:8.2f}"
        sDrf = sFor.format(dRef)
        sPrf = sFor.format(pRef)
        sDgo = sFor.format(dGOC)

#== Header Information ================================================

    EOL = "\n"
    COM = "--\n"
    
    sRule = "--============================================================\n"
    sLine = "-- Following Keywords & Data Must Be Moved to SOLUTION section\n"

    sUnit = "--    " + sDun + "      " + sPun + "      " + sDun + \
            "      " + sPun + "      " + sDun + "      " + sPun + "\n"
#                                          12345678  12345678
    sEqui = "    " + sDrf + "  " + sPrf + "    x.xx      x.xx    " + \
                     sDgo + "    x.xx    1       1   /\n"

#-- EQUIL card ------------------------------------------------------
    
    fSim.write(EOL)
    fSim.write(COM)
    fSim.write(sRule)
    fSim.write(sLine)
    fSim.write(sRule)
    fSim.write(COM)
    fSim.write(EOL)
    fSim.write("-- Note.  You will have to supply numerical values for Dwoc, Pwoc & Pgoc\n")
    fSim.write(EOL)
    fSim.write("EQUIL\n")
    fSim.write("--    Dref      Pref      Dwoc      Pwoc      Dgoc      Pgoc    PBVD    PDVD\n")
    fSim.write(sUnit)
    fSim.write(sEqui)
    fSim.write(EOL)

#-- PDVD & PBVD -----------------------------------------------------

    fSim.write("PDVD\n")
    fSim.write("--   TVDSS      Pdew\n")
    fSim.write("--    " + sDun + "      " + sPun + "\n")

    for pD in pDvD :
        dDep = clsUNI.I2X(pD[0],dUni)
        pDep = clsUNI.I2X(pD[1],pUni)
        sDep = sFor.format(dDep)
        sDew = sFor.format(pDep)
        fSim.write("    " + sDep + "  " + sDew + "\n")
    fSim.write("/\n")
    fSim.write(EOL)

    fSim.write("PBVD\n")
    fSim.write("--   TVDSS      Pbub\n")
    fSim.write("--    " + sDun + "      " + sPun + "\n")

    for pB in pBvD :
        dDep = clsUNI.I2X(pB[0],dUni)
        pDep = clsUNI.I2X(pB[1],pUni)
        sDep = sFor.format(dDep)
        sBub = sFor.format(pDep)
        fSim.write("    " + sDep + "  " + sBub + "\n")
    fSim.write("/\n")
    fSim.write(EOL)

#== Footer Information ================================================

    sLine = "-- Preceding Keywords & Data Must Be Moved to SOLUTION section\n"

    fSim.write(COM)
    fSim.write(sRule)
    fSim.write(sLine)
    fSim.write(sRule)
    fSim.write(COM)
    fSim.write(EOL)

#== No return value ===================================================

    return

#========================================================================
#  End of Module
#========================================================================

