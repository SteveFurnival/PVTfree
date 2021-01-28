
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Write Results from Regression Module
#  
#========================================================================

#!/usr/bin/python3

import numpy as NP

#========================================================================
#  Output Regression Information
#========================================================================

def outRegInfo(nReg,dSSQ,regO,regI,clsIO) :

    fReg = clsIO.fReg   #-- Regression Output File

#-- Header ----------------------------------------------------------

    ssqS = 50.0/dSSQ        #-- 0.5*100/SSQ to get normalised %Residual
    sLst = ""

#            12iReg12EXPx12OBSV12 Stage 12Weight12 Observed 12Calculated12 Residual 12 FracSSQ
    sLine = "  ----  ----  ----  -------  ------  ----------  ----------  ----------  ----------"
    sHead = "  iReg  EXP   Quan   Stage   Weight   Observed   Calculated   Residual    FracSSQ  "

    fReg.write(sLine+"\n")
    fReg.write(sHead+"\n")

#== Loop Over Elements (Residuals!) in the SSQ ========================    

    for iReg in range(nReg) :
            
        sReg  = "{:4d}".format(iReg+1)
            
        sExp  = regO[iReg][0]
        sQua  = regO[iReg][1].ljust(4)
        sStg  = regO[iReg][2]
            
        dWei  = regO[iReg][3]
        dObs  = regO[iReg][4]
        dCal  = regO[iReg][5]
        dRes  = regI[iReg]
        fSSQ  = ssqS*dRes*dRes

        sWei  =  "{:6.2f}".format(dWei)
        sObs  = "{:10.4f}".format(dObs)
        sCal  = "{:10.4f}".format(dCal)
        sRes  =  "{:8.4f}".format(dRes)
        sSSQ  = "{:10.4f}".format(fSSQ)

        if sExp != sLst : fReg.write(sLine+"\n")
        sLst = sExp
            
        sFReg = "  " + sReg + "  " + sExp + "  " + sQua + "  " + sStg + \
                "  " + sWei + "  " + sObs + "  " + sCal + "  " + sRes + "  " + sSSQ + "\n"

        fReg.write(sFReg)

    fReg.write(sLine+"\n")

#== No Return Value ===================================================

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

        sExt = regVarInfo(sTyp,clsREG,clsEOS,dicSAM)

        sLine = "  " + sVar + "  " + sNam + "  " + sMin + "  " + sVec + \
                "  " + sMax + "  " + sDlx + "  " + sBou + "  " + sExt + "\n"
        fReg.write(sLine)

    fReg.write(sHead)

    return

#========================================================================
#  Diagnostics on the Hessian Matrix for Output Purposes
#========================================================================

def diagnoseHess(nIter,detM,iOrd,eVal,eVec,dicREG,dicSAM,clsEOS,clsIO) :

    eTst = 1.0E-08      #-- Might want to change this value?
                        #-- Is consistent with value in calcYvec!

    fReg = clsIO.fReg
    nVar = len(eVal)

    sHead = "---------------------------------------\n"
    sLine = "  Regression Iteration (Diagnosis): " + str(nIter) + "\n"

    fReg.write("\n")
    fReg.write(sHead)
    fReg.write(sLine)
    fReg.write(sHead)

#-- Determinant Output ----------------------------------------------    

    sLabl = "  Determinant of the Hessian Matrix = {:12.5e}\n".format(detM)

    fReg.write("\n")
    fReg.write(sLabl)
    fReg.write("\n")

#== Eigenvalues =======================================================

    eMax = 0.0
    for iVar in range(nVar) : eMax = max(eMax,eVal[iVar])

#            12IV12vType12  xMin  12  xVec  12  xMax  12  delX  12  Bounded?  12Comps/Samps
    sLine = "  IV  vType  EigenValue  Remove?  Comps/Samp\n"
    sHead = "  --  -----  ----------  -------  ----------\n"
    fReg.write(sHead)
    fReg.write(sLine)
    fReg.write(sHead)

    for jVar in range(nVar) :

        iVar = iOrd[jVar]

        eRat = eVal[iVar]/eMax

        if eRat < eTst : sRem = "Yes".center(7)
        else           : sRem =  "No".center(7)

        clsREG = dicREG[iVar]

        sNam = clsREG.vNam.ljust(5)
        sTyp = clsREG.vTyp

        sNum =    "{:2d}".format(iVar+1)
        sVal = "{:10.3e}".format(eVal[iVar])

        sExt = regVarInfo(sTyp,clsREG,clsEOS,dicSAM)

        sLine = "  " + sNum + "  " + sNam + "  " + sVal + "  " + \
                       sRem + "  " + sExt + "\n"
        fReg.write(sLine)

    fReg.write("\n")

#== Eigenvectors =====================================================

    fReg.write("  Eigenvector Dependencies?\n")
    fReg.write("\n")

    sHead = "  ---  "
    sLine = "  VAR  "
    
    for iC in range(nVar-1) :
        sVar  = "{:3d}".format(iC+1)
        sHead = sHead + "---  "
        sLine = sLine + sVar + "  "

    sHead = sHead + "\n"
    sLine = sLine + "\n"

    fReg.write(sHead)
    fReg.write(sLine)
    fReg.write(sHead)
    fReg.write("\n")

    for iC in range(1,nVar) :

        iVar = iOrd[iC]
        
        sVar  = "{:3d}".format(iVar+1)
        sLine = "  " + sVar + "  "
        
        for jC in range(iC) :

            jVar = iOrd[jC]
            
            dotP = NP.inner(eVec[:,iVar],eVec[:,jVar])
            dTst = abs(dotP - 1.0)
            
            if dTst < 1.0E-05 : sDI = "Yes"
            else              : sDI = "No "
            
            sLine = sLine + sDI + "  "
            
        sLine = sLine + "\n"
        fReg.write(sLine)

    fReg.write("\n")

#========================================================================
#  End of Routine
#========================================================================

    return

#========================================================================
#  Information on a Regression Variable (for Output)
#========================================================================

def regVarInfo(sTyp,clsREG,clsEOS,dicSAM) :
    
    sExt = ""

    if   sTyp == "REG" or sTyp == "LBC" :   #-- Regular Variable
        iGrp = clsREG.vGrp
        for jGrp in range(len(iGrp)) :
            iOne = iGrp[jGrp]
            sThs = clsEOS.gNM(iOne)
            sExt = sExt + sThs + ", "
        sExt = sExt[:-2]
    elif sTyp == "KIJ" :                    #-- Binary
        iOne = clsREG.vOne
        sOne = clsEOS.gNM(iOne)
        sExt = sExt + sOne + "->"
        iGrp = clsREG.vGrp
        for jGrp in range(len(iGrp)) :
            iOne = iGrp[jGrp]
            sThs = clsEOS.gNM(iOne)
            sExt = sExt + sThs + ", "
        sExt = sExt[:-2]
    elif sTyp == "PLS" :                    #-- Plus Fraction Variable
        iSam = clsREG.vSam
        sSam = dicSAM[iSam].sNam
        sExt = sExt + sSam

#========================================================================
#  End of Routine
#========================================================================

    return sExt

#========================================================================
#  Write Header Before Initial Regression Residuals
#========================================================================

def writeInitResHead(clsIO) :

    fReg = clsIO.fReg   #-- Regression Output File

#-- Header Information ----------------------------------------------    

    fReg.write("\n")
    fReg.write("============================================================\n")
    fReg.write("  Regression Output:  Initial Residuals\n")
    fReg.write("============================================================\n")
    fReg.write("\n")

#== No Return Value ===================================================

    return

#========================================================================
#  Write Header Before Final Regression Residuals
#========================================================================

def writeFinalResHead(clsIO) :

    fReg = clsIO.fReg   #-- Regression Output File

#-- Header Information ----------------------------------------------    

    fReg.write("\n")
    fReg.write("============================================================\n")
    fReg.write("  Regression Output:  Final Residuals\n")
    fReg.write("============================================================\n")
    fReg.write("\n")

#== No Return Value ===================================================

    return

#========================================================================
#  Write a Vector
#========================================================================

def writeVector(sVec,vecT) :

    print(sVec," ",end="")

    for eVal in vecT :
        print(" {:10.3e}".format(eVal),end="")
    print(end="\n")

#== Return ============================================================        
        
    return

#========================================================================
#  Write a Matrix
#========================================================================

def writeMatrix(sMat,matX) :

    print(sMat)

    for dVec in matX :
        for dVal in dVec :
            print(" {:10.3e}".format(dVal),end="")
        print(end="\n")

#== Return ============================================================        
        
    return

#========================================================================
#  End of Module
#========================================================================

