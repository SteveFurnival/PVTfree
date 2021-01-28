#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Experiments Data Module
#  
#========================================================================

#!/usr/bin/python3

import numpy    as NP

import allData  as AD
import allProp  as AP
import calcExps as CX
import readGen  as RG

#========================================================================
#  Read Experiments Module
#  Read lines until we get to ENDEXP (or too many lines!)
#========================================================================

def readExps(clsIO,clsEOS,dicSAM,dicEXP,clsUNI) :

    iERR  = 0
    iLine = 0

    fInp = clsIO.fInp
    
#-- Experiment Names & Mnemonics ------------------------------------

    sNam = {'CCE'  :'CCE','CVD' :'CVD','DLE'  :'DLE','SEP' :'SEP',
            'FLASH':'FLS','PSAT':'SAT','SWELL':'SWL','GRAD':'GRD'}

#== Inactivate existing experiments ===================================

    for iExp in range(len(dicEXP)) : dicEXP[iExp].IsAct = False

#== Main Read Loop ====================================================    

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of EXP Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0                     : pass   #-- Blank line!
        elif tokS[0][:2]         == "--"   : pass   #-- Comment
        elif tokS[0][:3].upper() == "DEB"  : iERR = RG.readDebug(clsIO)
        elif tokS[0][:3].upper() == "OPT"  : iERR = RG.readOption(clsIO)

        elif tokS[0].upper() in sNam       :
            sHrt = sNam[tokS[0].upper()]
            iERR = readGen(sHrt,clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "PLOT" :
            if tokS[1].upper() == "NONE" :
                clsIO.Pall = False
        elif tokS[0][:4].upper() == "ENDE" :
            #print("ENDEXP k/w found")
            break

        if iERR < 0 : break

#======================================================================
#  Experiments Read OK => Run Them
#======================================================================

    if iERR == 0 : CX.calcExps("",clsEOS,dicSAM,dicEXP,clsUNI,clsIO)

    return iERR

#========================================================================
#  Read an Experiment
#========================================================================

def readGen(sTyp,clsIO,clsEOS,dicSAM,dicEXP,clsUNI) :

    nCom = clsEOS.nComp

    clsEXP = AD.classEXP(sTyp)  #-- Create a New Exp Class
    nX     = len(dicEXP)        #-- Current Size of Exp Dictionary

    #print("readGen: sTyp,nX ",sTyp,nX)

    dicEXP.update({nX:clsEXP})  #-- Add New Exp to the Dictionary

    iERR  = 0
    iLine = 0

    sName = None    #-- Sample Name for Exp
    sNInj = None    #-- Injection Sample Name (SWL only)
    Tres  = None    #-- Reservoir Temp
    sLCCE = None    #-- SL-definition for CCE
    
    colH  = []      #-- Column Headers
    colU  = []      #-- Column Units
    dTab  = []      #-- Data Table

    qData = False

#-- Independent and Observed Data for this Experiment ---------------

    sIndA = AP.classLIB().INDshrt.get(sTyp)
    sObsA = AP.classLIB().OBSshrt.get(sTyp)
    sCalA = AP.classLIB().CALshrt.get(sTyp)
    
    nIndA = len(sIndA)  #-- Num Allowed Independent Vars for this Exp
    nObsA = len(sObsA)  #-- Num Allowed Observed    Vars for this Exp
    nCalA = len(sCalA)  #-- Num Allowed Calculated  Vars for this Exp

    #print("sTyp,nIndA,nObsA,nCalA ",sTyp,nIndA,nObsA,nCalA)

    dWal = NP.ones(nObsA)      ; dWps = 1.0
    qPlt = NP.full(nObsA,True)

#== Read lines until blank found ======================================

    fInp = clsIO.fInp

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of Data for Experiment ",sTyp)
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if nTok == 0 : break
        else :

#-- Must be reading Column Headers ----------------------------------
            
            if   tokS[0][:2].upper() == "PR" or \
                 tokS[0][:2].upper() == "TE" or \
                 tokS[0][:2].upper() == "MO" or \
                 tokS[0][:2].upper() == "HE" :

                for iTok in range(nTok) : colH.append(tokS[iTok])                

                nHed = nTok
                nexL = next(fInp)               #-- Next line must define units

                tokS = nexL.split()
                nUni = len(tokS)

                if nHed != nUni :
                    print("Must have same number of Headers and Units for experiment ",sTyp," - Error")
                    return -1

                for iTok in range(nUni) : colU.append(tokS[iTok])

                qData = True
                
#-- Comment? --------------------------------------------------------                

            elif tokS[0][:2] == "--" : pass

#-- Weighting Factors -----------------------------------------------

            elif tokS[0][:4].upper() == "WEIG" :

                if sTyp == "CCE" or sTyp == "CVD" or \
                   sTyp == "DLE" or sTyp == "SEP" : qTyp = True
                else                              : qTyp = False

                iTok = 1                #-- Ignore 1st Token [WEIG]
                while iTok < nTok :  
                    sTok = tokS[iTok].upper()
                    if sTok == "PSAT" and qTyp :
                        iTok += 1
                        dTok  = tokS[iTok].upper()
                        dWps  = float(dTok)
                        iTok += 1
                    else :
                        for iObs in range(nObsA) :
                            sObs = sObsA[iObs]
                            if sTok == sObs :
                                iTok += 1
                                dTok  = tokS[iTok].upper()
                                dWal[iObs] = float(dTok)
                        iTok += 1

#-- Plot Oberseved Quantity -----------------------------------------

            elif tokS[0][:4].upper() == "PLOT" :

                iTok = 1                #-- Ignore 1st Token [PLOT]
                while iTok < nTok :  
                    sTok = tokS[iTok].upper()
                    if sTok == "NO" : qPlt[iTok-1] = False
                    else            : qPlt[iTok-1] = True
                    iTok += 1

#-- Must now be read the numeric data -------------------------------

            elif qData :

                dRow = NP.zeros(nHed)              #-- Store data on this row

                if   tokS[0][:1].upper() == "P" :  #-- Found PSAT
                    nLow =  1
                    isPD =  1
                elif tokS[0][:1].upper() == "D" or \
                     tokS[0][:1].upper() == "H"    :  #-- Found DREF (GRD-Only)
                    nLow =  1                         #--    or HREF
                    isPD = -1
                else :
                    nLow =  0
                    isPD =  0

                for iTok in range(nLow,nTok) :
                    dRow[iTok-nLow] = float(tokS[iTok])

                dTab.append(dRow)                  #-- Append Row to Table

                if   isPD ==  1 :
                    nRsat = len(dTab)
                    clsEXP.setPsatRow(nRsat-1)
                elif isPD == -1 :
                    nDref = len(dTab)
                    clsEXP.setDrefRow(nDref-1)

#-- Line following EXP type contains non-Table information ----------                

            else :  

                iTok = 0
                while iTok < nTok :
                    if   tokS[iTok].upper() == "SAMP" :
                        iTok += 1
                        sName = tokS[iTok]
                        #print("sName ",sName)
                    elif tokS[iTok].upper() == "TRES" :
                        iTok += 1
                        Tres = float(tokS[iTok])
                        iTok += 1
                        Tuni =       tokS[iTok]
                        #print("Tres,Tuni ",Tres,Tuni)
                    elif tokS[iTok].upper() == "SINJ" :
                        iTok += 1
                        sNInj = tokS[iTok]
                    elif tokS[iTok].upper() == "SLIQ" :
                        iTok += 1
                        sLCCE = tokS[iTok].upper()
                    iTok += 1
        
#========================================================================
#  Process the read data
#========================================================================

#-- Fluid Sample ----------------------------------------------------

    if sName == None :
        print("No sample defined (using SAMP=) for Experiment ",sTyp," - Error")
        return -1
    else :
        nSamp = len(dicSAM)
        qFoun = False
        for iS in range(nSamp) :
            if sName == dicSAM[iS].sNam :
                qFoun = True
                clsEXP.setSamp(iS,sName)
                break
        if not qFoun :
            print("Sample ",sName," specified for Experiment ",sTyp," not found - Error")
            return -1

#-- Injection Sample (SWL only) -------------------------------------

    if sTyp == 'SWL' :
        if sNInj == None :
            print("No Injection Sample defined (using SINJ) for Experiment SWL - Error")
            return -1
        else:
            nSamp = len(dicSAM)
            qFoun = False
            for iS in range(nSamp) :
                if sNInj == dicSAM[iS].sNam :
                    qFoun = True
                    clsEXP.setSinj(iS)
                    clsEXP.setSNin(sNInj)
                    break
            if not qFoun :
                print("Injection Sample ",sNInj," specified for Experiment SWL not found - Error")
                return -1
    else :
        if sNInj != None :
            print("Experiment ",sTyp," does not support SINJ= parameter [SWL only] - Error")
            return -1

#-- Reservoir Temperature -------------------------------------------

    if sTyp == 'CCE' or sTyp == 'CVD' or sTyp == 'DLE' or \
       sTyp == 'SWL' or sTyp == 'GRD' :
        if Tres == None :
            print("No Reservoir Temperature defined (using TRES=) for Experiment ",sTyp," - Error")
            return -1
        else :
            Tint = clsUNI.X2I(Tres,Tuni)
            if Tint == -999.9 : return -1
            else              : clsEXP.setTres(Tint,Tuni)

#-- CCE Liquid Saturation Definition (SL = VL/VTot or SL = VL/Vdew) -

    if sTyp == 'CCE' :
        if sLCCE == None :
            print("CCE Experiment has not set SLIQ-Type: SLiq = VLiq/Vdew assumed")
            clsEXP.setSLCCE("DEW")
        else:
            if   sLCCE == "DEW" or sLCCE == "TOT" :
                clsEXP.setSLCCE(sLCCE)
            else :
                print("SLIQ= argument of CCE Experiment = ",sLCCE,": Should be DEW or TOT - Error")
                return -1
    else :
        if sLCCE != None :
            print("Experiment ",sTyp," does not support SLIQ= parameter [CCE only] - Error")
            return -1
                
#== Process the Column Headers and Associated Units ===================

    nHead = len(colH)
    nUnit = len(colU)

    #print("colH ",colH)
    #print("colU ",colU)

    if nHead != nUnit :
        print("Number of Column Headers not equal to Column Units for Experiment ",sTyp," - Error")
        return -1

#-- Number of Rows of Data and Initialise Structures in classEXP ----

    nRowS = len(dTab)
    nObsL = nHead - nIndA

    if sTyp == "GRD" : nRowD = nRowS + 1
    else             : nRowD = nRowS

    clsEXP.createExpArrays(nIndA,nObsA,nCalA,nRowD)

    if sTyp == "GRD" : clsEXP.createCompArray(nRowD,nCom,"Z","C")   #-- zCal-Array

#-- SEP experiment; create Lsep & Vsep Arrays -----------------------    

    if sTyp == "SEP" :                  #-- Create Lsep/Vsep arrays
        clsEXP.createSepStages(nRowS)
        qLsep = False ; qVsep = False
        for i in range(nHead) :         #-- Test if they have been read
            if colH[i].upper() == "LSEP" : qLsep = True
            if colH[i].upper() == "VSEP" : qVsep = True
        if       qLsep and     qVsep :
            qSep = True
        elif     qLsep and not qVsep :
            print("Lsep provided but not Vsep - Error")
            return -1
        elif not qLsep and     qVsep :
            print("Vsep provided but not Lsep - Error")
            return -1
        else :
            qSep = False

#-- Calculated Data [Store Full Array Regardless of User Entry] -----

    nTemp = len(clsEXP.hCal)

    for i in range(nCalA) :
        clsEXP.hCal[i] =            sCalA[i]
        clsEXP.uCal[i] = clsUNI.DUN(sCalA[i])

#-- Independent Data ------------------------------------------------

    nFoun = 0
    for i in range(nIndA) :
        sInd = sIndA[i].upper()
        for j in range(nHead) :
            sCol = colH[j].upper()
            sUni = colU[j].upper()
            if sInd == sCol :
                nFoun += 1
                clsEXP.hInd[i] = sCol
                clsEXP.uInd[i] = sUni
                #print("i,sInd,j,sCol,sUni ",i,sInd,j,sCol,sUni)
                for k in range(nRowS) :
                    clsEXP.dInd[i][k] = clsUNI.X2I(dTab[k][j],sUni)
                    #print("k,dTab,cInd ",k,dTab[k][j],clsEXP.dInd[i][k])

    if nIndA != nFoun :
        print("Experiment ",sTyp," expects ",nIndA," columns of Independent data: Only ",nFoun," found - Error")
        return -1

#-- Independent/Observed Psat for selected Experiments --------------

    if sTyp == "CCE" or sTyp == "CVD" or sTyp == "DLE" or sTyp == "SEP" :
        nRsat = clsEXP.nRsat
        if nRsat < 0 :
            print("Experiment ",sTyp," must have one row starting with PSAT - Error")
            return -1
        else :
            PsatO = clsEXP.dInd[0][nRsat]
            PsatU = clsEXP.uInd[0]
            clsEXP.setPsatObs(PsatO)
            clsEXP.setPsatWei(dWps)
            clsEXP.setPsatUni(PsatU)

#-- Reference Depth and Pressure for GRD ----------------------------

    if sTyp == "GRD" :
        nDref = clsEXP.nDref
        if nDref < 0 :
            print("Experiment GRD must have one row starting with DREF - Error")
            return -1
        else :
            Dref = clsEXP.dInd[0][nDref]
            Pref = None
            for i in range(nHead) :
                sObs = colH[i].upper()
                sUni = colU[i].upper()
                if sObs[:2] == "PR" :
                    Pref = clsUNI.X2I(dTab[nDref][i],sUni)
                    break
            if Pref == None :
                print("Experiment GRD must a Pressure define on DREF row - Error")
                return -1

        clsEXP.setDref(Dref)
        clsEXP.setPref(Pref)
        
#-- Observed Data ---------------------------------------------------

    qObsA = [False for i in range(nObsA)]

    nFoun = 0
    iPlt  = 0
    for i in range(nHead) :
        sObs = colH[i].upper()
        sUni = colU[i].upper()
        for j in range(nObsA) :
            sCol = sObsA[j].upper()
            if sObs == sCol :
                clsEXP.hObs[j] = sCol
                clsEXP.uObs[j] = sUni
                clsEXP.uCal[j] = sUni  #-- Over-write default Cal-Unit with Obs-Unit
                clsEXP.qObs[j] = True
                clsEXP.qPlt[j] = qPlt[iPlt]
                iPlt         += 1
                for k  in range(nRowS) :
                    clsEXP.dObs[j][k] = clsUNI.X2I(dTab[k][i],sUni)
                    clsEXP.dWei[j][k] = dWal[j]
                nFoun += 1
        if sTyp == "SEP" and sObs == "LSEP" :
            for k in range(nRowS) : clsEXP.Lsep[k] = int(dTab[k][i])
        if sTyp == "SEP" and sObs == "VSEP" :
            for k in range(nRowS) : clsEXP.Vsep[k] = int(dTab[k][i])

    clsEXP.setUserObs(nFoun)

#-- If SEP experiment, default the Lsep/Vsep arrays if not set ------

    if sTyp == "SEP" and not qSep :
        for k in range(nRowS) :
            if   k == 0 :
                clsEXP.Lsep[k] =   1 ; clsEXP.Vsep[k] = 1
            elif k == nRowS - 1 :
                clsEXP.Lsep[k] =   0 ; clsEXP.Vsep[k] = 0
            else :
                clsEXP.Lsep[k] = k+1 ; clsEXP.Vsep[k] = 0

#== End of Routine ======================================================

    if iERR == 0 : clsEXP.isDef = True      #-- Exp is Defined!!

    return iERR
    
#========================================================================
#  End of Module
#========================================================================
