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

import calcExps as CX
import readGen  as RG

expDEPvars = { 'CCE':['PRES'],
               'CVD':['PRES'],
               'DLE':['PRES'],
               'SEP':['PRES','TEMP'],
               'FLS':['PRES','TEMP'],
               'SAT':['TEMP'],
               'SWL':['MOLE'],
               'GRD':['HEIG'] }

expOBSvars = { 'CCE':['RELV','SLIQ','ZFAC','DENO','VISG','VISO'],
               'CVD':['MREM','SLIQ','ZFAC'],
               'DLE':['BO'  ,'GOR' ,'DENO','BT'  ,'BG'  ,'ZFAC','GGRV','VISO'],
               'SEP':['BO'  ,'GOR' ,'DENO','GGRV'],
               'FLS':['ZFAC','DENO','VFRC'],
               'SAT':['PSAT','ZFAC','DENO'],
               'SWL':['PSAT','VSWL'],
               'GRD':['PRES','PSAT','DENS','ZC1' ,'ZC7+'] }

expOBSlong = { 'CCE':['Relative Volume',
                      'Liquid Dropout',
                      'Gas Z-Factor',
                      'Oil Density',
                      'Gas Viscosity',
                      'Oil Viscosity'],
               'CVD':['Moles Removed',
                      'Liquid Dropout',
                      'Gas Z-Factor'],
               'DLE':['Oil Formation Volume Factor'  ,
                      'Gas-Oil-Ratio' ,
                      'Oil Density',
                      'Total Formation Volume Factor'  ,
                      'Gas Formation Volume Factor'  ,
                      'Gas Z-Factor',
                      'Gas Gravity',
                      'Oil Viscosity'],
               'SEP':['Oil Formation Volume Factor'  ,
                      'Gas-Oil-Ratio' ,
                      'Oil Density',
                      'Gas Gravity'],
               'FLS':['Gas Z-Factor',
                      'Oil Density',
                      'Vapour Fraction'],
               'SAT':['Saturation Pressure',
                      'Gas Z-Factor',
                      'Oil Density'],
               'SWL':['Saturation Pressure',
                      'Swelling Factor'],
               'GRD':['Pressure',
                      'Saturation Pressure',
                      'Density',
                      'C1 Mole Fraction' ,
                      'C7+ Mole Fraction'] }

class classEXP :

#-- Creator ---------------------------------------------------------    

    def __init__(self,Name) :
        
        self.xName = Name
        self.nRsat = -1
        self.nDref = -1
        self.PsatO = -1.0
        self.IsAct = True

#-- Set and Get the Reservoir Temperature ---------------------------
        
    def setTres(self,Tres) : self.Tres = Tres

    def setTuni(self,Tuni) : self.Tuni = Tuni

    def getTres()          : return self.Tres

#-- Get and set the Sample Number -----------------------------------
    
    def setSamp(self,nSamp) : self.nSamp = nSamp

    def getSamp()           : return self.nSamp

#-- Get and set the Injection Sample Number (SWL only) --------------
    
    def setSinj(self,nSinj) : self.nSinj = nSinj

    def getSinj()           : return self.nSinj

#-- Set and Get the CCE Liquid Saturation Method (DEW or TOT) -------

    def setSLCCE(self,sLCCE) : self.sLCCE = sLCCE

    def getSLCCE()           : return self.sLCCE

#-- Observed/Calculated Psat [for relevant experiments] -------------

    def setPsatRow(self,nRsat) : self.nRsat = nRsat

    def setPsatObs(self,PsatO) : self.PsatO = PsatO

    def setPsatCal(self,PsatC) : self.PsatC = PsatC

    def setPsatWei(self,PsatW) : self.PsatW = PsatW

    def setPsatUni(self,PsatU) : self.PsatU = PsatU

#-- Set Dref Row for GRD Experiment ---------------------------------

    def setDrefRow(self,nDref) : self.nDref = nDref

    def setDref(self,Dref)     : self.Dref  = Dref
    def setPref(self,Pref)     : self.Pref  = Pref

#-- Create Arrays to Store Data -------------------------------------

    def createArrays(self,nDepA,nObsA,nRowS) :

        self.nDep = nDepA
        self.nObs = nObsA
        self.nRow = nRowS

        self.hDep = ["" for i in range(nDepA)]
        self.uDep = ["" for i in range(nDepA)]
        self.dDep = NP.zeros((nRowS,nDepA))

        self.hObs = ["" for i in range(nObsA)]
        self.uObs = ["" for i in range(nObsA)]
        self.qObs = NP.full(nObsA,False)
        self.qCal = NP.full(nObsA,True )
        self.qPlt = NP.full(nObsA,True )
        self.dObs = NP.zeros((nRowS,nObsA))

        self.dWei = NP.zeros((nRowS,nObsA))

        self.hCal = ["" for i in range(nObsA)]
        self.uCal = ["" for i in range(nObsA)]
        self.dCal = NP.zeros((nRowS,nObsA))

    def setUserObs(self,nObU) : self.nObU = nObU

    def createSepStages(self,nRowS) :

        self.Lsep = NP.zeros(nRowS,dtype=int)
        self.Vsep = NP.zeros(nRowS,dtype=int)

#========================================================================
#  Read Experiments Module
#  Read lines until we get to ENDEXP (or too many lines!)
#========================================================================

def readExps(clsIO,clsEOS,dicSAM,dicEXP,clsUNI) :

    iERR  = 0
    iLine = 0

    fInP = clsIO.fInP

#== Inactivate existing experiments ===================================

    nExp = len(dicEXP)
    for iExp in range(nExp) : dicEXP[iExp].IsAct = False

#== Main Read Loop ====================================================    

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of EXP Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0                    :  #-- Blank line!
            pass
        elif tokS[0][:2]         == "--"  :  #-- Comment
            pass
        elif tokS[0][:3].upper() == "DEB" :
            iERR = RG.readDebug(clsIO)
        elif tokS[0].upper()     == "CCE" :
            iERR = readGen("CCE",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "CVD" :
            iERR = readGen("CVD",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "DLE" :
            iERR = readGen("DLE",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "SEP" :
            iERR = readGen("SEP",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "FLASH" :
            iERR = readGen("FLS",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "PSAT" :
            iERR = readGen("SAT",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "SWELL" :
            iERR = readGen("SWL",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
        elif tokS[0].upper()     == "GRAD" :
            iERR = readGen("GRD",clsIO,clsEOS,dicSAM,dicEXP,clsUNI)
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

    if iERR == 0 :
        sExt = ""
        CX.calcExps(sExt,clsEOS,dicSAM,dicEXP,clsIO,clsUNI)

    return iERR,dicEXP

#========================================================================
#  Read an Experiment
#========================================================================

def readGen(sTyp,clsIO,clsEOS,dicSAM,dicEXP,clsUNI) :

    clsEXP = classEXP(sTyp)
    nX = len(dicEXP)

    #print("readGen: sTyp,nX ",sTyp,nX)

    dicEXP.update({nX:clsEXP})

    iERR  = 0
    iLine = 0

    sName = None
    sNInj = None
    Tres  = None
    sLCCE = None
    
    colH  = []
    colU  = []

    dTab  = []

    qData = False

#-- Dependent and Observed Data for this Experiment -----------------    

    sDepA = expDEPvars.get(sTyp)
    sObsA = expOBSvars.get(sTyp)
    
    nDepA = len(sDepA)
    nObsA = len(sObsA)

    dWal = [ 1.0 for i in range(nObsA)] ; dWps = 1.0
    qPlt = [True for i in range(nObsA)]

#== Read lines until blank found ======================================

    fInP = clsIO.fInP

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of Data for Experiment ",sTyp)
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if nTok == 0 :
            break
        else :

#-- Must be reading Column Headers ----------------------------------
            
            if   tokS[0][:2].upper() == "PR" or \
                 tokS[0][:2].upper() == "TE" or \
                 tokS[0][:2].upper() == "MO" or \
                 tokS[0][:2].upper() == "HE" :

                for iTok in range(nTok) :
                    colH.append(tokS[iTok])                

                nHed = nTok
                nexL = next(fInP)               #-- Next line must define units

                tokS = nexL.split()
                nUni = len(tokS)

                if nHed != nUni :
                    print("Must have same number of Headers and Units for experiment ",sTyp," - Error")
                    iERR = -1
                    return iERR

                for iTok in range(nUni) :
                    colU.append(tokS[iTok])

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
                    if sTok == "NO" :
                        qPlt[iTok-1] = False
                    else :
                        qPlt[iTok-1] = True
                    iTok += 1

#-- Must now be read the numeric data -------------------------------

            elif qData :

                dRow = [0.0 for i in range(nHed)]  #-- Store data on this row

                if   tokS[0][:1].upper() == "P" :  #-- Found PSAT
                    nLow =  1
                    isPD =  1
                elif tokS[0][:1].upper() == "D" :  #-- Found DREF (GRD-Only)
                    nLow =  1
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
        iERR = -1
        return iERR
    else :
        nSamp = clsEOS.NS
        qFoun = False
        for iS in range(nSamp) :
            if sName == dicSAM[iS].sName :
                qFoun = True
                clsEXP.setSamp(iS)
                break
        if not qFoun :
            print("Sample ",sName," specified for Experiment ",sTyp," not found - Error")
            iERR = -1
            return iERR

#-- Injection Sample (SWL only) -------------------------------------

    if sTyp == 'SWL' :
        if sNInj == None :
            print("No Injection Sample defined (using SINJ) for Experiment SWL - Error")
            iERR = -1
            return iERR
        else:
            nSamp = clsEOS.NS
            qFoun = False
            for iS in range(nSamp) :
                if sNInj == dicSAM[iS].sName :
                    qFoun = True
                    clsEXP.setSinj(iS)
                    break
            if not qFoun :
                print("Injection Sample ",sNInj," specified for Experiment SWL not found - Error")
                iERR = -1
                return iERR
    else :
        if sNInj != None :
            print("Experiment ",sTyp," does not support SINJ= parameter [SWL only] - Error")
            iERR = -1
            return iERR

#-- Reservoir Temperature -------------------------------------------

    if sTyp == 'CCE' or sTyp == 'CVD' or sTyp == 'DLE' or \
       sTyp == 'SWL' or sTyp == 'GRD' :
        if Tres == None :
            print("No Reservoir Temperature defined (using TRES=) for Experiment ",sTyp," - Error")
            iERR = -1
            return iERR
        else :
            Tint = clsUNI.X2I(Tres,Tuni)
            if Tint == -999.9 :
                iERR = -1
                return iERR
            else :
                clsEXP.setTres(Tint)
                clsEXP.setTuni(Tuni)

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
                iERR = -1
                return iERR
    else :
        if sLCCE != None :
            print("Experiment ",sTyp," does not support SLIQ= parameter [CCE only] - Error")
            iERR = -1
            return iERR
                
#== Process the Column Headers and Associated Units ===================

    nHead = len(colH)
    nUnit = len(colU)

    #print("colH ",colH)
    #print("colU ",colU)

    if nHead != nUnit :
        print("Number of Column Headers not equal to Column Units for Experiment ",sTyp," - Error")
        iERR = -1
        return iERR

#-- Number of Rows of Data and Initialise Structures in classEXP ----

    nRowS = len(dTab)
    nObsL = nHead - nDepA

    clsEXP.createArrays(nDepA,nObsA,nRowS)

#-- SEP experiment; create Lsep & Vsep Arrays -----------------------    

    if sTyp == "SEP" :                  #-- Create Lsep/Vsep arrays
        clsEXP.createSepStages(nRowS)
        qLsep = False ; qVsep = False
        for i in range(nHead) :         #-- Test if they have been read
            if colH[i].upper() == "LSEP" : qLsep = True
            if colH[i].upper() == "VSEP" : qVsep = True
        #print("qLsep,qVsep ",qLsep,qVsep)
        if       qLsep and     qVsep :
            qSep = True
        elif     qLsep and not qVsep :
            print("Lsep provided but not Vsep - Error")
            iERR = -1
            return iERR
        elif not qLsep and     qVsep :
            print("Vsep provided but not Lsep - Error")
            iERR = -1
            return iERR
        else :
            qSep = False

#-- Calculated Data [Store Full Array Regardless of User Entry] -----

    for i in range(nObsA) :
        clsEXP.hCal[i] =            sObsA[i]
        clsEXP.uCal[i] = clsUNI.DUN(sObsA[i])

#-- Dependent Data --------------------------------------------------

    nFoun = 0
    for i in range(nDepA) :
        sDep = sDepA[i].upper()
        for j in range(nHead) :
            sCol = colH[j].upper()
            sUni = colU[j].upper()
            if sDep == sCol :
                nFoun += 1
                clsEXP.hDep[i] = sCol
                clsEXP.uDep[i] = sUni
                #print("i,sDep,j,sCol,sUni ",i,sDep,j,sCol,sUni)
                for k in range(nRowS) :
                    clsEXP.dDep[k][i] = clsUNI.X2I(dTab[k][j],sUni)
                    #print("k,dTab,cDep ",k,dTab[k][j],clsEXP.dDep[k][i])

    if nDepA != nFoun :
        print("Experiment ",sTyp," expects ",nDepA," columns of Dependent data: Only ",nFoun," found - Error")
        iERR = -1
        return iERR

#-- Dependent/Observed Psat for selected Experiments ----------------

    if sTyp == "CCE" or sTyp == "CVD" or sTyp == "DLE" or sTyp == "SEP" :
        nRsat = clsEXP.nRsat
        if nRsat < 0 :
            print("Experiment ",sTyp," must have one row starting with PSAT - Error")
            iERR = -1
            return iERR
        else :
            PsatO = clsEXP.dDep[nRsat][0]
            PsatU = clsEXP.uDep[0]
            clsEXP.setPsatObs(PsatO)
            clsEXP.setPsatWei(dWps)
            clsEXP.setPsatUni(PsatU)

#-- Reference Depth and Pressure for GRD ----------------------------

    if sTyp == "GRD" :
        nDref = clsEXP.nDref
        if nDref < 0 :
            print("Experiment GRD must have one row starting with DREF - Error")
            iERR = -1
            return iERR
        else :
            Dref = clsEXP.dDep[nDref][0]
            Pref = None
            for i in range(nHead) :
                sObs = colH[i].upper()
                sUni = colU[i].upper()
                if sObs[:2] == "PR" :
                    Pref = clsUNI.X2I(dTab[nDref][i],sUni)
                    break
            if Pref == None :
                print("Experiment GRD must a Pressure define on DREF row - Error")
                iERR = -1
                return iERR

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
                    clsEXP.dObs[k][j] = clsUNI.X2I(dTab[k][i],sUni)
                    clsEXP.dWei[k][j] = dWal[j]
                nFoun += 1
        if sTyp == "SEP" and sObs == "LSEP" :
            for k in range(nRowS) :
                clsEXP.Lsep[k] = int(dTab[k][i])
        if sTyp == "SEP" and sObs == "VSEP" :
            for k in range(nRowS) :
                clsEXP.Vsep[k] = int(dTab[k][i])

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

    return iERR
    
#========================================================================
#  End of Module
#========================================================================
