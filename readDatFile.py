
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see mainPVTFreeWin.py for further details
#
#  Reads the Input DAT file and USR file
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP
import os
#import wxPython

import allData     as AD
import calcIntSamp as CI
import calcProps   as CP
import readGen     as RG
import readExps    as RX
import utilities   as UT
import writeOut    as WO

#========================================================================
#  Read DAT (Input) File
#========================================================================

def readDatFile(filN,clsEOS,dicSAM,dicEXP,clsUNI,clsIO) :

    iERR = 0

#-- Process the fileName which contains directory, rootname & extension (.dat)    

    path,root = UT.procFileName(filN)

    #print("readDatFile: path,root {:s} {:s}".format(path,root))

    clsIO = UT.updateIOFiles(path,root,clsIO)

#----------------------------------------------------------------------
#  Parse the Input File Line by Line
#----------------------------------------------------------------------

    fInp = clsIO.fInp
    fOut = clsIO.fInp

    iCnt = 0    #-- Count the number of lines read

    for curL in fInp :  #-- For Current line in the Input File

        iCnt += 1

#-- Split the Current Line into tokens using space ------------------

        tokenS = curL.split()
        nToken = len(tokenS)

#-- Blank Line or Comment? ------------------------------------------

        if   nToken        == 0    : pass
        elif tokenS[0][:2] == '--' : pass

#== INITial Section =================================================

        elif tokenS[0][:].upper() == 'INIT' :

            iERR,clsEOS = readINIT(dicSAM,clsUNI,clsIO,GUI=clsEOS)

            if iERR < 0 : break

#== DEFine Section =================================================

        elif tokenS[0][:3].upper() == 'DEF' :

            fDef = clsIO.fInp
            iERR,clsEOS,dicSAM = readDEF(fDef,dicSAM,clsIO)

            if iERR < 0 : break

#== EXP Section =======================================================

        elif tokenS[0][:3].upper() == 'EXP' :

            iERR = RX.readExps(clsIO,clsEOS,dicSAM,dicEXP,clsUNI)

            if iERR < 0 : break

#== Return ============================================================

    return iERR

#=======================================================================
#  Read INITial Section
#  If GUI mode, kwargs passes clsEOS, else created here after EOS k/w
#=======================================================================

def readINIT(dicSAM,clsUNI,clsIO,**kwargs) :

    iERR  = 0
    iCnt  = 0
    nSplt = 0

    fInp = clsIO.fInp   #-- Input  File (that we are reading)
    fOut = clsIO.fOut   #-- Output File

#----------------------------------------------------------------------
#  Parse Input File Line by Line
#----------------------------------------------------------------------

    for curL in fInp :
    
        iCnt += 1

#-- Current Line in curL; split into List tokenS ----------------------
    
        tokenS = curL.split()
        nToken = len(tokenS)

#== Process Options ===================================================

        if   nToken                == 0      : pass   #-- Blank Line
        elif tokenS[0][:2]         == "-- "  : pass   #-- Comment
        elif tokenS[0][:4].upper() == "ENDI" : break  #-- ENDINIT k/w => Exit

#-- EOS K/W Read => Create the clsEOS to hold data ------------------
        
        elif tokenS[0].upper()     == "EOS"  :

            EOS = tokenS[1].upper()

#-- If called from GUI, clsEOS passed-in ----------------------------
            
            qGUI = False
            
            for key,val in kwargs.items() :
                if key == 'GUI' :
                    qGUI   = True
                    clsEOS = val

            if not qGUI : clsEOS = AD.classEOS(EOS)

            clsEOS.setEOSparms(EOS)

            print("Equation of State Specified as " + EOS)

            sCom = "--"

            WO.outputHeader(fOut,sCom,clsIO)
            WO.outputEOS(fOut,sCom,clsIO,clsEOS)

#-- SPLIT k/w -------------------------------------------------------
                    
        elif tokenS[0][:3].upper() == "SPL"  :

            nSplt = int(tokenS[1])
            if nSplt < 2 or nSplt > 5 :
                print("SPLIT: nSplt " + str(nSplt) + " Out of Range, 2 =< Nsplt =< 5 - Error")
                iERR = -1
                break
            else:
                clsEOS.NumPsu = nSplt
                print("Plus Fraction will be split into " + str(nSplt) + " Pseudo-Components")

#-- SAMPLES k/w -----------------------------------------------------
                    
        elif tokenS[0][:4].upper() == "SAMP" :

            if nSplt == 0 :
                nSplt = 1
                clsEOS.NumPsu = nSplt
                
            iERR = readSAMP(tokenS,clsEOS,dicSAM,clsUNI,clsIO)
            if iERR < 0 : break

#-- DEBUG k/w -------------------------------------------------------

        elif tokenS[0][:3].upper() == "DEB"  :

            iERR = RG.readDebug(clsIO)
            if iERR < 0 : break

#-- OPTION k/w ------------------------------------------------------

        elif tokenS[0][:3].upper() == "OPT"  :

            iERR = RG.readOption(clsIO)
            if iERR < 0 : break

 #== Return values =====================================================

    return iERR,clsEOS

#========================================================================
#  Read SAMPle Keyword and Associated Data
#========================================================================

def readSAMP(tokenS,clsEOS,dicSAM,clsUNI,clsIO) :

    iERR = 0
    nUsr = 0
    #fLog = clsIO.fLog

    nSplt = clsEOS.NumPsu

    print("procSamp: Processing the SAMPLES keyword")

    cNam = []
    cCom = []

    qMolW = False
    qSpcG = False
    qAlfa = False

    nSamp = len(tokenS) - 1     #-- First token is SAMP keyword

    clsEOS.setNSamp(nSamp)

#----------------------------------------------------------------------
#  Process Samples
#----------------------------------------------------------------------

    for iS in range(nSamp) :
        
        sName  = tokenS[iS+1]       #-- Subsequent names are the Sample Short Names
        clsSAM = AD.classSample(sName)
        
        dicSAM[iS] = clsSAM         #-- Dictionary of Samples

#========================================================================
#  If more than one Sample, user must have entered the SPLIT keyword 1st
#========================================================================

    if nSamp > 1 and nSplt == 0 :

        print("More Than One Sample But No Precceding SPLIT keyword - Error")
        iERR = -1
        return iERR

#========================================================================
#  Read the following lines until a blank or a comment
#========================================================================

    fInp = clsIO.fInp   #-- Input File

    iLine = 0        

    for curL in fInp :

        #print("curL ", curL)

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines in SAMPLES keyword - Error")
            iERR = -1
            break

        tokenS = curL.split()                  #-- Split the line into its Tokens
        nToken = len(tokenS)

        if nToken == 0 : break                 #-- Found a Blank - Exit Loop
            

        if   tokenS[0]         == "--" : pass  #-- Eclipse-Style Comment => Break
        elif tokenS[0].upper() == "MW" :       #-- Mole Weight                

            for iS in range(nSamp) :
                mPlus = float(tokenS[iS+1])
                dicSAM[iS].setIntPlusMW(mPlus)
            qMolW = True
            
        elif tokenS[0].upper() == "SG" :       #-- Specific Gravity                

            for iS in range(nSamp) :
                sPlus = float(tokenS[iS+1])
                dicSAM[iS].setIntPlusSG(sPlus)
            qSpcG = True
            
        elif tokenS[0][:2].upper() == "AL" :   #-- Alpha-Coefficient                

            for iS in range(nSamp) :
                aPlus = float(tokenS[iS+1])
                dicSAM[iS].setIntPlusAL(aPlus)
            qAlfa = True
            
        else :

            if qMolW and qSpcG :
                pass
            else :
                print("Cannot read Compositions until at least MW & SG have been read - Error")
                break

            nUsr += 1
            cTmp = []

            cNam.append(tokenS[0])            #-- Load the Component Name

            for iS in range(nSamp) :          #-- Load Component Moles in Tmp
                cTmp.append(float(tokenS[iS+1]))

            cCom.append(cTmp)

#========================================================================
#  Load the Compositional Data into the Appropriate Structure
#========================================================================

    nCom = nUsr + nSplt - 1  #-- Post-Splitting, Fluid will have nCom-components
    comZ = NP.zeros(nUsr)

    for iS in range(nSamp) :

        clsSAM = dicSAM[iS]
        
        clsSAM.setIntComp(nCom,nUsr)  #-- Set N-Components and dimension Moles-Array
        
        sumT = 0.0
        for iC in range(nUsr) :
            comZ[iC] = cCom[iC][iS]
            sumT = sumT + comZ[iC]

#-- Ensure the Composition Sums to Unity or 100.0 -------------------

        #print("procSamp: iS,sumT ",iS,sumT)

        if sumT > 10.0 :
            if abs(sumT-100.0) > 0.5 :
                print("Sum of Moles for Sample ",iS+1," = [",dicSAM[iS].sName,"] Must Equal 100.0 - Error")
                iERR = -1
                return iERR
        else :
            if abs(sumT-1.0) > 0.005 :
                print("Sum of Moles for Sample ",iS+1," = [",dicSAM[iS].sName,"] Must Equal 1.0 - Error")
                iERR = -1
                return iERR            
        
#-- Force Normalisation if Slightly-Out -----------------------------        

        sumT = 1.0/sumT
        for iC in range(nUsr) :
            sNor = sumT*comZ[iC]
            if sNor < UT.macEPS : sNor = 1.0E-20   #-- If Z=0, set to small number
            clsSAM.sZI(iC,sNor)

#========================================================================
#  Default Component Properties
#========================================================================

    print("Have read Plus Fraction & Sample Compositions")

    clsEOS.setNComp(nCom)
    clsEOS.setNUser(nUsr)
    clsEOS.setNPseu(nSplt)
    clsEOS.setNSamp(nSamp)

#-- Load User Component Names into clsEOS class -------------------------

    for iC in range(nUsr) : clsEOS.sNM(iC,cNam[iC])

#== Assign and Calculate All Component Properties (including Splitting)        

    iERR = CP.allProps(clsEOS,dicSAM,clsUNI,clsIO)
    
#== Return Error Flag =================================================

    return iERR

#========================================================================
#  Open a User Defined Sample via File Dialog of USR format
#========================================================================

'''

def openUsrFile(self,clsEOS,dicSAM,clsUNI,clsIO) :

    iERR = 0

#== Use FileDialog to find a Suitable File ============================

    openFD = wx.FileDialog(self,"Open Input (.usr) File",os.getcwd(),"", \
                           "Input Files (*.usr)|*.usr", \
                           wx.FD_OPEN|wx.FD_FILE_MUST_EXIST)

    if openFD.ShowModal() == wx.ID_OK :          #-- Display The Dialog Pop-Up
        filU = openFD.GetPath()                  #-- Get Path/File Name
        iERR = readUsrFile(filU,dicSAM,clsEOS,clsUNI,clsIO)
    else :
        iERR = -1
        filU = "Not Defined"
        wx.MessageBox("Define USR Data or Cancel from 'Prepare New Sample Definition'", \
                      "No USR file Selected",wx.OK|wx.ICON_WARNING)

#== Return Error Flag =================================================

    return iERR

'''    

#========================================================================
#  Read Usr File
#========================================================================

def readUsrFile(filU,dicSAM,clsEOS,clsUNI,clsIO) :

    iERR = 0

    '''
    Looking for Short-Name, Long-Name, Component Names & Moles
    and Plus Fractions Properties (MW & SG)
    '''

    sNam = ''
    lNam = ''

    cNam = []
    cMol = []

    mPls = -1.0
    gPls = -1.0

#----------------------------------------------------------------------
#  Read USR File Line-By-Line
#----------------------------------------------------------------------

    fUsr = open(filU,'r')   #-- Open the Usr File
    iCnt = 0                #-- Count the number of lines read

    for curL in fUsr :  #-- For Current line in the USR File

        iCnt += 1

#-- Split the Current Line into tokens using space ------------------

        tokenS = curL.split()
        nToken = len(tokenS)

        #print("iCnt,nTok,curL ",iCnt,nToken,curL)
                
#-- Blank Line or Comment? ------------------------------------------

        if   nToken        == 0    : pass
        elif tokenS[0][:2] == '--' : pass

#-- Short-Name (used to 'index' the Sample) -------------------------

        elif tokenS[0][:4].upper() == "SNAM" : sNam = tokenS[1]

#-- Long-Name (for information only) --------------------------------

        elif tokenS[0][:4].upper() == "LNAM" :
            for iT in range(1,nToken) :         #-- 1st Token was K/W
                lNam = lNam + tokenS[iT] + " "

#-- Plus Fraction Mole Weight & Specific Gravity --------------------

        elif tokenS[0][:2].upper() == "MW" : mPls = float(tokenS[1])
        elif tokenS[0][:2].upper() == "SG" : gPls = float(tokenS[1])

#-- Else, Assumed to be a Component Name and #Moles ===================

        else :
            cNam.append(      tokenS[0] )
            cMol.append(float(tokenS[1]))

#== Define Sample Class and Assign Names, Moles, etc. =================

    if sNam == '' :
        print("No SNAME provided in USR File - Data Will Not Be Processed")
        iERR = -1
        return iERR

#-- Check New Sample Name is Unique ---------------------------------

    nOld = len(dicSAM)
    #print("readUsrFile: sNam,nOld ",sNam,nOld)

    qNew = True
    for clsSAM in dicSAM.values() :
        if clsSAM.sNam == sNam :
            qNew = False
            break

    if not qNew :
        print("SNAME provided in USR File already used - Data Will Not Be Processed")
        iERR = -1
        return iERR

#-- Plus Fraction MW & SG Physical? ===================================

    if mPls < 90.0 or mPls > 1000.0 :
        print("MWPLUS provided in USR File not valid - Data Will Not Be Processed")
        iERR = -1
        return iERR

    if gPls < 0.6 or gPls > 1.1 :
        print("SGPLUS provided in USR File not valid - Data Will Not Be Processed")
        iERR = -1
        return iERR

#== Create the New Sample Class =======================================

    clsSAM = AD.classSample(sNam)

    if lNam != '' : clsSAM.setLongName(lNam)

#-- Get the Carbon Number of the Plus Fraction (CXX+) ---------------

    sPls = cNam[-1]

    sPlt = sPls.split('C')  #-- Remove Leading  'C'
    sPls = sPlt[1]
    sPlt = sPls.split('+')  #-- Remove Trailing '+'
    nPls = int(sPlt[0])

    clsSAM.setUserPlus(nPls,mPls,gPls)

#== User Names and Moles ==============================================

    for iU in range(len(cNam)) :
        clsSAM.setUserName(cNam[iU])
        clsSAM.setUserComp(cMol[iU])

#== Create Internal Sample ============================================

    clsSAM.nCtot = len(cNam)

    CI.calcIntSamp(clsEOS,clsSAM)

    nOld = len(dicSAM)

    print("Number of User Comps: nUser ",clsEOS.nUser)

    if nOld == 0 :
        dicSAM[0] = clsSAM
        clsEOS.setNSamp(1)
        for iC in range(clsEOS.nUser) :
            print("iC,iNam ",iC,clsSAM.iNam[iC])
            clsEOS.sNM(iC,clsSAM.iNam[iC])
    else         :
        pass

#== And create properties =============================================

    CP.allProps(clsEOS,dicSAM,clsUNI,clsIO)

#== Return ============================================================

    return iERR

#=======================================================================
#  Pre-Process the SAV file to find the Required DEFI Section
#=======================================================================

def preProcDEF(clsEOS,dicSAM,clsIO) :

    iERR = 0
    iCnt = 0

    fDef = clsIO.fSav   #-- Required SAV file
    nSav = clsIO.nSav   #-- Get Latest DEF in SAV file

#-- Delete Existing EOS Class and Sample Dictionary -----------------    

    del clsEOS
    del dicSAM

#-- New Sample Dictionary -------------------------------------------

    dicSAM = {}

#----------------------------------------------------------------------
#  Loop Thru SAV File to Find Required DEFI section
#----------------------------------------------------------------------

    for curL in fDef :

        iCnt += 1

#-- Split Current Line into Tokens ----------------------------------

        tokenS = curL.split()
        nToken = len(tokenS)

#-- At least 2 Tokens?  First is DEF and 2nd is required number? ----

        if nToken >= 2 :

            if tokenS[0][:3].upper() == 'DEF' and int(tokenS[1]) == nSav :

                iERR,clsEOS,dicSAM = readDEF(fDef,dicSAM,clsIO)
                break

#== Return Values =====================================================

    return iERR,clsEOS,dicSAM

#=======================================================================
#  Read Define Section (Full Compositional Description)
#=======================================================================

def readDEF(fDef,dicSAM,clsIO) :

    iERR = 0
    iCnt = 0
    
    fOut = clsIO.fOut

#---------------------------------------------------------------------
#  Parse Input File Line by Line
#---------------------------------------------------------------------

    for curL in fDef :
    
        iCnt += 1

#-- Current Line in curL; split into List tokS
    
        tokS = curL.split( )
        nTok = len(tokS)
      
#-- Blank or Comment? -----------------------------------------------

        if nTok == 0 :  #-- Blank line
            pass
        else :

            if   tokS[0][:2] == "--"           : pass
            elif tokS[0][:4].upper() == "ENDD" : break
            elif tokS[0][:3].upper() == "DEB"  :
                iERR = RG.readDebug(clsIO)
            elif tokS[0][:3].upper() == "OPT"  :
                iERR = RG.readOption(clsIO)
            elif tokS[0].upper()     == "EOS"  :

#-- Equation of State -----------------------------------------------                    

                EOS = tokS[1].upper()

                clsEOS = AD.classEOS(EOS)

                clsEOS.setEOSparms(EOS)
                    
                print("Equation of State Specified as ",EOS)

                sCom = "--"

                WO.outputHeader(fOut,sCom,clsIO)
                WO.outputEOS(fOut,sCom,clsIO,clsEOS)

            elif tokS[0].upper() == "NCOMP" :

                nComp = int(tokS[1])
                clsEOS.setNComp(nComp)  #-- Dimension the Arrays
                    
            elif tokS[0].upper() == "NSAMP" :

                nSamp = int(tokS[1])
                clsEOS.nSamp = nSamp
                    
            elif tokS[0].upper() == "PROPS" :

                curL = next(fDef)  #-- Property Names
                curL = next(fDef)  #-- Property Units

                for iC in range(nComp) :
                    
                    curL = next(fDef)
                    tokS = curL.split( )
                    nTok = len(tokS)

                    if nTok != 17 :
                        print("DEFINE PROPS: Expecting 17 Columns, Only ",nTok," Read")
                        iERR = -1
                        break

                    sN =       tokS[ 0]
                    Mw = float(tokS[ 1])   ; Tc = float(tokS[ 2])
                    Pc = float(tokS[ 3])   ; Vc = float(tokS[ 4])
                    Zc = float(tokS[ 5])   ; AF = float(tokS[ 6])
                    Tb = float(tokS[ 7])   ; SG = float(tokS[ 8])
                    PA = float(tokS[ 9])   ; SS = float(tokS[10])
                    MA = float(tokS[11])   ; MB = float(tokS[12])
                    CA = float(tokS[13])   ; CB = float(tokS[14])
                    CC = float(tokS[15])   ; CD = float(tokS[16])

                    clsEOS.sNM(     iC,sN)
                    clsEOS.sPP("MW",iC,Mw) ; clsEOS.sPP("TC",iC,Tc)
                    clsEOS.sPP("PC",iC,Pc) ; clsEOS.sPP("VC",iC,Vc)
                    clsEOS.sPP("ZC",iC,Zc) ; clsEOS.sPP("AF",iC,AF)
                    clsEOS.sPP("TB",iC,Tb) ; clsEOS.sPP("SG",iC,SG)
                    clsEOS.sPP("PA",iC,PA) ; clsEOS.sPP("SS",iC,SS)
                    clsEOS.sPP("MA",iC,MA) ; clsEOS.sPP("MB",iC,MB)
                    clsEOS.sPP("CA",iC,CA) ; clsEOS.sPP("CB",iC,CB)
                    clsEOS.sPP("CC",iC,CC) ; clsEOS.sPP("CD",iC,CD)
                
            elif tokS[0].upper() == "BIP"   :

                curL = next(fDef)  #-- Component Names
                
                for iC in range(nComp) :
                    
                    curL = next(fDef)
                    tokS = curL.split( )
                    nTok = len(tokS)
                
                    if nTok != nComp + 1 :
                        print("DEFINE BIP: Expecting ",nComp+1," Columns, Only ",nTok," Read")
                        iERR = -1
                        break

                    iTok = 1
                    while iTok < nTok :
                        KIJ = float(tokS[iTok])
                        clsEOS.sIJ(iC,iTok-1,KIJ)
                        iTok += 1

            elif tokS[0].upper() == "SAMPLES" :
                
                curL = next(fDef)  #-- Sample Names
                tokS = curL.split()
                nTok = len(tokS)

                for iSamp in range(nSamp) :
                    sName = tokS[iSamp+1]
                    clsSAM = AD.classSample(sName)
                    dicSAM[iSamp] = clsSAM
                    dicSAM[iSamp].setIntComp(nComp,nComp)
                
                for iC in range(nComp) :
                    
                    curL = next(fDef)
                    tokS = curL.split( )

                    if nTok != nSamp + 1 :
                        print("DEFINE SAMPLES: Expecting ",nSamp+1," Columns, Only ",nTok," Read")
                        iERR = -1
                        break
                    
                    for iSamp in range(nSamp) :
                        ZI = float(tokS[iSamp+1])
                        if ZI < UT.macEPS : ZI = 1.0E-20  #-- Protect against Z = 0
                        dicSAM[iSamp].sZI(iC,ZI)

            else :
                pass

#========================================================================
#  Has any slop crept into the sample definitions?
#========================================================================

    for iSamp in range(nSamp) :
        sumT = 0.0
        for iC in range(nComp) :
            sumT = sumT + dicSAM[iSamp].gZI(iC)
        sumT = 1.0/sumT
        for iC in range(nComp) :
            zI = sumT*dicSAM[iSamp].gZI(iC)
            dicSAM[iSamp].sZI(iC,zI)

    clsEOS.nSamp = nSamp

#== Back-calculate the C7+ Properties =================================

    backCalcPlusFracProps(clsEOS,dicSAM)
            
#========================================================================
#  Do we need to sort the components?  Most to Least Volatility
#========================================================================

    clsEOS,dicSAM = CP.sortComponents(clsEOS,dicSAM)

#========================================================================
#  Output the data
#========================================================================

    sTit = "Initialisation from a saved DEFINE"

    WO.outputProps(sTit,clsEOS,dicSAM,clsIO)

#========================================================================
#  Write Fluid Description to the SAV file
#========================================================================

    WO.outputSave(sTit,clsEOS,dicSAM,clsIO)

#======================================================================
#  Generate (Approximate) Phase Plots
#======================================================================

    #CF.allSamplesPhasePlot(clsEOS,dicSAM,clsIO)

#== Return values =====================================================    

    return iERR,clsEOS,dicSAM

#========================================================================
#  From Full DEF, back-calculate the Plus Fraction Properties
#========================================================================

def backCalcPlusFracProps(clsEOS,dicSAM) :

    nCom = clsEOS.nComp
    nSam = len(dicSAM)

#== For each sample ===================================================    

    for iSam in range(nSam) :

        clsSAM = dicSAM[iSam]

        zPlus = 0.0 ; mPlus = 0.0 ; sPlus = 0.0

#-- For components with Mw > C(6.5) ---------------------------------        

        for iC in range(nCom) :

            Mw = clsEOS.gPP("MW",iC) ; SG = clsEOS.gPP("SG",iC)

            if Mw > 90.0 :

                ZI    = clsSAM.gZI(iC)
                zPlus = zPlus + ZI
                mPlus = mPlus + ZI*Mw
                sPlus = sPlus + ZI*Mw/SG

#-- SG by inverse weighting -----------------------------------------                

        sPlus = mPlus/sPlus ; mPlus = mPlus/zPlus

        clsSAM.iPlsMW = mPlus
        clsSAM.iPlsSG = sPlus
        clsSAM.iPlsAL = 1.0      #-- Default value for ALPHA (info lost)

#== No return values ==================================================

    return

#========================================================================
#  End of Module
#========================================================================

