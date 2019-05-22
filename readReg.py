#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Reads the Regression Data Module
#  
#========================================================================

#!/usr/bin/python3

import calcReg  as CR
import readGen  as RG
import readExps as RX

#========================================================================
#  Regression Variable Class
#========================================================================

class classREG :

    def __init__(self,vNam) :
        
        self.vNam = vNam

        if   vNam == "TCRIT" or vNam == "PCRIT" or vNam == "ACENF" or \
             vNam == "SHIFT" or vNam == "MULTA" or vNam == "MULTB"    :
            self.vTyp = "REG"
            self.vMin =  0.7
            self.vMax =  1.3
        elif vNam == "VCRIT" :
            self.vTyp = "LBC"
            self.vMin =  0.5
            self.vMax =  1.5
        elif vNam == "KIJ" :
            self.vTyp = "KIJ"
            self.vMin = -2.0
            self.vMax =  2.0
        elif vNam == "MPLUS" or vNam == "SPLUS" :
            self.vTyp = "PLS"
            self.vMin =  0.9
            self.vMax =  1.1
        elif vNam == "APLUS" :
            self.vTyp = "PLS"
            self.vMin =  0.5
            self.vMax =  2.5
        else :
            self.vTyp = "ERR"

    def setRegKIJone(self,vOne) : self.vOne = vOne

    def setRegSamp(self,vSam)   : self.vSam = vSam

    def setRegGroup(self,vGrp)  :
        self.vGrp = []
        self.vGrp = vGrp

#========================================================================
#  Read the Data between REG and ENDREG
#========================================================================

def readReg(clsIO,clsEOS,dicSAM,dicEXP,clsUNI) :

    nExp = len(dicEXP)

    print("readReg: Starting Regression with ",nExp," Experiments Defined")

    qExp = [False for i in range(nExp)]

#-- Process existing Experiments into a handy array -----------------

    rExp = procExistExp(dicEXP)

#-- Number of Iterations [clean-up later!] --------------------------

    mIter = None

#-- Dictionary of Regression Variables Held Locally -----------------

    dicREG = {}

#========================================================================
#  Read lines until we get to ENDREG (or too many lines!)
#========================================================================

    iERR  = 0
    iLine = 0

    fInP = clsIO.fInP

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of REG Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0 :  #-- Blank line!
            pass
        elif tokS[0].upper()     == "EXP"  :

            if nTok-1 > nExp :
                print("User appears to have declared more Experiments in Regression than previous defined - Error")
                iERR = -1
                break

            iTok = 1
            nTru = 0
            while iTok < nTok :
                sExp = tokS[iTok]
                #print("iTok,sExp ",iTok,sExp)
                for iExp in range(nExp) :
                    if sExp == rExp[iExp] :
                        qExp[iExp] = True
                        nTru      += 1
                iTok += 1

            if nTru < 1 :
                print("No Existing Experiments Selected within REG - Error")
                iERR = -1
                break

        elif tokS[0][:2].upper() == "NI"  :  #-- User Over-Writes Num Iterations

            mIter = int(tokS[1])

        elif tokS[0].upper()     == "VAR" :

            iERR = readRegVars(clsIO,dicREG,dicSAM,clsEOS)
            if iERR < 0 : break

        elif tokS[0][:3].upper() == "DEB" :

            iERR = RG.readDebug(clsIO)
            if iERR < 0 : break

        elif tokS[0][:4].upper() == "ENDR" :
            #print("ENDREG k/w found")
            break

        if iERR < 0 : break

    #print("readReg: Of which ",nTru," Experiments have been selected")
    #print("qExp ",qExp)

#======================================================================
#  Everything OK!  Run the Regression
#======================================================================

    if mIter == None : mIter = 5

    if iERR == 0 :

        nAct = 0
        for iExp in range(nExp) :
            if qExp[iExp] : nAct += 1

        print("readReg: Of Which ",nAct," Experiments Are Active")

        clsEOS,dicSAM,dicEXP = \
            CR.runRegression(mIter,clsIO,clsEOS,dicSAM,dicEXP,dicREG,qExp,clsUNI)

#== Return values =====================================================

    return iERR,clsEOS,dicSAM,dicEXP

#========================================================================
#  Read the Regression Variables
#========================================================================

def readRegVars(clsIO,dicREG,dicSAM,clsEOS) :

#========================================================================
#  Read lines until we get to ENDVAR (or too many lines!)
#========================================================================

    iERR  = 0
    iLine = 0
    qLBC  = False

    nVar  = 0

    fInP = clsIO.fInP

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of VAR [in REG] Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok                == 0       : pass  #-- Blank line!
        elif tokS[0][:2]         == "--"    : pass  #-- Comment!
        elif tokS[0][:4].upper() == "ENDV"  : break
        else :
            
            sTyp = tokS[0].upper()
            
            if   sTyp == "TCRIT" or sTyp == "PCRIT" or sTyp == "ACENF" or \
                 sTyp == "SHIFT" or sTyp == "MULTA" or sTyp == "MULTB" or \
                 sTyp == "VCRIT"                                          :

                iOff = 1
                qAll,iComs = checkComps4Reg(iOff,tokS,clsEOS)

                if not qAll :
                    print("Some invalid component names provided for ",sTyp," - Error")
                    iERR = -1
                    break
                else :
                    clsREG = classREG(sTyp)
                    clsREG.setRegGroup(iComs)
                
            elif sTyp == "KIJ" :

                sOne = tokS[1].upper()  #-- "From" Component
                iOne = checkComp4Reg(sOne,clsEOS)

                if iOne < 0 :
                    print("Some invalid component name ",sOne," provided for ",sTyp," - Error")
                    iERR = -1
                    break
                
                iOff = 2
                qAll,iComs = checkComps4Reg(iOff,tokS,clsEOS)
                
                if not qAll :
                    print("Some invalid component names provided for ",sTyp," - Error")
                    iERR = -1
                    break
                else :
                    clsREG = classREG(sTyp)
                    clsREG.setRegKIJone(iOne)
                    clsREG.setRegGroup(iComs)
                
            elif sTyp == "MPLUS" or sTyp == "SPLUS" or sTyp == "APLUS" :
                
                sName = tokS[1]                 #-- Sample Name
                iSamp = checkSamp4Reg(sName,dicSAM)
                
                if iSamp < 0 :
                    print(sTyp," Sample Name ",sName," Not Previously Defined - Error")
                    iERR = -1
                    break
                else :
                    clsREG = classREG(sTyp)
                    clsREG.setRegSamp(iSamp)

            elif sTyp[:3] == "LBC" :

                qLBC      = True
                iLow,iHig = setLBCgroups(clsEOS)

            else :
                print("Variable Type ",sTyp," Not Recognised - Error")
                iERR = -1
                break

#-- Increment the number of variables -------------------------------            

            if not qLBC : nVar += 1            

#== And add this variable to the dictionary of variables ==============            
            
        if not qLBC : dicREG[nVar-1] = clsREG

#== LBC Special =======================================================

    if qLBC :

        clsREG = classREG("VCRIT")
        clsREG.setRegGroup(iLow)
                    
        dicREG[nVar  ] = clsREG

        clsREG = classREG("VCRIT")
        clsREG.setRegGroup(iHig)
                    
        dicREG[nVar+1] = clsREG        

#----------------------------------------------------------------------
#  Test all variables are the same 'type'
#----------------------------------------------------------------------

    iERR,sReg = CR.regVarType(dicREG)

#== Return Error Flag =================================================                

    return iERR

#========================================================================
#  Process existing experiments into REG 'format'
#========================================================================

def procExistExp(dicEXP) :

    nExp = len(dicEXP)
    rExp = ["" for i in range(nExp)]

    nCCE = 0 ; nCVD = 0 ; nDLE = 0 ; nSEP = 0
    nFLS = 0 ; nSAT = 0 ; nSWL = 0 ; nGRD = 0

#== Do we have multiples of a given experiment type? ==================    
    
    for iExp in range(nExp) :
        
        clsEXP = dicEXP[iExp]
        
        xTyp  = clsEXP.xName

        rExp[iExp] = xTyp
        
        if   xTyp == "CCE" :
            nCCE += 1 ; rExp[iExp] = xTyp + ":" + str(nCCE)
        elif xTyp == "CVD" :
            nCVD += 1 ; rExp[iExp] = xTyp + ":" + str(nCVD)
        elif xTyp == "DLE" :
            nDLE += 1 ; rExp[iExp] = xTyp + ":" + str(nDLE)
        elif xTyp == "SEP" :
            nSEP += 1 ; rExp[iExp] = xTyp + ":" + str(nSEP)
        elif xTyp == "FLS" :
            nFLS += 1 ; rExp[iExp] = xTyp + ":" + str(nFLS)
        elif xTyp == "SAT" :
            nSAT += 1 ; rExp[iExp] = xTyp + ":" + str(nSAT)
        elif xTyp == "SWL" :
            nSWL += 1 ; rExp[iExp] = xTyp + ":" + str(nSWL)
        elif xTyp == "GRD" :
            nGRD += 1 ; rExp[iExp] = xTyp + ":" + str(nGRD)

#== Only one of a given type ==========================================            

    for iExp in range(nExp) :
        
        if   rExp[iExp][:3] == "CCE" and nCCE < 2 : rExp[iExp] = "CCE"
        elif rExp[iExp][:3] == "CVD" and nCVD < 2 : rExp[iExp] = "CVD"
        elif rExp[iExp][:3] == "DLE" and nDLE < 2 : rExp[iExp] = "DLE"
        elif rExp[iExp][:3] == "SEP" and nSEP < 2 : rExp[iExp] = "SEP"
        elif rExp[iExp][:3] == "FLS" and nFLS < 2 : rExp[iExp] = "FLS"
        elif rExp[iExp][:3] == "SAT" and nSAT < 2 : rExp[iExp] = "SAT"
        elif rExp[iExp][:3] == "SWL" and nSWL < 2 : rExp[iExp] = "SWL"
        elif rExp[iExp][:3] == "GRD" and nGRD < 2 : rExp[iExp] = "GRD"

#== Return Array of Experiment Types ==================================        

    return rExp

#========================================================================
#  Check a Sample Name is Valid
#========================================================================

def checkSamp4Reg(sName,dicSAM) :

    iSamp = -1
    nSamp = len(dicSAM)

    for iS in range(nSamp) :
        sThis = dicSAM[iS].sName
        if sThis == sName :
            iSamp = iS
            break

#== Return Sample Number ==============================================        

    return iSamp

#========================================================================
#  Check a Component Name is Valid
#========================================================================

def checkComp4Reg(sName,clsEOS) :

    iComp = -1
    nComp = clsEOS.NC

    for iC in range(nComp) :
        sThis = clsEOS.gPP("CN",iC)
        if sThis == sName :
            iComp = iC
            break

#== Return Component Number ===========================================        

    return iComp

#========================================================================
#  Check a List of Component Names is Valid
#========================================================================

def checkComps4Reg(nMin,comL,clsEOS) :

    nMax = len(comL)
    nLen = nMax - nMin
    qAll = True

    #print("checkComps4Reg: comL ",comL)

    iCom = [0 for i in range(nLen)]

    iAdd = 0
    for iL in range(nMin,nMax) :
        
        comT = comL[iL]
        iThs = checkComp4Reg(comT,clsEOS)
        
        if iThs < 0 :
            qAll = False
            #print("qAll = False: iL,comT ",iL,comT)
        else :
            iCom[iAdd] = iThs
            
        iAdd += 1

#== Return All-Component-Names-Found Flag and Array of Comp Nums ======

    return qAll,iCom

#========================================================================
#  LBC Special: Set Groups
#========================================================================

def setLBCgroups(clsEOS) :

    nCom = clsEOS.NC

    iLow = []
    iHig = []

    for iC in range(nCom) :

        Mw = clsEOS.gPP("MW",iC)

        if Mw < 60.0 :
            iLow.append(iC)
        else :
            iHig.append(iC)

#== Return Groups =====================================================

    return iLow,iHig

#========================================================================
#  End of Module
#========================================================================

