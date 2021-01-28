
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Defines the All the Classes used by PVTfree
#
#=======================================================================

#!/usr/bin/python3

import sys
import numpy as NP

from math import sqrt

#========================================================================
#  Equation of State Class
#========================================================================

class classEOS :

    NAM = []  #-- Component Names

    KIJ = []  #-- Binary Interaction Parameters

    PRP = []  #-- Component Names, Mw's, Tc's, etc.
    
    nProp = 16  #-- Number of Arrays in PRP
                #-- List of Arrays in pMap below

    #  pMap     Array   Array                   Unit
    #  Index    Code    Name                    Type
    #
    #   0       MW      Mole Weight             lb/lbmol
    #   1       TC      Critical Temperature    degR
    #   2       PC      Critical Pressure       psia
    #   3       VC      Critical Volume         ft3/lbmol
    #   4       ZC      Critical Z-Factor       dim
    #   5       AF      Acentric Factor         dim
    #   6       TB      Normal Boiling Point    degR
    #   7       SG      Specific Gravity        wrtAir=1
    #   8       PA      Parachor                dim
    #   9       SS      Volume Shift            dim
    #  10       MA      Omega-A Multiplier      dim
    #  11       MB      Omega-B Multiplier      dim
    #  12       CA      Ideal Gas Cp A-Coeff    Btu/lbmol.degR
    #  13       CB      Ideal Gas Cp B-Coeff    Btu/lbmol.degR^2
    #  14       CC      Ideal Gas Cp C-Coeff    Btu/lbmol.degR^3
    #  15       CD      Ideal Gas Cp D-Coeff    Btu/lbmol.degR^4
    
    pMap = {"MW": 0, "TC": 1, "PC": 2, "VC": 3, "ZC": 4, \
            "AF": 5, "TB": 6, "SG": 7, "PA": 8, "SS": 9, \
            "MA":10, "MB":11, "CA":12, "CB":13, "CC":14, \
            "CD":15                                      }

#-- Creator ---------------------------------------------------------
    
    def __init__(self,EOS) :

        self.Tref = 536.67  #-- Ref Temp (for Ideal Gas Enthalpy) [25 degC -> degR]

        self.EOS    = EOS
        self.nComp  = 0       #-- No components initially
        self.PlusCN = 7       #-- Plus Fraction Carbon Number
        self.NumPsu = 3       #-- Number of Pseudo-Components

        self.EtaL   = None    #-- Optionally store last cubic solve
        self.EtaV   = None

#-- Set the EoS-Specific Parameters ---------------------------------

    def setEOSparms(self,EOS) :

        self.EOS = EOS
        
        if EOS == "SRK" :              #-- User set Soave-Redlich-Kwong (SRK) EoS
            
            self.n1 = 1.0               #-- Coeffs in the Martin EoS
            self.n2 = 2.0
            
            self.OA = 0.4274802335403410        #-- OmegaA
            self.OB = 0.0866403499649577        #-- OmegaB

            self.w0 =  0.47979          #-- Coeffs in the alpha(AcenF)       
            self.w1 =  1.576
            self.w2 = -0.1925
            self.w3 =  0.025

        else :                          #-- By default, Peng-Robinson (PR) EoS
            
            sqr2    = 1.4142135623730951

            self.n1 = 2.0 + sqr2
            self.n2 = 2.0 - sqr2
            
            self.OA = 0.45723552892138225
            self.OB = 0.07779607390388849

            self.w0 =  0.379642
            self.w1 =  1.48503
            self.w2 = -0.164423
            self.w3 =  0.016666

        self.n3 = self.n1 + self.n2         #-- Useful Work Coeffs
        self.n4 = self.n1 * self.n2

        self.n5 = 1.0/(self.n1 - self.n2)

#-- Set the Number of Components ------------------------------------

    def setNComp(self,nComp) :

        self.nComp = nComp  #-- Now dimension the arrays

        self.PRP = NP.zeros((self.nProp,self.nComp))
        self.KIJ = NP.zeros((self.nComp,self.nComp))
        self.LIJ = NP.zeros((self.nComp,self.nComp))

        self.NAM = NP.empty(self.nComp,dtype="<U12")
        
    def incNComp(self,newNC) :

        nEXT = newNC - self.nComp              #-- Number of Elements to extend by
        
        aEXT = NP.zeros(nEXT)
        bEXT = NP.zeros(newNC)

        for iPRP in range(self.nProp) :       #-- Extend the Parameter Arrays
            self.PRP[iPRP].extend(aEXT)

        for iC in range(self.nComp) :          #-- Extend existings rows of KIJ
            self.KIJ[iC].extend(aEXT)
            self.LIJ[iC].extend(aEXT)
            
        for iC in range(self.nComp,newNC) :    #-- Add extra rows of KIJ
            self.KIJ.append(bEXT)
            self.LIJ.append(bEXT)
            self.NAM.append("")

        self.nComp = newNC                     #-- Update N-Components
            
#-- Set the Number of Samples ---------------------------------------

    def setNSamp(self,nSamp) : self.nSamp = nSamp
            
#-- Set the Number of User Components -------------------------------

    def setNUser(self,nUser) : self.nUser = nUser
            
#-- Set the Number of Pseudo-Components (in Splitting) --------------

    def setNPseu(self,nPseu) : self.nPseu = nPseu

#-- Set/Get a component Mole Weight, etc ----------------------------

    def sPP(self,typ,iC,valU) : self.PRP[self.pMap[typ]][iC] = valU

    def gPP(self,typ,iC) : return self.PRP[self.pMap[typ]][iC]

#-- Set/Get a component Name ----------------------------------------

    def sNM(self,iC,sNam) : self.NAM[iC] = sNam

    def gNM(self,iC)      : return self.NAM[iC]

#-- Set/Get a BIC ---------------------------------------------------

    def sIJ(self,iC,jC,KIJ) :

        self.KIJ[iC][jC] =       KIJ
        self.LIJ[iC][jC] = 1.0 - KIJ

    def gIJ(self,iC,jC)  : return self.KIJ[iC][jC]

    def gLIJ(self,iC,jC) : return self.LIJ[iC][jC]

#-- Construct the no-(P,T) EoS Coefficients -------------------------    

    def eosCoefsNoPT(self) :

        nCom = self.nComp

        self.aP = NP.zeros(nCom)
        self.aQ = NP.zeros(nCom)
        self.bI = NP.zeros(nCom)
        self.vS = NP.zeros(nCom)
        self.Mw = NP.zeros(nCom)
        self.Tc = NP.zeros(nCom)
        self.Pc = NP.zeros(nCom)
        self.Vc = NP.zeros(nCom)

        OA0 = self.OA
        OB0 = self.OB

        iTc = self.pMap["TC"]
        iPc = self.pMap["PC"]
        iVc = self.pMap["VC"]
        iAF = self.pMap["AF"]
        iMA = self.pMap["MA"]
        iMB = self.pMap["MB"]
        iMw = self.pMap["MW"]
        iSS = self.pMap["SS"]

        for iC in range(nCom) :
            
            Tc = self.PRP[iTc][iC]
            Pc = self.PRP[iPc][iC]
            AF = self.PRP[iAF][iC]
            MA = self.PRP[iMA][iC]
            MB = self.PRP[iMB][iC]

            OA = OA0*MA
            OB = OB0*MB

            mAcF = self.w0 + AF*(self.w1 + AF*(self.w2 + AF*self.w3))

            P = 1.0 + mAcF
            Q =     - mAcF/sqrt(Tc)

            self.bI[iC] = OB*Tc/Pc

            aSqr = Tc*sqrt(OA/Pc)

            self.aP[iC] = aSqr*P
            self.aQ[iC] = aSqr*Q

            self.Tc[iC] = Tc
            self.Pc[iC] = Pc
            self.Vc[iC] = self.PRP[iVc][iC]
            self.Mw[iC] = self.PRP[iMw][iC]
            self.vS[iC] = self.PRP[iSS][iC]

#-- Update the 1-Kij Array ------------------------------------------

        self.LIJ = NP.subtract(1.0,self.KIJ)

#-- Dummy Scalar, Vector and Matrix for EoS Coef Calculations -------

        self.dumS = 0.0
        self.dumV = NP.zeros( nCom )
        self.dumM = NP.zeros((nCom,nCom))

        return

#========================================================================
#  Class to Hold Sample Data
#========================================================================

class classSample :

#-- Short Sample Name will be used to distinguish this sample -------    

    def __init__(self,sNam) :
           
        self.sNam  = sNam
        self.qPhas = False
      
        self.uNam = []  #-- User Names & Composition
        self.uCom = []

        self.iNam = []  #-- Internal Names & Moles
        self.iCom = []

#-- Sample Long Name ------------------------------------------------
      
    def setLongName(self,lNam) : self.lNam = lNam

#-- User Sample Plus Fraction Carbon Number, Mole Weight and Specific Gravity   

    def setUserPlus(self,uPlsCN,uPlsMW,uPlsSG) :
      
        self.uPlsCN = uPlsCN
        self.uPlsMW = uPlsMW
        self.uPlsSG = uPlsSG

#-- Special Components ----------------------------------------------

    def setSpecial(self,spcH2,spcHE,spcN2,spcCO2,spcH2S,spcNeo) :

        self.spcH2  = spcH2
        self.spcHE  = spcHE
        self.spcN2  = spcN2
        self.spcCO2 = spcCO2
        self.spcH2S = spcH2S
        self.spcNeo = spcNeo

#-- Number of Components in SCN 7 thru 10 ---------------------------

    def setSCN(self,scnC7,scnC8,scnC9,scnC10) :

        self.scnC7  = scnC7
        self.scnC8  = scnC8
        self.scnC9  = scnC9
        self.scnC10 = scnC10

#-- Total Number of User Specified Components -----------------------

    def setTotCom(self,nCtot) : self.nCtot = nCtot

#-- Set User Names & Compositions -----------------------------------

    def setUserName(self,uNam) : self.uNam.append(uNam)
    def setUserComp(self,uCom) : self.uCom.append(uCom)

#-- Set Int Plus Fraction Properties --------------------------------

    def setIntPlusMW(self,MW) : self.iPlsMW = MW
    def setIntPlusSG(self,SG) : self.iPlsSG = SG
    def setIntPlusAL(self,AL) : self.iPlsAL = AL

#-- Number of Components and Composition ----------------------------

    def setIntComp(self,nCom,nUsr) :
        self.nCom = nCom
        self.nUsr = nUsr
        self.sCom = NP.zeros(nCom)

#-- Internal Name & Moles -------------------------------------------

    def setInternal(self,iNam,iMol) :
        self.iNam.append(iNam)
        self.iCom.append(iMol)

#-- Set and Get a Component's Moles ---------------------------------

    def sZI(self,iC,zI) : self.sCom[iC] = zI
    def gZI(self,iC)    : return self.sCom[iC]

    def gTot(self) :    #-- Return total Moles
        sTot = 0.0
        for iC in range(self.nCom) : sTot = sTot + self.sCom[iC]
        return sTot

#-- Values from the Michelsen Approx Phase Envelope Calculation -----    

    def setTempPsat(self,BubT,BubP,DewT,DewP,FlsT,FlsP,pCri,tCri) :

        self.qPhas = True
        
        self.BubT = BubT ; self.BubP = BubP
        self.DewT = DewT ; self.DewP = DewP
        self.FlsT = FlsT ; self.FlsP = FlsP

        nBub = len(self.BubT)
        nDew = len(self.DewT)
        nFls = len(self.FlsT)

#-- Bubble Point Line -----------------------------------------------        

        if nBub > 0 :
            TBmin = NP.amin(self.BubT)
            TBmax = NP.amax(self.BubT)
            PBmax = NP.amax(self.BubP)
            self.TBmin = TBmin
            self.TBmax = TBmax
            self.PBmax = PBmax
        else :
            self.TBmin = None
            self.TBmax = None
            self.PBmax = None

#-- Dew Point Line --------------------------------------------------            

        if nDew > 0 :
            TDmin = NP.amin(self.DewT)
            TDmax = NP.amax(self.DewT)
            PDmax = NP.amax(self.DewP)
            self.TDmin = TDmin
            self.TDmax = TDmax
            self.PDmax = PDmax
        else :
            self.TDmin = None
            self.TDmax = None
            self.PDmax = None

#-- Critical Point --------------------------------------------------            

        self.pCri = pCri
        self.tCri = tCri

#-- CricondenBar ----------------------------------------------------

        if   self.PBmax != None and self.PDmax != None :
            self.pBar = max(self.PBmax,self.PDmax)
        elif self.PBmax != None :
            self.pBar = self.PBmax
        elif self.PDmax != None :
            self.pBar = self.PDmax
        else :
            self.pBar = None

#-- CricondenTherm --------------------------------------------------

        if   self.TDmax != None : self.tHrm = self.TDmax
        elif self.TBmax != None : self.tHrm = self.TBmax
        else                    : self.tHrm = None

        #print("pBar,tHrm {:10.3f} {:10.3f}".format(self.pBar,self.tHrm))

#========================================================================
#  Experiment Class
#========================================================================

class classEXP :

#-- Creator ---------------------------------------------------------    

    def __init__(self,Name) :
        
        self.xName = Name
        self.nRsat = -1
        self.nDref = -1
        self.PsatO = -1.0
        self.IsAct = True
        self.isDef = False
        self.sLCCE = True
        self.xOut  = []

#-- Set the Reservoir Temperature and its Units ---------------------
        
    def setTres(self,Tres,Tuni) :
        self.Tres = Tres
        self.Tuni = Tuni

#-- Get and set the Sample Number -----------------------------------
    
    def setSamp(self,nSamp,sName) :
        self.nSamp = nSamp
        self.sName = sName

#-- Get the Injection Sample Number (SWL only) ----------------------
    
    def setSinj(self,nSinj) : self.nSinj = nSinj
    def setSNin(self,sNinj) : self.sNinj = sNinj

#-- Set the CCE Liquid Saturation Method (DEW or TOT) ---------------

    def setSLCCE(self,sLCCE) : self.sLCCE = sLCCE

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

    def createExpArrays(self,nIndA,nObsA,nCalA,nRowS) :

        self.nInd = nIndA       #-- Independent Data
        self.nObs = nObsA       #-- Observed    Data
        self.nCal = nCalA       #-- Calculated  Data
        self.nRow = nRowS       #-- #Rows of    Data

        self.hInd = NP.empty( nIndA,dtype="<U4" )
        self.uInd = NP.empty( nIndA,dtype="<U12")
        self.dInd = NP.zeros((nIndA,nRowS))     #-- Ind Data

        self.hObs = NP.empty(nObsA,dtype="<U4" )
        self.uObs = NP.empty(nObsA,dtype="<U12")
        self.qObs = NP.full(nObsA,False)
        self.qCal = NP.full(nCalA,True )
        self.qPlt = NP.full(nObsA,True )
        self.dObs = NP.zeros((nObsA,nRowS))
        self.dWei = NP.zeros((nObsA,nRowS))

        self.hCal = NP.empty( nCalA,dtype="<U4" )
        self.uCal = NP.empty( nCalA,dtype="<U12")
        self.dCal = NP.zeros((nCalA,nRowS))

    def setUserObs(self,nObU) : self.nObU = nObU

    def createCompArray(self,nRowS,nComS,YZ,OC) :

        if YZ == "Y" :
            if OC == "O" : self.yObs = NP.zeros((nRowS,nComS))
            else         : self.yCal = NP.zeros((nRowS,nComS))
        else :
            if OC == "O" : self.zObs = NP.zeros((nRowS,nComS))
            else         : self.zCal = NP.zeros((nRowS,nComS))

#-- Separator Liquid and Vapour Next Stage Routing ------------------    

    def createSepStages(self,nRowS) :

        self.Lsep = NP.zeros(nRowS,dtype=int)
        self.Vsep = NP.zeros(nRowS,dtype=int)

#-- Speed-Up Regression by storing/retreiving Psat information ------

    def setPsatInfo(self,qBub,pSat,logK) :
        self.qBub = qBub
        self.pSat = pSat
        self.logK = logK    #-- Vector [nComp]

    def getPsatInfo(self,qDif,nCom) :
        if not qDif : return None,None,NP.empty(nCom)
        else        : return self.qBub,self.pSat,self.logK

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
#  Used by Routines for Blackoil Output
#========================================================================

class classBlack :

#-- Creator ---------------------------------------------------------

    OutU = "FLD"
    tSim = "ECL"

    def __init__(self)         :
        self.xBlk  = "BLACK"
        self.qLsep = False
        self.qVsep = False
        self.qDep  = False

    def setSepPres(self,pSep)  : self.pSep  = pSep
    def setSepTemp(self,tSep)  : self.tSep  = tSep
    def setSepLsep(self,Lsep)  :
        self.Lsep  = Lsep
        self.qLsep = True
    def setSepVsep(self,Vsep)  :
        self.Vsep  = Vsep
        self.qVsep = True

    def setPresUni(self,pUni)  : self.pResU  = pUni

    def setPmax(self,Pmax)     : self.Pmax  = Pmax
    def setPmin(self,Pmin)     : self.Pmin  = Pmin
    def setPinc(self,Pinc)     : self.Pinc  = Pinc
    def setPsat(self,Psat)     : self.Psat  = Psat

    def setExpType(self,xTyp)  : self.xTyp  = xTyp

    def setTres(self,Tres,Tuni):
        self.Tres  = Tres
        self.tResU = Tuni

    def setSamp(self,sNum,sNam):
        self.sNum = sNum
        self.sNam = sNam

    def setDepUni(self,dUni)   : self.depU  = dUni

    def setDepRef(self,dRef,dUni) :
        self.dRef  = dRef
        self.dRuni = dUni

    def setPrsRef(self,pRef,pUni) :
        self.pRef  = pRef
        self.pRuni = pUni

    def setDepths(self,dMin,dMax,dInc,dUni) :
        self.dMin  = dMin
        self.dMax  = dMax
        self.dInc  = dInc
        self.dDuni = dUni

    def setRT(self,RT)         : self.RT    = RT

    def setDenSTO(self,dSTO)   : self.dSTO  = dSTO
    def setDenSTG(self,dSTG)   : self.dSTG  = dSTG
    def setDenSTW(self,dSTW)   : self.dSTW  = dSTW

    def setMwSTO(self,mSTO)    : self.mSTO  = mSTO
    def setMwSTO(self,mSTG)    : self.mSTG  = mSTG

    def setCo(self,Co)         : self.Co    = Co

    def setBoSlope(self,BoS)   : self.BoS   = BoS
    def setBoInter(self,BoI)   : self.BoI   = BoI

    def setBoSlope(self,UoS)   : self.UoS   = UoS
    def setBoInter(self,UoI)   : self.UoI   = UoI

    def setCoSlope(self,CoS)   : self.CoS   = CoS
    def setCoInter(self,CoI)   : self.CoI   = CoI

    def setOutUnit(self,OutU)  : self.OutU  = OutU

    def setSimType(self,tSim)  : self.tSim  = tSim
    def setOilKeyW(self,oTyp)  : self.oTyp  = oTyp
    def setGasKeyW(self,gTyp)  : self.gTyp  = gTyp

    def setPsepUni(self,pSepU) : self.pSepU = pSepU
    def setTsepUni(self,tSepU) : self.tSepU = tSepU

    def setPrefWat(self,pRefW) : self.pRefW = pRefW
    def setPrefUni(self,pRefU) : self.pRefU = pRefU
    def setSaltWat(self,bSalt) : self.bSalt = bSalt
    def setSaltUni(self,saltU) : self.saltU = saltU
    def setBrine(self,setBr)   : self.setBr = setBr

    def setBrefWat(self,bRefW) : self.bRefW = bRefW
    def setUrefWat(self,uRefW) : self.uRefW = uRefW
    def setCrefWat(self,cRefW) : self.cRefW = cRefW
    def setVrefWat(self,vRefW) : self.vRefW = vRefW

#-- These are use in the Lucas (Gas) Viscosity Correlation ----------    

    def setOilCrit(self,TcOil,PcOil,VcOil,ZcOil) :
        self.TcOil = TcOil
        self.PcOil = PcOil
        self.VcOil = VcOil
        self.ZcOil = ZcOil
        
    def setGasCrit(self,TcGas,PcGas,VcGas,ZcGas) :
        self.TcGas = TcGas
        self.PcGas = PcGas
        self.VcGas = VcGas
        self.ZcGas = ZcGas

#-- Base [EOS0] and Modifier [EOS1] holders for the EoS parameters --

    nEOS = 6

    EOS0 = {"aOil":1.0,"bOil":1.0,"sOil":1.0,"aGas":1.0,"bGas":1.0,"sGas":1.0}
    EOS1 = {"aOil":1.0,"bOil":1.0,"sOil":1.0,"aGas":1.0,"bGas":1.0,"sGas":1.0}

    VIS0 = {"rOil":1.0,"eOil":1.0,"uOil":1.0,"rGas":1.0,"eGas":1.0,"uGas":1.0}
    VIS1 = {"rOil":1.0,"eOil":1.0,"uOil":1.0,"rGas":1.0,"eGas":1.0,"uGas":1.0}

#-- Pointers for the dTab (Data Table) Storage -----------------------

    iPr =  0  #-- Pressure
    iRs =  1  #-- GOR
    iBo =  2  #-- Oil FVF
    iUo =  3  #-- Oil Viscosity
    iCo =  4  #-- Oil Compressibility
    iVo =  5  #-- Oil Viscosibility
    iRv =  6  #-- CGR
    iBg =  7  #-- Gas FVF
    iUg =  8  #-- Gas Viscosity
    iBd =  9  #-- Gas Dry FVF
    iUd = 10  #-- Gas Dry Viscosity

#-- Pointers for the cTab (Calc Table) Storage -----------------------

    iXo =  0  #-- Oil [in Liquid] Mole Fraction
    iYo =  1  #-- Oil [in Vapour] Mole Fraction
    iDo =  2  #-- Oil Density
    iDg =  3  #-- Gas Density
    iMo =  4  #-- Oil Molar Volume
    iMg =  5  #-- Gas Molar Volume

#========================================================================
#  Used to hold data needed by Compositional Output routines
#========================================================================

class classComp :

    def __init__(self)          :
        self.xCmp = "COMP"
        self.qDep = False

    def setOutUnit(self,OutU)   : self.OutU  = OutU
    def setSimType(self,tSim)   : self.tSim  = tSim

    def setTres(self,Tres,Tuni):
        self.Tres  = Tres
        self.tResU = Tuni

    def setSamp(self,sNum,sNam):
        self.sNum = sNum
        self.sNam = sNam

    def setDepUni(self,dUni)   : self.depU  = dUni

    def setDepRef(self,dRef,dUni) :
        self.dRef  = dRef
        self.dRuni = dUni

    def setDepGOC(self,dGOC)   : self.dGOC  = dGOC

    def setPrsRef(self,pRef,pUni) :
        self.pRef  = pRef
        self.pRuni = pUni

    def setDepths(self,dMin,dMax,dInc,dUni) :
        self.dMin  = dMin
        self.dMax  = dMax
        self.dInc  = dInc
        self.dDuni = dUni

    def setDenSTO(self,dSTO)   : self.dSTO  = dSTO
    def setDenSTG(self,dSTG)   : self.dSTG  = dSTG
    def setDenSTW(self,dSTW)   : self.dSTW  = dSTW

    def setPrefWat(self,pRefW) : self.pRefW = pRefW
    def setPrefUni(self,pRefU) : self.pRefU = pRefU
    def setSaltWat(self,bSalt) : self.bSalt = bSalt
    def setSaltUni(self,saltU) : self.saltU = saltU
    def setBrine(self,setBr)   : self.setBr = setBr

    def setBrefWat(self,bRefW) : self.bRefW = bRefW
    def setUrefWat(self,uRefW) : self.uRefW = uRefW
    def setCrefWat(self,cRefW) : self.cRefW = cRefW
    def setVrefWat(self,vRefW) : self.vRefW = vRefW

#========================================================================
#  This class holds data needed to perform a Grouping operation
#========================================================================

class classGRP :

    def __init__(self,newNam)    : self.newNam = newNam

    def setOldNames(self,oldNam) : self.oldNam = oldNam

    def setOldNums(self,oldNum)  : self.oldNum = oldNum

#========================================================================
#  The cInpOut class handles input and output, including graphics and debug
#========================================================================

class classInpOut :

    def __init__(self,Path) :

        self.Path = Path  #-- Directory Path

        self.Pall = True  #-- Plot All Experiments/Observations

        self.nSav = 0     #-- Number of times SAV file generated

#-- Debug Flags; 0 => All Off ---------------------------------------        

        self.Deb  = {"PSAT" : 0, "FLASH": 0, "STAB2": 0, "STAB1": 0, \
                     "RR"   : 0, "REG"  : 0, "BLACK": 0, "GROUP": 0, \
                     "CCE"  : 0, "CVD"  : 0, "DLE"  : 0, "SEP"  : 0, \
                     "SAT"  : 0, "FLS"  : 0, "SWL"  : 0, "GRD"  : 0, \
                     "ENV"  : 0, "TSAT" : 0}

#--------------------------------------------------------------------
#  Option Flags; STAB = 0 (SS/GDEM) or 1 (BFGS)
#--------------------------------------------------------------------

        self.Opt  = {"STAB" : 0}

        self.Tit1 = None
        self.Tit2 = None
        self.Tit3 = None

    def setFiles(self,rNam,fInp,fOut,fSav,fLog) :

        self.rNam = rNam  #-- Rootname for the Input File

        self.fInp = fInp  #-- Input  File
        self.fOut = fOut  #-- Output File
        self.fSav = fSav  #-- Save   File
        self.fLog = fLog  #-- Log    File       

        self.patR = self.Path + "\\" + rNam  #-- Path and Rootname

        #print("clsIO: patR ",self.patR)
        
    def setQXL( self,qXL ) : self.qXL  = qXL

    def setQREG(self,qReg) : self.qReg = qReg  #-- Regression File
    def setFREG(self,fReg) : self.fReg = fReg

    def setQDEB(self,qDeb) : self.qDeb = qDeb  #-- Debug File
    def setFDEB(self,fDeb) : self.fDeb = fDeb

    def setQ100(self,q100) : self.q100 = q100  #-- E100
    def setF100(self,f100) : self.f100 = f100

    def setQIMX(self,qIMX) : self.qIMX = qIMX  #-- IMEX
    def setFIMX(self,fIMX) : self.fIMX = fIMX

    def setQMOR(self,qMOR) : self.qMOR = qMOR  #-- Tempest-MORE
    def setFMOR(self,fMOR) : self.fMOR = fMOR

    def setQVIP(self,qVIP) : self.qVIP = qVIP  #-- VIP/Nexus
    def setFVIP(self,fVIP) : self.fVIP = fVIP

    def setQ300(self,q300) : self.q300 = q300  #-- E300
    def setF300(self,f300) : self.f300 = f300

    def setQGEM(self,qGEM) : self.qGEM = qGEM  #-- GEM
    def setFGEM(self,fGEM) : self.fGEM = fGEM

    def setDocStr(self,sDS) : self.sDocStr = sDS  #-- DocString

#-- Re-direct 'print' if in interactive mode ------------------------

    def redirectPrint(self) :

        #self.Orig  = sys.stdout     #-- Save Standard Output
        #sys.stdout = self.fLog      #-- Re-direct to the fLog File

        pass

#========================================================================
#  End of Module
#========================================================================
