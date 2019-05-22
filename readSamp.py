
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Samples data for Fluid Initialisation
#
#=======================================================================

#!/usr/bin/python3

import numpy     as NP

import constants as CO
import calcProps as CP

#== One Sample ========================================================

class classSample :

    def __init__(self,name) :
        self.sName = name
        self.qPhas = False

    def setPlsMw(self,Mw)   : self.mPlus = Mw
    def setPlsSG(self,Sg)   : self.sPlus = Sg
    def setPlsAl(self,Al)   : self.aPlus = Al

    def setNComp(self,Nc)   :
        self.nComp = Nc
        self.sComp = NP.zeros(Nc)

    def setNSamp(self,Ns)   : self.NS = Ns

    def setNUser(self,Nu)   : self.NU = Nu

    def incNComp(self,newNC) :

        nEXT = newNC - self.nComp
        
        aEXT = [0.0 for i in range(nEXT)]     #-- Dummy vector for Props

        self.sComp.extend(aEXT)
        
        self.nComp = newNC

    def sZI(self,iC,zC) : self.sComp[iC] = zC

    def gZI(self,iC)    : return self.sComp[iC]

    def gTot(self) :
        sTot = 0.0
        for iC in range(self.nComp) :
            sTot = sTot + self.sComp[iC]
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

        if   self.TDmax != None :
            self.tHrm = self.TDmax
        elif self.TBmax != None :
            self.tHrm = self.TBmax
        else :
            self.tHrm = None

        #print("pBar,tHrm {:10.3f} {:10.3f}".format(self.pBar,self.tHrm))

#========================================================================
#  clsIO  is the Input/Output Class
#  tokU is a List containing the User Defined Samples
#========================================================================

def readSamp(clsIO,tokS,clsEOS,dicSAM,clsUNI) :

    iERR = 0
    nUsr = 0

    nSplt = clsEOS.NP

    print("procSamp: Processing the SAMPLES keyword")

    cNam = []
    cCom = []

    qMolW = False
    qSpcG = False
    qAlfa = False

    nSamp = len(tokS) - 1                     #-- First token is SAMP keyword

    clsEOS.setNSamp(nSamp)

#----------------------------------------------------------------------
#  Process Samples
#----------------------------------------------------------------------

    for iS in range(nSamp) :
        
        sName = tokS[iS+1]

        clsSAM = classSample(sName)

        dicSAM[iS] = clsSAM

        #print("iS,dS ",dicSAM[iS].sName)

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

    iLine = 0        

    for curL in clsIO.fInP :

        #print("curL ", curL)

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines in SAMPLES keyword - Error")
            iERR = -1
            break

        tokS = curL.split()                  #-- Split the line into its Tokens
        nTok = len(tokS)

        if nTok == 0 :                       #-- Found a Blank - Exit Loop
            break

        if   tokS[0] == "--" :               #-- Eclipse-Style Comment => Break
            pass
        elif tokS[0].upper() == "MW" :       #-- Mole Weight                

            for iS in range(nSamp) :
                mPlus = float(tokS[iS+1])
                dicSAM[iS].setPlsMw(mPlus)

            qMolW = True
            
        elif tokS[0].upper() == "SG" :       #-- Specific Gravity                

            for iS in range(nSamp) :
                sPlus = float(tokS[iS+1])
                dicSAM[iS].setPlsSG(sPlus)
                
            qSpcG = True
            
        elif tokS[0][:2].upper() == "AL" :   #-- Alpha-Coefficient                

            for iS in range(nSamp) :
                aPlus = float(tokS[iS+1])
                dicSAM[iS].setPlsAl(aPlus)

            qAlfa = True
            
        else :

            if qMolW and qSpcG :
                pass
                #print("Plus Fraction MW & SG have been read - OK")
            else :
                print("Cannot read Compositions until at least MW & SG have been read - Error")
                break

            nUsr += 1
            cTmp = []

            cNam.append(tokS[0])              #-- Load the Component Name

            for iS in range(nSamp) :          #-- Load Component Moles in Tmp
                cTmp.append(float(tokS[iS+1]))

            cCom.append(cTmp)

#========================================================================
#  Load the Compositional Data into the Appropriate Structure
#========================================================================

    nCom = nUsr + nSplt - 1  #-- Post-Splitting, Fluid will have nCom-components

    comZ = NP.zeros(nUsr)

    for iS in range(nSamp) :

        clsSAM = dicSAM[iS]
        
        clsSAM.setNComp(nCom)  #-- Set N-Components and dimension Moles-Array
        
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
            if sNor < CO.macEPS : sNor = 1.0E-20   #-- If Z=0, set to small number
            clsSAM.sZI(iC,sNor)

    print("Have read Plus Fraction & Sample Compositions")

#========================================================================
#  Default Component Properties
#========================================================================

    clsEOS.setNComp(nCom)
    clsEOS.setNUser(nUsr)

#-- Load User Component Names into clsEOS class -------------------------

    for iC in range(nUsr) :
        clsEOS.sPP("CN",iC,cNam[iC])

#== Assign and Calculate All Component Properties (including Splitting)        

    iERR = CP.allProps(clsEOS,dicSAM,clsIO,clsUNI)
    
    return iERR,dicSAM

