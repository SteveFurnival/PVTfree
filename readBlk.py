
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Data for Blackoil Table Generation Module
#
#=======================================================================

#!/usr/bin/python3

import blkDrive as BD
import readGen  as RG

class classBLK :

#-- Creator ---------------------------------------------------------

    OutU = "FLD"
    tSim = "ECL"

    def __init__(self)         :
        self.xBlk  = "BLACK"
        self.qLsep = False
        self.qVsep = False

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

    def setTres(self,Tres)     : self.Tres  = Tres
    def setTuni(self,Tuni)     : self.tResU = Tuni

    def setSample(self,sNum)   : self.sNum  = sNum

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
    #def setCMGtype(self,cmgT)  : self.cmgT  = cmgT

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
#  Read the BLACKOIL data for Input File
#========================================================================

def readBlack(clsIO,clsEOS,dicSAM,clsUNI) :

    iERR  = 0
    iLine = 0

    #print("readBlack: Entered")

    pSep = [] ; tSep = []
    Lsep = [] ; Vsep = []
    
    pMax = None
    pMin = None
    pInc = None

    sSam = None
    Tres = None

    xTyp = None

    pRef = None
    Salt = None

    fInP = clsIO.fInP

#== Create the clsBLK class to hold relevant data ======================    

    clsBLK = classBLK()

#----------------------------------------------------------------------
#  Loop over lines until ENDB[LACK] found
#----------------------------------------------------------------------

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of BLACK Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0 :                        #-- Blank line!
            pass
        elif tokS[0][:2]         == "--"   :    #-- Comment
            pass
        elif tokS[0][:2].upper() == "PS"   :    #-- Separator Train Pressures

            pSuni = tokS[1].lower()             #-- 1st Argument must be units
            qPuni = pSuni in clsUNI.sPres

            if qPuni :
                clsBLK.setPsepUni(pSuni.upper())
            else :
                print("First argument after PSEP must be valid Pressure Units [PSIA, PSIG, BARA, BARG] - Error")
                iERR = -1
                break

            iTok = 2
            while iTok < nTok :
                pVal = float(tokS[iTok])
                pSep.append(clsUNI.X2I(pVal,pSuni))
                iTok += 1

        elif tokS[0][:2].upper() == "TS"   :    #-- Separator Train Temps

            tSuni = tokS[1].lower()             #-- 1st Argument must be units
            tPuni = tSuni in clsUNI.sTemp

            if tPuni :
                clsBLK.setTsepUni(tSuni.upper())
            else :
                print("First argument after TSEP must be valid Tempearture Units [degR, degF, degC, Kelv] - Error")
                iERR = -1
                break

            iTok = 2
            while iTok < nTok :
                tVal = float(tokS[iTok])
                tSep.append(clsUNI.X2I(tVal,tSuni))
                iTok += 1

        elif tokS[0][:2].upper() == "LS" :      #-- Next Stage for Sep-Liq

            iTok = 1
            while iTok < nTok :
                lSep = int(tokS[iTok])
                Lsep.append(lSep)
                iTok += 1
            
        elif tokS[0][:2].upper() == "VS" :      #-- Next Stage for Sep-Vap

            iTok = 1
            while iTok < nTok :
                vSep = int(tokS[iTok])
                Vsep.append(vSep)
                iTok += 1
                       
        elif tokS[0][:2].upper() == "PR" :      #-- Pressure-RANGE (PRange)

            pRuni = tokS[1].lower()             #-- 1st Argument must be units
            qPuni = pRuni in clsUNI.sPres

            if qPuni :
                clsBLK.setPresUni(pSuni.upper())
            else :
                print("First argument after PRANGE must be valid Pressure Units [PSIA, PSIG, BARA, BARG] - Error")
                iERR = -1
                break

            iTok = 2
            while iTok < nTok :
                sVal = tokS[iTok].upper()
                if   sVal == "PMAX" :
                    iTok += 1
                    pMax  = float(tokS[iTok])
                elif sVal == "PMIN" :
                    iTok += 1
                    pMin  = float(tokS[iTok])
                elif sVal == "PINC" :
                    iTok += 1
                    pInc  = float(tokS[iTok])
                iTok += 1

        elif tokS[0][:2].upper() == "UN"   :    #-- UNITS keyword [FIELD/METRIC]

            uTyp = tokS[1].upper()

            if   uTyp[:1] == "M" : clsBLK.setOutUnit("MET")
            elif uTyp[:1] == "S" : clsBLK.setOutUnit("SI")      #-- CMG-Only!!
            else                 : clsBLK.setOutUnit("FLD")

        elif tokS[0][:2].upper() == "SI"   :    #-- Simulator Type [CMG/ECL/MOR/VIP or NEX]

            uTyp = tokS[1].upper()

            if   uTyp[:1] == "M" :                      #-- Tempest-MORE
                oTyp,gTyp = argsMORE(tokS)
                clsBLK.setSimType("MOR")
                clsBLK.setOilKeyW(oTyp)
                clsBLK.setGasKeyW(gTyp)
            elif uTyp[:1] == "C" :                      #-- CMG-IMEX
                oTyp      = argsIMEX(tokS)
                clsBLK.setSimType("CMG")
                clsBLK.setOilKeyW(oTyp)
            elif uTyp[:1] == "V" or uTyp[:1] == "N" :   #-- VIP or Nexus
                oTyp,gTyp = argsVIP(tokS)
                clsBLK.setSimType("VIP")
                clsBLK.setOilKeyW(oTyp)
                clsBLK.setGasKeyW(gTyp)
            elif uTyp[:1] == "E" :                      #-- Eclipse-100
                oTyp,gTyp = argsE100(tokS)
                clsBLK.setSimType("ECL")
                clsBLK.setOilKeyW(oTyp)
                clsBLK.setGasKeyW(gTyp)
            else :
                iERR = -1
                print("Error - User has not set the Simulator Type for B/O Table Output")
                break

        elif tokS[0][:1].upper() == "W" :       #-- Properties to define Water

            iTok = 1
            while iTok < nTok:
                if   tokS[iTok][:1].upper() == "P" :    #-- Reference Pressure
                    iTok += 1
                    pRef  = float(tokS[iTok])
                    iTok += 1
                    wUni  =       tokS[iTok]
                elif tokS[iTok][:1].upper() == "S" :    #-- Salinity
                    iTok += 1
                    Salt  = float(tokS[iTok])
                    iTok += 1
                    sUni  =       tokS[iTok]
                iTok += 1

        elif tokS[0][:3].upper() == "DEB" :

            iERR = RG.readDebug(clsIO)
            if iERR < 0 : break

        elif tokS[0][:4].upper() == "ENDB" :
            #print("readBlack: ENDBLACK k/w read")
            break

        elif tokS[0][:3].upper() == "CCE" or \
             tokS[0][:3].upper() == "CVD" or \
             tokS[0][:3].upper() == "DLE"    :

            xTyp = tokS[0][:3].upper()
        
#-- Line must contains non-Table information ------------------------                

        else :  

            iTok = 0
            while iTok < nTok :
                if   tokS[iTok].upper() == "SAMP" :
                    iTok += 1
                    sSam = tokS[iTok]
                elif tokS[iTok].upper() == "TRES" :
                    iTok += 1
                    Tres = float(tokS[iTok])
                    iTok += 1
                    Tuni =       tokS[iTok]
                iTok += 1
        
#----------------------------------------------------------------------
#  Process user data: Starting with Separator Train P's and T's
#----------------------------------------------------------------------

    nPsep = len(pSep)
    nTsep = len(tSep)
    nLsep = len(Lsep)
    nVsep = len(Vsep)

    if nPsep == 0 :
        print("No Separator Pressure Stages Defined via PSEP - Error")
        iERR = -1
        return iERR

    if nTsep == 0 :
        print("No Separator Temperature Stages Defined via TSEP - Error")
        iERR = -1
        return iERR

    if nPsep != nTsep :
        print("Number of PSEP Stages = ",nPsep," Must Equal Number of TSEP Stages = ",nTsep," - Error")
        iERR = -1
        return iERR

    if nLsep > 0 :
        if nLsep != nTsep :
            print("Number of LSEP Stages = ",nLsep," Must Equal Number of TSEP Stages = ",nTsep," - Error")
            iERR = -1
            return iERR
        clsBLK.setSepLsep(Lsep)
    else :
        Lsep = [0 for i in range(nTsep)]
        for iSep in range(nTsep-1) : Lsep[iSep] = iSep + 1
        clsBLK.setSepLsep(Lsep)

    if nVsep > 0 :
        if nVsep != nTsep :
            print("Number of VSEP Stages = ",nLsep," Must Equal Number of TSEP Stages = ",nTsep," - Error")
            iERR = -1
            return iERR
        clsBLK.setSepVsep(Vsep)
    else :
        Vsep = [0 for i in range(nTsep)]        
        clsBLK.setSepVsep(Vsep)

    clsBLK.setSepPres(pSep)
    clsBLK.setSepTemp(tSep)

#-- Pressure range (Max, Min: Inc) for table generation -------------    

    if pMax == None :
        print("Maximum Pressure (PMAX) Has Not Been Set With PRANGE - Error")
        iERR = -1
        return iERR
    else :
        clsBLK.setPmax(clsUNI.X2I(pMax,pRuni))
        
    if pMin == None :
        print("Minimum Pressure (PMIN) Has Not Been Set With PRANGE - Error")
        iERR = -1
        return iERR
    else :
        clsBLK.setPmin(clsUNI.X2I(pMin,pRuni))
        
    if pInc == None :
        print("Pressure Increment (PINC) Has Not Been Set With PRANGE - Error")
        iERR = -1
        return iERR
    else :
        clsBLK.setPinc(clsUNI.X2I(pInc,pRuni))

#-- Valid Reservoir Sample to work on? ------------------------------        

    if sSam == None :
        print("No Sample Name Has Been Read via SAMP argument - Error")
        iERR = -1
        return iERR
    else :
        nSam = len(dicSAM)
        qFnd = False
        for iSam in range(nSam) :
            sArr = dicSAM[iSam].sName
            if sSam.upper() == sArr.upper() :
                clsBLK.setSample(iSam)
                qFnd = True
                break
        if not qFnd :
            print("Did Not Find Sample ",sSam," In Previously Defined Samples - Error")
            iERR = -1
            return iERR

#-- Reservoir Temperature and its Units -----------------------------        

    if Tres == None :
        print("No Reservoir Temperature Has Been Read via TRES argument - Error")
        iERR = -1
        return iERR
    else            :
        clsBLK.setTuni(Tuni)
        clsBLK.setTres(clsUNI.X2I(Tres,Tuni))

#-- Experiment Type -------------------------------------------------        

    if xTyp == None :
        print("Experiment Type (CCE, CVD or DLE) Not Specified - Error")
        iERR = -1
        return iERR
    else            :
        clsBLK.setExpType(xTyp)

#-- Have we read the WATER card with pRef and Salt? -----------------

    if pRef == None and Salt == None :

        clsBLK.setBrine(False)
        clsBLK.setPrefWat(14.6959)   #-- Standard Pressure [psia]
        clsBLK.setPrefUni("psia")
        clsBLK.setSaltWat(0.0)       #-- Pure Water
        clsBLK.setSaltUni("mfrac")

    else :

        clsBLK.setBrine(True)

        if Salt != None :
            clsBLK.setSaltUni(sUni)
            clsBLK.setSaltWat(clsUNI.X2I(Salt,sUni))
        else :
            clsBLK.setSaltWat(0.0)       #-- Pure Water
            clsBLK.setSaltUni("mfrac")
        
        if pRef != None :
            clsBLK.setPrefUni(wUni)
            clsBLK.setPrefWat(clsUNI.X2I(pRef,wUni))
        else :
            clsBLK.setPrefWat(14.6959)   #-- Standard Pressure [psia]
            clsBLK.setPrefUni("psia")
            

#========================================================================
#  Data All Valid, Proceed to Calculate Tables
#========================================================================

    if iERR == 0 :

        print("readBlack: All BlackOil Data Read OK")

        BD.calcBlack(clsIO,clsBLK,clsEOS,dicSAM,clsUNI)

#========================================================================
#  End of Routine
#========================================================================

    return iERR

#========================================================================
#  E100 Simulator Arguments
#========================================================================

def argsE100(tokS) :

#-- Default Types for Oil & Gas Tables ------------------------------    

    oTyp = None
    gTyp = None

    nTok = len(tokS)    #-- First two args must be SIM & ECL

#== Addtional arguments? ==============================================

    for iArg in range(2,nTok) :
        
        sArg = tokS[iArg].upper()

        if   sArg == "PVCO" : oTyp = "PVCO"
        elif sArg == "PVDO" : oTyp = "PVDO"
        elif sArg == "PVTO" : oTyp = "PVTO"

        if   sArg == "PVDG" : gTyp = "PVDG"
        elif sArg == "PVTG" : gTyp = "PVTG"

#== Not set?  Defaults are PVTO & PVTG ================================

    if oTyp == None :
        print("User Has Not Set the Preferred E100 OIL Keyword - PVTO Assumed")
        oTyp = "PVTO"

    if gTyp == None :
        print("User Has Not Set the Preferred E100 GAS Keyword - PVTG Assumed")
        gTyp = "PVTG"

#== Return values =====================================================

    return oTyp,gTyp

#========================================================================
#  Tempest-MORE Simulator Arguments
#========================================================================

def argsMORE(tokS) :

#-- Default Types for Oil & Gas Tables ------------------------------    

    oTyp = None
    gTyp = None

    nTok = len(tokS)    #-- First two args must be SIM & MOR

#== Addtional arguments? ==============================================

    for iArg in range(2,nTok) :
        
        sArg = tokS[iArg].upper()

        if   sArg == "OPVT" : oTyp = "OPVT"
        elif sArg == "OPVD" : oTyp = "OPVD"
        elif sArg == "PVTO" : oTyp = "PVTO"

        if   sArg == "GPVT" : gTyp = "GPVT"
        elif sArg == "PVTG" : gTyp = "PVTG"

#== Not set?  Defaults are OPVT & GPVT ================================

    if oTyp == None :
        print("User Has Not Set the Preferred MORE OIL Keyword - OPVT Assumed")
        oTyp = "OPVT"

    if gTyp == None :
        print("User Has Not Set the Preferred MORE GAS Keyword - GPVT Assumed")
        gTyp = "GPVT"

#== Return values =====================================================

    return oTyp,gTyp

#========================================================================
#  CMG-IMEX Simulator Arguments
#========================================================================

def argsIMEX(tokS) :

#-- IMEX only takes one Table Type ----------------------------------

    Typ = None

    nTok = len(tokS)    #-- First two args must be SIM & CMG

#== Addtional arguments? ==============================================

    for iArg in range(2,nTok) :
        
        sArg = tokS[iArg].upper()

        if   sArg == "PVTVO"   : Typ = "PVTVO"
        elif sArg == "PVTG"    : Typ = "PVTG"
        elif sArg == "PVT"     : Typ = "PVT"
        elif sArg == "PVTCOND" : Typ = "PVTCOND"

#== Not set?  Defaults are PVTO & PVTG ================================

    if Typ == None :
        print("User Has Not Set the Preferred IMEX PVT Keyword - PVTVO Assumed")
        Typ = "PVTVO"

#== Return values =====================================================

    return Typ

#========================================================================
#  VIP/Nexus Simulator Arguments
#========================================================================

def argsVIP(tokS) :

#-- Default Types for Oil & Gas Tables ------------------------------    

    oTyp = None
    gTyp = None

    nTok = len(tokS)    #-- First two args must be SIM & VIP or NEX

#== Addtional arguments? ==============================================

    for iArg in range(2,nTok) :
        
        sArg = tokS[iArg].upper()

        if   sArg == "BOOTAB" : oTyp = "BOOTAB"
        elif sArg == "BODTAB" : oTyp = "BODTAB"

        if   sArg == "BOGTAB" : gTyp = "BOGTAB"
        elif sArg == "BDGTAB" : gTyp = "BDGTAB"

#== Not set?  Defaults are OPVT & GPVT ================================

    if oTyp == None :
        print("User Has Not Set the Preferred VIP/Nexus OIL Keyword - BOOTAB Assumed")
        oTyp = "BOOTAB"

    if gTyp == None :
        print("User Has Not Set the Preferred VIP/Nexus GAS Keyword - BOGTAB Assumed")
        gTyp = "BOGTAB"

#== Return values =====================================================

    return oTyp,gTyp

#========================================================================
#  End of Module
#========================================================================
