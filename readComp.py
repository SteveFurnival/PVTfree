
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read Data for Compositional Data Output
#
#=======================================================================

#!/usr/bin/python3

import datetime

import numpy     as NP

import constants as CO
import calcWater as CW
import writeOut  as WO

class classCMP :

    def __init__(self)          : self.xCmp  = "COMP"

    def setOutUnit(self,OutU)   : self.OutU  = OutU
    def setSimType(self,tSim)   : self.tSim  = tSim

    def setTres(self,Tres)     : self.Tres  = Tres
    def setTuni(self,Tuni)     : self.tResU = Tuni

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
#  Read the Data
#========================================================================

def readComp(clsEOS,dicSAM,clsIO,clsUNI) :

    iERR  = 0
    iLine = 0

    tRes  = None
    tUni  = None

    pRef = None
    Salt = None

    fInP = clsIO.fInP

#== Create the clsCMP class to hold relevant data ======================    

    clsCMP = classCMP()

#----------------------------------------------------------------------
#  Loop over lines until ENDC[OMP] found
#----------------------------------------------------------------------

    for curL in fInP :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of COMP Data - Error")
            iERR = -1
            break

        tokS = curL.split()
        nTok = len(tokS)

        if   nTok == 0 :                        #-- Blank line!
            pass
        elif tokS[0][:2]         == "--"   :    #-- Comment
            pass
        elif tokS[0][:4].upper() == "ENDC" :    #-- ENDCOMP read
            break

        elif tokS[0][:2].upper() == "UN"   :    #-- UNITS keyword [FIELD/METRIC]

            uTyp = tokS[1].upper()

            if uTyp[:1] == "M" : clsCMP.setOutUnit("MET")
            else               : clsCMP.setOutUnit("FLD")

        elif tokS[0][:2].upper() == "SI"   :    #-- Simulator Type [ECL/MOR]

            uTyp = tokS[1].upper()

            if   uTyp[:1] == "M" : clsCMP.setSimType("MOR")
            elif uTyp[:1] == "E" : clsCMP.setSimType("ECL")
            elif uTyp[:1] == "C" or \
                 uTyp[:1] == "G" : clsCMP.setSimType("CMG")
            elif uTyp[:1] == "V" or \
                 uTyp[:1] == "N" : clsCMP.setSimType("VIP")

        elif tokS[0][:2].upper() == "TR"   :    #-- Reservoir Temperature

            tRes = float(tokS[1])
            tUni =       tokS[2]

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

#----------------------------------------------------------------------
#  Check Data Has Been Defined
#----------------------------------------------------------------------

    if tRes == None :
        print("No Reservoir Temperature Has Been Read via TRES argument - Error")
        iERR = -1
        return iERR
    else            :
        clsCMP.setTuni(tUni)
        clsCMP.setTres(clsUNI.X2I(tRes,tUni))

#-- Have we read the WATER card with pRef and Salt? -----------------

    if pRef == None and Salt == None :

        clsCMP.setBrine(False)
        clsCMP.setPrefWat(14.6959)   #-- Standard Pressure [psia]
        clsCMP.setPrefUni("psia")
        clsCMP.setSaltWat(0.0)       #-- Pure Water
        clsCMP.setSaltUni("mfrac")

    else :

        clsCMP.setBrine(True)

        if Salt != None :
            clsCMP.setSaltUni(sUni)
            clsCMP.setSaltWat(clsUNI.X2I(Salt,sUni))
        else :
            clsCMP.setSaltWat(0.0)       #-- Pure Water
            clsCMP.setSaltUni("mfrac")
        
        if pRef != None :
            clsCMP.setPrefUni(wUni)
            clsCMP.setPrefWat(clsUNI.X2I(pRef,wUni))
        else :
            clsCMP.setPrefWat(14.6959)   #-- Standard Pressure [psia]
            clsCMP.setPrefUni("psia")
            
#----------------------------------------------------------------------
#  Data Read OK: Write Output
#----------------------------------------------------------------------

    if clsCMP.tSim == "ECL" : writeE300(clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "MOR" : writeMORE(clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "VIP" : writeVIP( clsEOS,dicSAM,clsCMP,clsIO,clsUNI)
    if clsCMP.tSim == "CMG" : writeGEM( clsEOS,dicSAM,clsCMP,clsIO,clsUNI)

#========================================================================
#  End of Routine
#========================================================================

    return iERR

#========================================================================
#  Write GEM Output
#========================================================================

def writeGEM(clsEOS,dicSAM,clsCMP,clsIO,clsUNI) :

#-- Initialisation --------------------------------------------------    

    sCom = "**"
    nCom = clsEOS.NC

    dVec = NP.zeros(nCom)

#-- Open the GEM File, if not already open -------------------------

    if not clsIO.qGEM :
        pathGEM = clsIO.patR + ".gem"
        fGEM    = open(pathGEM,'w')
        clsIO.setQGEM(True)
        clsIO.setFGEM(fGEM)
    else :
        fGEM = clsIO.fGEM

#== Headers ===========================================================        

    print("Writing Compositional Description for GEM")

    fSim = clsIO.fGEM
    sInp = clsIO.fInP.name

    OutU = clsCMP.OutU

    if OutU == "MET" : sUni = "METRIC"
    else             : sUni = "FIELD "

#-- Write Header Information to Simulator File ----------------------        

    WO.outputHeader(fSim,sCom,clsIO)
    WO.outputEOS(fSim,sCom,clsIO,clsEOS)

#-- Header ----------------------------------------------------------    

#               123456789012345678901234567890123456789012345678901234567890
    sHead = "**==========================================================\n"
    fSim.write(sHead)
    sLabl = "**  CMG-GEM EOS Model generated by PVTfree Program\n"
    fSim.write(sLabl)
    sLabl = "**  From dataset " + sInp + "\n"
    fSim.write(sLabl)
    sLabl = "**  User specified " + sUni + " Units\n"
    fSim.write(sLabl)
    fSim.write(sHead)
    fSim.write("\n")

#-- Equation of State -----------------------------------------------    

    sLabl = "** Equation of State: Peng-Robinson (PR) or Soave-Redlich-Kwong (SRK)\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sEoS = clsEOS.EOS
    sOut = "*MODEL  *" + sEoS + "\n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Number of Components --------------------------------------------

    sLabl = "** Number of Components; All assumed to be USER components\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sOut = "*NC  " + str(nCom) + "  " + str(nCom) + "\n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Component Names -------------------------------------------------    

    sLabl = "** Component Names\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sOut = "*COMPNAME  "
    for iC in range(nCom) :
        sNam = clsEOS.gPP("CN",iC)
        sOut = sOut + "'" + sNam + "'  "
    sOut = sOut + "\n"
    fSim.write(sOut)
    fSim.write("\n")

    fSim.write("*EOSSET  1\n")
    fSim.write("\n")

#-- Mole Weights ----------------------------------------------------    

    sUni = ": [gm/gmol]"
    sCom = "**  Molecular Weights " + sUni + "\n"
    sKey = "*MW"
    nLen = 7
    sFor = "{:8.3f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MW",iC)
    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Specific Gravity ------------------------------------------------

    sCom = "**  Specific Gravity\n"
    sKey = "*SG"
    nLen = 7
    sFor = "{:8.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("SG",iC)
    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Boiling Point Temperature ---------------------------------------

    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("TB",iC)
    if OutU == "MET" :
        sUni = ": [degC]"
        for iC in range(nCom) : dVec[iC] = clsUNI.I2X(dVec[iC],"degc")
    else :
        sUni = ": [degF]"
        for iC in range(nCom) : dVec[iC] = clsUNI.I2X(dVec[iC],"degf")
    sCom = "**  Boiling Point Temperature " + sUni + "\n"
    sKey = "*TB"
    nLen = 7
    sFor = "{:8.3f}"

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Temperatures -------------------------------------------
    
    sUni = ": [Kelvin]"
    for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("TC",iC),"kelv")
    sCom = "**  Critical Temperatures " + sUni + "\n"
    sKey = "*TCRIT"
    nLen = 7
    sFor = "{:8.3f}"

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Pressures ----------------------------------------------    

    sUni = ": [atm]"
    sFor = "{:8.4f}"
    for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("PC",iC),"atm")
    sCom = "**  Critical Pressures " + sUni + "\n"
    sKey = "*PCRIT"
    nLen = 7

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Volumes ------------------------------------------------

    sUni = ": [m3/kgmol]"
    sFor = "{:8.3f}"
    for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("VC",iC),"m3/kgmol")
    sCom = "**  Critical Volumes " + sUni + "\n"
    sKey = "*VCRIT"
    nLen = 7

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Z-Factors ----------------------------------------------

    sCom = "**  Critical Z-Factors " + "\n"
    sKey = "*ZCRIT"
    nLen = 7
    sFor = "{:8.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("ZC",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Acentric Factors ------------------------------------------------

    sCom = "**  Acentric Factors " + "\n"
    sKey = "*AC"
    nLen = 7
    sFor = "{:8.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("AF",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Parachors -------------------------------------------------------    

    sCom = "**  Parachors " + "\n"
    sKey = "*PCHOR"
    nLen = 7
    sFor = "{:8.3f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("PA",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)
    
#-- Omega-A's -------------------------------------------------------

    sCom = "**  Omega-A Values " + "\n"
    sKey = "*OMEGA"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MA",iC)*clsEOS.OA

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Omega-B's -------------------------------------------------------

    sCom = "**  Omega-B Values " + "\n"
    sKey = "*OMEGB"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MB",iC)*clsEOS.OB

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Volume Shifts ---------------------------------------------------

    sCom = "**  Volume Shifts " + "\n"
    sKey = "*VSHIFT"
    nLen = 7
    sFor = "{:9.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("SS",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Binary Interaction Parameters -----------------------------------

    sCom = "** Binary Interaction Parameters\n"
    fSim.write(sCom)
    fSim.write("\n")
    sOut = "*BIN  "
    
    for iC in range(1,nCom) :
        for jC in range(iC) :
            dVal = clsEOS.gIJ(iC,jC) ; sVal = "{:9.6f}  ".format(dVal)
            sOut = sOut + sVal
        sOut = sOut + "\n"
        fSim.write(sOut)
        sOut = "      "

    fSim.write("\n")

#-- Reservoir Temperature -------------------------------------------

    tRes = clsCMP.Tres

    if OutU == "MET" :
        sUni = ": [degC]"
        tRes = clsUNI.I2X(tRes,"degC")
    else             :
        sUni = ": [degF]"
        tRes = clsUNI.I2X(tRes,"degF")
    sTr = "{:10.3f}\n".format(tRes)
    sLabl = "** Reservoir Temperature " + sUni + "\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sOut = "*TRES  " + sTr
    fSim.write(sOut)
    fSim.write("\n")

#-- Brine Density at Standard (Stock Tank) Conditions ---------------    

    CW.calcPVTW(clsCMP,clsIO,clsUNI)  #-- Brine Properties

    dSTW = clsCMP.dSTW

#-- Aqueous Phase Properties ----------------------------------------    

    sLabl = "** Aqueous Phase Properties\n"
    fSim.write(sLabl)
    fSim.write("\n")

    pRefW = clsCMP.pRefW    #-- Ref Pres
    bRefW = clsCMP.bRefW    #-- Bw = Vres/Vsur = RhoSur/RhoRes
    cRefW = clsCMP.cRefW    #-- Compressibility
    bSalt = clsCMP.bSalt

    dRes  = dSTW/bRefW      #-- Reservoir Density

    if OutU == "MET" :
        dSTW  = clsUNI.I2X(dSTW ,"kg/m3")
        dRes  = clsUNI.I2X(dRes ,"kg/m3")
        pRefW = clsUNI.I2X(pRefW,"kpa")
        cRefW = clsUNI.I2X(cRefW,"1/kpa")

    sDsur = "{:10.3f}\n".format(dSTW)
    sDres = "{:10.3f}\n".format(dRes)
    sPref = "{:10.3f}\n".format(pRefW)
    sCref = "{:10.3e}\n".format(cRefW)
    sSalt = "{:10.5f}\n".format(bSalt)

    fSim.write("*DENW   " + sDsur)
    fSim.write("*DENWS  " + sDres)
    fSim.write("*CW     " + sCref)
    fSim.write("*REFPW  " + sPref)
    fSim.write("\n")

    fSim.write("*AQUEOUS-VISCOSITY *KESTIN\n")
    fSim.write("*SALINITY          *WTFRAC  " + sSalt)
    fSim.write("\n")

#== No return value ===================================================

    return

#========================================================================
#  Write VIP-COMP or Nexus-COMP Output
#========================================================================

def writeVIP(clsEOS,dicSAM,clsCMP,clsIO,clsUNI) :

#-- Initialisation --------------------------------------------------    

    sCom = "C "
    nCom = clsEOS.NC

    dVec = NP.zeros(nCom)

#-- Open the VIP File, if not already open -------------------------

    if not clsIO.qVIP :
        pathVIP = clsIO.patR + ".vip"
        fVIP    = open(pathVIP,'w')
        clsIO.setQVIP(True)
        clsIO.setFVIP(fVIP)
    else :
        fVIP = clsIO.fVIP

#== Headers ===========================================================        

    print("Writing Compositional Description for Tempest-MORE")

    fSim = clsIO.fVIP
    sInp = clsIO.fInP.name

    OutU = clsCMP.OutU

    if OutU == "MET" : sUni = "METRIC"
    else             : sUni = "FIELD "

#-- Write Header Information to Simulator File ----------------------        

    WO.outputHeader(fSim,sCom,clsIO)
    WO.outputEOS(fSim,sCom,clsIO,clsEOS)

#-- Header ----------------------------------------------------------    

#               123456789012345678901234567890123456789012345678901234567890
    sHead = "C ==========================================================\n"
    fSim.write(sHead)
    sLabl = "C   VIP/Nexus EOS Model generated by PVTfree Program\n"
    fSim.write(sLabl)
    sLabl = "C   From dataset " + sInp + "\n"
    fSim.write(sLabl)
    sLabl = "C   User specified " + sUni + " Units\n"
    fSim.write(sLabl)
    fSim.write(sHead)
    fSim.write("\n")

#-- Equation of State -----------------------------------------------    

    sLabl = "C  Equation of State: Peng-Robinson (PR) or Soave-Redlich-Kwong (SRK)\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sEoS = clsEOS.EOS
    sOut = "EOS  " + sEoS + "  1  \n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Components ------------------------------------------------------

    sLabl = "C  Component Names\n"
    fSim.write(sLabl)
    fSim.write("\n")
    fSim.write("COMPONENTS\n")
    sOut = "  "
    for iC in range(nCom) :
        sNam = clsEOS.gPP("CN",iC)
        sOut = sOut + sNam + "  "
    sOut = sOut + "\n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Main Properties Table -------------------------------------------

    sLabl = "C  Main Fluid Properties\n"
    fSim.write(sLabl)
    fSim.write("\n")

    if OutU == "MET" :
        sOut = "PROPERTIES  C  KPA\n"
    else :
        sOut = "PROPERTIES  F  PSIA\n"
    fSim.write(sOut)

    sOut = "COMP       MW          TC          PC        ZC       ACENTRIC      OMEGAA      OMEGAB      VSHIFT      PCHOR\n"
    fSim.write(sOut)
    
    for iC in range(nCom) :
        sN = clsEOS.gPP("CN",iC) + "  "
        Mw = clsEOS.gPP("MW",iC) ; sMw = "{:10.3f}  ".format(Mw)
        Tc = clsEOS.gPP("TC",iC)
        Pc = clsEOS.gPP("PC",iC)
        if OutU == "MET" :
            Tc = clsUNI.I2X(Tc,"degc")
            Pc = clsUNI.I2X(Pc,"kpa" )
        else :
            Tc = clsUNI.I2X(Tc,"degf")
        sTc = "{:10.3f}  ".format(Tc)
        sPc = "{:10.3f}  ".format(Pc)
        Zc = clsEOS.gPP("ZC",iC) ; sZc = "{:10.5f}  ".format(Zc)
        AF = clsEOS.gPP("AF",iC) ; sAF = "{:10.5f}  ".format(AF)
        MA = clsEOS.gPP("MA",iC) ; sOA = "{:10.7f}  ".format(MA*clsEOS.OA)
        MB = clsEOS.gPP("MB",iC) ; sOB = "{:10.7f}  ".format(MB*clsEOS.OB)
        SS = clsEOS.gPP("SS",iC) ; sSS = "{:10.7f}  ".format(SS)
        PA = clsEOS.gPP("PA",iC) ; sPA = "{:10.3f}  ".format(PA)
        sOut = sN + sMw + sTc + sPc + sZc + sAF + sOA + sOB + sSS + sPA + "\n"
        fSim.write(sOut)

    fSim.write("\n")

#-- BInary Interaction Parameters -----------------------------------    

    sLabl = "C  Binary Interaction Parameters\n"
    fSim.write(sLabl)
    fSim.write("\n")

    for iC in range(1,nCom) :

        sCom = clsEOS.gPP("CN",iC)

        sOut = "DJK  " + sCom + "\n"
        fSim.write(sOut)

        for jC in range(iC) :
            sCom = clsEOS.gPP("CN",jC)
            dVal = clsEOS.gIJ(iC,jC)
            sVal = "  {:10.5f}\n".format(dVal)
            sOut = sCom + sVal
            fSim.write(sOut)

    fSim.write("\n")

#-- End of EOS Section ----------------------------------------------

    sOut = "C  End of EOS Section\n"
    fSim.write(sOut)
    fSim.write("\n")
    
    fSim.write(sHead)
    fSim.write("ENDEOS\n")
    fSim.write(sHead)
    fSim.write("\n")

#== Water Properties ==================================================

    fSim.write("C\n")
    fSim.write("C  Water Properties\n")
    fSim.write("C\n")
    fSim.write("\n")

    CW.calcPVTW(clsCMP,clsIO,clsUNI)  #-- Brine Properties

    pRefW = clsCMP.pRefW    #-- Ref Pres
    dSTW  = clsCMP.dSTW     #-- Stock Tank Density
    bRefW = clsCMP.bRefW    #-- Ref Bw
    cRefW = clsCMP.cRefW    #-- Ref Comp
    uRefW = clsCMP.uRefW    #-- Ref Visc
    vRefW = clsCMP.vRefW    #-- Viscosibility

    if OutU == "MET" :
        pRefW = clsUNI.I2X(pRefW,"kpa")
        dSTW  = clsUNI.I2X(dSTW ,"kg/m3")
        cRefW = clsUNI.I2X(cRefW,"1/kpa")
        vRefW = clsUNI.I2X(vRefW,"1/kpa")

    vRefW = uRefW*vRefW     #-- d[Visc]/dp = Visc*Viscosibility

    sPref = "{:10.3f}  ".format(pRefW)
    sDsur = "{:10.3f}  ".format(dSTW)
    sBref = "{:10.5f}  ".format(bRefW)
    sCref = "{:10.3e}  ".format(cRefW)
    sUref = "{:10.5f}  ".format(uRefW)
    sVref = "{:10.3e}  ".format(vRefW)

#-- Write Water Properties ------------------------------------------

    fSim.write("PVTW  IPVTW  PBASEW       DWB          BWI        CW           VW        VWP\n")
    sOut = "      1      " + sPref + sDsur + sBref + sCref + sUref + sVref + "\n"
    fSim.write(sOut)
    fSim.write("\n")

#== No return values ==================================================

    return

#========================================================================
#  Write Tempest-MORE Output
#========================================================================

def writeMORE(clsEOS,dicSAM,clsCMP,clsIO,clsUNI) :

#-- Initialisation --------------------------------------------------    

    sCom = "--"
    nCom = clsEOS.NC

    dVec = NP.zeros(nCom)

#-- Open the MORE File, if not already open ------------------------

    if not clsIO.qMOR :
        pathMOR = clsIO.patR + ".mor"
        fMOR    = open(pathMOR,'w')
        clsIO.setQMOR(True)
        clsIO.setFMOR(fMOR)
    else :
        fMOR = clsIO.fMOR

#== Headers ===========================================================        

    print("Writing Compositional Description for Tempest-MORE")

    fSim = clsIO.fMOR
    sInp = clsIO.fInP.name

    OutU = clsCMP.OutU

    if OutU == "MET" : sUni = "METRIC"
    else             : sUni = "FIELD "

#-- Write Header Information to Simulator File ----------------------        

    WO.outputHeader(fSim,sCom,clsIO)
    WO.outputEOS(fSim,sCom,clsIO,clsEOS)

#-- Header ----------------------------------------------------------    

#               123456789012345678901234567890123456789012345678901234567890
    sHead = "--==========================================================\n"
    fSim.write(sHead)
    sLabl = "--  Tempest-MORE EOS Model generated by PVTfree Program\n"
    fSim.write(sLabl)
    sLabl = "--  From dataset " + sInp + "\n"
    fSim.write(sLabl)
    sLabl = "--  User specified " + sUni + " Units\n"
    fSim.write(sLabl)
    fSim.write(sHead)
    fSim.write("\n")

#----------------------------------------------------------------------
#  MORE-Specific Data
#----------------------------------------------------------------------

    sLabl = "-- Component Names (CNAM) are specified in the INPU section\n"
    fSim.write(sLabl)
    fSim.write("\n")

    sOut = "CNAM  "
    for iC in range(nCom) :
        sOut = sOut + clsEOS.gPP("CN",iC) + "  "
    sOut = sOut + "WATR\n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Composition -----------------------------------------------------

    sOut = "-- Compositions\n"
    fSim.write(sOut)
    fSim.write("\n")

    nSam = len(dicSAM)

    for iSam in range(nSam) :

        clsSAM = dicSAM[iSam]

        sNam = clsSAM.sName

        sCom = ""
        sKey = "SCMP  " + sNam + "\n"
        nLen = 6
        sFor = "{:10.7f}"
        for iC in range(nCom) : dVec[iC] = clsSAM.gZI(iC)

        writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Fluid Section Header --------------------------------------------    

    fSim.write(sHead)
    sLabl = "FLUI EOS\n"
    fSim.write(sLabl)
    fSim.write(sHead)
    fSim.write("\n")

    sLabl = "-- Equation of State\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sEoS = clsEOS.EOS
    sOut = "EQUA  " + sEoS + "\n"
    fSim.write(sOut)
    fSim.write("\n")

#-- Reservoir Temperature -------------------------------------------

    tRes = clsCMP.Tres

    if OutU == "MET" :
        sUni = ": [degC]"
        tRes = clsUNI.I2X(tRes,"degC")
    else             :
        sUni = ": [degF]"
        tRes = clsUNI.I2X(tRes,"degF")
    sTr = "{:10.3f}\n".format(tRes)
    sLabl = "-- Reservoir Temperature " + sUni + "\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sOut = "TEMP  " + sTr
    fSim.write(sOut)
    fSim.write("\n")

#-- Main Properties Table -------------------------------------------

    sLabl = "-- Main Fluid Properties\n"
    fSim.write(sLabl)
    fSim.write("\n")

    sOut = "PROP  CNam       Mw          Tc          Pc       AcF         Zc          SG         Para\n"
    fSim.write(sOut)

    if OutU == "MET" :
        sOut = "--                          Kelv        bara\n"
    else :
        sOut = "--                          degR        psia\n"
    fSim.write(sOut)
    
    for iC in range(nCom) :
        sN = clsEOS.gPP("CN",iC) + "  "
        Mw = clsEOS.gPP("MW",iC) ; sMw = "{:10.3f}  ".format(Mw)
        Tc = clsEOS.gPP("TC",iC)
        Pc = clsEOS.gPP("PC",iC)
        if OutU == "MET" :
            Tc = clsUNI.I2X(Tc,"kelv")
            Pc = clsUNI.I2X(Pc,"bara") ; sPc = "{:10.4f}".format(Pc)
        else :
            sPc = "{:10.3f}  ".format(Pc)
        sTc = "{:10.3f}  ".format(Tc)
        AF = clsEOS.gPP("AF",iC) ; sAF = "{:10.5f}  ".format(AF)
        Zc = clsEOS.gPP("ZC",iC) ; sZc = "{:10.5f}  ".format(Zc)
        SG = clsEOS.gPP("SG",iC) ; sSG = "{:10.5f}  ".format(SG)
        PA = clsEOS.gPP("PA",iC) ; sPA = "{:10.3f}  ".format(PA)
        sOut = "      " + sN + sMw + sTc + sPc + sAF + sZc + sSG + sPA + "\n"
        fSim.write(sOut)

    fSim.write("/\n")
    fSim.write("\n")

#-- Omega-A & Omega-B Multiplers -----------------------------------    

    sCom = "--  Omega-A Multipliers " + "\n"
    sKey = "OMGA  MULT\n"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MA",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

    sCom = "--  Omega-B Multipliers " + "\n"
    sKey = "OMGB  MULT\n"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MB",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Volume Shift Parameters -----------------------------------------    

    sCom = "--  Volume Shifts " + "\n"
    sKey = "VOLU\n"
    nLen = 6
    sFor = "{:11.8f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("SS",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- BIPs ------------------------------------------------------------    

    fSim.write("-- Binary Iteraction Parameters (Lower Triangle)\n")
    fSim.write("\n")
    fSim.write("INTE\n")

    sCom = ""
    sKey = "ROW"
    nLen = 7
    sFor = "{:9.6f}"

    for iC in range(1,nCom) :
        bVec = NP.zeros(iC)
        nCom = iC
        for jC in range(iC) :
            bVec[jC] = clsEOS.gIJ(iC,jC)
        writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,bVec)

    fSim.write("/\n")
    fSim.write("\n")

#-- Brine Density at Standard (Stock Tank) Conditions ---------------    

    mFrc = clsCMP.bSalt #-- Mass Fraction of Salt in the Brine

    dSTW,comW = CW.calcRoweChouDen(mFrc,CO.tStand,CO.pStand,clsUNI)  #-- Stock Tank Brine Density

    clsCMP.setDenSTW(dSTW)

    clsCMP.setDenSTO(-1.0)  #-- In EoS Mode, E300 calculates STO Density
    clsCMP.setDenSTG(-1.0)  #-- Ditto

#-- Brine Properties at Reference Pressure & Tres (for PVTW keyword) 

    CW.calcPVTW(clsCMP,clsIO,clsUNI)

#== Write DENSITY and PVTW keywords ===================================

    WO.outputECLDensity(OutU,fSim,clsCMP,clsUNI)
    WO.outputECLPVTW(OutU,fSim,clsCMP,clsUNI)

#== No return value ===================================================

    return

#========================================================================
#  Write E300 Output
#========================================================================

def writeE300(clsEOS,dicSAM,clsCMP,clsIO,clsUNI) :

#-- Open the E300 File, if not already open ------------------------

    if not clsIO.q300 :
        path300 = clsIO.patR + ".e300"
        f300    = open(path300,'w')
        clsIO.setQ300(True)
        clsIO.setF300(f300)
    else :
        f300 = clsIO.f300

    sCom = "--"

#== Headers ===========================================================        

    print("Writing Compositional Description for E300")

    fSim = clsIO.f300
    sInp = clsIO.rNam

    OutU = clsCMP.OutU

    if OutU == "MET" : sUni = "METRIC"
    else             : sUni = "FIELD "

#-- Write Header Information to Simulator File ----------------------        

    WO.outputHeader(fSim,sCom,clsIO)
    WO.outputEOS(fSim,sCom,clsIO,clsEOS)

#-- Header ----------------------------------------------------------    

#               123456789012345678901234567890123456789012345678901234567890
    sHead = "--==========================================================\n"
    fSim.write(sHead)
    sLabl = "--  E300 EOS Model generated by PVTfree Program\n"
    fSim.write(sLabl)
    sLabl = "--  From dataset " + sInp + "\n"
    fSim.write(sLabl)
    sLabl = "--  User specified " + sUni + " Units\n"
    fSim.write(sLabl)
    fSim.write(sHead)
    fSim.write("\n")

#-- Number of Components --------------------------------------------
    
    nCom = clsEOS.NC
    
    sLabl = "--  Number of Components\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sLabl = "NCOMPS\n"
    fSim.write(sLabl)
    sLabl = "  {:2d}  /\n".format(nCom)
    fSim.write(sLabl)
    fSim.write("\n")

#-- Equation of State -----------------------------------------------    

    EOS = clsEOS.EOS
    if EOS == "SRK" : sExt = "Soave-Redlich-Kwong (SRK)"
    else            : sExt = "Peng-Robinson (PR)"
    sLabl = "--  Equation of State: " + sExt + "\n"
    fSim.write(sLabl)
    fSim.write("\n")
    sLabl = "EOS\n"
    fSim.write(sLabl)
    sLabl = "  " + str(EOS) + "  /\n"
    fSim.write(sLabl)
    fSim.write("\n")

    if EOS == "PR" :
        sLabl = "--  Modified Form of the Peng-Robinson EOS\n"
        fSim.write(sLabl)
        fSim.write("\n")
        sLabl = "PRCORR\n"
        fSim.write(sLabl)
        fSim.write("\n")

#== Component Properties ==============================================

    sVec = [ "" for i in range(nCom)]
    dVec = NP.zeros(nCom)

#-- Component Names -------------------------------------------------    

    sCom = "--  Component Names " + "\n"
    sKey = "CNAMES\n"
    nLen = 7
    nFor = 8
    for iC in range(nCom) : sVec[iC] = clsEOS.gPP("CN",iC)
    writeVectorS(fSim,sCom,sKey,nCom,nLen,nFor,sVec)

#-- Mole Weights ----------------------------------------------------    

    if OutU == "MET" : sUni = ": [kg/kgmol]"
    else             : sUni = ": [lb/lbmol]"
    sCom = "--  Molecular Weights " + sUni + "\n"
    sKey = "MW\n"
    nLen = 7
    sFor = "{:8.3f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MW",iC)
    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Temperatures -------------------------------------------
    
    if OutU == "MET" :
        sUni = ": [Kelvin]"
        for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("TC",iC),"kelv")
    else             :
        sUni = ": [degrees Rankine]"
        for iC in range(nCom) : dVec[iC] = clsEOS.gPP("TC",iC)
    sCom = "--  Critical Temperatures " + sUni + "\n"
    sKey = "TCRIT\n"
    nLen = 7
    sFor = "{:8.3f}"

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Pressures ----------------------------------------------    

    if OutU == "MET" :
        sUni = ": [barsa]"
        sFor = "{:8.4f}"
        for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("PC",iC),"bara")
    else             :
        sUni = ": [psia]"
        sFor = "{:8.3f}"
        for iC in range(nCom) : dVec[iC] = clsEOS.gPP("PC",iC)
    sCom = "--  Critical Pressures " + sUni + "\n"
    sKey = "PCRIT\n"
    nLen = 7

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Volumes ------------------------------------------------

    if OutU == "MET" :
        sUni = ": [m3/kgmol]"
        sFor = "{:8.3f}"
        for iC in range(nCom) : dVec[iC] = clsUNI.I2X(clsEOS.gPP("VC",iC),"m3/kgmol")
    else             :
        sUni = ": [ft3/lbmol]"
        sFor = "{:8.4f}"
        for iC in range(nCom) : dVec[iC] = clsEOS.gPP("VC",iC)
    sCom = "--  Critical Volumes " + sUni + "\n"
    sKey = "VCRIT\n"
    nLen = 7

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Critical Z-Factors ----------------------------------------------

    sCom = "--  Critical Z-Factors " + "\n"
    sKey = "ZCRIT\n"
    nLen = 7
    sFor = "{:8.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("ZC",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Acentric Factors ------------------------------------------------

    sCom = "--  Acentric Factors " + "\n"
    sKey = "ACF\n"
    nLen = 7
    sFor = "{:8.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("AF",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Omega-A's -------------------------------------------------------

    sCom = "--  Omega-A Values " + "\n"
    sKey = "OMEGAA\n"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MA",iC)*clsEOS.OA

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Omega-B's -------------------------------------------------------

    sCom = "--  Omega-B Values " + "\n"
    sKey = "OMEGAB\n"
    nLen = 6
    sFor = "{:11.9f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("MB",iC)*clsEOS.OB

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Parachors -------------------------------------------------------    

    sCom = "--  Parachors " + "\n"
    sKey = "PARACHOR\n"
    nLen = 7
    sFor = "{:8.3f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("PA",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)
    
#-- Volume Shifts ---------------------------------------------------

    sCom = "--  Volume Shifts " + "\n"
    sKey = "SSHIFT\n"
    nLen = 7
    sFor = "{:9.6f}"
    for iC in range(nCom) : dVec[iC] = clsEOS.gPP("SS",iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#== Binary Interation Coefficients ====================================

    fSim.write("-- Binary Iteraction Parameters\n")
    fSim.write("\n")
    fSim.write("BIC\n")

    sCom = ""
    sKey = ""
    nLen = 7
    sFor = "{:9.6f}"

    for iC in range(1,nCom) :
        bVec = NP.zeros(iC)
        kCom = iC
        for jC in range(iC) :
            bVec[jC] = clsEOS.gIJ(iC,jC)
        writeVectorD(fSim,sCom,sKey,kCom,nLen,sFor,bVec)

    fSim.write("/\n")
    fSim.write("\n")

#-- Reservoir Temperature -------------------------------------------    

    tRes = clsCMP.Tres

    if OutU == "MET" :
        sUni = ": [degC]"
        tRes = clsUNI.I2X(tRes,"degC")
    else             :
        sUni = ": [degF]"
        tRes = clsUNI.I2X(tRes,"degF")
    sLabl = "-- Reservoir Temperature " + sUni + "\n"
    fSim.write(sLabl)
    fSim.write("\n")
    fSim.write("RTEMP\n")
    sLabl = "  {:7.3f}".format(tRes) +  "  /\n"
    fSim.write(sLabl)
    fSim.write("\n")

#-- Composition -----------------------------------------------------

    sCom = "--  Composition\n"
    sKey = "ZI\n"
    nLen = 6
    sFor = "{:10.7f}"
    for iC in range(nCom) : dVec[iC] = dicSAM[0].gZI(iC)

    writeVectorD(fSim,sCom,sKey,nCom,nLen,sFor,dVec)

#-- Brine Density at Standard (Stock Tank) Conditions ---------------    

    mFrc = clsCMP.bSalt #-- Mass Fraction of Salt in the Brine

    dSTW,comW = CW.calcRoweChouDen(mFrc,CO.tStand,CO.pStand,clsUNI)  #-- Stock Tank Brine Density

    clsCMP.setDenSTW(dSTW)

    clsCMP.setDenSTO(-1.0)  #-- In EoS Mode, E300 calculates STO Density
    clsCMP.setDenSTG(-1.0)  #-- Ditto

#-- Brine Properties at Reference Pressure & Tres (for PVTW keyword) 

    CW.calcPVTW(clsCMP,clsIO,clsUNI)

#== Write DENSITY and PVTW keywords ===================================

    WO.outputECLDensity(OutU,fSim,clsCMP,clsUNI)

    WO.outputECLPVTW(OutU,fSim,clsCMP,clsUNI)

#========================================================================
#  End of Routine
#========================================================================

    return

#========================================================================
#  Write a Vector of Float (Double) Values
#========================================================================

def writeVectorD(fSim,sComm,sKW,nCom,nLen,sForm,dVec) :

    nSETs = 0
    n1SET = nLen
    nSumm = nCom

    nSETa = []

    while (nSumm-n1SET) > 0 :
        nSETs += 1
        nSETa.append(n1SET)
        nSumm = nSumm - n1SET

    nSETs += 1
    nSETa.append(nSumm)

    iCom = 0

    if sKW != "" :
        if sKW == "ROW" :
            fSim.write("  " + sKW)
        else :
            fSim.write(sComm)
            fSim.write("\n")
            fSim.write(sKW)

    for iS in range(nSETs) :

        sLine = "  "
        nThis = nSETa[iS]

        for iC in range(nThis) :
            dThis = dVec[iCom]
            sThis = sForm.format(dThis)
            sLine = sLine + sThis + "  "
            iCom += 1

        if sKW != "" :
            if sKW == "ROW" or sKW[:1] == "*" :
                sLine = sLine + "\n"
            else :
                if   iCom == nCom : sLine = sLine + "  /\n"
                else              : sLine = sLine + "\n"
        else :
            sLine = sLine + "\n"

        fSim.write(sLine)

    if sKW != "" and sKW != "ROW" : fSim.write("\n")

#========================================================================
#  End of Routine
#========================================================================

    return

#========================================================================
#  Write a Vector of String Values
#========================================================================

def writeVectorS(fSim,sComm,sKW,nCom,nLen,nForm,sVec) :

    nSETs = 0
    n1SET = nLen
    nSumm = nCom

    nSETa = []

    while (nSumm-n1SET) > 0 :
        nSETs += 1
        nSETa.append(n1SET)
        nSumm = nSumm - n1SET

    nSETs += 1
    nSETa.append(nSumm)

    iCom = 0

    fSim.write(sComm)
    fSim.write("\n")
    fSim.write(sKW)

    for iS in range(nSETs) :

        sLine = "  "
        nThis = nSETa[iS]

        for iC in range(nThis) :
            sThis = sVec[iCom]
            sThis = sThis.center(nForm)
            sLine = sLine + sThis + "  "
            iCom += 1

        if iCom == nCom : sLine = sLine + "  /\n"
        else            : sLine = sLine + "\n"

        fSim.write(sLine)

    fSim.write("\n")

#========================================================================
#  End of Routine
#========================================================================

    return

#========================================================================
#  End of Module
#========================================================================

