#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Component Properties Module
#  
#========================================================================

#!/usr/bin/python3

from   copy     import deepcopy
from   math     import log, exp, sqrt

import numpy         as NP
import scipy.special as SS

import calcEOS    as CE
import constants  as CO
import calcPhsPlt as CP
import stabTest   as ST
import writeOut   as WO

#== Component Types: Library (L) or Single Carbon Number [SCN] (L) ======

typeC = {'N2' :'L','CO2':'L','H2S':'L','C1' :'L','C2' :'L','C3' :'L',\
         'IC4':'L','NC4':'L','IC5':'L','NC5':'L',\
         'C6' :'S','C7' :'S','C8' :'S','C9' :'S','C10':'S',\
         'C11':'S','C12':'S','C13':'S','C14':'S','C15':'S',\
         'C16':'S','C17':'S','C18':'S','C19':'S','C20':'S',\
         'C21':'S','C22':'S','C23':'S','C24':'S','C25':'S',\
         'C26':'S','C27':'S','C28':'S','C29':'S','C30':'S',\
         'C31':'S','C32':'S','C33':'S','C34':'S','C35':'S',\
         'C36':'S','C37':'S','C38':'S','C39':'S','C40':'S',\
         'C41':'S','C42':'S','C43':'S','C44':'S','C45':'S' }

#----------------------------------------------------------------------
#  Library Components (N2, CO2, H2S and C1 thru NC5)
#                         Units
#  Mole Weights           lb/lbmol
#  Critical Temperatures  degF
#  Critical Pressures     psia
#  Acentric Factors       dim
#  Critical Volumes       ft3/lbmol
#  Specific Gravities     dim
#  Boiling Points         degF
#  Parachors              dim
#  Ideal Gas Specific Heat Coefficients
#  CpA                    Btu/lbmol.degF
#  CpB                    Btu/lbmol.degF^2
#  CpC                    Btu/lbmol.degF^3
#  CpD                    Btu/lbmol.degF^4
#----------------------------------------------------------------------

lMolW = {'N2' :  28.013,'CO2':  44.010,'H2S':  34.076,\
         'C1' :  16.043,'C2' :  30.070,'C3' :  44.097,\
         'IC4':  58.124,'NC4':  58.124,'IC5':  72.151,'NC5': 72.151}
lTcrt = {'N2' :-232.400,'CO2':  87.900,'H2S': 212.700,\
         'C1' :-116.630,'C2' :  90.090,'C3' : 206.010,\
         'IC4': 274.980,'NC4': 305.650,'IC5': 369.100,'NC5':385.700}
lPcrt = {'N2' : 493.000,'CO2':1070.600,'H2S':1306.000,\
         'C1' : 667.800,'C2' : 707.800,'C3' : 616.300,\
         'IC4': 529.100,'NC4': 550.700,'IC5': 490.400,'NC5':488.600}
lAcnF = {'N2' :   0.045,'CO2':   0.231,'H2S':   0.100,\
         'C1' :   0.012,'C2' :   0.091,'C3' :   0.145,\
         'IC4':   0.176,'NC4':   0.193,'IC5':   0.227,'NC5':  0.251}
lVcrt = {'N2' :   1.443,'CO2':   1.505,'H2S':   1.564,\
         'C1' :   1.590,'C2' :   2.370,'C3' :   3.250,\
         'IC4':   4.208,'NC4':   4.080,'IC5':   4.899,'NC5':  4.870}
lSpcG = {'N2' :   0.470,'CO2':   0.507,'H2S':   0.500,\
         'C1' :   0.330,'C2' :   0.450,'C3' :   0.508,\
         'IC4':   0.563,'NC4':   0.584,'IC5':   0.625,'NC5':  0.631}
lTboi = {'N2' :-320.400,'CO2':-109.300,'H2S': -76.600,\
         'C1' :-258.690,'C2' :-127.480,'C3' : -43.670,\
         'IC4':  10.900,'NC4':  31.100,'IC5':  82.120,'NC5': 96.920}
lPara = {'N2' :  41.000,'CO2':  70.000,'H2S':  41.000,\
         'C1' :  77.000,'C2' : 108.000,'C3' : 150.300,\
         'IC4': 181.500,'NC4': 189.900,'IC5': 225.000,'NC5':231.500}

lCpAI = {'N2' : 7.45E+00,'CO2': 4.73E+00,'H2S': 7.63E+00,\
         'C1' : 4.60E+00,'C2' : 1.29E+00,'C3' :-1.01E+00,\
         'IC4':-3.32E-01,'NC4': 2.27E+00,'IC5':-2.28E+00,'NC5':-8.67E-01}
lCpBI = {'N2' :-1.80E-03,'CO2': 9.75E-03,'H2S': 1.91E-04,\
         'C1' : 6.92E-03,'C2' : 2.37E-02,'C3' : 4.07E-02,\
         'IC4': 5.11E-02,'NC4': 4.40E-02,'IC5': 6.73E-02,'NC5': 6.47E-02}
lCpCI = {'N2' : 1.98E-06,'CO2':-4.13E-06,'H2S': 1.79E-06,\
         'C1' : 8.83E-07,'C2' :-5.12E-06,'C3' :-1.17E-05,\
         'IC4':-1.36E-05,'NC4':-8.17E-06,'IC5':-2.01E-05,'NC5':-1.90E-05}
lCpDI = {'N2' :-4.79E-10,'CO2': 7.01E-10,'H2S':-4.82E-10,\
         'C1' :-4.64E-10,'C2' : 3.57E-10,'C3' : 1.32E-09,\
         'IC4': 1.19E-09,'NC4':-1.16E-10,'IC5': 2.35E-09,'NC5': 2.17E-09}

#-- Volume Shifts: Peng-Robinson (PR) and (Soave) Redlich-Kwong (RK)

lVSPR = {'N2' :-0.1927,'CO2':-0.0817,'H2S':-0.1288,\
         'C1' :-0.1595,'C2' :-0.1134,'C3' :-0.0863,\
         'IC4':-0.0844,'NC4':-0.0675,'IC5':-0.0608,'NC5':-0.0390}
lVSRK = {'N2' :-0.0079,'CO2': 0.0833,'H2S': 0.0466,\
         'C1' : 0.0234,'C2' : 0.0605,'C3' : 0.0825,\
         'IC4': 0.0830,'NC4': 0.0975,'IC5': 0.1022,'NC5': 0.1209}

#-- BICs for N2-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

lKN2P = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.130,\
         'C1' : 0.025,'C2' : 0.010,'C3' : 0.090,\
         'IC4': 0.095,'NC4': 0.095,'IC5': 0.100,'NC5': 0.100}
lKN2S = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.120,\
         'C1' : 0.020,'C2' : 0.060,'C3' : 0.080,\
         'IC4': 0.080,'NC4': 0.080,'IC5': 0.080,'NC5': 0.080}

#-- BICs for CO2-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

lKCOP = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.135,\
         'C1' : 0.105,'C2' : 0.130,'C3' : 0.125,\
         'IC4': 0.120,'NC4': 0.115,'IC5': 0.115,'NC5': 0.115}
lKCOS = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.120,\
         'C1' : 0.120,'C2' : 0.150,'C3' : 0.150,\
         'IC4': 0.150,'NC4': 0.150,'IC5': 0.150,'NC5': 0.150}

#-- BICs for H2S-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

lKHSP = {'N2' : 0.130,'CO2': 0.135,'H2S': 0.000,\
         'C1' : 0.070,'C2' : 0.085,'C3' : 0.080,\
         'IC4': 0.075,'NC4': 0.075,'IC5': 0.070,'NC5': 0.070}
lKHSS = {'N2' : 0.120,'CO2': 0.120,'H2S': 0.000,\
         'C1' : 0.080,'C2' : 0.070,'C3' : 0.070,\
         'IC4': 0.060,'NC4': 0.060,'IC5': 0.060,'NC5': 0.060}

#-- BICs for C1-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

lKC1P = {'N2' : 0.025,'CO2': 0.105,'H2S': 0.070,\
         'C1' : 0.000,'C2' : 0.000,'C3' : 0.000,\
         'IC4': 0.000,'NC4': 0.000,'IC5': 0.000,'NC5': 0.000}
lKC1S = {'N2' : 0.020,'CO2': 0.120,'H2S': 0.080,\
         'C1' : 0.000,'C2' : 0.000,'C3' : 0.000,\
         'IC4': 0.000,'NC4': 0.000,'IC5': 0.000,'NC5': 0.000}
#----------------------------------------------------------------------
#  SCN Mole Weights [SCN's being C6 thru C45]
#----------------------------------------------------------------------

scnMW = {'C6' : 84.0,'C7' : 96.0,'C8' :107.0,'C9' :121.0,'C10':134.0,\
         'C11':147.0,'C12':161.0,'C13':175.0,'C14':190.0,'C15':206.0,\
         'C16':222.0,'C17':237.0,'C18':251.0,'C19':263.0,'C20':275.0,\
         'C21':291.0,'C22':300.0,'C23':312.0,'C24':324.0,'C25':337.0,\
         'C26':349.0,'C27':360.0,'C28':372.0,'C29':382.0,'C30':394.0,\
         'C31':404.0,'C32':415.0,'C33':426.0,'C34':437.0,'C35':445.0,\
         'C36':456.0,'C37':464.0,'C38':475.0,'C39':484.0,'C40':495.0,\
         'C41':502.0,'C42':512.0,'C43':521.0,'C44':531.0,'C45':539.0 }

#========================================================================
#  Define Properties of the User Defined Components
#========================================================================

def allProps(clsEOS,dicSAM,clsIO,clsUNI) :

    iERR = 0

    nComp = clsEOS.NC
    nUser = clsEOS.NU
    nSplt = clsEOS.NP
    nSamp = clsEOS.NS

    print("allProps: Attempting to assign Component Properties")

#-- Find the type (Library or SCN) of first Nc-1 components ---------    

    for iC in range(nUser-1) :

        cName = clsEOS.gPP("CN",iC)
        
        typC  = typeC.get(cName)

        if typC == None :
            print("allProps: Component ",iC+1,cName," is not Library or SCN - Error")
            iERR = -1
            return iERR

        #print("iC,cNam,typC ",iC,cName,typC)

#== Library Components ================================================        

        if   typC == 'L' :

            isLibComp(iC,cName,clsEOS)

#== SCN Components ====================================================            

        elif typC == 'S' :  #-- SCN     Component

            molWt = scnMW[cName]

            heavyComp(molWt,-1.0,iC,clsEOS)

            print("Component ({:2d}) = {:4s} is assumed to be a SCN - Properties Assigned/Calculated".format(iC+1,cName))

#== Not Library or SCN ================================================
            
        else :              

            print("Component Name ",cName," Not Recognised - Error")
            iERR = -1            

#=======================================================================
#  Process the Plus Fraction(s)
#  If more than one sample, split using Whitson Quadrature Method
#=======================================================================

#-- Mole Weight of Penultimate Component: 'Last' component = nComp-1 --

    cPlus = clsEOS.gPP("CN",nUser-1)     #-- Name of the User Plus Fraction

    if nSamp > 1 or nSplt > 1 :
        
        print("Component ({:2d}) = {:4s} is Plus Fraction to be Split into {:1d} Pseudos".format(nUser,cPlus,nSplt))

        iERR = splitPlus(dicSAM,clsEOS)

        if iERR < 0 :
            return iERR

        print("Plus Fraction Splitting Completed")

    else :

        molWt = dicSAM[0].mPlus
        specG = dicSAM[0].sPlus

        print("Component ({:2d}) = {:4s} is Plus Fraction with No Split Requested - Properties Calculated".format(nUser+1,cPlus))

        heavyComp(molWt,specG,nUser-1,clsEOS)

#=======================================================================
#  Sort the Binary Interaction Parameters
#=======================================================================

    nComp = clsEOS.NC

    #print("compProps: Set BICs, N = ",nComp)

    assignBIPs(clsEOS)

    print("Binary Interaction Parameters (BIPs) Assigned/Calculated")

#========================================================================
#  Do we need to sort the components?  Most to Least Volatility
#========================================================================

    clsEOS,dicSAM = sortComponents(clsEOS,dicSAM)

#========================================================================
#  Output the data
#========================================================================

    sTit = "Initialisation"

    WO.outputProps(clsIO,clsEOS,dicSAM,sTit)

#========================================================================
#  Write Fluid Description to the SAV file
#========================================================================

    WO.outputSave(sTit,clsEOS,dicSAM,clsIO)

#======================================================================
#  Generate (Approximate) Phase Plots
#======================================================================

    CP.allSamplesPhasePlot(clsEOS,dicSAM,clsIO)

#======================================================================
#  End of Routine
#======================================================================

    return iERR

#======================================================================
#  Is a Library Component
#======================================================================

def isLibComp(iC,cName,clsEOS) :

    Pcrit = lPcrt[cName]
    Vcrit = lVcrt[cName]
    Tcrit = lTcrt[cName] + CO.dF2dR            #-- Want this stored in degR

    Zcrit = Pcrit*Vcrit/(CO.gasCon*Tcrit)

    clsEOS.sPP("CN",iC,cName)

    clsEOS.sPP("MW",iC,lMolW[cName])
    clsEOS.sPP("TC",iC,Tcrit)
    clsEOS.sPP("PC",iC,Pcrit)
    clsEOS.sPP("VC",iC,Vcrit)
    clsEOS.sPP("ZC",iC,Zcrit)
    clsEOS.sPP("AF",iC,lAcnF[cName])
    clsEOS.sPP("SG",iC,lSpcG[cName])
    clsEOS.sPP("TB",iC,lTboi[cName] + CO.dF2dR)  #-- Want this stored in degR
    clsEOS.sPP("PA",iC,lPara[cName])
    clsEOS.sPP("MA",iC,1.0)
    clsEOS.sPP("MB",iC,1.0)

    if clsEOS.EOS == "SRK" :
        clsEOS.sPP("SS",iC,lVSRK[cName])
    else :
        clsEOS.sPP("SS",iC,lVSPR[cName])

    clsEOS.sPP("CA",iC,lCpAI[cName])
    clsEOS.sPP("CB",iC,lCpBI[cName])
    clsEOS.sPP("CC",iC,lCpCI[cName])
    clsEOS.sPP("CD",iC,lCpDI[cName])

    print("Component ({:2d}) = {:4s} is in Internal Library - Properties Assigned".format(iC+1,cName))

#======================================================================
#  End of Routine
#======================================================================

    return

#=======================================================================
#  Split the Plus Fraction into Pseudo-Components
#  Reference: Ingolf Soreide PhD Thesis
#  http://www.ipt.ntnu.no/~curtis/courses/Theses/Ingolf-Soreide-NTH-PhD-Thesis.pdf
#  Section 3.2.6
#=======================================================================

def splitPlus(dicSAM,clsEOS) :

    iERR = 0

    nUser = clsEOS.NU  #-- Component Number of the User Plus Fraction
    nComp = clsEOS.NC  #-- Component Number of Heaviest Pseudo
    nSamp = clsEOS.NS  #-- Number of Samples
    nSplt = clsEOS.NP  #-- Number of Pseudo-Components (=NC-NU+1)

    #print("splitPlus: nUser,nComp,nSamp,nSplt ",nUser,nComp,nSamp,nSplt)

    cPlus = clsEOS.gPP("CN",nUser-1)                 #-- Name of the User Plus Fraction
    cPlus = procPlusFracName(nUser-1,cPlus,clsEOS)   #-- If already split, process-name

    #print("splitPlus: cPlus ",cPlus)

#----------------------------------------------------------------------
#  Gauss-Laguerre Zeros (xZ) and Weights (wT) depending on Num Pseudos
#  Now using Scipy Special Function to get Zeros & Weights
#----------------------------------------------------------------------
    
    xZ,wT = SS.roots_laguerre(nSplt)

#-- Whitson Rule: Heaviest Pseudo = MIN(500,2.5*max-Plus-Frac) ------            

    maxMw = 0.0
    for iS in range(nSamp) :
        maxMw = max(maxMw,dicSAM[iS].mPlus)

    maxMw = min(2.5*maxMw,500.0)

#-- eta = Min Mole Wt in Plus Frac, i.e. if C7+ => Mw(C6+1/2) -------
        
    eta   = clsEOS.gPP("MW",nUser-2) + 6.0  #-- Penultimate Comp + 6.0 A.U.
    beta0 = (maxMw - eta)/xZ[nSplt-1]       #-- Eqn.(3.33)

#--------------------------------------------------------------------
#  alfa - User alpha-Coeff (defaulted to 1.0 if not sepcified)
#  mObs - User Plus Fraction Mole Weight
#  delt - Eqn.(3.34)
#  gama - The Gamma Function - now use Scipy Special Function
#--------------------------------------------------------------------

    alfa = NP.zeros(nSamp)
    mObs = NP.zeros(nSamp)
    delt = NP.zeros(nSamp)
    gama = NP.zeros(nSamp)

    for iSam in range(nSamp) :
        clsSAM = dicSAM[iSam]
        if clsSAM.aPlus < 1.0E-06 : clsSAM.aPlus = 1.0
        alfa[iSam] = clsSAM.aPlus
        mObs[iSam] = clsSAM.mPlus
        delt[iSam] = exp(alfa[iSam]*beta0/(mObs[iSam] - eta) - 1.0)
        gama[iSam] = SS.gamma(alfa[iSam])

#--------------------------------------------------------------------
#  mSpl - Pseudo-Component Mole Weights, Eqn.(3.32)
#  fCon - "Constant" part (no delt-dependence) of Eqn.(3.30)
#--------------------------------------------------------------------

    mSpl = NP.zeros( nSplt)
    fCon = NP.zeros((nSamp,nSplt))

    for iSpl in range(nSplt) :
        mSpl[iSpl] = eta + beta0*xZ[iSpl]
        for iSam in range(nSamp) :
            expO = (alfa[iSam]-1.0)/gama[iSam]
            fCon[iSam][iSpl] = pow(xZ[iSpl],expO)

    #print("mSpl ",mSpl)
    #print("fCon ",fCon)

#-- Adjust the delta's to match the Measured Mole Weight by Sample --

    delt,zSpl = matchAllDelta(alfa,delt,mObs,mSpl,wT,xZ,fCon)

#-- Ensure the sample plus fraction compositions normalise ----------

    sumZ = NP.zeros(nSamp)

    for iSam in range(nSamp) :
        iOffs = nUser - 1
        for iSp in range(nSplt) :
            sumZ[iSam] = sumZ[iSam] + zSpl[iSam][iSp]
        sPls = dicSAM[iSam].sComp[nUser-1]
        #print("sumZ,sPls ",sumZ[iSa],sPls)
        sPls = sPls/sumZ[iSam]
        for iSp in range(nSplt) :
            dicSAM[iSam].sZI(iOffs,sPls*zSpl[iSam][iSp])
            iOffs += 1
            #print("iSa,iSp,sTst ",iSa,iSp,sTst)

#== Now, check the specific gravities against Soreide correlation =====

    nPlus = 0
    fPlus = 0.0
    zPlus = NP.zeros(nSplt)
    gPlus = NP.zeros(nSplt)
        
    for iSam in range(nSamp) :
        sPlus = dicSAM[iSam].sPlus
        if sPlus > 0.0 :
            for iSpl in range(nSplt) :
                zPlus[iSpl] = zSpl[iSam][iSpl]
            fCThs = adjustFC(nSplt,zPlus,mSpl,sPlus)
            nPlus = nPlus + 1
            fPlus = fPlus + fCThs

    if nPlus > 0 :
        fCavg = fPlus/float(nPlus)
        for iSpl in range(nSplt) :
            gPlus[iSpl] = soreideSpcG(mSpl[iSpl],fCavg)
        #print("gPlus ",gPlus)
    else :
        print("No Plus Fraction Specific Gravities for Any of the Samples - Error")
        iERR = -1
        return iERR

#== Now can define the Pseudo-Component Properties ====================

    iC = nUser - 1
    for iSpl in range(nSplt) :
        molWt =  mSpl[iSpl]
        SpecG = gPlus[iSpl]
        heavyComp(molWt,SpecG,iC,clsEOS)
        cName = cPlus + "P" + str(iSpl+1)
        clsEOS.sPP("CN",iC,cName)
        iC += 1

#=======================================================================
#  End of Routine
#=======================================================================

    return iERR
            
#=======================================================================
#  Assign Binary Interaction Parameters (BIPs)
#=======================================================================

def assignBIPs(clsEOS) :

    nComp = clsEOS.NC

#-- Initially, set all values to Zero -------------------------------
    
    for iC in range(nComp) :
        for jC in range(nComp) :
            clsEOS.sIJ(iC,jC,0.0)

#-- Then set Library/Other Components -------------------------------

    for iC in range(nComp) :

        cINam = clsEOS.gPP("CN",iC)

        #print("iC,cINam ",iC,cINam)

        if   cINam.upper() == "N2" :
            
            for jC in range(iC+1,nComp) :
                
                cJNam = clsEOS.gPP("CN",jC)
                typJC = typeC.get(cJNam)  #-- Other Component Library or Not?

                if typJC == "L" :
                    if clsEOS.EOS == "SRK" : KIJ = lKN2S[cJNam]
                    else                   : KIJ = lKN2P[cJNam]
                else :
                    if clsEOS.EOS == "SRK" : KIJ = 0.080
                    else                   : KIJ = 0.115
                        
                #print("N2 :jC,cJNam,typJC,KIJ ",jC,cJNam,typJC,KIJ)
                
                clsEOS.sIJ(iC,jC,KIJ)
                clsEOS.sIJ(jC,iC,KIJ)
            
        elif cINam.upper() == "CO2" :
            
            for jC in range(iC+1,nComp) :
                
                cJNam = clsEOS.gPP("CN",jC)
                typJC = typeC.get(cJNam)  #-- Other Component Library or Not?
                
                if typJC == "L" :
                    if clsEOS.EOS == "SRK" : KIJ = lKCOS[cJNam]
                    else                   : KIJ = lKCOP[cJNam]
                else :
                    if clsEOS.EOS == "SRK" : KIJ = 0.150
                    else                   : KIJ = 0.115
                        
                #print("CO2:jC,cJNam,typJC,KIJ ",jC,cJNam,typJC,KIJ)
                
                clsEOS.sIJ(iC,jC,KIJ)
                clsEOS.sIJ(jC,iC,KIJ)
            
        elif cINam.upper() == "H2S" :
            
            for jC in range(iC+1,nComp) :
                
                cJNam = clsEOS.gPP("CN",jC)
                typJC = typeC.get(cJNam)  #-- Other Component Library or Not?
                
                if typJC == "L" :
                    if clsEOS.EOS == "SRK" : KIJ = lKHSS[cJNam]
                    else                   : KIJ = lKHSP[cJNam]
                else :
                    if clsEOS.EOS == "SRK" : KIJ = 0.030
                    else                   : KIJ = 0.055
                        
                #print("H2S:jC,cJNam,typJC,KIJ ",jC,cJNam,typJC,KIJ)
                
                clsEOS.sIJ(iC,jC,KIJ)
                clsEOS.sIJ(jC,iC,KIJ)
            
        elif cINam.upper() == "C1" :
            
            for jC in range(iC+1,nComp) :
                
                cJNam = clsEOS.gPP("CN",jC)
                typJC = typeC.get(cJNam)  #-- Other Component Library or Not?

                if typJC == "L" :
                    if clsEOS.EOS == "SRK" : KIJ = lKC1S[cJNam]
                    else                   : KIJ = lKC1P[cJNam]
                else :
                    if clsEOS.EOS == "SRK" : KIJ = 0.000
                    else :
                        mWJ = clsEOS.gPP("MW",jC)
                        if mWJ > 90.0 : KIJ = 0.14*clsEOS.gPP("SG",jC) - 0.0668
                        else          : KIJ = 0.000
                            
                #print("C1 :iC,jC,cJNam,typJC,KIJ ",iC,jC,cJNam,typJC,KIJ)

                clsEOS.sIJ(iC,jC,KIJ)
                clsEOS.sIJ(jC,iC,KIJ)

        else :
            pass

#=======================================================================
#  End of Routine
#=======================================================================

    return

#=======================================================================
#  Calculates/Assigns Properties for an SCN or Plus Fraction Component
#=======================================================================

def heavyComp(molWt,SpecG,iC,clsEOS) :

    if SpecG < 0.0 :
        cF    = 0.29                    #-- 0.29 < cF < 0.31: Take Average!
        SpecG = soreideSpcG(molWt,cF)

    Tboil = soreideTboil(molWt,SpecG)

    Tcrit =    keslerLeeTcrit(SpecG,Tboil)
    Pcrit =    keslerLeePcrit(SpecG,Tboil)
    Vcrit = riaziDaubertVcrit(SpecG,Tboil)

    AcenF = LeeKeslerACF(Tboil,SpecG,Tcrit,Pcrit)
            
    Zcrit = Pcrit*Vcrit/(CO.gasCon*Tcrit)

    ParaC = Parachor(molWt)

    IdCpA = IdealGasCpA(molWt) ; IdCpB = IdealGasCpB(molWt)
    IdCpC = IdealGasCpC(molWt) ; IdCpD = IdealGasCpD(molWt)

    clsEOS.sPP("MW",iC,molWt)
    clsEOS.sPP("TC",iC,Tcrit) ; clsEOS.sPP("PC",iC,Pcrit)
    clsEOS.sPP("VC",iC,Vcrit) ; clsEOS.sPP("ZC",iC,Zcrit)
    clsEOS.sPP("AF",iC,AcenF) ; clsEOS.sPP("SG",iC,SpecG)
    clsEOS.sPP("TB",iC,Tboil) ; clsEOS.sPP("PA",iC,ParaC)
    clsEOS.sPP("MA",iC,1.0)   ; clsEOS.sPP("MB",iC,1.0)

    clsEOS.sPP("CA",iC,IdCpA) ; clsEOS.sPP("CB",iC,IdCpB)
    clsEOS.sPP("CC",iC,IdCpC) ; clsEOS.sPP("CD",iC,IdCpD)
    
    VolSh = calcVolShift(clsEOS,iC)  #-- Have to store (Pc,Tc) before this call

    clsEOS.sPP("SS",iC,VolSh)
            
#-- Adjust Acentric Factor to match Tboil ---------------------------

    AcenF = adjustAcentFac(clsEOS,iC)

    clsEOS.sPP("AF",iC,AcenF)

    return 

#=======================================================================
#  Soreide Specific Gravity Correlation
#  Whitson & Brule, Eqn.(5.44)
#  Note, split into function (this routine) and derivative (follows)
#=======================================================================

def soreideSpcG(molWt,fCons) :

    aCons = 0.28554

    specG = aCons + fCons*soreideDerv(molWt)

    return specG

def soreideDerv(molWt) :

    bCons = 65.94
    cCons =  0.129969

    sDerv = pow((molWt-bCons),cCons)

    return sDerv

#=======================================================================
#  Soreide Normal Boiling Point Temperature Correlation [in degR]
#  Whitson & Brule, Eqn.(5.45)
#=======================================================================

def soreideTboil(molWt,SpecG) :

    mWexp = -0.03522
    sGexp =  3.26600

    exArg = -4.922E-03*molWt - 4.7685*SpecG + 3.462E-03*molWt*SpecG

    Tboil = 1928.3 - 1.695E+05*pow(molWt,mWexp)*pow(SpecG,sGexp)*exp( exArg )

    return Tboil

#=======================================================================
#  Kesler-Lee Critical Temperature Correlation [in degR]
#  Whitson & Brule, Eqn.(5.45)
#=======================================================================

def keslerLeeTcrit(SpecG,Tboil) :

    third =         (0.4244 + 0.1174*SpecG)*Tboil
    fourt = 1.0E+05*(0.4669 - 3.2623*SpecG)/Tboil
    
    Tcrit = 341.5 + 811.0*SpecG + third + fourt

    return Tcrit
            
#=======================================================================
#  Kesler-Lee Critical Temperature Correlation [in psia]
#  Whitson & Brule, Eqn.(5.45)
#=======================================================================

def keslerLeePcrit(SpecG,Tboil) :

    rSG1 =  1.0/SpecG
    rSG2 = rSG1*rSG1

    third = -1.0E-03*(0.24244 + 2.2898*rSG1 + 0.11857*rSG2)
    fourt =  1.0E-07*(1.46850 + 3.6480*rSG1 + 0.47227*rSG2)
    fifth = -1.0E-10*(0.42019               + 1.69770*rSG2)

    lnPc = 8.3634 - 0.0566*rSG1 + Tboil*(third + Tboil*(fourt + Tboil*fifth))

    Pcrit = exp(lnPc)

    return Pcrit

#=======================================================================
#  Riazi-Daubert Critical Volume Correlation [in ft3/lbmol]
#  Whitson & Brule, Eqn.(5.64)
#=======================================================================

def riaziDaubertVcrit(SpecG,Tboil) :

    Vcrit = 7.0434E-07*pow(Tboil,2.3829)*pow(SpecG,-1.6830)

    return Vcrit

#=======================================================================
#  Lee-Kesler correlation for acentric factor based on normal boiling
#  point and specific gravity.
#  Lee, B.I. and Kesler, M.G.: AIChE J. 21 (1975) pg. 510.
#=======================================================================

def LeeKeslerACF(Tb,SG,Tc,Pc) :

    Third = 1.0/3.0
    Tbr   = Tb/Tc

    if Tb > 600 and Tbr >= 0.8 :
        watK = pow(Tb,Third)/SG
        AcF = -7.904 + 0.1352*watK - 0.007465*watK*watK + \
               8.359*Tbr + (1.408 - 0.01063*watK)/Tbr
    else :
        Tbr6 = pow(Tbr,6)
        lTbr = log(Tbr)
        lPbr=  log(CO.pStand/Pc)
        AcF  = (lPbr - 5.92714 + 6.09648/Tbr + 1.28862*lTbr - 0.169347*Tbr6) \
             / (15.2518 - 15.6875/Tbr - 13.4721*lTbr + 0.43577*Tbr6)

    return AcF

#========================================================================
#  Heavy Component/Plus Fraction Parachors
#========================================================================

def Parachor(molWt) :

    PC = 11.4 + molWt * (3.23 - 0.0022 * molWt)

    return PC

#========================================================================
#  SCN/Plus Fraction Ideal Gas Specific Heat Coefficients
#========================================================================

def IdealGasCpA(Mw) :

    CpA = -1.39E-02 * Mw + 1.408E-01

    return CpA
    
def IdealGasCpB(Mw) :

    CpB = 9.48E-04 * Mw - 6.450E-03

    return CpB
    
def IdealGasCpC(Mw) :

    CpC = -2.79E-07 * Mw + 8.300E-07

    return CpC
    
def IdealGasCpD(Mw) :

    CpD = ICpD =  3.46E-10 * Mw - 3.140E-10

    return CpD

#========================================================================
#  Soreide Specific Gravity Correlation
#========================================================================

def adjustFC(nP,zP,mP,sP) :

    EPS = 1.0E-08

    fC= 0.29
    dF = 1.0
    iC = 0

    while abs(dF) > EPS :

        sM = 0.0
        sG = 0.0
        sD = 0.0

        for iP in range(nP) :

            zM = zP[iP]*mP[iP]
            dG = soreideDerv(mP[iP])
            gI = soreideSpcG(mP[iP],fC)

            sM = sM + zM
            sG = sG + zM/gI
            sD = sD + zM*dG/(gI*gI)

        sC = sM/sG
        dV = sC*sD/sG

        fU = sC - sP

        dF = - fU/dV

        iC += 1

        fC = fC + dF

    return fC

#========================================================================
#  Whitson Splitting Algorithm
#  Adjust delta parameter by sample to ensure calculated plus fraction
#  mole weight matches the laboratory measured value
#========================================================================

def matchAllDelta(alfa,dOld,mObs,mSpl,wT,xZ,fCon) :

    nSam = len(dOld)
    nSpl = len(xZ)

    dNew = NP.zeros( nSam)
    zSpl = NP.zeros((nSam,nSpl))

    for iSam in range(nSam) :
        dNew[iSam] = dOld[iSam]

#== Iterate sample by sample ==========================================

    for iSam in range(nSam) :

        zSam = NP.zeros(nSpl)
        fSam = NP.zeros(nSpl)

        for iSpl in range(nSpl) :
            fSam[iSpl] = fCon[iSam][iSpl]

        aSam = alfa[iSam]
        dSam = dOld[iSam]
        mSam = mObs[iSam]

        dNew,zS = matchOneDelta(aSam,dSam,mSam,xZ,wT,fSam,mSpl,zSam)

        for iSpl in range(nSpl) :
            zSpl[iSam][iSpl] = zSam[iSpl]

    return dNew,zSpl

#========================================================================
#  One Sample Iteration of delta to match plus fraction mole weight
#========================================================================

def matchOneDelta(alfa,dOld,mObs,xZ,wT,fC,mSpl,zS) :

    nSplt = len(wT)

    EPS  = 1.0E-12
    dDel = 1.0
    dNew = dOld
    iCnt = 1

    while abs(dDel) > EPS :

        zSum = 0.0
        mSum = 0.0
        zDer = 0.0
        mDer = 0.0

        for iSp in range(nSplt) :

            funC = quadFunc(xZ[iSp],fC[iSp],alfa,dNew)

            zS[iSp] = wT[iSp]*funC

            zSum = zSum + zS[iSp]
            mSum = mSum + zS[iSp]*mSpl[iSp]

            derV = wT[iSp]*quadDerv(xZ[iSp],funC,alfa,dNew)

            zDer = zDer + derV
            mDer = mDer + derV*mSpl[iSp]

        derV = mDer/mSum - zDer/zSum
        mSum = mSum/zSum
        if derV == 0.0 : break
        derV = mSum*derV

        resI = mSum - mObs
        dDel = resI / derV

        dNew = dNew - dDel

        iCnt += 1

        if iCnt > 10 :
            print("matchOneDelta: Not Converging - Error")
            break

    return dNew,zS

#========================================================================
#  Whitson Splitting Quadrature Function
#========================================================================
    
def quadFunc(x0,fC,alfa,delt) :

    qF = fC*pow(1.0+delt,alfa)/pow(delt,x0)

    return qF

#========================================================================
#  Exact Derivative (wrt delta) of Whitson Quadrature Function
#========================================================================

def quadDerv(x0,fC,alfa,delt) :

    dF = (fC/delt)*(alfa/(1.0+log(delt)) - x0)
    
    return dF

#========================================================================
#  Calculates the Volume Shift Coefficient for an SCN or Pseudo-Component
#========================================================================

def calcVolShift(clsEOS,iC) :

#-- Get component Tc, Pc, AcF, Mw & SG --    

    Pcrit = clsEOS.gPP("PC",iC)
    Tcrit = clsEOS.gPP("TC",iC)
    AcenF = clsEOS.gPP("AF",iC)
    MolWt = clsEOS.gPP("MW",iC)
    SpecG = clsEOS.gPP("SG",iC)

#-- Reduced Pressure and Temperature (of Standard Conditions) --    

    pRed = CO.pStand/Pcrit
    tRed = CO.tStand/Tcrit

    pByT = pRed/tRed
    alsq = CE.sqrtAlpha(tRed,AcenF,clsEOS)

#-- OmegaA & OmegaB Coefficients --    

    omgA = clsEOS.OA*clsEOS.gPP("MA",iC)
    omgB = clsEOS.OB*clsEOS.gPP("MB",iC)

#-- EoS (A,B) Coefficients --    

    AEoS = omgA*(pByT/tRed)*pow(alsq,2)
    BEoS = omgB* pByT

#== Setup and Solve the Cubic EoS =====================================    

    iLiq = 1    #-- iLiq = 1 => Liquid, -1 => Vapour, 0 => GFE-Test

    eta,dumS = CE.setupSolveCubicD(iLiq,AEoS,BEoS,clsEOS)

    #print("calcVolShift: AEoS,BEoS,eta ",AEoS,BEoS,eta)

#-- 2-Parameter Z-Factor, EoS and "Measured" Molar Volumes --

    ZLiq = eta + BEoS

    VEoS = ZLiq*CO.gasCon*CO.tStand/CO.pStand
    VObs = MolWt/(CO.denWat*SpecG)     #-- Factor converts SG to lb/ft3

    vShf = VEoS - VObs

    bEoS = omgB*CO.gasCon*Tcrit/Pcrit

    SS = vShf/bEoS

    #print("test: VEoS,VObs,vShf,bEoS,SS {:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(VEoS,VObs,vShf,bEoS,SS))

    return SS

#========================================================================
#  Soreide Scheme to Adjust Acentric Factor to Match Tboil
#========================================================================

def adjustAcentFac(clsEOS,iC) :

    pRes = CO.pStand
    tRes = clsEOS.gPP("TB",iC)

    #print("adjustAF: iC,TB {:2d} {:10.3f}".format(iC,tRes))

    tCrt = clsEOS.gPP("TC",iC)
    pCrt = clsEOS.gPP("PC",iC)
    AcnF = clsEOS.gPP("AF",iC)
    OmgA = clsEOS.gPP("MA",iC)*clsEOS.OA
    OmgB = clsEOS.gPP("MB",iC)*clsEOS.OB

    tRed = tRes/tCrt
    pRed = pRes/pCrt

    BEoS = OmgB*pRed/ tRed           #-- B-coeff not a function of AF
    ACon = OmgA*pRed/(tRed*tRed)     #-- None-AF part of A-coeff

#== Iterate in AF ===================================================

    dw = 1.0
    iW = 1

    while abs(dw) > 1.0E-8 :

        E2,E1,eLiq,eVap,AEoS,dAdw = eosCoefs1C(ACon,BEoS,tRed,AcnF,clsEOS)

#-- log(Fug-Coeff) and its AF-Derivative ----------------------------    

        logFL,derFL = logPhiDerv1C(E2,E1,eLiq,AEoS,BEoS,dAdw,clsEOS)
        logFV,derFV = logPhiDerv1C(E2,E1,eVap,AEoS,BEoS,dAdw,clsEOS)

        #print("adjustAF: FL,FV,dL,dV {:10.3e} {:10.3e} {:10.3e} {:10.3e}".format(logFL,logFV,derFL,derFV))

        res = logFL - logFV
        jac = derFL - derFV

        dw = - res/jac

        #print("adjustAF: iC,iW,AF,res,jac,dw {:2d} {:2d} {:8.5f} {:10.3e} {:10.3e} {:10.3e}".format(iC+1,iW,AcnF,res,jac,dw))

        AcnF = AcnF + dw

        iW += 1

        if iW > 10 : break

    return AcnF

#========================================================================
#  EoS Coefficients for Single Component
#========================================================================

def eosCoefs1C(ACon,BEoS,tRed,AcnF,clsEOS) :

    sqAl = CE.sqrtAlpha(tRed,AcnF,clsEOS)
    dmdw = CE.calcdMdw(AcnF,clsEOS)

    AEoS = ACon*sqAl*sqAl
    dAdw = 2.0*AEoS*(1.0-sqrt(tRed))*dmdw/sqAl

#-- Cubic Coefficients of f(x) = E3*x^3 + E2*x^2 + E1*x + E0 = 0 ----   

    E0 = -clsEOS.n4*BEoS*BEoS
    E2 =  clsEOS.n3*BEoS - 1.0
    E1 = -E0 - E2 - 1.0 + AEoS

    nRoot,eLiq,eVap = CE.solveCubic(BEoS,E2,E1,E0)

    #print("adjustAF: nR,eS,eL {:1d} {:8.5f} {:8.5f}".format(nRoot,eLiq,eVap))

    return E2,E1,eLiq,eVap,AEoS,dAdw

#========================================================================
#  Fugacity Coefficient (and its AF-derivative) of a Single Component
#========================================================================

def logPhiDerv1C(E2,E1,eta,AEoS,BEoS,dAdw,clsEOS) :

    dfde = E1 + eta*(2.0*E2 + eta*3.0)
    
    dedw = -eta*dAdw/dfde

    sig1 = eta + clsEOS.n1*BEoS
    sig2 = eta + clsEOS.n2*BEoS

    C1 = - log(eta)
    C2 =       eta + BEoS - 1.0
    C3 = clsEOS.n5*AEoS/BEoS
    C4 =   log(sig2/sig1)

    logF = C1 + C2 + C3*C4

    dC1dw = - dedw/eta
    dC2dw =   dedw
    dC3dw = clsEOS.n5*dAdw/BEoS
    dC4dw = (1.0/sig2 - 1.0/sig1)*dedw

    derF = dC1dw + dC2dw + C3*dC4dw + C4*dC3dw

    return logF,derF

#========================================================================
#  Plus Fraction Name Manipulation
#  If Already split, i.e. now in regression, remove all after '+'
#========================================================================

def procPlusFracName(nOld,cOld,clsEOS) :

    cSpl = cOld.split('+')      #-- Split the Name using Plus-Char

    cNew = cSpl[0] + "+"

    #print("procPlusFracName: cOld,cNew ",cOld,cNew)

    clsEOS.sPP("CN",nOld,cNew)

    return cNew

#========================================================================
#  Sort Components from Most to Least Volatile - Requirement of RR Solver
#========================================================================

def sortComponents(clsEOS,dicSAM) :

    nCom = clsEOS.NC
    nPrp = clsEOS.nPR
    nSam = clsEOS.NS

#-- Use Wilson K-Values as Measured of Volatility -------------------    

    wilK = ST.wilsonK(CO.pStand,CO.tStand,clsEOS)

    iSor = NP.argsort(wilK)[::-1]  #-- [::-1] => High to Low Ordering

    #print("sortComponents: iSor ",iSor)

#== Take Deep Copy of the EOS & Samples Classes =======================

    clsEOSC = deepcopy(clsEOS)
    dicSAMC = deepcopy(dicSAM)

#== Update the EOS Class and Dictionary of Samples ====================

    for iPrp in range(nPrp) :

        for iC in range(nCom) :
            clsEOS.PRP[iPrp][iSor[iC]] = clsEOSC.PRP[iPrp][iC]

    for iC in range(nCom) :

        for jC in range(nCom) :
            clsEOS.KIJ[iSor[iC]][iSor[jC]] = clsEOSC.KIJ[iC][jC]

    for iSam in range(nSam) :

        for iC in range(nCom) :
            dicSAM[iSam].sComp[iSor[iC]] = dicSAMC[iSam].sComp[iC]
            
#== Return updated EOS Class and Dictionary of Samples ================

    return clsEOS,dicSAM

#========================================================================
#  End of Module
#========================================================================
