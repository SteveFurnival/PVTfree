
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Calculates an Internal Fluid Definition from a User Definition
#  
#========================================================================

from   math          import exp
from   scipy.special import gamma

import allData   as AD
import allProp   as AP
import calcProps as CP
import utilities as UT

#========================================================================
#  Calculates Internal Sample Description from User Supplied Data
#========================================================================

def calcIntSamp(clsEOS,dicSAM,clsUNI,clsIO) :

   nSamp = len(dicSAM)
   nComp = clsEOS.nComp

#== Process Existing Samples ==========================================

   iPlus = 0
   for iC in range(nComp) :
      if clsEOS.gPP("MW",iC) > 90.0 :
         iPlus = iC
         break

   for iS in range(nSamp-1) :
      clsSAM = dicSAM[iS]
      clsSAM.sCom[iPlus] = UT.moleFracC7P(clsSAM.sCom,clsEOS)
      for iC in range(iPlus+1,nComp) : clsSAM.sCom[iC] = 0.0

#== New Sample ========================================================

   clsSAM = dicSAM[nSamp-1]   #-- Get Working Sample-Class

#-- Extract eXact (full supplied) string or Part of string ----------

   X = 'X'
   P = 'P'

#-- Total Number of Components and Plus Fraction Properties ---------   

   nCom = clsSAM.nCtot

   nPlus = clsSAM.uPlsCN
   mPlus = clsSAM.uPlsMW
   sPlus = clsSAM.uPlsSG

#-- Inorganics and Neo-Pentane --------------------------------------   

   iH2  = findComp("H2" ,X,clsSAM)
   iHE  = findComp("HE" ,X,clsSAM)
   iN2  = findComp("N2" ,X,clsSAM)
   iCO2 = findComp("CO2",X,clsSAM)
   iH2S = findComp("H2S",X,clsSAM)
   iNeo = findComp("Neo",P,clsSAM)

#-- First Few Hydrocarbons ------------------------------------------   

   iC1  = findComp("C1" ,X,clsSAM)
   iC2  = findComp("C2" ,X,clsSAM)
   iC3  = findComp("C3" ,X,clsSAM)

   iIC4 = findComp("iC4",X,clsSAM)
   iNC4 = findComp("nC4",X,clsSAM)
   iIC5 = findComp("iC5",X,clsSAM)
   iNC5 = findComp("nC5",X,clsSAM)

#-- SCN's C6 to C11 -------------------------------------------------   

   iC6  = findComp("C6" ,X,clsSAM)
   iC7  = findComp("C7" ,X,clsSAM)
   iC8  = findComp("C8" ,X,clsSAM)
   iC9  = findComp("C9" ,X,clsSAM)
   iC10 = findComp("C10",X,clsSAM)
   iC11 = findComp("C11",X,clsSAM)

#== Process ===========================================================

#-- N2 [Nitrogen]: Also H2 & HE (will be zero or trace) -------------

   zN2 = 0.0
   
   if iN2  >= 0 : zN2  = zN2  + clsSAM.uCom[iN2]
   if iH2  >= 0 : zN2  = zN2  + clsSAM.uCom[iH2]
   if iHE  >= 0 : zN2  = zN2  + clsSAM.uCom[iHE]

   if zN2  > 0.0 : clsSAM.setInternal("N2",zN2)

#-- CO2 and H2S -----------------------------------------------------

   zCO2 = 0.0 ; zH2S = 0.0

   if iCO2 >= 0 : zCO2 = zCO2 + clsSAM.uCom[iCO2]
   if iH2S >= 0 : zH2S = zH2S + clsSAM.uCom[iH2S]

   if zCO2 > 0.0 : clsSAM.setInternal("CO2",zCO2)
   if zH2S > 0.0 : clsSAM.setInternal("H2S",zH2S)

#-- C1 -> nC5 -------------------------------------------------------

   zC1  = 0.0 ; zC2  = 0.0 ; zC3  = 0.0
   zIC4 = 0.0 ; zNC4 = 0.0 ; zIC5 = 0.0 ; zNC5 = 0.0

   if iC1  >= 0 : zC1  = zC1  + clsSAM.uCom[iC1]
   if iC2  >= 0 : zC2  = zC2  + clsSAM.uCom[iC2]
   if iC3  >= 0 : zC3  = zC3  + clsSAM.uCom[iC3]
   if iIC4 >= 0 : zIC4 = zIC4 + clsSAM.uCom[iIC4]
   if iNC4 >= 0 : zNC4 = zNC4 + clsSAM.uCom[iNC4]
   if iIC5 >= 0 : zIC5 = zIC5 + clsSAM.uCom[iIC5]
   if iNC5 >= 0 : zNC5 = zNC5 + clsSAM.uCom[iNC5]

   if iNeo >= 0 : zIC5 = zIC5 + clsSAM.uCom[iNeo]

   if zC1  > 0.0 : clsSAM.setInternal("C1" ,zC1 )
   if zC2  > 0.0 : clsSAM.setInternal("C2" ,zC2 )
   if zC3  > 0.0 : clsSAM.setInternal("C3" ,zC3 )
   if zIC4 > 0.0 : clsSAM.setInternal("IC4",zIC4)
   if zNC4 > 0.0 : clsSAM.setInternal("NC4",zNC4)
   if zIC5 > 0.0 : clsSAM.setInternal("IC5",zIC5)
   if zNC5 > 0.0 : clsSAM.setInternal("NC5",zNC5)

#-- SCN's C6 -> C10 -------------------------------------------------

   uC7P = []

   if iC6  >= 0   : zC6  = clsSAM.uCom[iC6]
   if zC6  >  0.0 : clsSAM.setInternal("C6",zC6)

   if nPlus > 7  : uC7P = sumSCN(iC6,iC7,uC7P,clsSAM)
   if nPlus > 8  : uC7P = sumSCN(iC7,iC8,uC7P,clsSAM)
   if nPlus > 9  : uC7P = sumSCN(iC8,iC9,uC7P,clsSAM)
   
   if nPlus > 10 : #-- C10
      
      uC7P = sumSCN(iC9,iC10,uC7P,clsSAM)

      iCN = 11
      for iC in range(iC10+1,nCom-1) :
         zNam = 'C' + str(iCN)
         uC7P.append(clsSAM.uCom[iC])
         iCN = iCN + 1

   nC7P = len(uC7P)

   #print("nC7P ",nC7P)

#-- User Plus Fraction ----------------------------------------------

   zPlus = "C" + str(nPlus) + "+"

#== Internal [C7+] Plus Fraction Properties ===========================   
   
   clsLIB = AP.classLIB()

   zC7P = 0.0 ; mC7P = 0.0 ; sC7P = 0.0 ; jC7P = 7

   for iC7P in range(nC7P) :
      zNam = 'C'  + str(jC7P)
      zC7P = zC7P + uC7P[iC7P]
      mC7P = mC7P + uC7P[iC7P]*clsLIB.scnMW[zNam]
      sC7P = sC7P + uC7P[iC7P]*clsLIB.scnMW[zNam]/clsLIB.scnSG[zNam]
      jC7P = jC7P + 1
      #print("zNam,z,Mw,SG ",zNam,uC7P[iC7P],clsLIB.scnMW[zNam],clsLIB.scnSG[zNam])

#-- Add User Plus Fraction ------------------------------------------

   zPls = clsSAM.uCom[nCom-1]

   #print("zPls,mPls,sPls ",zPls,clsSAM.uPlsMW,clsSAM.uPlsSG)
      
   zC7P = zC7P + zPls
   mC7P = mC7P + zPls*clsSAM.uPlsMW
   sC7P = sC7P + zPls*clsSAM.uPlsMW/clsSAM.uPlsSG

#-- Renormalise the MW & SG -----------------------------------------      

   zInt = 'C7+'
   sC7P = mC7P/sC7P
   mC7P = mC7P/zC7P
      
   #print("zNam,zC7P,mC7P,sC7P {:s} {:7.5f} {:7.3f} {:7.5f}".format(zInt,zC7P,mC7P,sC7P))

#== Calculate Whitson-Alpha Parameter, if Appropriate =================
   
   if nC7P > 5 :

#-- Minimum & Maximum Values ----------------------------------------
      
      aMin = 0.5 ; fMin = fitAlpha(aMin,zC7P,mC7P,uC7P)
      aMax = 2.5 ; fMax = fitAlpha(aMax,zC7P,mC7P,uC7P)

#== Brent's Method to minimise fSSQ and hence optimal-Alpha ===========      

      fMid = 1.0E+6

#-- Mid-Value less than Min or Max? ---------------------------------      

      while fMid > fMin or fMid > fMax :

         aMid = 0.5*(aMin+aMax)
         fMid = fitAlpha(aMid,zC7P,mC7P,uC7P)

         if fMid > fMax : aMin = aMid ; fMin = fMid
         if fMid > fMin : aMax = aMid ; fMax = fMid

#== Call Brent ========================================================

      alfa = brentAlpha(aMin,aMid,aMax,zC7P,mC7P,uC7P)

   else :

      alfa = 1.0  #-- Default Value if Fit Not Possible

#== Add Values to the Sample Class ====================================

   clsSAM.setInternal(zInt,zC7P)
   
   clsSAM.setIntPlusMW(mC7P)
   clsSAM.setIntPlusSG(sC7P)
   clsSAM.setIntPlusAL(alfa)

   nUser = len(clsSAM.iCom)
   nPseu = clsEOS.NumPsu      #-- #Pseudo's to split from C7+

   nCint = nUser + nPseu - 1

   clsSAM.setIntComp(nCint,nUser)

   clsEOS.setNComp(nCint)
   clsEOS.setNUser(nUser)
   clsEOS.setNPseu(nPseu)
   clsEOS.setNSamp(nSamp)

#-- And define component names in EOS-class -------------------------

   for iC in range(nUser) : clsEOS.sNM(iC,clsSAM.iNam[iC])

#-- And copy 'internal' composition to 'working' --------------------

   for iC in range(nUser) : clsSAM.sZI(iC,clsSAM.iCom[iC])

#== Calculate Properties ==============================================   

   iERR = CP.allProps(clsEOS,dicSAM,clsUNI,clsIO)
    
#== Return ============================================================      

   return

#========================================================================
#  Assemble C7-C10 SCN's from PNA splits
#========================================================================

def sumSCN(IL,IU,uC7P,clsSAM) :

   zSum = 0.0

   for iC in range(IL+1,IU+1) :
      zSum = zSum + clsSAM.uCom[iC]

#== Return appended uC7P list =========================================

   uC7P.append(zSum)

   return uC7P

#========================================================================
#  Find a Special Component in User Set
#========================================================================

def findComp(zCom,XorP,clsSAM) :

   jC = -1
 
   if XorP == 'X' :                    #-- Exact Match Required
      for iC in range(clsSAM.nCtot) :
         if zCom.upper() == clsSAM.uNam[iC].upper() :
            jC = iC
            break
   else :                              #-- Partial Match Required
      for iC in range(clsSAM.nCtot) :
         if zCom.upper() in clsSAM.uNam[iC].upper() :
            jC = iC
            break

#== Return position in User Array =====================================   

   return jC

#========================================================================
#  Brent's Method (BM) to Minimise Whitson PDF for Alpha
#  BM - see, "Numerical Recipes", Chapter 9
#========================================================================

def brentAlpha(ax,bx,cx,z7,M7,zO) :

   itMax = 100
   cGold = 0.381966
   zEps  = 1.0E-10

   #print("brentMin: ax,bx,cx {:10.3e} {:10.3e} {:10.3e}".format(ax,bx,cx))

#-- a & b must be in ascending order --------------------------------

   if ax < cx : a = ax
   else       : a = cx

   if ax > cx : b = ax
   else       : b = cx

#-- Initialisation --------------------------------------------------

   tol = 1.0E-03     #-- Don't need the value too accurately

   v = bx
   w = v
   x = v
   e = 0.0           #-- This will be distance moved on step before last

   fx = fitAlpha(x,z7,M7,zO)

   fv = fx
   fw = fx

#----------------------------------------------------------------------
#  Main Loop
#----------------------------------------------------------------------

   for iT in range(itMax) :

      xm = 0.5*(a+b)
      Q1 = True

      tol1 = tol*abs(x) + zEps
      tol2 = 2.0*tol1

#== Convergence Test ==================================================

      if abs(x-xm) <= (tol2-0.5*(b-a)) : break

#== Trial parabolic fit ===============================================

      if abs(e) > tol1 :

         r = (x-w)*(fx-fv)
         q = (x-v)*(fx-fw)
         p = (x-v)*q - (x-w)*r
         q = 2.0*(q-r)

         if q > 0.0 : p = -p

         q     = abs(q)
         eTemp = e
         e     = d

         if abs(p) >= abs(0.5*q*eTemp) or p <= q*(a-x) or p>= q*(b-x) :
            Q1 = True
         else :
            Q1 = False
            d  = p/q
            u  = x+d
            if u-a < tol2 or b-u < tol2 :
               if xm-x > 0.0 : d =  tol1
               else          : d = -tol1

#== Golden Section Step ===============================================

      if Q1 :
         if x >= xm : e = a-x
         else       : e = b-x
         d = cGold*e

#== Arrive with d from Parabolic Fit or Golden Search =================

      if abs(d) >= tol1 :
         u = x+d
      else :
         if d > 0.0 : u = x+tol1
         else       : u = x-tol1

#== New Function Call =================================================

      fu = fitAlpha(u,z7,M7,zO)

      #print("iT,u,fu {:2d} {:10.3e} {:10.3e}".format(iT,u,fu))

#== Housekeeping ======================================================

      if fu <= fx :

         if u >= x : a = x
         else      : b = x

         v = w ; fv = fw
         w = x ; fw = fx
         x = u ; fx = fu

      else :

         if u < x : a = u
         else     : b = u

         if fu <= fw or w == x :
            v = w ; fv = fx
            w = u ; fw = fu
         else :
            v = u ; fv = fu

#== Return value ======================================================

   return x

#========================================================================
#  Generate the Whitson PDF for a Given Alpha-Parameter
#========================================================================

def fitAlpha(alfa,zC7P,mC7P,uC7P) :

#-- Initialisation --------------------------------------------------

   eta = 90.0

   beta = (mC7P - eta)/alfa
   gama = gamma(alfa)
   MwBU = eta

   #print("alfa,zC7P,mC7P,beta,gama {:3.1f} {:6.4f} {:6.2f} {:5.3f}".format(alfa,zC7P,mC7P,gama))

   sum1 = 0.0
   sum2 = 0.0
   fSSQ = 0.0

#== Loop over supplied mole fractions =================================

   for iC7P in range(len(uC7P)) :

      MwBL = MwBU
      MwBU = MwBL + 14.0

      if iC7P == len(uC7P)-1 : MwBU = 10000.0

      P0L,P1L = P0P1(alfa,beta,eta,gama,MwBL)
      P0U,P1U = P0P1(alfa,beta,eta,gama,MwBU)

      zPls = zC7P*(P0U-P0L)

      sum1 = sum1 + zPls

      MwAv = eta + beta*(P1U-P1L)/(P0U-P0L)

      sum2 = sum2 + zPls*MwAv

      delZ = (uC7P[iC7P]-zPls)/zPls

      if iC7P == len(uC7P)-1 : delZ = 0.0

      fSSQ = fSSQ + delZ*delZ

      #print("iC7P,uC7P,zPls,MwAv,delZ,fSSQ {:2d} {:6.4f} {:6.4f} {:6.2f} {:7.4f} {:7.4f}".format(iC7P,uC7P[iC7P],zPls,MwAv,delZ,fSSQ))

   sum2 = sum2/sum1

#== Return SSQ ========================================================

   return fSSQ

#========================================================================
#  Whitson & Brule Table A.4 & Section 5.3.2
#========================================================================

def P0P1(alfa,beta,eta,gama,Mw) :

#-- Initialisation --------------------------------------------------

   P0 = 0.0
   P1 = 0.0

   if Mw <= eta : return P0,P1

   y = (Mw - eta)/beta
   q = exp(-y)*y**alfa/gama

   T = 1.0/alfa
   C = 1.0
   S = T

#== Build Eqn.(5.21) ==================================================

   while abs(T) > 1.0E-8 :

      T = T*y/(alfa+C)

      S = S + T
      C = C + 1.0

#== Return P0 & P1 ====================================================

   P0 = q*S
   P1 = P0 - q/alfa

   return P0,P1

#========================================================================
#  End of Module
#========================================================================
