#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Units Module
# 
#========================================================================

#=======================================================================
#  If new Unit-Types added, insert array in (alphabetic) position below
#  and define the appropriate number of entries, then define a new
#  dimension variable and ensure it is added to mUnit
#  If new units added to existing arrays, append them to the end of
#  the existing list(s)
#  Make sure the unit conversion factors are updated accordingly
#=======================================================================

class classUnits :

    def __init__(self) :

        sArea = ["ft2","m2","km2","acre"]
        sCGRz = ["scf/scf","stb/scf","stb/mscf","sm3/sm3","sm3/ksm3"]
        sComp = ["1/psi","1/bar","1/kpa"]
        sDens = ["lb/ft3","kg/m3","gm/cm3","degapi"]
        sDimZ = ["dim"]
        sGfvf = ["rcf/scf","rb/scf","rb/mscf","rm3/sm3","rm3/ksm3"]
        sGORz = ["scf/scf","scf/stb","mscf/stb","sm3/sm3","ksm3/sm3"]
        sHght = ["ft","m","in","mm"]
        sKHzz = ["md.ft","md.m"]
        sMole = ["frac","perc"]
        sMolV = ["ft3/lbmol","m3/kgmol","cm3/gmol"]
        sMolW = ["lb/lbmol","kg/kgmol","gm/gmol"]
        sNull = [""]
        sOfvf = ["rcf/scf","rb/stb","rm3/sm3"]
        sRatG = ["scf/d","mscf/d","mmscf/d","sm3/d","ksm3/d","msm3/d"]
        sRatL = ["scf/d","stb/d","mstb/d","mmstb/d","sm3/d","ksm3/d","msm3/d"]
        sSalt = ["mfrac","ppm","wt%"]
        sVolG = ["scf","mscf","mmscf","bscf","sm3","ksm3","msm3"]
        sVolL = ["scf","stb","mstb","mmstb","sm3","ksm3","msm3"]
        sVisc = ["cp","mpa.s"]

        self.sPres = ["psia","psig","bara","barg","kpa","atm"]
        self.sTemp = ["degr","degf","degc","kelv"]

        mArea = len(sArea)
        mCGRz = len(sCGRz)
        mComp = len(sComp)
        mDens = len(sDens)
        mDimZ = len(sDimZ)
        mGfvf = len(sGfvf)
        mGORz = len(sGORz)
        mHght = len(sHght)
        mKHzz = len(sKHzz)
        mMole = len(sMole)
        mMolV = len(sMolV)
        mMolW = len(sMolW)
        mNull = len(sNull)
        mOfvf = len(sOfvf)
        mPres = len(self.sPres)
        mRatG = len(sRatG)
        mRatL = len(sRatL)
        mSalt = len(sSalt)
        mTemp = len(self.sTemp)
        mVolG = len(sVolG)
        mVolL = len(sVolL)
        mVisc = len(sVisc)

        mTotU = mArea + mCGRz + mComp + mDens + mDimZ + mGfvf + mGORz + \
                mHght + mKHzz + mMole + mMolV + mMolW + mNull + mOfvf + \
                mPres + mRatG + mRatL + mSalt + mTemp + mVolG + mVolL + \
                mVisc

        mR = 0          ; tArea = [mR,mR,mR,mR]
        mR = mR + mArea ; tCGRz = [mR,mR,mR,mR,mR]
        mR = mR + mCGRz ; tComp = [mR,mR,mR]
        mR = mR + mComp ; tDens = [mR,mR,mR,mR]
        mR = mR + mDens ; tDimZ = [mR]
        mR = mR + mDimZ ; tGfvf = [mR,mR,mR,mR,mR]
        mR = mR + mGfvf ; tGORz = [mR,mR,mR,mR,mR]
        mR = mR + mGORz ; tHght = [mR,mR,mR,mR]
        mR = mR + mHght ; tKHzz = [mR,mR]
        mR = mR + mKHzz ; tMole = [mR,mR,mR]
        mR = mR + mMole ; tMolV = [mR,mR,mR]
        mR = mR + mMolV ; tMolW = [mR,mR,mR]
        mR = mR + mMolW ; tNull = [mR]
        mR = mR + mNull ; tOfvf = [mR,mR,mR]
        mR = mR + mOfvf ; tPres = [mR,mR,mR,mR,mR,mR]
        mR = mR + mPres ; tRatG = [mR,mR,mR,mR,mR,mR]
        mR = mR + mRatG ; tRatL = [mR,mR,mR,mR,mR,mR,mR]
        mR = mR + mRatL ; tSalt = [mR,mR,mR]
        mR = mR + mSalt ; tTemp = [mR,mR,mR,mR]
        mR = mR + mTemp ; tVolG = [mR,mR,mR,mR,mR,mR,mR]
        mR = mR + mVolG ; tVolL = [mR,mR,mR,mR,mR,mR,mR]
        mR = mR + mVolL ; tVisc = [mR,mR]

#== Assemble all individual units into one set ==========================

        self.sUnit = []                  #-- Create null list

        self.sUnit = self.sUnit + sArea  #-- And add each list in turn ...
        self.sUnit = self.sUnit + sCGRz
        self.sUnit = self.sUnit + sComp
        self.sUnit = self.sUnit + sDens
        self.sUnit = self.sUnit + sDimZ
        self.sUnit = self.sUnit + sGfvf
        self.sUnit = self.sUnit + sGORz
        self.sUnit = self.sUnit + sHght
        self.sUnit = self.sUnit + sKHzz
        self.sUnit = self.sUnit + sMole
        self.sUnit = self.sUnit + sMolV
        self.sUnit = self.sUnit + sMolW
        self.sUnit = self.sUnit + sNull
        self.sUnit = self.sUnit + sOfvf
        self.sUnit = self.sUnit + self.sPres
        self.sUnit = self.sUnit + sRatG
        self.sUnit = self.sUnit + sRatL
        self.sUnit = self.sUnit + sSalt
        self.sUnit = self.sUnit + self.sTemp
        self.sUnit = self.sUnit + sVolG
        self.sUnit = self.sUnit + sVolL
        self.sUnit = self.sUnit + sVisc

        self.uType = []

        self.uType = self.uType + tArea
        self.uType = self.uType + tCGRz
        self.uType = self.uType + tComp
        self.uType = self.uType + tDens
        self.uType = self.uType + tDimZ
        self.uType = self.uType + tGfvf
        self.uType = self.uType + tGORz
        self.uType = self.uType + tHght
        self.uType = self.uType + tKHzz
        self.uType = self.uType + tMole
        self.uType = self.uType + tMolV
        self.uType = self.uType + tMolW
        self.uType = self.uType + tNull
        self.uType = self.uType + tOfvf
        self.uType = self.uType + tPres
        self.uType = self.uType + tRatG
        self.uType = self.uType + tRatL
        self.uType = self.uType + tSalt
        self.uType = self.uType + tTemp
        self.uType = self.uType + tVolG
        self.uType = self.uType + tVolL
        self.uType = self.uType + tVisc

        mTest = len(self.sUnit)

        #print("mTotU,mTest ",mTotU,mTest)

#=======================================================================
#  Unit Conversion Internal <-> External
#  Formula is xInt = [C(1)*xExt + C(2)]/[C(3)*xExt + C(4)]
#  Each variable has 4-rows of C(1)'s, C(2)'s, C(3)'s and C(4)'s
#=======================================================================

        dArea = [ [1.0      , 0.0, 0.0, 1.0],
                  [10.7639  , 0.0, 0.0, 1.0],
                  [10.7639E6, 0.0, 0.0, 1.0],
                  [43560.0  , 0.0, 0.0, 1.0] ]

        dCGRz = [ [ 1.0       , 0.0, 0.0, 1.0 ],
                  [ 5.61457   , 0.0, 0.0, 1.0 ],
                  [ 0.00561457, 0.0, 0.0, 1.0 ],
                  [ 1.0       , 0.0, 0.0, 1.0 ],
                  [ 0.001     , 0.0, 0.0, 1.0 ] ]

        dComp = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                  [ 0.0689474, 0.0, 0.0, 1.0 ],
                  [ 6.89474  , 0.0, 0.0, 1.0 ] ]

        dDens = [ [ 1.0     ,    0.0 , 0.0,   1.0],
                  [ 0.062428,    0.0 , 0.0,   1.0],
                  [62.428   ,    0.0 , 0.0,   1.0],
                  [ 0.0     , 8833.56, 1.0, 131.5] ]

        dDimZ = [ [1.0, 0.0, 0.0, 1.0 ] ]

        dGfvf = [ [1.0       , 0.0, 0.0, 1.0],
                  [5.61457   , 0.0, 0.0, 1.0],
                  [0.00561457, 0.0, 0.0, 1.0],
                  [1.0       , 0.0, 0.0, 1.0],
                  [0.001     , 0.0, 0.0, 1.0] ]

        dGORz = [ [   1.0     , 0.0, 0.0, 1.0],
                  [   0.178108, 0.0, 0.0, 1.0],
                  [ 178.108   , 0.0, 0.0, 1.0],
                  [   1.0     , 0.0, 0.0, 1.0],
                  [1000.0     , 0.0, 0.0, 1.0] ]

        dKHzz = [ [1.0    , 0.0, 0.0, 1.0],
                  [3.28084, 0.0, 0.0, 1.0] ]

        dHght = [ [1.0      , 0.0, 0.0, 1.0],
                  [3.28084  , 0.0, 0.0, 1.0],
                  [0.0833333, 0.0, 0.0, 1.0],
                  [3.28084E3, 0.0, 0.0, 1.0] ]

        dMole = [ [1.0 , 0.0, 0.0, 1.0],
                  [0.01, 0.0, 0.0, 1.0] ]

        dMolV = [ [ 1.0     , 0.0, 0.0, 1.0],
                  [ 0.062428, 0.0, 0.0, 1.0],
                  [62.428   , 0.0, 0.0, 1.0] ]

        dMolW = [ [1.0, 0.0, 0.0, 1.0],
                  [1.0, 0.0, 0.0, 1.0],
                  [1.0, 0.0, 0.0, 1.0] ]

        dNull = [ [1.0, 0.0, 0.0, 1.0 ] ]

        dOfvf = [ [1.0, 0.0, 0.0, 1.0 ],
                  [1.0, 0.0, 0.0, 1.0 ],
                  [1.0, 0.0, 0.0, 1.0 ] ]

        dPres = [ [ 1.0     ,  0.0   , 0.0, 1.0],
                  [ 1.0     , 14.6959, 0.0, 1.0],
                  [14.5038  ,  0.0   , 0.0, 1.0],
                  [14.5038  , 14.6959, 0.0, 1.0],
                  [ 0.145038,  0.0   , 0.0, 1.0],
                  [14.6959  ,  0.0   , 0.0, 1.0] ]

        dRatG = [ [ 1.0E0   , 0.0, 0.0, 1.0],
                  [ 1.0E3   , 0.0, 0.0, 1.0],
                  [ 1.0E6   , 0.0, 0.0, 1.0],
                  [35.3147E0, 0.0, 0.0, 1.0],
                  [35.3147E3, 0.0, 0.0, 1.0],
                  [35.3147E6, 0.0, 0.0, 1.0] ]

        dRatL = [ [ 1.0E0    , 0.0, 0.0, 1.0],
                  [ 5.61458E0, 0.0, 0.0, 1.0],
                  [ 5.61458E3, 0.0, 0.0, 1.0],
                  [ 5.61458E6, 0.0, 0.0, 1.0],
                  [35.3147E0 , 0.0, 0.0, 1.0],
                  [35.3147E3 , 0.0, 0.0, 1.0],
                  [35.3147E6 , 0.0, 0.0, 1.0] ]

        dSalt = [ [1.0    , 0.0, 0.0, 1.0],
                  [1.0E-06, 0.0, 0.0, 1.0],
                  [1.0E-04, 0.0, 0.0, 1.0]  ]

        dTemp = [ [1.0,   0.0 , 0.0, 1.0],
                  [1.0, 459.67, 0.0, 1.0],
                  [1.8, 491.67, 0.0, 1.0],
                  [1.8,   0.0 , 0.0, 1.0] ]

        dVolG = [ [ 1.0E0   , 0.0, 0.0, 1.0],
                  [ 1.0E3   , 0.0, 0.0, 1.0],
                  [ 1.0E6   , 0.0, 0.0, 1.0],
                  [ 1.0E9   , 0.0, 0.0, 1.0],
                  [35.3147E0, 0.0, 0.0, 1.0],
                  [35.3147E3, 0.0, 0.0, 1.0],
                  [35.3147E6, 0.0, 0.0, 1.0] ]

        dVolL = [ [ 1.0E0    , 0.0, 0.0, 1.0],
                  [ 5.61458E0, 0.0, 0.0, 1.0],
                  [ 5.61458E3, 0.0, 0.0, 1.0],
                  [ 5.61458E6, 0.0, 0.0, 1.0],
                  [35.3147E0 , 0.0, 0.0, 1.0],
                  [35.3147E3 , 0.0, 0.0, 1.0],
                  [35.3147E6 , 0.0, 0.0, 1.0] ]

        dVisc = [ [1.0, 0.0, 0.0, 1.0 ],
                  [1.0, 0.0, 0.0, 1.0 ] ]

#========================================================================
#  Assemble all individual units into one set
#========================================================================

        self.dUnit = []                     #-- Empty list

        self.dUnit = self.dUnit + dArea     #-- Add each array in turn ...
        self.dUnit = self.dUnit + dCGRz
        self.dUnit = self.dUnit + dComp
        self.dUnit = self.dUnit + dDens
        self.dUnit = self.dUnit + dDimZ
        self.dUnit = self.dUnit + dGfvf
        self.dUnit = self.dUnit + dGORz
        self.dUnit = self.dUnit + dHght
        self.dUnit = self.dUnit + dKHzz
        self.dUnit = self.dUnit + dMole
        self.dUnit = self.dUnit + dMolV
        self.dUnit = self.dUnit + dMolW
        self.dUnit = self.dUnit + dNull
        self.dUnit = self.dUnit + dOfvf
        self.dUnit = self.dUnit + dPres
        self.dUnit = self.dUnit + dRatG
        self.dUnit = self.dUnit + dRatL
        self.dUnit = self.dUnit + dSalt
        self.dUnit = self.dUnit + dTemp
        self.dUnit = self.dUnit + dVolG
        self.dUnit = self.dUnit + dVolL
        self.dUnit = self.dUnit + dVisc

#=======================================================================
#  Default Units for a Given Quantity Type
#=======================================================================

        self.defU = {'AREA':'ft2'    ,'BG'  :'rcf/scf','BO'  :'rcf/scf',
                     'BT'  :'rcf/scf','DENO':'lb/ft3' ,'DENS':'lb/ft3' ,
                     'GGRV':'dim'    ,'GOR' :'scf/scf',
                     'HEIG':'ft'     ,'KH'  :'md.ft'  ,
                     'MOLE':'frac'   ,'MREM':'frac'   ,
                     'PRES':'psia'   ,'PSAT':'psia'   ,
                     'RELV':'dim'    ,
                     'SKIN':'dim'    ,'SLIQ':'frac'   ,
                     'TEMP':'degr'   ,
                     'VISG':'cp'     ,'VISO':'cp'     ,'VOLG':'scf'    ,
                     'VOLL':'scf'    ,'VSWL':'dim'    ,'VFRC':'dim'    ,
                     'ZC1' :'frac'   ,'ZC7+':'frac'   ,'ZFAC':'dim'    }

#========================================================================
#  Convert eXternal (User) Units to Internal Units
#========================================================================

    def X2I(self,extVal,usrUni) :

        iA = None
        try :
            iA = self.sUnit.index( usrUni.lower() )
        except ValueError :
            iA = -1
            print("X2I: User specified unit ",usrUni," not recognised - Error")

        if iA > -1 :
            intVal = (self.dUnit[iA][0]*extVal + self.dUnit[iA][1])/ \
                     (self.dUnit[iA][2]*extVal + self.dUnit[iA][3])
        else :
            intVal = -999.9

        return intVal

#========================================================================
#  Convert Internal Units to eXternal (User) Units
#========================================================================

    def I2X(self,intVal,usrUni) :

        iA = None
        try :
            iA = self.sUnit.index( usrUni.lower() )
        except ValueError :
            iA = -1
            print("I2X: User specified unit ",usrUni," not recognised - Error")

        extVal = -(self.dUnit[iA][3]*intVal - self.dUnit[iA][1])/ \
                  (self.dUnit[iA][2]*intVal - self.dUnit[iA][0])

        return extVal

#========================================================================
#  For a given User Unit, Find the Type and Hence Default (Internal) Unit
#========================================================================

    def IUN(self,usrUni) :

        iA = None
        try :
            iA = self.sUnit.index( usrUni.lower() )
        except ValueError :
            iA = -1
            print("IUN: User specified unit ",usrUni," not recognised - Error")

        iT = self.uType[iA]

        intUni = self.sUnit[iT]

        return intUni

#========================================================================
#  For a given User Type, Find the Default (Internal) Unit
#========================================================================

    def DUN(self,usrUni) :

        return self.defU.get(usrUni)

#========================================================================
#  End of Module
#========================================================================


