#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Units Module
# 
#========================================================================

#======================================================================
#  If new units are required, add the formatted (fTy) and non-formatted
#  or lowercase (sTy) unit strings, with suitable I/O formats (i.e.
#  float or exponential) and the sets of conversions factors - see below.
#======================================================================

class classUnits :

    def __init__(self) :

        self.fTy = {}     #-- Formatted (not for comparison purposes)
        self.sTy = {}     #-- Lowercase (    for comparison purposes)
        self.Con = {}
        self.For = {}     #-- Formats for I/O

#=======================================================================
#  Unit Conversion Internal <-> External
#  Formula is xInt = [C(1)*xExt + C(2)]/[C(3)*xExt + C(4)]
#  Each variable has 4-rows of C(1)'s, C(2)'s, C(3)'s and C(4)'s
#=======================================================================

        self.fTy["area"] = ["ft2","m2","km2","acre"]
        self.sTy["area"] = ["ft2","m2","km2","acre"]
        self.For["area"] = ["{:10.3e}","{:10.3e}","{:8.3f}","{:8.3f}"]
        self.Con["area"] = [ [1.0      , 0.0, 0.0, 1.0],
                             [10.7639  , 0.0, 0.0, 1.0],
                             [10.7639E6, 0.0, 0.0, 1.0],
                             [43560.0  , 0.0, 0.0, 1.0] ]

        self.fTy["cpa" ] = ["Btu/lbmol.R" ,"J/gmol.K" ,"kJ/kgmol.K" ,"cal/gmol.K", "kcal/kgmol.K"]
        self.sTy["cpa" ] = ["btu/lbmol.r" ,"j/gmol.k" ,"kj/kgmol.k" ,"cal/gmol.k", "kcal/kgmol.k"]
        self.For["cpa" ] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["cpa" ] = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                             [ 0.23885  , 0.0, 0.0, 1.0 ],
                             [ 0.23885  , 0.0, 0.0, 1.0 ],
                             [ 3.23784  , 0.0, 0.0, 1.0 ],
                             [ 3.23784  , 0.0, 0.0, 1.0 ] ]

        self.fTy["cpb" ] = ["Btu/lbmol.R2","J/gmol.K2","kJ/kgmol.K2","cal/gmol.K2","kcal/kgmol.K2"]
        self.sTy["cpb" ] = ["btu/lbmol.r2","j/gmol.k2","kj/kgmol.k2","cal/gmol.k2","kcal/kgmol.k2"]
        self.For["cpb" ] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["cpb" ] = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                             [ 0.42993  , 0.0, 0.0, 1.0 ],
                             [ 0.42993  , 0.0, 0.0, 1.0 ],
                             [ 5.82811  , 0.0, 0.0, 1.0 ],
                             [ 5.82811  , 0.0, 0.0, 1.0 ] ]

        self.fTy["cpc" ] = ["Btu/lbmol.R3","J/gmol.K3","kJ/kgmol.K3","cal/gmol.K3","kcal/kgmol.K3"]
        self.sTy["cpc" ] = ["btu/lbmol.r3","j/gmol.k3","kj/kgmol.k3","cal/gmol.k3","kcal/kgmol.k3"]
        self.For["cpc" ] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["cpc" ] = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                             [ 0.77387  , 0.0, 0.0, 1.0 ],
                             [ 0.77387  , 0.0, 0.0, 1.0 ],
                             [10.4906   , 0.0, 0.0, 1.0 ],
                             [10.4906   , 0.0, 0.0, 1.0 ] ]

        self.fTy["cpd" ] = ["Btu/lbmol.R4","J/gmol.K4","kJ/kgmol.K4","cal/gmol.K4","kcal/kgmol.K4"]
        self.sTy["cpd" ] = ["btu/lbmol.r4","j/gmol.k4","kj/kgmol.k4","cal/gmol.k4","kcal/kgmol.k4"]
        self.For["cpd" ] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["cpd" ] = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                             [ 1.39297  , 0.0, 0.0, 1.0 ],
                             [ 1.39297  , 0.0, 0.0, 1.0 ],
                             [18.8831   , 0.0, 0.0, 1.0 ],
                             [18.8831   , 0.0, 0.0, 1.0 ] ]

        self.fTy["cgr" ] = ["scf/scf","stb/scf","stb/Mscf","sm3/sm3","sm3/ksm3"]
        self.sTy["cgr" ] = ["scf/scf","stb/scf","stb/mscf","sm3/sm3","sm3/ksm3"]
        self.For["cgr" ] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["cgr" ] = [ [ 1.0       , 0.0, 0.0, 1.0 ],
                             [ 5.61457   , 0.0, 0.0, 1.0 ],
                             [ 0.00561457, 0.0, 0.0, 1.0 ],
                             [ 1.0       , 0.0, 0.0, 1.0 ],
                             [ 0.001     , 0.0, 0.0, 1.0 ] ]

        self.fTy["comp"] = ["1/psi","1/bar","1/kPa"]
        self.sTy["comp"] = ["1/psi","1/bar","1/kpa"]
        self.For["comp"] = ["{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["comp"] = [ [ 1.0      , 0.0, 0.0, 1.0 ],
                             [ 0.0689474, 0.0, 0.0, 1.0 ],
                             [ 6.89474  , 0.0, 0.0, 1.0 ] ]

        self.fTy["dens"] = ["lb/ft3","kg/m3","gm/cm3","degAPI"]
        self.sTy["dens"] = ["lb/ft3","kg/m3","gm/cm3","degapi"]
        self.For["dens"] = ["{:8.3f}","{:8.3f}","{:8.5f}","{:8.3f}"]
        self.Con["dens"] = [ [ 1.0     ,    0.0 , 0.0,   1.0],
                             [ 0.062428,    0.0 , 0.0,   1.0],
                             [62.428   ,    0.0 , 0.0,   1.0],
                             [ 0.0     , 8833.56, 1.0, 131.5] ]

        self.fTy["dim" ] = ["dim"]
        self.sTy["dim" ] = ["dim"]
        self.For["dim" ] = ["{:8.5f}"]
        self.Con["dim" ] = [ [1.0, 0.0, 0.0, 1.0 ] ]

        self.fTy["enth"] = ["Btu/lbmol" ,"J/gmol" ,"kJ/kgmol" ,"cal/gmol", "kcal/kgmol"]
        self.sTy["enth"] = ["btu/lbmol" ,"j/gmol" ,"kj/kgmol" ,"cal/gmol", "kcal/kgmol"]
        self.For["enth"] = ["{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}","{:10.3e}"]
        self.Con["enth"] = [ [ 1.0     , 0.0, 0.0, 1.0 ],
                             [ 2.32601 , 0.0, 0.0, 1.0 ],
                             [ 2.32601 , 0.0, 0.0, 1.0 ],
                             [ 0.55556 , 0.0, 0.0, 1.0 ],
                             [ 0.55556 , 0.0, 0.0, 1.0 ] ]

        self.fTy["gfvf"] = ["rcf/scf","rb/scf","rb/Mscf","rm3/sm3","rm3/ksm3"]
        self.sTy["gfvf"] = ["rcf/scf","rb/scf","rb/mscf","rm3/sm3","rm3/ksm3"]
        self.For["gfvf"] = ["{:10.3e}","{:10.3e}","{:8.3f}","{:10.3e}","{:8.3f}"]
        self.Con["gfvf"] = [ [1.0       , 0.0, 0.0, 1.0],
                             [5.61457   , 0.0, 0.0, 1.0],
                             [0.00561457, 0.0, 0.0, 1.0],
                             [1.0       , 0.0, 0.0, 1.0],
                             [0.001     , 0.0, 0.0, 1.0] ]

        self.fTy["gor" ] = ["scf/scf","scf/stb","Mscf/stb","sm3/sm3","ksm3/sm3"]
        self.sTy["gor" ] = ["scf/scf","scf/stb","mscf/stb","sm3/sm3","ksm3/sm3"]
        self.For["gor" ] = ["{:8.3f}","{:8.3f}","{:8.5f}","{:8.3f}","{:8.5f}"]
        self.Con["gor" ] = [ [   1.0     , 0.0, 0.0, 1.0],
                             [   0.178108, 0.0, 0.0, 1.0],
                             [ 178.108   , 0.0, 0.0, 1.0],
                             [   1.0     , 0.0, 0.0, 1.0],
                             [1000.0     , 0.0, 0.0, 1.0] ]

        self.fTy["leng"] = ["ft","m","in","mm"]
        self.sTy["leng"] = ["ft","m","in","mm"]
        self.For["leng"] = ["{:8.3f}","{:8.3f}","{:10.3e}","{:10.3e}"]
        self.Con["leng"] = [ [1.0      , 0.0, 0.0, 1.0],
                             [3.28084  , 0.0, 0.0, 1.0],
                             [0.0833333, 0.0, 0.0, 1.0],
                             [3.28084E3, 0.0, 0.0, 1.0] ]

        self.fTy["ift" ] = ["dyne/cm","mN/m"]
        self.sTy["ift" ] = ["dyne/cm","mn/m"]
        self.For["ift" ] = ["{:8.5f}","{:8.5f}"]
        self.Con["ift" ] = [ [ 1.0, 0.0, 0.0, 1.0],
                             [ 1.0, 0.0, 0.0, 1.0] ]

        self.fTy["inte"] = ["None"]
        self.sTy["inte"] = ["none"]
        self.For["inte"] = ["{:1.0f}"]
        self.Con["inte"] = [ [ 1, 0, 0, 1] ]

        self.fTy["kh"  ] = ["mD.ft","mD.m"]
        self.sTy["kh"  ] = ["md.ft","md.m"]
        self.For["kh"  ] = ["{:10.3e}","{:10.3e}"]
        self.Con["kh"  ] = [ [1.0    , 0.0, 0.0, 1.0],
                             [3.28084, 0.0, 0.0, 1.0] ]

        self.fTy["mole"] = ["Frac","Perc"]
        self.sTy["mole"] = ["frac","perc"]
        self.For["mole"] = ["{:8.5f}","{:8.3f}"]
        self.Con["mole"] = [ [1.0 , 0.0, 0.0, 1.0],
                             [0.01, 0.0, 0.0, 1.0] ]

        self.fTy["molv"] = ["ft3/lbmol","m3/kgmol","cm3/gmol"]
        self.sTy["molv"] = ["ft3/lbmol","m3/kgmol","cm3/gmol"]
        self.For["molv"] = ["{:8.3f}","{:8.3f}","{:8.5f}"]
        self.Con["molv"] = [ [ 1.0     , 0.0, 0.0, 1.0],
                             [ 0.062428, 0.0, 0.0, 1.0],
                             [62.428   , 0.0, 0.0, 1.0] ]

        self.fTy["molw"] = ["lb/lbmol","kg/kgmol","gm/gmol"]
        self.sTy["molw"] = ["lb/lbmol","kg/kgmol","gm/gmol"]
        self.For["molw"] = ["{:8.3f}","{:8.3f}","{:8.3f}"]
        self.Con["molw"] = [ [1.0, 0.0, 0.0, 1.0],
                             [1.0, 0.0, 0.0, 1.0],
                             [1.0, 0.0, 0.0, 1.0] ]

        self.fTy["null"] = [""]
        self.sTy["null"] = [""]
        self.For["null"] = ["{:8.5f}"]
        self.Con["null"] = [ [1.0, 0.0, 0.0, 1.0 ] ]

        self.fTy["ofvf"] = ["rcf/scf","rb/stb","rm3/sm3"]
        self.sTy["ofvf"] = ["rcf/scf","rb/stb","rm3/sm3"]
        self.For["ofvf"] = ["{:8.5f}","{:8.5f}","{:8.5f}"]
        self.Con["ofvf"] = [ [1.0, 0.0, 0.0, 1.0 ],
                             [1.0, 0.0, 0.0, 1.0 ],
                             [1.0, 0.0, 0.0, 1.0 ] ]

        self.fTy["pres"] = ["psia","psig","bara","barg","kPa","atm"]
        self.sTy["pres"] = ["psia","psig","bara","barg","kpa","atm"]
        self.For["pres"] = ["{:8.3f}","{:8.3f}","{:8.4f}","{:8.4f}","{:8.2f}","{:8.34}"]
        self.Con["pres"] = [ [ 1.0     ,  0.0   , 0.0, 1.0],
                             [ 1.0     , 14.6959, 0.0, 1.0],
                             [14.5038  ,  0.0   , 0.0, 1.0],
                             [14.5038  , 14.6959, 0.0, 1.0],
                             [ 0.145038,  0.0   , 0.0, 1.0],
                             [14.6959  ,  0.0   , 0.0, 1.0] ]

        self.fTy["ratg"] = ["scf/d","Mscf/d","MMscf/d","sm3/d","ksm3/d","Msm3/d"]
        self.sTy["ratg"] = ["scf/d","mscf/d","mmscf/d","sm3/d","ksm3/d","msm3/d"]
        self.For["ratg"] = ["{:10.3e}","{:8.3f}","{:8.5f}","{:10.3e}","{:8.3f}","{:8.5f}"]
        self.Con["ratg"] = [ [ 1.0E0   , 0.0, 0.0, 1.0],
                             [ 1.0E3   , 0.0, 0.0, 1.0],
                             [ 1.0E6   , 0.0, 0.0, 1.0],
                             [35.3147E0, 0.0, 0.0, 1.0],
                             [35.3147E3, 0.0, 0.0, 1.0],
                             [35.3147E6, 0.0, 0.0, 1.0] ]

        self.fTy["ratl"] = ["scf/d","stb/d" ,"Mstb/d" ,"MMstb/d","sm3/d" ,"ksm3/d","Msm3/d"]
        self.sTy["ratl"] = ["scf/d","stb/d" ,"mstb/d" ,"mmstb/d","sm3/d" ,"ksm3/d","msm3/d"]
        self.For["ratl"] = ["{:10.3e}","{:8.3f}","{:8.5f}","{:10.3e}","{:8.3f}","{:8.5f}"]
        self.Con["ratl"] = [ [ 1.0E0    , 0.0, 0.0, 1.0],
                             [ 5.61458E0, 0.0, 0.0, 1.0],
                             [ 5.61458E3, 0.0, 0.0, 1.0],
                             [ 5.61458E6, 0.0, 0.0, 1.0],
                             [35.3147E0 , 0.0, 0.0, 1.0],
                             [35.3147E3 , 0.0, 0.0, 1.0],
                             [35.3147E6 , 0.0, 0.0, 1.0] ]

        self.fTy["salt"] = ["mFrac","PPM","WT%"]
        self.sTy["salt"] = ["mfrac","ppm","wt%"]
        self.For["salt"] = ["{:8.5f}","{:8.1f}","{:8.3f}"]
        self.Con["salt"] = [ [1.0    , 0.0, 0.0, 1.0],
                             [1.0E-06, 0.0, 0.0, 1.0],
                             [1.0E-04, 0.0, 0.0, 1.0]  ]

        self.fTy["temp"] = ["degR","degF","degC","Kelv"]
        self.sTy["temp"] = ["degr","degf","degc","kelv"]
        self.For["temp"] = ["{:8.3f}","{:8.3f}","{:8.3f}","{:8.3f}"]
        self.Con["temp"] = [ [1.0,   0.0 , 0.0, 1.0],
                             [1.0, 459.67, 0.0, 1.0],
                             [1.8, 491.67, 0.0, 1.0],
                             [1.8,   0.0 , 0.0, 1.0] ]

        self.fTy["volg"] = ["scf","Mscf","MMscf","Bscf","sm3","ksm3","Msm3"]
        self.sTy["volg"] = ["scf","mscf","mmscf","bscf","sm3","ksm3","msm3"]
        self.For["volg"] = ["{:10.3e}","{:8.1f}","{:8.3f}","{:8.5f}","{:10.3e}","{:8.1f}","{:8.3f}"]
        self.Con["volg"] = [ [ 1.0E0   , 0.0, 0.0, 1.0],
                             [ 1.0E3   , 0.0, 0.0, 1.0],
                             [ 1.0E6   , 0.0, 0.0, 1.0],
                             [ 1.0E9   , 0.0, 0.0, 1.0],
                             [35.3147E0, 0.0, 0.0, 1.0],
                             [35.3147E3, 0.0, 0.0, 1.0],
                             [35.3147E6, 0.0, 0.0, 1.0] ]
    
        self.fTy["voll"] = ["scf","stb","Mstb","MMstb","sm3","ksm3","Msm3"]
        self.sTy["voll"] = ["scf","stb","mstb","mmstb","sm3","ksm3","msm3"]
        self.For["voll"] = ["{:10.3e}","{:8.1f}","{:8.3f}","{:8.5f}","{:10.3e}","{:8.1f}","{:8.3f}"]
        self.Con["voll"] = [ [ 1.0E0    , 0.0, 0.0, 1.0],
                             [ 5.61458E0, 0.0, 0.0, 1.0],
                             [ 5.61458E3, 0.0, 0.0, 1.0],
                             [ 5.61458E6, 0.0, 0.0, 1.0],
                             [35.3147E0 , 0.0, 0.0, 1.0],
                             [35.3147E3 , 0.0, 0.0, 1.0],
                             [35.3147E6 , 0.0, 0.0, 1.0] ]

        self.fTy["visc"] = ["cp","mPa.s"]
        self.sTy["visc"] = ["cp","mpa.s"]
        self.For["visc"] = ["{:8.3f}","{:8.3f}"]
        self.Con["visc"] = [ [1.0, 0.0, 0.0, 1.0 ],
                             [1.0, 0.0, 0.0, 1.0 ] ]

#=======================================================================
#  Default Units for a Given Quantity Type
#=======================================================================

        self.defU = {'AREA':'ft2'     ,
                     'BG'  :'rcf/scf' ,'BO'  :'rcf/scf'  ,'BT'  :'rcf/scf'  ,
                     'CPG' :'Btu/lbmol.R','CPO' :'Btu/lbmol.R',
                     'DENG':'lb/ft3'  ,'DENO':'lb/ft3'   ,'DENS':'lb/ft3'   ,
                     'GGRV':'dim'     ,'GOR' :'scf/scf'  ,
                     'HEIG':'ft'      ,'HG'  :'Btu/lbmol','HO'  :'Btu/lbmol',
                     'IFT' :'dyne/cm' ,'KH'  :'md.ft'   ,
                     'MOLE':'frac'    ,'MREM':'frac'     ,
                     'MWO' :'lb/lbmol','MWG' :'lb/lbmol' ,
                     'PRES':'psia'    ,'PSAT':'psia'     ,
                     'RELV':'dim'     ,
                     'SKIN':'dim'     ,'SLIQ':'frac'     ,
                     'TEMP':'degr'    ,
                     'VISG':'cp'      ,'VISO':'cp'       ,'VOLG':'scf'      ,
                     'VOLL':'scf'     ,'VSWL':'dim'      ,'VFRC':'dim'      ,
                     'ZC1' :'frac'    ,'ZC7+':'frac'     ,'ZFAC':'dim'      }

#========================================================================
#  Min/Max by Unit Type in Default Units
#========================================================================

        self.Lim = {}

        self.Lim["area"] = [   0.0   ,  30.0E+9]  #-- ft2
        self.Lim["dens"] = [   0.01  ,  99.0   ]  #-- lb/ft3
        self.Lim["dim" ] = [   0.01  ,  99.0   ]  #-- dim
        self.Lim["frac"] = [   0.0   ,   1.0   ]  #-- frac (of 1)
        self.Lim["gfvf"] = [   1.0E-6,  10.0   ]  #-- rcf/scf
        self.Lim["ggrv"] = [   0.5   ,   2.0   ]  #-- wrt Air=1
        self.Lim["gor" ] = [   0.0   ,   1.0E+6]  #-- scf/scf
        self.Lim["kh"  ] = [   0.0   ,   1.0E+8]  #-- mD.ft
        self.Lim["leng"] = [  -5.0E+3,  30.0E+3]  #-- ft
        self.Lim["ofvf"] = [   0.9   ,   5.0   ]  #-- rcf/scf
        self.Lim["pres"] = [   0.1   ,  30.0E+3]  #-- psia
        self.Lim["temp"] = [ 450.0   , 999.0   ]  #-- degR
        self.Lim["visc"] = [   0.001 , 100.0E+6]  #-- cP

#========================================================================
#  Convert eXternal (User) Units to Internal Units
#========================================================================

    def X2I(self,xVal,usrU) :

        iVal = -999.0
        usrL = usrU.lower()

        for key,val in self.sTy.items() :
            if usrL in val :
                iA    = val.index(usrL)
                dCoef = self.Con[key][iA]
                iVal  = (dCoef[0]*xVal + dCoef[1])/(dCoef[2]*xVal + dCoef[3])
                break

        if iVal < -900.0 :
            print("X2I: usrU,usrL ",usrU,usrL)
            crash = 1.0/0.0

        return iVal

#========================================================================
#  Convert Internal Units to eXternal (User) Units
#========================================================================

    def I2X(self,iVal,usrU) :

        xVal = -999.0
        usrL = usrU.lower()

        for key,val in self.sTy.items() :
            if usrL in val :
                iA    = val.index(usrL)
                dCoef = self.Con[key][iA]
                xVal  = -(dCoef[3]*iVal - dCoef[1])/(dCoef[2]*iVal - dCoef[0])
                break

        if xVal < -900.0 :
            print("I2X: iVal,usrU,usrL ",iVal,usrU,usrL)
            crash = 1.0/0.0

        return xVal

#========================================================================
#  For a given User Type, Find the Default (Internal) Unit
#========================================================================

    def DUN(self,usrUni) : return self.defU.get(usrUni)

#========================================================================
#  Check that a quantity is valid,i.e. uMin < uInt < uMax
#========================================================================

    def CHU(self,dExt,usrU ) :

        qChk = False
        usrL = usrU.lower()
        dInt = self.X2I(dExt,usrL)  #-- Convert to Internal Units

        for key,val in self.sTy.items() :
            if usrL in val :
                if dInt > self.Lim[key][0] and \
                   dInt < self.Lim[key][1]     : qChk = True
                break

        return qChk

#========================================================================
#  Convert Type-A unit to Type-B unit via Internal unit
#========================================================================

    def A2B(self,aVal,aUni,bUni) :

        iVal = self.X2I(aVal,aUni)
        bVal = self.I2X(iVal,bUni)

        return bVal

#========================================================================
#  End of Module
#========================================================================
