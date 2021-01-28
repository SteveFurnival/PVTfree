
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Coeffients (and derivatives) of 2-Component EoS/Viscosity for
#  Blackoil Table Generation
#
#=======================================================================

#!/usr/bin/python3

from math import sqrt

#========================================================================
#  EOS Phase A-Coefficient
#========================================================================

def Acoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    aOil0 = clsBLK.EOS0["aOil"]
    aOilM = clsBLK.EOS1["aOil"]

    aGas0 = clsBLK.EOS0["aGas"]
    aGasM = clsBLK.EOS1["aGas"]

    aOil = aOil0*aOilM
    aGas = aGas0*aGasM

    aCof = xOil*xOil*aOil + 2.0*xOil*xGas*sqrt(aOil*aGas) + xGas*xGas*aGas

    return aCof
    
#========================================================================
#  d[EOS Phase A-Coefficient]/dAoil
#========================================================================

def dAdAo(xOil,clsBLK) :

    xGas = 1.0 - xOil

    aOil0 = clsBLK.EOS0["aOil"]
    aOilM = clsBLK.EOS1["aOil"]

    aGas0 = clsBLK.EOS0["aGas"]
    aGasM = clsBLK.EOS1["aGas"]

    aOil = aOil0*aOilM
    aGas = aGas0*aGasM

    dAdAo = xOil*xOil*aOil0 + xOil*xGas*sqrt(aOil*aGas)/aOilM

    return dAdAo
    
#========================================================================
#  d[EOS Phase A-Coefficient]/dAgas
#========================================================================

def dAdAg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    aOil0 = clsBLK.EOS0["aOil"]
    aOilM = clsBLK.EOS1["aOil"]

    aGas0 = clsBLK.EOS0["aGas"]
    aGasM = clsBLK.EOS1["aGas"]

    aOil = aOil0*aOilM
    aGas = aGas0*aGasM

#== srfhere: 1st term was xGas*xGas*aGas = Error? =====================    

    dAdAg = xGas*xGas*aGas0 + xOil*xGas*sqrt(aOil*aGas)/aGasM

    return dAdAg
    
#========================================================================
#  EOS Phase B-Coefficient
#========================================================================

def Bcoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    bOil0 = clsBLK.EOS0["bOil"]
    bOilM = clsBLK.EOS1["bOil"]

    bGas0 = clsBLK.EOS0["bGas"]
    bGasM = clsBLK.EOS1["bGas"]

    bOil = bOil0*bOilM
    bGas = bGas0*bGasM

    bCof = xOil*bOil + xGas*bGas

    return bCof

#========================================================================
#  d[EOS Phase B-Coefficient]/dBoil
#========================================================================

def dBdBo(xOil,clsBLK) :

    bOil0 = clsBLK.EOS0["bOil"]

    dBdBo = xOil*bOil0

    return dBdBo

#========================================================================
#  d[EOS Phase B-Coefficient]/dBgas
#========================================================================

def dBdBg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    bGas0 = clsBLK.EOS0["bGas"]

    dBdBg = xGas*bGas0

    return dBdBg

#========================================================================
#  EOS Phase C-Coefficient
#========================================================================

def Ccoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    bOil0 = clsBLK.EOS0["bOil"]
    sOilM = clsBLK.EOS1["sOil"]

    bGas0 = clsBLK.EOS0["bGas"]
    sGasM = clsBLK.EOS1["sGas"]

    cOil = bOil0*sOilM
    cGas = bGas0*sGasM

    cCof = xOil*cOil + xGas*cGas

    return cCof

#========================================================================
#  d[EOS Phase C-Coefficient]/dBoil
#========================================================================

def dCdBo(xOil,clsBLK) :

    bOil0 = clsBLK.EOS0["bOil"]
    sOilM = clsBLK.EOS1["sOil"]

    dCdBo = xOil*bOil0*sOilM

    return dCdBo

#========================================================================
#  d[EOS Phase C-Coefficient]/dBgas
#========================================================================

def dCdBg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    bGas0 = clsBLK.EOS0["bGas"]
    sGasM = clsBLK.EOS1["sGas"]

    dCdBg = xGas*bGas0*sGasM

    return dCdBg

#========================================================================
#  d[EOS Phase C-Coefficient]/dSoil
#========================================================================

def dCdSo(xOil,clsBLK) :

    bOil0 = clsBLK.EOS0["bOil"]
    bOilM = clsBLK.EOS1["bOil"]

    dCdSo = xOil*bOil0

    return dCdSo

#========================================================================
#  d[EOS Phase C-Coefficient]/dSgas
#========================================================================

def dCdSg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    bGas0 = clsBLK.EOS0["bGas"]
    bGasM = clsBLK.EOS1["bGas"]

    dCdSg = xGas*bGas0

    return dCdSg

#========================================================================
#  EOS Phase U-Coefficient
#========================================================================

def Ucoef(vOil,bOil,cOil) :

    uCof = vOil + cOil - bOil

    return uCof

#========================================================================
#  EOS Phase W-Coefficient
#========================================================================

def Wcoef(vOil,bOil,cOil) :

    vShf = vOil + cOil

    wCof = vShf*vShf + 2.0*bOil*vShf - bOil*bOil

    return wCof

#========================================================================
#  Viscosity Phase R-Coefficient
#========================================================================

def Rcoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    rOil0 = clsBLK.VIS0["rOil"]
    rOilM = clsBLK.VIS1["rOil"]

    rGas0 = clsBLK.VIS0["rGas"]
    rGasM = clsBLK.VIS1["rGas"]

    rOil = rOil0*rOilM
    rGas = rGas0*rGasM

    rCof = xOil*rOil + xGas*rGas

    return rCof

#========================================================================
#  d[Viscosity Phase R-Coefficient]/dRoil
#========================================================================

def dRdRo(xOil,clsBLK) :

    rOil0 = clsBLK.VIS0["rOil"]

    dRdRo = xOil*rOil0

    return dRdRo

#========================================================================
#  d[Viscosity Phase R-Coefficient]/dRgas
#========================================================================

def dRdRg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    rGas0 = clsBLK.VIS0["rGas"]

    dRdRg = xGas*rGas0

    return dRdRg

#========================================================================
#  Viscosity Phase Eta-Coefficient
#========================================================================

def Ecoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    eOil0 = clsBLK.VIS0["eOil"]
    eOilM = clsBLK.VIS1["eOil"]

    eGas0 = clsBLK.VIS0["eGas"]
    eGasM = clsBLK.VIS1["eGas"]

    eOil = eOil0*eOilM
    eGas = eGas0*eGasM

    eCof = xOil*eOil + xGas*eGas

    return eCof

#========================================================================
#  d[Viscosity Phase Eta-Coefficient]/dEoil
#========================================================================

def dEdEo(xOil,clsBLK) :

    eOil0 = clsBLK.VIS0["eOil"]

    dEdEo = xOil*eOil0

    return dEdEo

#========================================================================
#  d[Viscosity Phase Eta-Coefficient]/dRgas
#========================================================================

def dEdEg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    eGas0 = clsBLK.VIS0["eGas"]

    dEdEg = xGas*eGas0

    return dEdEg

#========================================================================
#  Viscosity Phase Mu-Coefficient
#========================================================================

def Mcoef(xOil,clsBLK) :

    xGas = 1.0 - xOil

    uOil0 = clsBLK.VIS0["uOil"]
    uOilM = clsBLK.VIS1["uOil"]

    uGas0 = clsBLK.VIS0["uGas"]
    uGasM = clsBLK.VIS1["uGas"]

    uOil = uOil0*uOilM
    uGas = uGas0*uGasM

    uCof = xOil*uOil + xGas*uGas

    return uCof

#========================================================================
#  d[Viscosity Phase Mu-Coefficient]/dUoil
#========================================================================

def dMdUo(xOil,clsBLK) :

    uOil0 = clsBLK.VIS0["uOil"]

    dMdUo = xOil*uOil0

    return dMdUo

#========================================================================
#  d[Viscosity Phase Mu-Coefficient]/dUgas
#========================================================================

def dMdUg(xOil,clsBLK) :

    xGas = 1.0 - xOil

    uGas0 = clsBLK.VIS0["uGas"]

    dMdUg = xGas*uGas0

    return dMdUg

#========================================================================
#  End of Module
#========================================================================
