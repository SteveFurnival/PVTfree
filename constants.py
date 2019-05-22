
#========================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Define "Constants" to use Globally in PVTfree
#
#========================================================================

#!/usr/bin/python3

macEPS = 2.22E-16   #-- 'Double Precision' Minimum Number
mGDEM1 = 4          #-- SS Iteration/GDEM Step for Psat & Stab
mGDEM2 = 5          #-- SS Iteration/GDEM Step for Flash

Third  = 1.0/3.0    #-- Used in cubicSolver

btuCon =   1.98588  #-- Gas Constant          [Btu/(lbmol.degR)]
denAir =   0.0765   #-- Density of Air at STC [lb/ft3]
denWat =  62.428    #-- Density of Pure Water [lb/ft3]
dF2dR  = 459.67     #-- 0 degF                [degR]
gasCon =  10.7315   #-- Gas Constant          [psia.ft3/(lbmol.degR)]
molAir =  28.97     #-- Mole Weight of Air    [lb/lbmol]
pStand =  14.6959   #-- Standard Pressure     [psia]
tStand = 519.67     #-- Std Temp (60 degF)    [degR]
volMol = 379.483    #-- Vol 1 lbmol Ideal Gas [scf]

lbcMwX = -0.5       #-- LBC Viscosity Mw-Exp
lbcPcX = -2.0/3.0   #-- LBC Viscosity Pc-Exp
lbcTcX =  1.0/6.0   #-- LBC Viscosity Tc-Exp

#========================================================================
#  End of Module
#========================================================================

