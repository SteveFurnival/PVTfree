
#========================================================================
#  Hassanzadeh, Pooladi-Darvish, Elsharkawy, Keith & Leonenko
#  "Predicting PVT data for CO2-brine mixtures for black-oil simulation
#   of CO2 geological storage"
#  Int J of Greenhouse Gas Control, 2, (2008), 65-77
#  https://scholar.harvard.edu/files/davidkeith/files/81.hassanzadeh.pvtCO2blackoil.e.pdf
#
#  Also Spycher & Pruess
#  http://escholarship.org/uc/item/0kp5q3vn
#========================================================================

#!/usr/bin/python3

from math import acos,copysign,cos,exp,log,sqrt

#== Constants [and Global Variables - make class vars?] ===============

dC2dK = 273.15          #-- degC to Kelvin

MwS = 58.44             #-- Mole Weight Salt (NaCl)
MwW = 18.0152           #-- Mole Weight Water
MwC = 44.01             #-- Mole Weight CO2

aCO20 =  7.54E7         #-- aCO2 = aCO20 + tKel*aCO21
aCO21 = -4.13E4

aH2OCO2 = 7.89E7        #-- bar.cm3/sqrt(K)/gmol2

bCO2  = 27.80           #-- cm3/gmol
bH2O  = 18.18

vCO2  = 32.6            #-- cm3/gmol
vH2O  = 18.1

Rcon  = 83.1447         #-- Gas Constant [bar.cm3/gmol/K]
Mcon  = 55.508          #-- NaCl Molality Constant [1000/MwS gmol/kg]

Third = 1.0/3.0

pStnd =   1.01325       #-- Standard Pressure    [bar]
tStnd = 288.71          #-- Standard Temperature [Kelvin]

#========================================================================
#  Main Program
#========================================================================

def runMain() :

#-- Set some test input values --------------------------------------

    pBar = [ 50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]
    tCel = [ 30.0, 60.0, 90.0]
    
    PPM  = [ 0.0, 55249.0,104712.0,189573.0]  #-- Molality [0,1,2,4 gmol/kg]

#----------------------------------------------------------------------
#  Test Range of (salt,T,p)
#----------------------------------------------------------------------

    for salt in PPM :

#-- Molality [gmol/m3] ----------------------------------------------

        molaL = PPM2Molal(salt)

        print("PPM,molaL {:9.1f} {:7.4f}".format(salt,molaL))

#-- Temperature -----------------------------------------------------        

        for temp in tCel :

            tKel = temp + dC2dK

            gStar = lambZeta(tKel,molaL)

#-- Pressure --------------------------------------------------------            

            for pRes in pBar :

                calcProps(pRes,tKel,molaL,gStar)

#== No return values ==================================================

    return

#========================================================================
#  Calculate Properties
#========================================================================

def calcProps(pBar,tKel,molaL,gStar) :

    mol2 = 2.0*molaL    #-- '2' = Stoichiometry of Na+Cl

    m2pM = mol2 + Mcon  #-- Mcon = 1000/Mw(water) = 55.508

    tK12 = sqrt(tKel)
    tCel =      tKel - dC2dK

    RT  =  Rcon*tKel
    RTp  =   RT/pBar

#== Redlich-Kwong EoS for Gas [CO2] Phase =============================    

    aCO2 = aCO2temp(tKel)

    aGas = aCO2
    bGas = bCO2

#-- Cubic Coeffs ----------------------------------------------------

    apT12 = aGas/(pBar*tK12)

    E2 = - RT/pBar
    E1 =   E2*bGas + apT12 - bGas*bGas
    E0 = -    bGas * apT12

#== Solve Cubic =======================================================

    Vgas = solveCubic(E2,E1,E0)     #-- S&P Eqn.A2

#----------------------------------------------------------------------
#  Fugacity Coefficients, see Spycher & Pruess Appendix A
#-----------------------------------------------------------------------

    Vmb = Vgas - bGas
    Vpb = Vgas + bGas

    Zm1  = Vgas/RTp - 1.0       #-- pV/RT - 1

    cCon = - log(Vmb/RTp)
    aCon = aGas/(RT*tK12*bGas)*log(Vgas/Vpb)
    aMul =    2.0*aCon /aGas
    bMul = (Zm1 - aCon)/bGas

    phiCO2 = exp(cCon + aCO2   *aMul + bCO2*bMul)   #-- S&P Eqn.A8
    phiH2O = exp(cCon + aH2OCO2*aMul + bH2O*bMul)

#== (A,B) Coefficients ================================================

    pCon = (pBar - 1.0)/RT
    
    KCO2 = K0CO2(tCel)*exp(pCon*vCO2)   #-- S&P Eqn.05
    KH2O = K0H2O(tCel)*exp(pCon*vH2O)

    rA = pBar*phiH2O/ KH2O              #-- S&P Eqn.10 (actually 1/A)
    BP = pBar*phiCO2/(KCO2*Mcon*gStar)  #-- S&P Eqn.17

    yH2O = Mcon*(1.0 - BP)/((rA - BP)*m2pM + mol2*BP)  #-- S&P Eqn.B7
    xCO2 =   BP*(1.0 - yH2O)                           #-- S&P Eqn.B2

    mCO2 = xCO2*m2pM/(1.0 - xCO2)       #-- S&P Eqn.B6

    xSal = mol2/(m2pM + mCO2)           #-- S&P Eqn.B3

    xH2O = 1.0 - xCO2 - xSal            #-- S&P Wqn.B5
    yCO2 = 1.0 - yH2O

    #print("p,T,yH2O,xCO2,xH2O,xSal,mCO2 {:8.3f} {:8.3f} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f}".format(pBar,tCel,yH2O,xCO2,xH2O,xSal,mCO2))

#== Brine Density [Rowe & Chou Correlation] ===========================

    denB,comB = calcRoweChouDen(molaL,tKel,pBar)

#== Brine Viscosity [cP] and Viscosibility [cP/bar] ===================

    visB,visC = calcKestinVisc(molaL,tKel,pBar)

#-- Brine Mole Weight -----------------------------------------------

    xBri = xH2O + xSal

    MwB = (xH2O*MwW + xSal*MwS)/xBri

    Vphi = partMolVolCO2(tCel)  #-- CO2 partial molar volume

#== Saturated Brine Density Hass Eqn.20 ===============================

    rSat = 1000.0*(xBri*MwB + xCO2*MwC)/(1000.0*xBri*MwB/denB+xCO2*Vphi)

    print("p,T,denB,visB,MwB,Vphi,rSat {:8.3f} {:8.3f} {:8.3f} {:8.5f} {:8.3f} {:8.3f} {:8.3f}".format(pBar,tKel,denB,visB,MwB,Vphi,rSat))

#== No return value ===================================================

    return

#========================================================================
#  Solve the Cubic x^3 + E2 x^2 + E1 x + E0 = 0
#  See Numerical Recipes Section 5.6 [E2 = a, E1 = b, E0 = c]
#  http://www.aip.de/groups/soe/local/numres/bookfpdf/f5-6.pdf
#========================================================================

def solveCubic(E2,E1,E0) :

    F2 = E2*Third                     #--  a/3
    F1 = E1*Third                     #--  b/3

    F22 = F2*F2                       #-- (a/3)^2

    Q =    F22 - F1                   #-- Eqn.(5.6.10)
    R = F2*F22 - F1*F2*1.5 + 0.5*E0   #-- Eqn.(5.6.10)

#-- Discriminant ----------------------------------------------------    

    D = R*R - Q*Q*Q

#== One Real Root =====================================================
    
    if   D > 0.0 :
        
        argA = abs(R) + sqrt(D)        #-- Arg of Eqn.(5.6.13)
        A13  = argA**Third
        A    = -copysign(A13,R)        #-- Eqn.(5.6.15)
        
        if A == 0.0 : B = 0.0          #-- Eqn.(5.6.16)
        else        : B = Q/A
        
        x1 = A + B                     #-- Eqn.(5.6.17)
        x3 = x1

        nR = 1

#== Three Real Roots: Ignore Central Root as Unphysical in EoS ========
        
    elif D < 0.0 :
        
        sqrtQ = sqrt(Q)
        theta = acos(R/(Q*sqrtQ))      #-- Eqn.(5.6.11)
        cTh13 =  cos(Third*theta)
        
        x1 = -2.0*sqrtQ*cTh13          #-- Eqn.(5.6.12-x1)
        
#-- Next expression have used cos(A-B) = cosA.cosB + sinA.sinB where
#-- cos(-2.PI/3) = -1/2, sin(-2.PI/3) = -sqrt(3)/2 and sinA = sqrt(1-cosA^2)

        x3 = sqrtQ*(cTh13+NP.sqrt(3.0*(1.0-cTh13*cTh13)))

        nR = 3

#== Error in NR Text: see http://numerical.recipes/forum/showthread.php?t=865
        
    else :
        
        sqrtQ = sqrt(Q)
        
        if   R > 0.0 :
            x1 = -2.0*sqrtQ
            x3 =      sqrtQ
        elif R < 0.0 :
            x1 =     -sqrtQ
            x3 =  2.0*sqrtQ
        else :
            x1 = 0.0
            x3 = 0.0

        nR = 3

#== Return Values =====================================================

    vHig = x3 - F2

    if nR > 1 : print("3-Roots Detected")

    return vHig

#========================================================================
#  RK-EOS a-Coeff of CO2 is Temperature-Dependent
#========================================================================

def aCO2temp(tKel) : return aCO20 + aCO21*tKel

#========================================================================
#  log kappa0-CO2(gas): Eqn.( 9)
#========================================================================

def K0CO2(tCel) :
    return 10.0**(1.189+tCel*(1.304E-2-tCel*5.446E-5))

#========================================================================
#  log kappa0-H2O: Eqn.(11)
#========================================================================

def K0H2O(tC) :
    return 10.0**(-2.209+tC*(3.097E-2+tC*(-1.098E-4+tC*2.048E-7)))

#========================================================================
#  log kappa0-CO2(liq): Eqn.(10)  Only if T<Tcrit and P<Psat
#  Shouldn't be required in a saline aquifer; keep just in case ...
#========================================================================

def logK0CO2Liq(tCel) :
    return 1.169+tCel*(1.368E-2-tCel*5.380E-5)

#========================================================================
#  Partial Molar Volume of CO2 [cm3/gmol]; Hassanzadeh Eqn.(21)
#========================================================================

def partMolVolCO2(tC) :
    return 37.51+tC*(-9.585E-2+tC*(8.740E-4+tC*-5.044E-7*tC))

#========================================================================
#  Lambda/Zeta Spycher & Pruess (no Pressure dependency!)
#========================================================================

def lambZeta(tKel,molaL) :

    rTk = 1.0/tKel

    lamb = 2.217E-4*tKel + rTk*(  1.074 + rTk*2648.0)
    zeta = 1.300E-5*tKel + rTk*(-20.120 + rTk*5259.0)

    gStar = exp(molaL*(2.0*lamb + molaL*zeta))

    return gStar

#========================================================================
#  PPM (gm/kg) to Molality (gmol/m3)
#  Whitson & Brule Eqn.(9.2c); assumes brine is NaCl
#========================================================================

def PPM2Molal(PPM) :

    if PPM > 0.0 :
        molaL = 17.1/(1.0E+6/PPM-1.0)
    else :
        molaL = 0.0

    return molaL

#========================================================================
#  Salt (NaCl) Mass Fraction from Molality (reverse of above)
#========================================================================

def molal2MF(molaL) :

    if molaL > 0.0 :
        mF = 1.0/(17.1/molaL+1.0)
    else :
        mF = 0.0

    return mF        

#========================================================================
#  Rowe and Chou Density Correlation
#========================================================================

def calcRoweChouDen(molaL,tKel,pBar) :

    pCon = 0.9806657         #-- Pressure in kgf/cm2 (from Google!!)

#-- Molality [gmol/kg] => mass fraction -----------------------------

    mFrc = molal2MF(molaL)

    mFr2 = mFrc*mFrc

    pkgF = pBar/pCon    
    pkg2 = pkgF*pkgF

    T1 = tKel
    
    T2 = T1*T1
    T3 = T1*T2

    a1,da1dT = rc_A1(T1,T2,T3)
    a2,da2dT = rc_A2(T1,T2,T3)
    a3,da3dT = rc_A3(T1,T2,T3)
    a4,da4dT = rc_A4(T1,T2,T3)
    a5,da5dT = rc_A5(T1,T2,T3)
    a6,da6dT = rc_A6(T1,T2,T3)
    a7,da7dT = rc_A7(T1,T2,T3)
    a8,da8dT = rc_A8(T1,T2,T3)

#== Density [kg/m3] ===================================================

    dRec = a1 - a2*pkgF - a3*pkg2 \
              + a4*mFrc + a5*mFr2 \
              - a6*pkgF*mFrc - a7*pkgF*mFr2 - 0.5*a8*pkg2*mFrc

    denW = 1000.0/dRec  #-- Density [kg/m3]

#== Derivatives [wrt Temp,Pres & Salinity] ============================

    dadT = da1dT - da2dT*pkgF - da3dT*pkg2 \
                 + da4dT*mFrc + da5dT*mFr2 \
                 - da6dT*pkgF*mFrc - da7dT*pkgF*mFr2 - 0.5*da8dT*pkg2*mFrc

    dDdT = - denW*dadT/dRec

    dDdP = - denW*(- a2 - 2.0*a3*pkgF - a6*mFrc -     a7*mFr2      -     a8*pkgF*mFrc)
    dDdS = - denW*(  a4 + 2.0*a5*mFrc - a6*pkgF - 2.0*a7*pkgF*mFrc - 0.5*a8*pkg2)

#== Return Values =====================================================    

    dDdP = dDdP*pCon                #-- 1/kgf/cm2 -> 1/bar (inverse of way-in!)

    comP = dDdP/denW                #-- Compressibility [1/bar]

    return denW,comP

#========================================================================
#  Rowe and Chou A1-Coefficient
#========================================================================

def rc_A1(T1,T2,T3) :

    a1 =  5.916365E+00
    a2 = -1.035794E-02
    a3 =  9.270048E-06
    a4 = -1.127522E+03
    a5 =  1.006741E+05

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A2-Coefficient
#========================================================================

def rc_A2(T1,T2,T3) :

    a1 =  5.2049100E-03
    a2 = -1.0482101E-05
    a3 =  8.3285320E-09
    a4 = -1.1702939E+00
    a5 =  1.0227830E+02

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A3-Coefficient
#========================================================================

def rc_A3(T1,T2,T3) :

    a1 =  1.185470E-08
    a2 = -6.599143E-11
    a3 =  0.0
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A4-Coefficient
#========================================================================

def rc_A4(T1,T2,T3) :

    a1 = -2.51660E+00
    a2 =  1.11766E-02
    a3 = -1.70522E-05
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A5-Coefficient
#========================================================================

def rc_A5(T1,T2,T3) :

    a1 =  2.84851E+00
    a2 = -1.54305E-02
    a3 =  2.23982E-05
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A6-Coefficient
#========================================================================

def rc_A6(T1,T2,T3) :

    a1 = -1.48140E-03
    a2 =  8.29639E-06
    a3 = -1.24690E-08
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A7-Coefficient
#========================================================================

def rc_A7(T1,T2,T3) :

    a1 =  2.7141E-03
    a2 = -1.5391E-05
    a3 =  2.2655E-08
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Rowe and Chou A8-Coefficient
#========================================================================

def rc_A8(T1,T2,T3) :

    a1 =  6.2158E-07
    a2 = -4.0075E-09
    a3 =  6.5972E-12
    a4 =  0.0
    a5 =  0.0

    a    = a1 + a2*T1 +     a3*T2 + a4/T1 +     a5/T2
    dadT =      a2    + 2.0*a3*T1 - a4/T2 - 2.0*a5/T3

    return a,dadT

#========================================================================
#  Kestin Viscosity Correlation
#========================================================================

def calcKestinVisc(molaL,tKel,pBar) :

    pMPa = 0.1*pBar
    tCel =     tKel - dC2dK

    mu0w20 = 1.002  #-- Pure water at 20 degC

    a = kestinA(molaL)
    b = kestinB(molaL)
    c = kestinC(tCel)

    d = a + b*c

    mu0w = mu0w20*10.0**c

    muws = mu0w  *10.0**d

    pCor = 0.001*(0.8 + 0.01*(tCel - 90.0)*exp(-0.25*molaL))

    visW = muws*(1.0 + pCor*pMPa)   #--   Viscosity    [cP]

    dvdM = muws*pCor                #-- d[Viscosity/dp [cP/MPa]
    dvdp = dvdM/0.1                 #--                [cP/bar]

    visC = dvdp/visW                #-- Viscosibility [= (1/mu)(dmu/dp) 1/bar]

    return visW,visC

#========================================================================
#  Kestin A-Sum
#========================================================================

def kestinA(M) :

    a1 =  3.324E-02
    a2 =  3.624E-03
    a3 = -1.879E-04

    A = M*(a1 + M*(a2 + M*a3))

    return A

#========================================================================
#  Kestin B-Sum
#========================================================================

def kestinB(M) :

    b1 = -3.96E-02
    b2 =  1.02E-02
    b3 = -7.02E-04

    B = M*(b1 + M*(b2 + M*b3))

    return B

#========================================================================
#  Kestin C-Sum
#========================================================================

def kestinC(degC) :

    c1 =  1.2378E+00
    c2 = -1.3030E-03
    c3 =  3.0600E-06
    c4 =  2.5500E-08

    T = 20.0 - degC

    sumC = T*(c1 + T*(c2 + T*(c3 + T*c4)))

    C = sumC/(96.0 + degC)

    return C

#========================================================================
#  Vesovic et al CO2 Viscosity,
#  Hassanzadeh et al, Int.J.of Greenhouse Gas Control, 2, (2008), pp.65-77
#  https://scholar.harvard.edu/files/davidkeith/files/81.hassanzadeh.pvtco2blackoil.e.pdf
#========================================================================

def viscCO2(tKel,kgm3) :

    tStar = tKel/251.196
    
    lnT = log(tStar)
    tS2 = tStar*tStar
    tS3 = tStar*tS2

#-- Eqn.(C6) --------------------------------------------------------

    lnPsi = 0.235156 + lnT*(-0.491266     + lnT*(5.211155E-02 \
                     + lnT*( 5.347906E-02 - lnT* 1.537102E-02)))

    eta0 = 1.00697*sqrt(tKel)/exp(lnPsi)

#-- Excess Viscosity ------------------------------------------------

    rho2 = kgm3*kgm3
    rho4 = rho2*rho2
    rho6 = rho2*rho4
    rho8 = rho2*rho6

    etaX = 0.4071119E-02*rhoV     + 0.7198037E-04*rho2 \
         + 0.2411697E-16*rho6/tS3 + 0.2971072E-22*rho8 \
         - 0.1627888E-22*rho8/tStar

#-- Critical Component ----------------------------------------------

    #etaC = kgm3*(5.5934E-03 + kgm3*(6.1757E-05 + rho2*2.6430E-11))

#== Total Viscosity [uPa.s] ===========================================

    etaU = eta0 + etaX #+etaC

#  NB. Hassanzedah is wrong with 'critical' term - it's way more complicated
#      in 2nd Vesovic & Wakeham paper.  But it's also very small effect,
#      even near (Pc,Tc) so best to ignore

    #print("eta0,etaX,etaC {:12.5e} {:12.5e} {:12.5e}".format(eta0,etaX,etaC))

#== Return value [cP] =================================================

    return 1.0E-3*etaU

#========================================================================
#  Run Main
#========================================================================

runMain()

#========================================================================
#  End of Module
#========================================================================
