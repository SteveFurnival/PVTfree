
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTFree.py for further details
#
#  Defines the All the (Library) Properies used by PVTfree
#
#=======================================================================

#!/usr/bin/python3

#========================================================================
#  Library of Component and Experiment Properties 
#========================================================================

class classLIB :

    def __init__(self) :

#== Single Carbon Number [SCN] Mole Weight & Specific Gravity [C6, C7, ..., C45]        
           
        self.scnMW = {'C6' : 84.0,'C7' : 96.0,'C8' :107.0,'C9' :121.0,'C10':134.0,\
                      'C11':147.0,'C12':161.0,'C13':175.0,'C14':190.0,'C15':206.0,\
                      'C16':222.0,'C17':237.0,'C18':251.0,'C19':263.0,'C20':275.0,\
                      'C21':291.0,'C22':300.0,'C23':312.0,'C24':324.0,'C25':337.0,\
                      'C26':349.0,'C27':360.0,'C28':372.0,'C29':382.0,'C30':394.0,\
                      'C31':404.0,'C32':415.0,'C33':426.0,'C34':437.0,'C35':445.0,\
                      'C36':456.0,'C37':464.0,'C38':475.0,'C39':484.0,'C40':495.0,\
                      'C41':502.0,'C42':512.0,'C43':521.0,'C44':531.0,'C45':539.0 }

        self.scnSG = {'C6' :0.690,'C7' :0.727,'C8' :0.749,'C9' :0.768,'C10':0.782,\
                      'C11':0.793,'C12':0.804,'C13':0.815,'C14':0.826,'C15':0.836,\
                      'C16':0.843,'C17':0.851,'C18':0.856,'C19':0.861,'C20':0.866,\
                      'C21':0.871,'C22':0.876,'C23':0.881,'C24':0.885,'C25':0.888,\
                      'C26':0.892,'C27':0.896,'C28':0.899,'C29':0.902,'C30':0.905,\
                      'C31':0.909,'C32':0.912,'C33':0.915,'C34':0.917,'C35':0.920,\
                      'C36':0.922,'C37':0.925,'C38':0.927,'C39':0.929,'C40':0.931,\
                      'C41':0.933,'C42':0.934,'C43':0.936,'C44':0.938,'C45':0.940 }

#== Library (N2,CO2,H2S,C1 ... NC5) or SCN Component? =================

        self.typeC = {'N2' :'L','CO2':'L','H2S':'L',\
                      'C1' :'L','C2' :'L','C3' :'L',\
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

        self.lMolW = {'N2' :  28.013,'CO2':  44.010,'H2S':  34.076,\
                      'C1' :  16.043,'C2' :  30.070,'C3' :  44.097,\
                      'IC4':  58.124,'NC4':  58.124,'IC5':  72.151,'NC5': 72.151}
        self.lTcrt = {'N2' :-232.400,'CO2':  87.900,'H2S': 212.700,\
                      'C1' :-116.630,'C2' :  90.090,'C3' : 206.010,\
                      'IC4': 274.980,'NC4': 305.650,'IC5': 369.100,'NC5':385.700}
        self.lPcrt = {'N2' : 493.000,'CO2':1070.600,'H2S':1306.000,\
                      'C1' : 667.800,'C2' : 707.800,'C3' : 616.300,\
                      'IC4': 529.100,'NC4': 550.700,'IC5': 490.400,'NC5':488.600}
        self.lAcnF = {'N2' :   0.045,'CO2':   0.231,'H2S':   0.100,\
                      'C1' :   0.012,'C2' :   0.091,'C3' :   0.145,\
                      'IC4':   0.176,'NC4':   0.193,'IC5':   0.227,'NC5':  0.251}
        self.lVcrt = {'N2' :   1.443,'CO2':   1.505,'H2S':   1.564,\
                      'C1' :   1.590,'C2' :   2.370,'C3' :   3.250,\
                      'IC4':   4.208,'NC4':   4.080,'IC5':   4.899,'NC5':  4.870}
        self.lSpcG = {'N2' :   0.470,'CO2':   0.507,'H2S':   0.500,\
                      'C1' :   0.330,'C2' :   0.450,'C3' :   0.508,\
                      'IC4':   0.563,'NC4':   0.584,'IC5':   0.625,'NC5':  0.631}
        self.lTboi = {'N2' :-320.400,'CO2':-109.300,'H2S': -76.600,\
                      'C1' :-258.690,'C2' :-127.480,'C3' : -43.670,\
                      'IC4':  10.900,'NC4':  31.100,'IC5':  82.120,'NC5': 96.920}
        self.lPara = {'N2' :  41.000,'CO2':  70.000,'H2S':  41.000,\
                      'C1' :  77.000,'C2' : 108.000,'C3' : 150.300,\
                      'IC4': 181.500,'NC4': 189.900,'IC5': 225.000,'NC5':231.500}

        self.lCpAI = {'N2' : 7.45E+00,'CO2': 4.73E+00,'H2S': 7.63E+00,\
                      'C1' : 4.60E+00,'C2' : 1.29E+00,'C3' :-1.01E+00,\
                      'IC4':-3.32E-01,'NC4': 2.27E+00,'IC5':-2.28E+00,'NC5':-8.67E-01}
        self.lCpBI = {'N2' :-1.80E-03,'CO2': 9.75E-03,'H2S': 1.91E-04,\
                      'C1' : 6.92E-03,'C2' : 2.37E-02,'C3' : 4.07E-02,\
                      'IC4': 5.11E-02,'NC4': 4.40E-02,'IC5': 6.73E-02,'NC5': 6.47E-02}
        self.lCpCI = {'N2' : 1.98E-06,'CO2':-4.13E-06,'H2S': 1.79E-06,\
                      'C1' : 8.83E-07,'C2' :-5.12E-06,'C3' :-1.17E-05,\
                      'IC4':-1.36E-05,'NC4':-8.17E-06,'IC5':-2.01E-05,'NC5':-1.90E-05}
        self.lCpDI = {'N2' :-4.79E-10,'CO2': 7.01E-10,'H2S':-4.82E-10,\
                      'C1' :-4.64E-10,'C2' : 3.57E-10,'C3' : 1.32E-09,\
                      'IC4': 1.19E-09,'NC4':-1.16E-10,'IC5': 2.35E-09,'NC5': 2.17E-09}

#-- Volume Shifts: Peng-Robinson (PR) and (Soave) Redlich-Kwong (RK)

        self.lVSPR = {'N2' :-0.1927,'CO2':-0.0817,'H2S':-0.1288,\
                      'C1' :-0.1595,'C2' :-0.1134,'C3' :-0.0863,\
                      'IC4':-0.0844,'NC4':-0.0675,'IC5':-0.0608,'NC5':-0.0390}
        self.lVSRK = {'N2' :-0.0079,'CO2': 0.0833,'H2S': 0.0466,\
                      'C1' : 0.0234,'C2' : 0.0605,'C3' : 0.0825,\
                      'IC4': 0.0830,'NC4': 0.0975,'IC5': 0.1022,'NC5': 0.1209}

#-- BICs for N2-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

        self.lKN2P = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.130,\
                      'C1' : 0.025,'C2' : 0.010,'C3' : 0.090,\
                      'IC4': 0.095,'NC4': 0.095,'IC5': 0.100,'NC5': 0.100}
        self.lKN2S = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.120,\
                      'C1' : 0.020,'C2' : 0.060,'C3' : 0.080,\
                      'IC4': 0.080,'NC4': 0.080,'IC5': 0.080,'NC5': 0.080}

#-- BICs for CO2-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

        self.lKCOP = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.135,\
                      'C1' : 0.105,'C2' : 0.130,'C3' : 0.125,\
                      'IC4': 0.120,'NC4': 0.115,'IC5': 0.115,'NC5': 0.115}
        self.lKCOS = {'N2' : 0.000,'CO2': 0.000,'H2S': 0.120,\
                      'C1' : 0.120,'C2' : 0.150,'C3' : 0.150,\
                      'IC4': 0.150,'NC4': 0.150,'IC5': 0.150,'NC5': 0.150}

#-- BICs for H2S-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

        self.lKHSP = {'N2' : 0.130,'CO2': 0.135,'H2S': 0.000,\
                      'C1' : 0.070,'C2' : 0.085,'C3' : 0.080,\
                      'IC4': 0.075,'NC4': 0.075,'IC5': 0.070,'NC5': 0.070}
        self.lKHSS = {'N2' : 0.120,'CO2': 0.120,'H2S': 0.000,\
                      'C1' : 0.080,'C2' : 0.070,'C3' : 0.070,\
                      'IC4': 0.060,'NC4': 0.060,'IC5': 0.060,'NC5': 0.060}

#-- BICs for C1-Other for Peng-Robinson (P) and Soave-Redlich-Kwong (S)

        self.lKC1P = {'N2' : 0.025,'CO2': 0.105,'H2S': 0.070,\
                      'C1' : 0.000,'C2' : 0.000,'C3' : 0.000,\
                      'IC4': 0.000,'NC4': 0.000,'IC5': 0.000,'NC5': 0.000}
        self.lKC1S = {'N2' : 0.020,'CO2': 0.120,'H2S': 0.080,\
                      'C1' : 0.000,'C2' : 0.000,'C3' : 0.000,\
                      'IC4': 0.000,'NC4': 0.000,'IC5': 0.000,'NC5': 0.000}
       
#----------------------------------------------------------------------
#  List of Allowed Observations/Calculations by Experiment Type
#----------------------------------------------------------------------

#-- Independent Variables by Exp ------------------------------------        
           
        self.INDshrt = { 'CCE':['PRES'],
                         'CVD':['PRES'],
                         'DLE':['PRES'],
                         'SEP':['PRES','TEMP'],
                         'FLS':['PRES','TEMP'],
                         'SAT':['TEMP'],
                         'SWL':['MOLE'],
                         'GRD':['HEIG'] }

        self.INDlong = { 'CCE':['Pressure'],
                         'CVD':['Pressure'],
                         'DLE':['Pressure'],
                         'SEP':['Pressure','Temperature'],
                         'FLS':['Pressure','Temperature'],
                         'SAT':['Temperature'],
                         'SWL':['Moles'],
                         'GRD':['Height'] }

#-- Observed Variables by Exp [Short Mnemonics & Long Descriptions]        

        self.OBSshrt = { 'CCE':['RELV','SLIQ','ZFAC','DENO','VISO','VISG'],
                         'CVD':['MREM','SLIQ','ZFAC'],
                         'DLE':['BO'  ,'GOR' ,'DENO','BT'  ,'BG'  ,'ZFAC','GGRV','VISO'],
                         'SEP':['BO'  ,'GOR' ,'DENO','GGRV'],
                         'FLS':['ZFAC','DENO','VFRC'],
                         'SAT':['PSAT','ZFAC','DENO'],
                         'SWL':['PSAT','VSWL'],
                         'GRD':['PRES','PSAT','DENS','ZC1' ,'ZC7+'] }

        self.OBSlong = { 'CCE':['Relative Volume','Liquid Dropout',
                                'Gas Z-Factor'   ,'Oil Density',
                                'Oil Viscosity'  ,'Gas Viscosity'],
                         'CVD':['Moles Removed'  ,'Liquid Dropout',
                                'Gas Z-Factor'],
                         'DLE':['Oil Form Vol Factor'  ,
                                'Gas-Oil-Ratio'  ,'Oil Density',
                                'Tot Form Vol Factor'  ,
                                'Gas Form Vol Factor'  ,
                                'Gas Z-Factor'   ,'Gas Gravity',
                                'Oil Viscosity'],
                         'SEP':['Oil Formation Volume Factor'  ,
                                'Gas-Oil-Ratio'  ,'Oil Density',
                                'Gas Gravity'],
                         'FLS':['Gas Z-Factor'   ,'Oil Density',
                                'Vapour Fraction'],
                         'SAT':['Saturation Pressure',
                                'Gas Z-Factor'   ,'Oil Density'],
                         'SWL':['Saturation Pressure',
                                'Swelling Factor'],
                         'GRD':['Pressure'       ,'Saturation Pressure',
                                'Density'        ,'C1 Mole Fraction' ,
                                'C7+ Mole Fraction'] }

#-- Calculated Variables by Exp -------------------------------------

        self.CALshrt = { 'CCE':['RELV','SLIQ','ZFAC','DENO','VISO','VISG',
                                'MWO' ,'MWG' ,'DENG',
                                'IFT' ,'HO'  ,'HG'  ,'CPO' ,'CPG' ],
                         'CVD':['MREM','SLIQ','ZFAC',
                                'MWO' ,'MWG' ,'DENO','DENG','VISO','VISG'],
                         'DLE':['BO'  ,'GOR' ,'DENO','BT'  ,'BG'  ,'ZFAC','GGRV','VISO',
                                'MWO' ,'MWG' ,'DENG','VISO','VISG'],
                         'SEP':['BO'  ,'GOR' ,'DENO','GGRV',
                                'MWO' ,'MWG' ,'DENG','VISO','VISG'],
                         'FLS':['ZFAC','DENO','VFRC',
                                'MWO' ,'MWG' ,'DENG','VISO','VISG'],
                         'SAT':['PSAT','ZFAC','DENO',
                                'MWO' ,'MWG' ,'DENG','VISO','VISG'],
                         'SWL':['PSAT','VSWL',
                                'MWO' ,'MWG' ,'DENO','DENG','VISO','VISG'],
                         'GRD':['PRES','PSAT','DENS','ZC1' ,'ZC7+',
                                'MWO' ,'MWG' ,'DENO','DENG','VISO','VISG'] }

        self.CALlong = { 'CCE':['Relative Volume','Liquid Dropout' ,
                                'Gas Z-Factor'   ,'Oil Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ,
                                'Oil Mole Weight','Gas Mole Weight',
                                                  'Gas Density'    ,
                                'Surface Tension',
                                'Oil Enthalpy'   ,'Gas Enthalpy'   ,
                                'Oil Spec Heat'  ,'Gas Spec Heat'  ],
                         'CVD':['Moles Removed'  ,'Liquid Dropout' ,
                                'Gas Z-Factor'   ,
                                'Oil Mole Weight','Gas Mole Weight',
                                'Oil Density'    ,'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ],
                         'DLE':['Oil FVF'        ,'Gas-Oil-Ratio'  ,
                                'Oil Density'    ,'Tot FVF'        ,
                                'Gas FVF'        ,'Gas Z-Factor'   ,
                                'Gas Gravity'    ,'Oil Viscosity'
                                'Oil Mole Weight','Gas Mole Weight',
                                'Gas Density'    ,'Gas Viscosity'  ],
                         'SEP':['Oil FVF'        ,'Gas-Oil-Ratio'  ,
                                'Oil Density'    ,'Gas Gravity'    ,
                                'Oil Mole Weight','Gas Mole Weight',
                                                  'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ],
                         'FLS':['Gas Z-Factor'   ,'Oil Density',
                                'Vapour Frac'    ,
                                'Oil Mole Weight','Gas Mole Weight',
                                                  'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ],
                         'SAT':['Sat Pressure'   ,
                                'Gas Z-Factor'   ,'Oil Density'    ,
                                'Oil Mole Weight','Gas Mole Weight',
                                                  'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ],
                         'SWL':['Sat Pressure'   ,'Swelling Factor',
                                'Oil Mole Weight','Gas Mole Weight',
                                'Oil Density'    ,'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ],
                         'GRD':['Pressure'       ,'Sat Pressure'   ,
                                'Density'        ,'C1 Mole Frac'   ,
                                'C7+ Mole Frac'  ,
                                'Oil Mole Weight','Gas Mole Weight',
                                'Oil Density'    ,'Gas Density'    ,
                                'Oil Viscosity'  ,'Gas Viscosity'  ] }

#-- Units Type per Mnemonic -----------------------------------------        

        self.EXPunit = { 'BG'  :'gfvf','BO'  :'ofvf','BT'  :'ofvf',
                         'CPG' :'cpa' ,'CPO' :'cpa' ,'HG'  :'enth','HO'  :'enth',
                         'DENG':'dens','DENO':'dens','DENS':'dens',
                         'GGRV':'dim' ,'GOR' :'gor' ,'HEIG':'leng','IFT' :'ift' ,
                         'MOLE':'mole','MREM':'mole','MWO' :'molw','MWG' :'molw',
                         'PRES':'pres','PSAT':'pres',
                         'RELV':'dim' ,'SLIQ':'mole','TEMP':'temp',
                         'VFRC':'mole','VISG':'visc','VISO':'visc','VSWL':'dim' ,
                         'ZC1 ':'mole','ZC7+':'mole','ZFAC':'dim'               }

#========================================================================
#  End of Module
#========================================================================
