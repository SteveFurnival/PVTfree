
#=======================================================================
#
#  Module in PVTfree, Copyright (C) Steve Furnival, HoBoil Limited
#  see PVTfree.py for further details
#
#  Read General Data, i.e. Debug, Options, ...
#
#=======================================================================

#!/usr/bin/python3

def readDebug(clsIO) :

#-- Open the Debug File, if not already open ------------------------

    if not clsIO.qDeb :
        pathDeb = clsIO.patR + ".deb"
        fDeb    = open(pathDeb,'w')
        clsIO.setQDEB(True)
        clsIO.setFDEB(fDeb)
    else :
        fDeb = clsIO.fDeb

#-- Initialisation --------------------------------------------------    

    print("Reading Debug Instructions")

    iERR  = 0
    iLine = 0

#== Input File ========================================================

    fInp = clsIO.fInp

#----------------------------------------------------------------------
#  Loop until ENDDEB read
#----------------------------------------------------------------------

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of DEBUG Data - Error")
            iERR = -1
            break

#-- Split current line into Tokens for processing -------------------

        tokS = curL.split()
        nTok = len(tokS)

#== Start processing tokens ===========================================

        if   nTok                == 0        : pass    #-- Blank line!
        elif tokS[0][:2]         == "--"     : pass    #-- Comment
        elif tokS[0][:6].upper() == "ENDDEB" : break   #-- ENDDEBUG => Exit

        elif tokS[0][:2].upper() == "PS"    : clsIO.Deb["PSAT"]  = int(tokS[1])
        elif tokS[0][:3].upper() == "FLA"   : clsIO.Deb["FLASH"] = int(tokS[1])
        elif tokS[0][:5].upper() == "STAB2" : clsIO.Deb["STAB2"] = int(tokS[1])
        elif tokS[0][:5].upper() == "STAB1" : clsIO.Deb["STAB1"] = int(tokS[1])
        elif tokS[0][:2].upper() == "RR"    : clsIO.Deb["RR"]    = int(tokS[1])
        elif tokS[0][:3].upper() == "REG"   : clsIO.Deb["REG"]   = int(tokS[1])
        elif tokS[0][:2].upper() == "BL"    : clsIO.Deb["BLACK"] = int(tokS[1])
        elif tokS[0][:3].upper() == "GRO"   : clsIO.Deb["GROUP"] = int(tokS[1])
        elif tokS[0][:3].upper() == "ENV"   : clsIO.Deb["ENV"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "CCE"   : clsIO.Deb["CCE"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "CVD"   : clsIO.Deb["CVD"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "DLE"   : clsIO.Deb["DLE"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "SEP"   : clsIO.Deb["SEP"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "SAT"   : clsIO.Deb["SAT"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "FLS"   : clsIO.Deb["FLS"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "SWL"   : clsIO.Deb["SWL"]   = int(tokS[1])
        elif tokS[0][:3].upper() == "GRD"   : clsIO.Deb["GRD"]   = int(tokS[1])
        elif tokS[0][:2].upper() == "TS"    : clsIO.Deb["TSAT"]  = int(tokS[1])

#== Return values =====================================================

    #print("clsIO.Deb ",clsIO.Deb)

    return iERR

#========================================================================
#  Read Options
#========================================================================

def readOption(clsIO) :

#-- Open the Debug File, if not already open ------------------------

    if not clsIO.qDeb :
        pathDeb = clsIO.patR + ".deb"
        fDeb    = open(pathDeb,'w')
        clsIO.setQDEB(True)
        clsIO.setFDEB(fDeb)
    else :
        fDeb = clsIO.fDeb

#-- Initialisation --------------------------------------------------    

    print("Reading Option Instructions")

    iERR  = 0
    iLine = 0

#== Input File ========================================================

    fInp = clsIO.fInp

#----------------------------------------------------------------------
#  Loop until ENDOPT read
#----------------------------------------------------------------------

    for curL in fInp :

        iLine += 1

        if iLine > 100 :
            print("Too Many Lines of OPTION Data - Error")
            iERR = -1
            break

#-- Split current line into Tokens for processing -------------------

        tokS = curL.split()
        nTok = len(tokS)

#== Start processing tokens ===========================================

        if   nTok                == 0        : pass    #-- Blank line!
        elif tokS[0][:2]         == "--"     : pass    #-- Comment
        elif tokS[0][:6].upper() == "ENDOPT" : break   #-- ENDOPTION => Exit

        elif tokS[0][:4].upper() == "STAB"   : clsIO.Opt["STAB"] = int(tokS[1])

#== Return values =====================================================

    return iERR

#========================================================================
#  End of Module
#========================================================================
