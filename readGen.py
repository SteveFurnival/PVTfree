
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

    fInP = clsIO.fInP

#----------------------------------------------------------------------
#  Loop until ENDDEB read
#----------------------------------------------------------------------

    for curL in fInP :

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

        elif tokS[0][:2].upper() == "PS"   :

            nVal = int(tokS[1])
            clsIO.Deb["PSAT"] = nVal

        elif tokS[0][:3].upper() == "FLA"   :

            nVal = int(tokS[1])
            clsIO.Deb["FLASH"] = nVal

        elif tokS[0][:5].upper() == "STAB2"   :

            nVal = int(tokS[1])
            clsIO.Deb["STAB2"] = nVal

        elif tokS[0][:5].upper() == "STAB1"   :

            nVal = int(tokS[1])
            clsIO.Deb["STAB1"] = nVal

        elif tokS[0][:2].upper() == "RR"   :

            nVal = int(tokS[1])
            clsIO.Deb["RR"] = nVal

        elif tokS[0][:3].upper() == "REG"   :

            nVal = int(tokS[1])
            clsIO.Deb["REG"] = nVal

        elif tokS[0][:2].upper() == "BL"   :

            nVal = int(tokS[1])
            clsIO.Deb["BLACK"] = nVal

        elif tokS[0][:3].upper() == "GRO"   :

            nVal = int(tokS[1])
            clsIO.Deb["GROUP"] = nVal

        elif tokS[0][:3].upper() == "ENV"   :

            nVal = int(tokS[1])
            clsIO.Deb["ENV"]   = nVal

        elif tokS[0][:3].upper() == "CCE"   :

            nVal = int(tokS[1])
            clsIO.Deb["CCE"] = nVal

        elif tokS[0][:3].upper() == "CVD"   :

            nVal = int(tokS[1])
            clsIO.Deb["CVD"] = nVal

        elif tokS[0][:3].upper() == "DLE"   :

            nVal = int(tokS[1])
            clsIO.Deb["DLE"] = nVal

        elif tokS[0][:3].upper() == "SEP"   :

            nVal = int(tokS[1])
            clsIO.Deb["SEP"] = nVal

        elif tokS[0][:3].upper() == "SAT"   :

            nVal = int(tokS[1])
            clsIO.Deb["SAT"] = nVal

        elif tokS[0][:3].upper() == "FLS"   :

            nVal = int(tokS[1])
            clsIO.Deb["FLS"] = nVal

        elif tokS[0][:3].upper() == "SWL"   :

            nVal = int(tokS[1])
            clsIO.Deb["SWL"] = nVal

        elif tokS[0][:3].upper() == "GRD"   :

            nVal = int(tokS[1])
            clsIO.Deb["GRD"] = nVal

        elif tokS[0][:2].upper() == "TS"   :

            nVal = int(tokS[1])
            clsIO.Deb["TSAT"] = nVal

#== Return values =====================================================

    #print("clsIO.Deb ",clsIO.Deb)

    return iERR

#========================================================================
#  End of Module
#========================================================================
