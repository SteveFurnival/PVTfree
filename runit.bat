@echo off
rem
rem =======================================================
rem
rem  Command Prompt (Windows) Script to run PVTfree
rem
rem  The path to the location of the Python code may
rem  have to be changed by the user
rem
rem  The script will look for a single argument (%1)
rem  which is the rootname (no .dat) of the datafile
rem  that the user wishes to run
rem
rem  Copyright (C) Steve Furnival, HoBoil Limited, 2019.
rem  see PVTfree.py for further details
rem
rem =======================================================
rem
cd GraphOut
del *.png
cd ..
py D:\steve\PVTfree\PVTfree.py %1
