# Microsoft Developer Studio Project File - Name="LSD_updLSC" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=LSD_updLSC - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "LSD_updLSC.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "LSD_updLSC.mak" CFG="LSD_updLSC - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "LSD_updLSC - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "LSD_updLSC - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "LSD_updLSC - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x419 /d "NDEBUG"
# ADD RSC /l 0x419 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "LSD_updLSC - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /check:bounds /compile_only /dbglibs /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x419 /d "_DEBUG"
# ADD RSC /l 0x419 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "LSD_updLSC - Win32 Release"
# Name "LSD_updLSC - Win32 Debug"
# Begin Source File

SOURCE=.\approx.f90
NODEP_F90_APPRO=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\classes.f90
# End Source File
# Begin Source File

SOURCE=.\correction.f90
NODEP_F90_CORRE=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\golden_section.f90
DEP_F90_GOLDE=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\group_lines.f90
NODEP_F90_GROUP=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\initial_guess_gauss.f90
NODEP_F90_INITI=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\input.f90
DEP_F90_INPUT=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\line_strength_correction.f90
DEP_F90_LINE_=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\LSD_updLSC.f90
DEP_F90_LSD_U=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\map1.f90
# End Source File
# Begin Source File

SOURCE=.\mini.f90
# End Source File
# Begin Source File

SOURCE=.\mini_reg.f90
# End Source File
# Begin Source File

SOURCE=.\model_spectrum.f90
NODEP_F90_MODEL=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\modules.f90
DEP_F90_MODUL=\
	".\Debug\gaussians.mod"\
	".\Debug\regions.mod"\
	".\Debug\spectra.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\obs_LSD.f90
NODEP_F90_OBS_L=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=.\output.f90
NODEP_F90_OUTPU=\
	".\Debug\common_vars.mod"\
	
# End Source File
# Begin Source File

SOURCE=..\..\Compac_66\X86\DF\CXML\LIB\CXML.LIB
# End Source File
# End Target
# End Project
