# This File is part of the Freemol package
# (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
# Please read the full copyright statment in the COPYING file
########################################################################
# Main Makefile for Freemol by F.Mariotti: (c) F.Mariotti
########################################################################
#
#	Makefile for FREEMOL
#
#	Fabio Mariotti
#	University of Fribourg
#
########################################################################
# $Id: Makefile,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
########################################################################
#
#
include Makeflags
#
#	TARGETS
#
# Default goes to HELP
#
   all : help


#
# Freemol Program
#
   freemol : printme binaries helpdocs



.PHONY: all others includes utilities moduledata modules programs
.PHONY: cleanall cleanwork cleanobj cleanlib cleanmid
.PHONY: help freemol
#
   printme:
	echo " "
	echo "    Make is running on Freemol"
	echo "    Version: " `./config/printversion`
	echo " "

   binaries: libraries 

   libraries : others includes utilities moduledata modules programs

#

   others:
	$(MAKE) -C others

   includes: others
	$(MAKE) -C $(FMINCLUDES)

   utilities: includes
	$(MAKE) -C $(FMUTILS)

   moduledata: includes utilities
	$(MAKE) -C $(FMMODULEDATA)

   modules: includes moduledata utilities
	$(MAKE) -C $(FMMODULES)

   programs: includes moduledata modules utilities
	$(MAKE) -C $(FMPROGRAMS)

   helpdocs: 
	$(MAKE) -C $(FMHELPDOCS)

   cleanall: cleanprog cleanwork cleanobjs cleanlibs cleanmid cleanmod cleandoc

#
#   cleanprog:
#	rm -f $(FMPROG)
#

   cleanprog:
	$(MAKE) -C $(FMPROGRAMS) clean
	rm -f bin/*.exe

   cleanwork:
	$(FMCLP) $(FMHOME)

   cleanobjs:
	find $(FMHOME) \( -name '*.lst' -o -name '*.o' -o -name 'test' \) -print -exec rm -f \{\} \;
	find $(FMHOME) -name '*timestamp*' -print -exec rm -f \{\} \;

   cleanlibs:
	find $(FMHOME) -name '*.a' -print -exec rm -f \{\} \;

   cleanmid:
	echo "CleanMid removed"
#	find $(FMHOME) -name 'fm*.F90' -print -exec rm -f \{\} \;

   cleanmod:
	find $(FMHOME) -name '*.mod' -print -exec rm -f \{\} \;
	find $(FMHOME) -name '*.d' -print -exec rm -f \{\} \;
#	rm -f $(FMHOME)/MODFILES/*

   cleandoc:
	find $(FMHOME) -name 'fundoc.txt' -print -exec rm -f \{\} \;
	rm -f ./docs/*.txt
	find $(FMHOME)/help -name '*.dvi' -print -exec rm -f \{\} \;
	find $(FMHOME)/help -name '*.log' -print -exec rm -f \{\} \;
	find $(FMHOME)/help -name '*.aux' -print -exec rm -f \{\} \;
	find $(FMHOME)/help -name '*.toc' -print -exec rm -f \{\} \;
	find $(FMHOME)/help -name '*.lot' -print -exec rm -f \{\} \;

   cleanconfig:
	rm -f $(FMHOME)/Makefile $(FMHOME)/Makeflags $(FMHOME)/Make.machine
	find $(FMHOME) -name 'Make.machine' -print -exec rm -f \{\} \;
	find $(FMHOME) -name 'Makefile.here' -print -exec rm -f \{\} \;
	rm -f $(FMHOME)/utils/Makefile
	rm -f $(FMHOME)/includes/Makefile
	rm -f $(FMHOME)/modules/Makefile
	rm -f $(FMHOME)/moduledata/Makefile
	rm -f $(FMHOME)/programs/Makefile
	rm -f $(FMHOME)/programs/adfrom/Makefile
	rm -f $(FMHOME)/programs/adfpsi/Makefile
	rm -f $(FMHOME)/programs/pdrivers/Makefile
	rm -f $(FMHOME)/programs/CSMG/Makefile
	rm -f $(FMHOME)/programs/Freemol/Makefile
	rm -f $(FMHOME)/programs/fit1Dpol/Makefile
	rm -rf $(FMHOME)/FFILES
	rm -rf $(FMHOME)/MODFILES
#	rm -rf $(FMHOME)/bin
#
   help:
	echo " "
	echo " Makefile for Freemol by F.Mariotti"
	echo " "
	echo " Available targets:"
	echo " "
	echo " binaries             Builds binaries and libraries"
	echo " libraries             \" libraries"
	echo " includes              \" Only includes library"
	echo " utilities             \" Only utility library"
	echo " moduledata            \" Only Data modules library"
	echo " modules               \" Only modules library"
	echo " programs              \" Only programs"
	echo " helpdocs              \" help documents tree"
	echo " all                  Print this help message"
	echo " freemol               Make the program!"
	echo " "
	echo " Clean targets:"
	echo " "
	echo " cleanprog            Cleanup the Freemol executable"
	echo " cleanwork             \" backup files tmp files etc..."
	echo " cleanobjs             \" all the object files"
	echo " cleanlibs             \" all library files"
	echo " cleanmid              \" mid files fm*.F90"
	echo " cleanmod              \" all the *.mod files"
	echo " cleandoc              \" all created documents"
	echo " cleanall              \" all targets above"
	echo " cleanconfig           \" all files created by configuration"
	echo " "
	echo " help                 This target"
	echo " "
	echo " "
	echo " "

#
  tarit:
	tar -C .. -czvf ../`basename \`pwd\``.`date +%d%b%Y`.tar.gz `basename \`pwd\``

#

