########################################################################
#
#       Makefile for FRIMOL external packages
#
########################################################################
# Frimol by F.Mariotti: (c) F.Mariotti
########################################################################
# $Id: Makefile,v 1.1.1.1 2009/01/12 16:56:08 mariotti Exp $
########################################################################
#
#
include ../Makeflags
#
.PHONY: all printme blas lapack clean
#
#       TARGETS
#
#    today : 
#	echo "Today do nothing"

      all : printme blas lapack
#
   printme:
	echo " "
	echo "    Make is running on Frimol other packages"
	echo "    Version: " `../config/printversion`
	echo " "
#
      blas:
	$(MAKE) -C blas

#
    lapack:
	$(MAKE) -C lapack

#
    clean:
	$(MAKE) -C blas clean
	$(MAKE) -C lapack clean

#
