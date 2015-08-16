#! /bin/sh
#
##############################################################################
# This File is part of the Freemol package
# Copyright (C) 2006 Fabio Mariotti, <fabio.mariotti@scriptsforscience.com>
# Please read the full copyright statment in the COPYING file
# This file goes under the GPL section of the COPYING document
##############################################################################
#
# We create a distribution directory with source code which includes all the
# files with a known license and builts libraries and/or programs.
#
# At present we should get a distribution directory with all the tools to
# work with the freemol environment.
#
##############################################################################
#
# This script should run as:
# builddist.sh <DIR>
#
# This script should create a dist directory and copy all needed files in the
# root/dist directory and run builddist.sh in any subdirectory. Because it
# does reside in the config/ directory it should build the config/ directory
# structure itself.
#
##############################################################################
# BASIC CONFIGURATION DATA
##############################################################################
# FILES
fm_files=''
# DIRS
fm_dirs=''
#
##############################################################################
# check command line dir
##############################################################################
#
# We need one parameter at least
if [ "a$1" = "a" ];
then
  echo "Please give a name on the command line."
  exit;
fi;
#
# save the build dir name and shift pars
bdir=$1
shift
#
# do checks on the dir
#
if [ -d $bdir ];
then
  echo "The build Directory exists. At present we do not accept upgrtades. Remove it and run me again."
  exit;
else
  echo mkdir $bdir
  echo "Dist Dir created."
fi;
#

