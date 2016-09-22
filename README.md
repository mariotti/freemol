# freemol
freemol 2003

# Copyright notice
Please note that some code lines might be a copy of other sources. Free source code, including minuit. As I just started I might have missed the reference (if it exists anymore). I will do an effort to update the relative documentation.

Minuit in particular, I started from an original F77 (free) code and translated into F90. I contacted CERN for a confirmation.

# Intro

I want to revive an old project. At present I start with some "stubs" routines.

The target are molecules, the name actually changed from frimol to freemol. Was FRIbourg MOLecules now it is FREE (beer) MOLecules.

# Stubs description

# Copyright and Free code list

This is one of the original copyright (copyleft?) statment for the freeread routines.

    module linetools
    !
    !H
    !H--------------------------------------------------------------
    !H      Some Routine are adapted from FREEREAD package taken
    !H      from CCL.
    !H      CCL: Computational Chemistry List http://www.ccl.org/
    !H      A statment from Jan Labanowski:
    !H      This software was taken from anon. ftp on
    !H      rani.chem.yale.edu. It is a set of routines which allow 
    !H      format free input in FORTRAN
    !H
    !H      Jan Labanowski
    !H      jkl@osc.edu
    !H--------------------------------------------------------------
    !

