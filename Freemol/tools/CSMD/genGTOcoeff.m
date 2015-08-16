function [gnorm, zeta] = genGTOcoeff(nc,ac,ar)
#
# Produce the exponent and normalization
# for a given input of nuclear charge
# atomic charge and atomic radii
#
# This File is part of the Freemol package
# (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
# Please read the full copyright statment in the COPYING file
#
# nc = input("Insert the Nuclei Charge: ");
# ac = input("Insert the Atomic Charge: ");
# ar = input(" Insert the Atomic Radii: ");
#
# From the atomic radii we generate the zeta
#
zeta = - log(0.5)/(ar*ar);
#
# From the integration to fit the number of electrons
# we get the gto coefficient. numelec = nc - ac
#
ne = nc - ac;
global gto_zeta;
gto_zeta = zeta;
#[val ier nfun err] = quad("gtorz",0,Inf);
val = quad("gtorz",0,Inf);
gnorm = ne/val;
#
endfunction

