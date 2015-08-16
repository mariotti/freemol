function y=gtorz(x)
# This File is part of the Freemol package
# (C) 2003 Fabio Mariotti, <fabio.mariotti@scriptsforscience.org>
# Please read the full copyright statment in the COPYING file
  global gto_zeta
     y=exp(- gto_zeta .* x .* x ) .* x .* x .* 4.0 .* pi;
endfunction
