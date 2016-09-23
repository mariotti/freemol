# Freemol package

This package is in Alpha stage.

The code itself has been tested for years, since the '90.

# Current Status

 - Testing configuration
   - ./config/autoconfig and ./config/configure, are working on FEDORA 24 with gfortran
   - local blas and lapack are compiled into library. see ./lib
   - On OS X it compiles

# Libs features

  - Read/Write files in .INI (windows ini files) or .mld (molden files) formats,
    a want to be easy to read replacement for XML formats.
  - 

# Programs

## adfrom
   A code to work with ADF (Amsterdam Density Functional)
   
## ch4sym2cart
   Convert coordinates from symmetric to cartesian for CH4 (Methane) like molecules

## CSMG
   Compute the CSM "Continuous Symmetry Measure" of weighted point objects using gaussians
   functions as weight and pure spacial overlap/superposition.

## fit1Dpol
   Fit about anything in 1 dimension. Anything might mean you need to add your functions,
   but it has a nice polynomial setup.
   It has nice features to work with multi-columns ASCII files.
   You can script which data to fit with ranged commands like
   copy, add, sub, shift, scale etc..

## Freemol
   Just a placeholder at present

## XY4Coord
   Coordinates transformation

## XY4PolySphere
   Coordinates transformation, the polyspherical version
