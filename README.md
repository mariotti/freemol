# freemol
freemol 2003

# Copyright notice/Licences: Please read this
Please note that some code lines might be a copy of other sources,
usually free source code, as this project started many years ago.

Minuit in particular: I started it from an original F77 (free from the 90s) code
and translated into F90. I contacted CERN (current minuit copyright holder) for a confirmation but I could not get a clear answer. I think it is obvious as it is a really old code which was free on the net many years ago.

There is also some free code from CCL (Computational Chemistry List), it is mentioned as comment in the code.

If code, packages or else might not fit with this "main" repository licence they are present in the folder "others/". Please read each package or code licence. If a package/code is there, then the redistribution is granted but under the given package licence which you should read.

 ## Why so complicated?
 
 There is a different solution: Get the required packages from other/external sources.
 
 But I want to create a "self-consistend" package which will run without external dependencies. For these reasons:
 
    - It might need to compile for systems without a network connection.
    - If an hardware system (for example) is not supported by the main distribution we can create patches which apply only to the given distribution.
    - Security: It is better to: download, check, compile, run, check, put in "others/", check, compile, run, check then a simple download from the net.
    - Then the obvious: I promote fortran as the most portable source code. I need all code in a unique package.
    
Warning: For example blas and lapack are definitly more optimised if you use your system ones.

# News

Added a branch to add travis CI services, to test.

# Intro

The framework is at ./Freemol

# Framework

## controlled configure and autoconfigure

## automatic makefile creation

## controlled help comment extraction

## section for manuals

# Programs

## adfrom
   A code to work with ADF (Amsterdam Density Functional)
   
## ch4sym2cart
   Convert coordinates from symmetric to Cartesian for CH4 (Methane) alike molecules

## CSMG
   Compute the CSM "Continuous Symmetry Measure" of weighted point objects using Gaussian
   functions. See this paper for details:
   
   [https://fabiomariotti.files.wordpress.com/2013/07/csm_paper.pdf](https://fabiomariotti.files.wordpress.com/2013/07/csm_paper.pdf)
   
   And eventually this short blog:

   [https://blog.techottis.ch/2016/04/01/continuous-symmetry-measure-an-old-work/](https://blog.techottis.ch/2016/04/01/continuous-symmetry-measure-an-old-work/)

   This is an example of input

    # from symmetry_mod file
    [molecule] nrec=3 format=nscxyz 
    o  1 8 1.000000 2.000000 3.000000
    h  2 1 2.000000 2.000000 3.000000
    h  3 1 1.000000 3.000000 3.000000
    [x-csmg-cgauss] 
    1 2.883843 1.591247
    2 3.838395 7.701635
    3 3.838395 7.701635
    [x-csmg-symop] 
     E    1.0000  0.0000  0.0000    0.0  0.0000  0.0000  0.0000  F 
     C2_2   0.707107 0.707107 0.000000 180.000000  0.0000  0.0000  0.0000  F 
     S1   -0.707107 0.707107 0.000000 180.0000  0.0000  0.0000  0.0000  I 
     S2   0.000000 0.000000 1.000000 180.0000  0.0000  0.0000  0.0000  I 

### Dependencies
See _Libraries_

    #
    USEDLIBS = -L$(FMLIBDIR) -lmodules -lmoduledata -lutils -lincludes \
               -lfmlapack -lfmblas
    #


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

# Libraries
  These libraries and tools are used in the programs listed above. For a working example of their
  usage please check the programs. In order to link to these libraries or tools please check also the
  generated makefiles.

## Includes

### linetools: format free input

### messages: a logging facility

### pcmdline: to handle the command line

### strtools: string manipulation routines

### vartypes: standardize variables type within codes

### chemconst: constants defined at different precision

## Utils

### mathtools: few simple and common matrix tools

### qnumbers: a small tool for quantum numbers

## Modules

### baseio: a file manager

### extio: extensions to baseio for easy file manipulation

### osec, sections: a tool to work with sectioned files like molden format

### minuit: a multidimensional minimization tool

## Moduledata

### molecules: a molden based structure and tools to store molecules data and read in "free" format

## Tools

### CSMD supporting tools

## others

### local blas and lapack
