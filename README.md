# QuICPack
A Quantum Information and Computation Library

The package contains routines developed in FORTRAN77/ FORTRAN90 for computations involved in Quantum Information and Computation, occasional applications in Quantum Many-Body Physics, and related areas. Please acknowledge if you happen to find something useful.

quicpack is *free*, and can be redistributed/modified according to requirement. Please note that using quicpack requires basic knowledge of FORTRAN77/ 90, and it uses

lapack (linear algebra pack available at http://www.netlib.org/lapack/)
blas (basic linear algebra subprograms available at http://www.netlib.org/blas/)
nlopt (non-linear optimization routine available at http://ab-initio.mit.edu/wiki/index.php/NLopt)
arpack (https://people.sc.fsu.edu/~jburkardt/f_src/arpack/arpack.html)
subroutines and functions available in  Numerical Recipes in FORTRAN77/FORTRAN90
Disclaimer: quicpack is custom-built, does not work in all cases, and is always under severe update.

Using the routines:

(1) Unbundle all the files to the main program folder. 
(2) Compile your code main.f90 (the one which uses quicpack routines) using 
            $ gfortran -lblas -llapack -I/usr/local/include -lm -lnlopt *.f90 main.f90 -o out.out
      Note here that *.f90 includes all the FORTRAN scripts required in main.f90, along with the quicpack routines.
(3) To run the executable file out.out, use
            $ ./out.out

Alternatively, create a shell script.
