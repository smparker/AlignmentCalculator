# Alignment Calculator v 1.0.0 #

Alignment Calculator simulates the adiabatic and nonadiabatic alignment of rigid molecules subjected to intense linearly polarized light.

### Required software ###

* [Intel MKL](https://software.intel.com/en-us/intel-mkl) or other BLAS/LAPACK installation
* [BOOST](http://www.boost.org)
* [Sundials / CVode](http://computation.llnl.gov/projects/sundials-suite-nonlinear-differential-algebraic-equation-solvers/sundials-software) (will be optional in future update)
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)

### Building ###

To build, simply run the configure script and link any libraries not in your default path, then make.

* ./configure CPPFLAGS=-I/opt/local/include -I/usr/local/include LDFLAGS=-L/opt/local/lib -lsundials_cvode -lsundials_nvecserial
* make

If problems arise, the configure file may be rebuilt using automake/autoconf. Use the following: 

* aclocal
* autoheader
* automake -a
* autoconf

Then follow the first set of instructions. 

### Running the Code  ###

* Run the executable followed by the name of the input file (i.e. "./alignmentcalculator inputs.json")
* Example input files are provided for convenience 


### Questions or Comments? ###

Message me and I'll respond as soon as I'm able