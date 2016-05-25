/**
 * \file "adiabaticPropagator.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_NONADIABATICPROP
#define ALIGNMENTCALCULATOR_NONADIABATICPROP

#include "propagatorBase.hpp"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <cvode/cvode_band.h>       /* prototype for CVBand */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

/**
 * @brief Nonadiabatic Calculation Propagator
 * @details Class for managing data and outputs for aligning molecules in a nonadiabatic field (laser pulse)
 *
 */
class nonadiabaticPropagator : public propagatorBase
{
public:
  bool firstRun_; ///< Flag for propagator to identify if the calculation has run previously
  int noutputs_; ///< Number of output times
  int index_flag_; ///< Flag for passing information to the cvode propagator
  double t0_; ///< Initial time
  double dt_; ///< Time Step
  double time_; ///< Current Time
  double tFinal_; ///< Final Time
  double atol_,rtol_; ///< absolute and relative error tolerances
  std::string output_file_name_; ///< Name of the output file
  std::ofstream out_file_; ///< Output file stream
  std::vector<pulse> pulses_; ///< vector containing all pulse objects (for calculating pulse trains)
  std::vector<N_Vector> atols_,ys_; ///< CVode tolerance and RHS vector storage
  std::shared_ptr<matrices> scratch_matrices_,scratch_ydot_; ///< Useful scratch space to avoid reallocation of memory
  std::vector<void*> cvode_managers_; ///< CVode objects
  nonadiabaticPropagator(inputParameters &IP); ///< Constructor
  void initializeCVODE(inputParameters &IP); ///< Initialize all differential equation solvers
  void initializeOutputs(inputParameters &IP); ///< Initialize output streams
  static int evalRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data); ///< Evaluate the right hand side function (Liouville von Neumann equation)
  void step(); ///< Step density matrices by time step
  void run(); ///< Run full simulation
  void printOutputs(); ///< Print all observables to output file
};

/// Function for checking the proper return of CVode functions
int check_flag(void *flagvalue, char *funcname, int opt);

/// row helper function for one-to-one correspondence between matrix coordinate and upper triangular storage scheme
inline int row_index(int i, int N )
{
    double row = (-2*N - 1 + sqrt( (4*N*(N+1) - 8*(double)i - 7) )) / -2;
    if( row == (double)(int) row ) row -= 1;
    return row;
}

/// column helper function for one-to-one correspondence between matrix coordinate and upper triangular storage scheme
inline int column_index( int i, int N )
{
    int row = row_index(i, N);
    return  i - N * row + row*(row+1) / 2;
}


#endif