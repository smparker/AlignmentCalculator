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


class nonadiabaticPropagator : public propagatorBase
{
public:
  bool firstRun_;
  int index_flag_;
  double t0_, dt_, time_;
  double tFinal_;
  double noutputs_;
  double atol_,rtol_;
  std::vector<pulse> pulses_;
  std::vector<N_Vector> atols_,ys_;
  std::vector<void*> cvode_managers_;
  nonadiabaticPropagator(inputParameters &IP);
  void initializeCVODE();
  static int evalRHS(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  void step();
  void run();
};


int check_flag(void *flagvalue, char *funcname, int opt);

void PrintFinalStats(void *cvode_mem);

#endif