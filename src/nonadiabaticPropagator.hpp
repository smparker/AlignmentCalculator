/**
 * \file "adiabaticPropagator.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_NONADIABATICPROP
#define ALIGNMENTCALCULATOR_NONADIABATICPROP

#include "propagatorBase.hpp"

class nonadiabaticPropagator : public propagatorBase
{
public:
  double t0_;
  double tFinal_;
  double noutputs_;
  std::vector<pulse> pulses_;

  nonadiabaticPropagator(inputParameters &IP);
  void initializeCVODE();
  void evalRHS();
  void step();
  void run();
};

#endif