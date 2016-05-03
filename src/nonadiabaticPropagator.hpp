/**
 * \file "adiabaticPropagator.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_NONADIABATICPROP
#define ALIGNMENTCALCULATOR_NONADIABATICPROP

#include "molecules.hpp"
#include "matrix.hpp"
#include "inputs.hpp"
#include "outputs.hpp"
#include <memory>

class nonadiabaticPropagator
{
public:
  double partitionFxn_;
  double t0_;
  double tFinal_;
  double noutputs_;
  double temperature_;
  std::shared_ptr<moleculeBase> molecule_;
  std::shared_ptr<basisSubset> basisSets_;
  std::shared_ptr<matrices> fieldFreeHamiltonians_;
  std::shared_ptr<matrices> intHamiltonians_;
  std::shared_ptr<matrices> densities_;
  std::shared_ptr<arrays> populations_;
  std::vector<observable> observables_;
  std::vector<pulse> pulses_;

  nonadiabaticPropagator();
  nonadiabaticPropagator(inputParameters &IP);
  void determineSymmetry();
  void setupOutputs();
  void initializeCVODE();
  void evalRHS();
  void step();
  void run();
};

#endif