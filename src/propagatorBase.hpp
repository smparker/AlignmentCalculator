/**
 * \file "propagatorBase.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_PROPAGATORBASE
#define ALIGNMENTCALCULATOR_PROPAGATORBASE

#include "molecules.hpp"
#include "matrix.hpp"
#include "inputs.hpp"
#include "outputs.hpp"
 #include <memory>

class propagatorBase
{
public:
  double partition_function_;
  double temperature_;
  MOLSYM symmetry_;
  std::shared_ptr<moleculeBase> molecule_;
  std::shared_ptr<basisSubsets> basisSets_;
  std::shared_ptr<matrices> fieldFreeHamiltonians_;
  std::shared_ptr<matrices> intHamiltonians_;
  std::shared_ptr<matrices> densities_;
  std::shared_ptr<arrays> populations_;
  std::vector<observable> observables_;

  propagatorBase();
  propagatorBase(inputParameters &IP);
  MOLSYM determineSymmetry(inputParameters &IP);
  void setupOutputs();
};

#endif