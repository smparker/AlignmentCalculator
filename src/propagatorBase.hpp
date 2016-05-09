/**
 * \file "propagatorBase.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_PROPAGATORBASE
#define ALIGNMENTCALCULATOR_PROPAGATORBASE

#include "outputs.hpp"

class propagatorBase
{
public:
  double partition_function_, temperature_;
  MOLSYM symmetry_;
  std::shared_ptr<moleculeBase> molecule_;
  std::shared_ptr<basisSubsets> basisSets_;
  std::shared_ptr<matrices> fieldFreeHamiltonians_;
  std::shared_ptr<matrices> intHamiltonians_;
  std::shared_ptr<matrices> densities_;
  std::shared_ptr<arrays> populations_;
  std::vector<std::shared_ptr<observable>> observables_;

  propagatorBase();
  propagatorBase(inputParameters &IP);
  MOLSYM determineSymmetry(inputParameters &IP);
  virtual void initializeOutputs(inputParameters &IP) = 0;
  virtual void printOutputs() = 0;
};

#endif