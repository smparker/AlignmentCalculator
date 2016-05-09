/**
 * \file "outputs.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_OUTPUTS
#define ALIGNMENTCALCULATOR_OUTPUTS

#include "molecules.hpp"
#include <fstream>

class observable
{
public:
  std::string id_tag_;
  observable(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  virtual void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) = 0;
  virtual double evaluate_(std::shared_ptr<matrices> densities_) = 0;
};


class obsCosTheta3D : public observable
{
public:
  obsCosTheta3D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  double evaluate_(std::shared_ptr<matrices> densities_);
};

// class obsCosTheta2D : public observable
// {

// };

// class obsEnergy : public observable
// {

// };

// class obsCosChi : public observable
// {

// };

// class obsJ : public observable
// {

// };

// class obsK : public observable
// {

// };

// class obsM : public observable
// {

// };

#endif