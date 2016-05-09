/**
 * \file "outputs.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_OUTPUTS
#define ALIGNMENTCALCULATOR_OUTPUTS

#include "molecules.hpp"
#include <fstream>

/**
 * @brief      Base class for observables that are obtained from a density matrix as Tr(O*rho). Other outputs are obtained from other members of the propagator class, such as the wavefunction density or the list of basis states used in the calculation.
 */
class observable
{
public:
  std::string id_tag_;
  std::shared_ptr<matrices> operator_matrix_;
  observable(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  virtual void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) = 0;
  virtual double evaluate_(std::shared_ptr<matrices> densities_);
};


class obsCosTheta3D : public observable
{
public:
  obsCosTheta3D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

// class obsCosTheta2D : public observable
// {

// };

// class obsEnergy : public observable
// {

// };

class obsCosChi : public observable
{
public:
  obsCosChi(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

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