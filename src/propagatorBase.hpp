/**
 * \file "propagatorBase.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_PROPAGATORBASE
#define ALIGNMENTCALCULATOR_PROPAGATORBASE

#include "outputs.hpp"

/**
 * @brief Base class for adiabatic and nonadiabatic calculations
 * @details Class for managing setup of calculations and outputting data containing functionality universal to all jobtypes
 *
 */
class propagatorBase
{
public:
  double partition_function_; ///< Partition function for the initial state
  double temperature_; ///< Initial temperature
  MOLSYM symmetry_; ///< Symmetry of the molecule determining the coordinate system
  std::shared_ptr<moleculeBase> molecule_; ///< molecule object
  std::shared_ptr<basisSubsets> basisSets_; ///< |JKM> states arranged in symmetry coupled subsets
  std::shared_ptr<matrices> fieldFreeHamiltonians_; ///< Rigid rotor Hamiltonian matrices
  std::shared_ptr<matrices> intHamiltonians_; ///< Interaction Hamiltonian prefactors (i.e. not including field strength)
  std::shared_ptr<matrices> densities_; ///< Density matrix storage space
  std::shared_ptr<arrays> populations_; ///< Boltzmann population for corresponding to the basis states
  std::vector<std::shared_ptr<observable>> observables_; ///< Observables to be calculated during propagation

  propagatorBase(); ///< Default constructor, no real functionality
  propagatorBase(inputParameters &IP); ///< Recommended constructor

  /**
   * @brief Outputs basis information to file
   * @details Prints file at end of propagation including a list of all basis states, energies, and thermal populations
   */
  void outputBasisStats();
  MOLSYM determineSymmetry(inputParameters &IP);
  void initialize_();
  virtual void initializeOutputs(inputParameters &IP) = 0;
  virtual void printOutputs() = 0;
  virtual void transformObservables();
  void removeSmallPopulations();
};

#endif