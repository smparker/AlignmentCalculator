/**
 * \file "molecules.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_MOLECULES
#define ALIGNMENTCALCULATOR_MOLECULES

#include "inputs.hpp"
#include "matrix.hpp"

/// Container class for diagonal polarizability components
struct polarizability
{
  double aXX_;
  double aYY_;
  double aZZ_;
};

/// Container for rotational constants
struct rotationalConstants
{
  double Ae_;
  double Be_;
  double Ce_;
};

/// Storage for the basis function |JKM> quantum numbers
struct basis
{
  int J; ///< Angular momentum quantum number
  int K; ///< Projection of angular momentum onto body-fixed z axis
  int M; ///< Projection of angular momentum onto space-fixed z axis
  basis(int, int, int);
};

typedef std::vector<basis> basisSubset; ///< vector of basis set objects
typedef std::vector<std::shared_ptr<basisSubset>> basisSubsets; ///< vector of pointers to basis subsets
typedef std::vector<std::shared_ptr<matrixComp>> matrices; ///< vector of pointers to operator matrices
typedef std::vector<std::shared_ptr<std::vector<double>>> arrays; ///< vector of pointers to data arrays


/**
 * @brief Molecule base class
 * @details Class for storing and generating molecule data common to all rigid molecule symmetries
 *
 */
class moleculeBase
{
public:
  double even_j_degen_; ///< partition function factor to account for nuclear spin degeneracy of even J states in linear molecules
  double odd_j_degen_; ///< partition function factor to account for nuclear spin degeneracy of odd J states in linear molecules
  double partition_function_; ///< Rotational partition function
  polarizability pol_; ///< Polarizability components
  rotationalConstants rot_; ///< Rotational constants
  std::shared_ptr<matrices> Us_,invUs_; ///< Transformation matrices used if field-free Hamiltonian is not diagonal

  moleculeBase(); ///< Empty constructor
  moleculeBase(inputParameters &); ///< Constructor based on input parameters

  /**
   * @brief create basis sets
   * @details Creates the basis subsets for the molecule based on selection rules
   *
   * @param JMAX maximum value for J
   * @return pointer to basis subset array
   */
  virtual std::shared_ptr<basisSubsets> createBasisSets(int JMAX) = 0;

  /// Creates field-free rigid rotor Hamiltonians for the basis subsets provided
  virtual std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;

  /// Calculates thermal populations and partition function based on the field free Hamiltonian and basis set information
  virtual std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double) = 0;

  /// Creates a set of matrices to be used as density matrix storage or as scratch space with population data placed on the diagonals
  virtual std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>) = 0;

  /// Creates the interaction Hamiltonian prefactors (i.e. all terms except the field intensity)
  virtual std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;
};

/**
 * @brief Linear Molecules
 * @details Class derived from molecule base for linear molecules
 *
 * @param IP Input parameters
 */
class linearMolecule : public moleculeBase
{
public:
  linearMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double);
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
};

/**
 * @brief Symmetric Top Molecules
 * @details Class derived from molecule base for symmetric top molecules
 *
 * @param IP Input parameters
 */
class symmetricTopMolecule : public moleculeBase
{
public:
  MOLSYM symmetry_; ///< Used to differentiate oblate from prolate tops
  symmetricTopMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double);
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
};

/**
 * @brief Asymmetric Top Molecules
 * @details Class derived from molecule base for asymmetric top molecules
 *
 * @param IP Input parameters
 */
class asymmetricTopMolecule : public moleculeBase
{
public:
  double Xe_,Ye_,Ze_; ///< Placeholders to track coordinate system
  asymmetricTopMolecule(inputParameters &IP);
  void constructTransformationMatrices(std::shared_ptr<matrices>); ///< Create U and invU for basis set transformations
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double);
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
};

/**
 * @brief Field Matter Interaction Matrix Element
 * @details Calculates the coupling between two |JKM> states in an off resonance field
 *
 * @param J J of State 1
 * @param K K of State 1
 * @param M M of State 1
 * @param Q Interaction Quantum Number
 * @param S Other Interaction Quantum Number
 * @param j J of State 2
 * @param k K of State 2
 * @param m M of State 2
 * @return Coupling Matrix Element
 */
double FMIME (int J, int K, int M, int Q, int S, int j, int k, int m);


#endif