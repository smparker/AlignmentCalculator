/**
 * \file "molecules.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_MOLECULES
#define ALIGNMENTCALCULATOR_MOLECULES

#include "inputs.hpp"
#include "matrix.hpp"

struct polarizability
{
  double aXX_;
  double aYY_;
  double aZZ_;
};

struct rotationalConstants
{
  double Ae_;
  double Be_;
  double Ce_;
};

struct basis
{
  int J;
  int K;
  int M;
  basis(int, int, int);
};

typedef std::vector<basis> basisSubset;
typedef std::vector<std::shared_ptr<basisSubset>> basisSubsets;
typedef std::vector<std::shared_ptr<matrixComp>> matrices;
typedef std::vector<std::vector<double>*> arrays;

class moleculeBase
{
public:
  polarizability pol_;
  rotationalConstants rot_;

  moleculeBase();
  moleculeBase(inputParameters &);
  virtual std::shared_ptr<basisSubsets> createBasisSets(int JMAX) = 0;
  virtual std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;
  virtual std::shared_ptr<arrays> initializePopulations() = 0;
  virtual std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>) = 0;
  virtual std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;

  virtual double calculatePartitionFxn() = 0;
};


class linearMolecule : public moleculeBase
{
public:
  linearMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations();
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
  double calculatePartitionFxn();
};

class symmetricTopMolecule : public moleculeBase
{
public:
  symmetricTopMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations();
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
  double calculatePartitionFxn();
};

class asymmetricTopMolecule : public moleculeBase
{
public:
  std::shared_ptr<matrices> U,invU; ///< Transformation matrices
  asymmetricTopMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations();
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
  double calculatePartitionFxn();
};


#endif