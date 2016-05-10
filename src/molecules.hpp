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
typedef std::vector<std::shared_ptr<std::vector<double>>> arrays;

class moleculeBase
{
public:
  double even_j_degen_, odd_j_degen_;
  double partition_function_;
  polarizability pol_;
  rotationalConstants rot_;

  moleculeBase();
  moleculeBase(inputParameters &);
  virtual std::shared_ptr<basisSubsets> createBasisSets(int JMAX) = 0;
  virtual std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;
  virtual std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double) = 0;
  virtual std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>) = 0;
  virtual std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets) = 0;

};


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

class symmetricTopMolecule : public moleculeBase
{
public:
  MOLSYM symmetry_;
  symmetricTopMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double);
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
};

class asymmetricTopMolecule : public moleculeBase
{
public:
  std::shared_ptr<matrices> U,invU; ///< Transformation matrices
  asymmetricTopMolecule(inputParameters &IP);
  std::shared_ptr<basisSubsets> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets);
  std::shared_ptr<arrays> initializePopulations(std::shared_ptr<basisSubsets>,std::shared_ptr<matrices>,double);
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets);
};

double FMIME (int J, int K, int M, int Q, int S, int j, int k, int m);


#endif