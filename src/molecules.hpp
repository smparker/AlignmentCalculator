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
  double aXX;
  double aYY;
  double aZZ;
};

struct rotationalConstants
{
  double Ae;
  double Be;
  double Ce;
};

struct basis
{
  int J;
  int K;
  int M;
};

typedef std::vector<std::shared_ptr<std::vector<basis>>> basisSubset;
typedef std::vector<std::shared_ptr<matrixComp>> matrices;
typedef std::vector<std::vector<double>*> arrays;

class moleculeBase
{
public:
  polarizability pol_;
  rotationalConstants rot_;

  moleculeBase();
  moleculeBase(inputParameters &);
  // virtual std::shared_ptr<basisSubset> createBasisSets(int JMAX) = 0;
  // virtual std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubset> sets) = 0;
  // virtual std::shared_ptr<arrays> initializePopulations() = 0;
  // virtual std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>) = 0;
  // virtual std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubset> sets) = 0;

  // virtual double calculatePartitionFxn() = 0;
};


class linearMolecule : public moleculeBase
{

};

class symmetricTopMolecule : public moleculeBase
{

};

class asymmetricTopMolecule : public moleculeBase
{

};


#endif