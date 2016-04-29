/**
 * \file "molecules.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_MOLECULES
#define ALIGNMENTCALCULATOR_MOLECULES

typedef std::vector<std::shared_ptr<std::vector<basis>>> basisSubset;
typedef std::vector<std::shared_ptr<matrixComp>> matrices;
typedef std::vector<std::vector<double>*> arrays;

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

class moleculeBase
{
  polarizability pol_;
  rotationalConstants rot_;

  std::shared_ptr<basisSubset> createBasisSets(int JMAX);
  std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubset> sets);
  std::shared_ptr<arrays> initializePopulations();
  std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>);
  std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubset> sets);

  std::shared_ptr<arrays> populations_;

  double calculatePartitionFxn();

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