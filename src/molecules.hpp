/**
 * \file "molecules.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_MOLECULES
#define ALIGNMENTCALCULATOR_MOLECULES

struct polarizability
{
  double aXX;
  double aYY;
  double aZZ;
}

struct rotationalConstants
{
  double Ae;
  double Be;
  double Ce;
}

struct basis
{
  int J;
  int K;
  int M;
}

class moleculeBase
{
  polarizability pol_;
  rotationalConstants rot_;


  std::vector<std::vector<basis>*> createBasisSets(int JMAX);
  double calculatePartitionFxn();

}


class linearMolecule : public moleculeBase
{

}

class symmetricTopMolecule : public moleculeBase
{

}

class asymmetricTopMolecule : public moleculeBase
{

}


#endif