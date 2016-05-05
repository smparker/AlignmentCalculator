#include "molecules.hpp"
#include "constants.hpp"
#include <map>
#include <gsl/gsl_sf_coupling.h> // Wigner 3-j symbols

basis::basis(int j, int k, int m) : J(j), K(k), M(m) {}

/**************************
  Molecule Base Class
***************************/

moleculeBase::moleculeBase(inputParameters &IP) :
  even_j_degen_(IP.even_j_degeneracy_),
  odd_j_degen_(IP.odd_j_degeneracy_)
{}

/**************************
  Linear Molecules
***************************/

linearMolecule::linearMolecule(inputParameters &IP) :
  moleculeBase(IP)
{
  rot_.Be_  = IP.rotational_constants_[0];
  pol_.aZZ_ = IP.polarizabilities_[0];
  pol_.aXX_ = pol_.aYY_ = IP.polarizabilities_[1];
}

std::shared_ptr<basisSubsets> linearMolecule::createBasisSets(int JMAX)
{
  // Makes basis sets assuming M is conserved and J couples only to other J of the same parity
  auto basis_set = std::make_shared<basisSubsets>();
  std::map<int,std::shared_ptr<basisSubset>> oddBasisSets, evenBasisSets;
  // Even J symmetry sets
  for (int jj = 0; jj <= JMAX; jj+=2)
  {
    for (int mm = -1*jj; mm <= jj; mm++)
    {
      if (evenBasisSets.find(mm) == evenBasisSets.end())
        evenBasisSets[mm] = std::make_shared<basisSubset>();
      evenBasisSets[mm]->push_back( basis(jj,0,mm) );
    }
  }
  // Odd J symmetry sets
  for (int jj = 1; jj <= JMAX; jj+=2)
  {
    for (int mm = -1*jj; mm <= jj; mm++)
    {
      if (evenBasisSets.find(mm) == evenBasisSets.end())
        evenBasisSets[mm] = std::make_shared<basisSubset>();
      evenBasisSets[mm]->push_back( basis(jj,0,mm) );
    }
  }

  // Append all sets to a single list
  for (auto set : evenBasisSets)
    basis_set->push_back(set.second);
  for (auto set : oddBasisSets)
    basis_set->push_back(set.second);

  return basis_set;

}

std::shared_ptr<matrices> linearMolecule::createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets)
{
  auto ffHams = std::make_shared<matrices>();
  for (auto &set : *sets)
  {
    int N = set->size();
    ffHams->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      int J = set->at(ii).J;
      ffHams->back()->element(ii,ii) = rot_.Be_*J*(J+1);
    }
  }
  return ffHams;
}

std::shared_ptr<matrices> linearMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{
  auto intHams = std::make_shared<matrices>();
  double coeff = (-1.0/4.0)*abs(pol_.aZZ_ - pol_.aXX_);
  for (auto &set : *sets)
  {
    int N = set->size();
    intHams->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = 0; jj < N; jj++)
      {
        double coupling = FMIME( set->at(ii).J, set->at(ii).K, set->at(ii).M,0,0,set->at(jj).J, set->at(jj).K ,set->at(jj).M);
          intHams->back()->element(ii,jj) += (2.0/3.0)*coeff*coupling;
      }
    }
  }
  return intHams;
}

std::shared_ptr<arrays> linearMolecule::initializePopulations(std::shared_ptr<basisSubsets> sets, std::shared_ptr<matrices> ffHam, double temperature)
{
  double temp;
  temperature == 0.0 ? temp = 1.0e-30 : temp = temperature; ///< Avoids divide by zero error
  auto pops = std::make_shared<arrays>();
  partition_function_ = 0.0;
  for (int ii = 0; ii < sets->size(); ii++)
  {
    auto basisSet = sets->at(ii);
    auto Hamiltonian = ffHam->at(ii);
    int N = basisSet->size();

    pops->push_back(std::make_shared<std::vector<double>>(N,0.0));

    for (int jj = 0; jj < N; jj++)
    {
      double popTemp = exp(-1.0*Hamiltonian->element(jj,jj).real()/(temp*CONSTANTS::BOLTZ));
      (basisSet->at(jj).J % 2 == 0) ? popTemp *= even_j_degen_ : popTemp *= odd_j_degen_;
      pops->back()->at(jj) = popTemp;
      partition_function_ += popTemp;
    }
  }
  // Scale everything by the partition function
  for (auto &p : *pops)
  {
    std::transform(p->begin(), p->end(), p->begin(), [&](double a){return a/partition_function_;});
  }
  return pops;
}

std::shared_ptr<matrices> linearMolecule::initializeDensities(std::shared_ptr<arrays>)
{}


#if 0
  hamiltonians_.clear();
  cosines_.clear();
  partitions_.clear();
  eigenvectorMatrices_.clear();
  eigenvalueArrays_.clear();
  basisSets_.clear();

  //temp correction if 0 K
  if (temp_ == 0.0) temp_ = 1.0e-100;


  /// Make Basis Sets


  // Construct field free Hamiltonians (diagonal), cosine matrices and partition functions
  for (auto basisSet : basisSets_)
  {
    int N = basisSet->size();
    hamiltonians_.push_back(make_shared<matrixReal>(N,N));
    cosines_.push_back(make_shared<matrixReal>(N,N));
    partitions_.push_back(make_shared<vector<double>>());
    for (int ii = 0; ii < basisSet->size(); ii++)
    {
      int J = (*basisSet)[ii][0];
      double energy = Be_*J*(J+1);
      hamiltonians_.back()->element(ii,ii) = energy;
      double nuclear_degeneracy = 1.0;
      if ( (J % 2 == 0) && useNuclearDegen_) nuclear_degeneracy = 3.0;
      partitions_.back()->push_back(nuclear_degeneracy*exp(-1.0*energy/(temp_*3.1669e-6)));
    }
  }
  // hamiltonians_[0]->printMem(); ///< Check how much memory we're using

  // Fill the cosine matrices
  for (int basisNumber = 0; basisNumber < basisSets_.size(); basisNumber++)
  {
    shared_ptr<basisStates> bSet = basisSets_[basisNumber];
    for (int ii = 0; ii < bSet->size(); ii++)
    {
      for (int jj = 0; jj < bSet -> size(); jj++)
      {
        double interaction = FMIME( (*bSet)[ii][0], 0, (*bSet)[ii][2], 0, 0, (*bSet)[jj][0], 0, (*bSet)[jj][2]);
        (ii == jj) ? (cosines_[basisNumber]->element(ii,jj) = (1.0/3.0)*(1.0+2.0*interaction)) :( cosines_[basisNumber]->element(ii,jj) = (2.0/3.0)*interaction) ;
      }
    }
  }

  // Sum and scale all of the partition functions
  double partitionSum = 0.0;
  for (auto partition : partitions_)
    partitionSum += accumulate(partition->data(),partition->data()+partition->size(),0.0);
  cout << partitionSum << endl;
  for (auto partition : partitions_)
    transform(partition->begin(),partition->end(),partition->begin(),[&](double a){return a/partitionSum;});

  // Allocate calculation space for diagonalizations
  for (auto basisSet : basisSets_)
  {
    int N = basisSet->size();
    eigenvectorMatrices_.push_back(make_shared<matrixReal>(N,N));
    eigenvalueArrays_.push_back(make_shared<vector<double>>(N,0.0));
  }
  // eigenvectorMatrices_[0]->printMem();
#endif



/**************************
  Symmetric Top Molecules
***************************/

symmetricTopMolecule::symmetricTopMolecule(inputParameters &IP) :
  moleculeBase(IP)
{

}

std::shared_ptr<basisSubsets> symmetricTopMolecule::createBasisSets(int JMAX)
{}

std::shared_ptr<matrices> symmetricTopMolecule::createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets)
{}

std::shared_ptr<arrays> symmetricTopMolecule::initializePopulations(std::shared_ptr<basisSubsets> sets, std::shared_ptr<matrices> ffHam, double temperature)
{}

std::shared_ptr<matrices> symmetricTopMolecule::initializeDensities(std::shared_ptr<arrays>)
{}

std::shared_ptr<matrices> symmetricTopMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{}


/**************************
  Asymmetric Top Molecules
***************************/

asymmetricTopMolecule::asymmetricTopMolecule(inputParameters &IP) :
  moleculeBase(IP)
{

}

std::shared_ptr<basisSubsets> asymmetricTopMolecule::createBasisSets(int JMAX)
{

}

std::shared_ptr<matrices> asymmetricTopMolecule::createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets)
{

}

std::shared_ptr<arrays> asymmetricTopMolecule::initializePopulations(std::shared_ptr<basisSubsets> sets, std::shared_ptr<matrices> ffHam, double temperature)
{

}

std::shared_ptr<matrices> asymmetricTopMolecule::initializeDensities(std::shared_ptr<arrays>)
{

}

std::shared_ptr<matrices> asymmetricTopMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{

}

/**************************
  Molecule helper functions
***************************/

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
double FMIME (int J, int K, int M, int Q, int S, int j, int k, int m)
{
  double coeff = pow(-1.0,k+m)*sqrt((2.0*J + 1.0)*(2.0*j+1.0));
  double J1    = gsl_sf_coupling_3j(2*J, 4, 2*j, 2*M, 2*Q, -2*m);
  double J2    = gsl_sf_coupling_3j(2*J, 4, 2*j, 2*K, 2*S, -2*k);
  return coeff*J2*J1;
}
