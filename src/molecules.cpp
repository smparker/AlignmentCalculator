#include "molecules.hpp"
#include <map>

basis::basis(int j, int k, int m) : J(j), K(k), M(m) {}

/**************************
  Molecule Base Class
***************************/

moleculeBase::moleculeBase(inputParameters &IP)
{

}

/**************************
  Linear Molecules
***************************/

linearMolecule::linearMolecule(inputParameters &IP) :
  moleculeBase(IP)
{
  createBasisSets(IP.max_j);
  // Be_ = 
}

std::shared_ptr<basisSubsets> linearMolecule::createBasisSets(int JMAX)
{
  auto basis_set = std::make_shared<basisSubsets>();
  std::map<int,std::shared_ptr<basisSubset>> oddBasisSets, evenBasisSets;
  for (int jj = 0; jj <= JMAX; jj++)
  {
    if (jj % 2 == 0)
    {
      for (int mm = -1*jj; mm <= jj; mm++)
      {
        if (evenBasisSets.find(mm) == evenBasisSets.end())
          evenBasisSets[mm] = std::make_shared<basisSubset>();
        evenBasisSets[mm]->push_back( basis(jj,0,mm) );
      }
    }
    else
    {
      for (int mm = -1*jj; mm <= jj; mm++)
      {
        if (oddBasisSets.find(mm) == oddBasisSets.end())
          oddBasisSets[mm] = std::make_shared<basisSubset>();
        oddBasisSets[mm]->push_back( basis(jj,0,mm) );
      }
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
{}

std::shared_ptr<arrays> linearMolecule::initializePopulations()
{}

std::shared_ptr<matrices> linearMolecule::initializeDensities(std::shared_ptr<arrays>)
{}

std::shared_ptr<matrices> linearMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{}

double linearMolecule::calculatePartitionFxn()
{}

#if 0
hamiltnians_.clear();
  cosines_.clear();
  partitions_.clear();
  eigenvectorMatrices_.clear();
  eigenvalueArrays_.clear();
  basisSets_.clear();

  //temp correction if 0 K
  if (temp_ == 0.0) temp_ = 1.0e-100;

  // Construct the basis sets
  map<int,shared_ptr<basisStates>> oddBasisSets;
  map<int,shared_ptr<basisStates>> evenBasisSets;
  for (int jj = 0; jj <= jStates_; jj++)
  {
    if (jj % 2 == 0)
    {
      for (int mm = -1*jj; mm <= jj; mm++)
      {
        if (evenBasisSets.find(mm) == evenBasisSets.end()) evenBasisSets[mm] = make_shared<basisStates>();
        evenBasisSets[mm]->push_back({jj,0,mm});
      }
    }
    else
    {
      for (int mm = -1*jj; mm <= jj; mm++)
      {
        if (oddBasisSets.find(mm) == oddBasisSets.end()) oddBasisSets[mm] = make_shared<basisStates>();
        oddBasisSets[mm]->push_back({jj,0,mm});
      }
    }
  }
  // Append all sets to a single list
  for (auto set : evenBasisSets)
    basisSets_.push_back(set.second);
  for (auto set : oddBasisSets)
    basisSets_.push_back(set.second);

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

std::shared_ptr<arrays> symmetricTopMolecule::initializePopulations()
{}

std::shared_ptr<matrices> symmetricTopMolecule::initializeDensities(std::shared_ptr<arrays>)
{}

std::shared_ptr<matrices> symmetricTopMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{}

double symmetricTopMolecule::calculatePartitionFxn()
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

std::shared_ptr<arrays> asymmetricTopMolecule::initializePopulations()
{

}

std::shared_ptr<matrices> asymmetricTopMolecule::initializeDensities(std::shared_ptr<arrays>)
{

}

std::shared_ptr<matrices> asymmetricTopMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{

}

double asymmetricTopMolecule::calculatePartitionFxn()
{

}





