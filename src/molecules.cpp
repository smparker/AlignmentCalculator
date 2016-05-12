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
  pol_.aZZ_ = IP.polarizabilities_[2];
  pol_.aXX_ = pol_.aYY_ = IP.polarizabilities_[0];
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
  double coeff = (-1.0/4.0)*std::abs(pol_.aZZ_ - pol_.aXX_);
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
    std::transform(p->begin(), p->end(), p->begin(), [&](double a){return a/partition_function_;});

  return pops;
}

std::shared_ptr<matrices> linearMolecule::initializeDensities(std::shared_ptr<arrays> pops)
{
  auto DMs = std::make_shared<matrices>();
  for (auto &p : *pops)
  {
    int N = p->size();
    DMs->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      DMs->back()->element(ii,ii) = p->at(ii);
  }
  return DMs;
}


/**************************
  Symmetric Top Molecules
***************************/

symmetricTopMolecule::symmetricTopMolecule(inputParameters &IP) :
  moleculeBase(IP)
{
  rot_.Ae_ = IP.rotational_constants_[2];
  pol_.aXX_ = IP.polarizabilities_[2];
  rot_.Be_ = IP.rotational_constants_[1];
  pol_.aYY_ = IP.polarizabilities_[1];
  rot_.Ce_ = IP.rotational_constants_[0];
  pol_.aZZ_ = IP.polarizabilities_[0];
  if (rot_.Ce_ == rot_.Be_)
    symmetry_ = MOLSYM::SYMMETRIC_PROLATE;
  else
    symmetry_ = MOLSYM::SYMMETRIC_OBLATE;
}

std::shared_ptr<basisSubsets> symmetricTopMolecule::createBasisSets(int JMAX)
{
  auto basis_set = std::make_shared<basisSubsets>();
  std::map<int,std::shared_ptr<basisSubset>> oKoJBasisSets, oKeJBasisSets, eKoJBasisSets, eKeJBasisSets;
  for (int jj = 0; jj <= JMAX; jj++)
  {
    for (int kk = -1*jj; kk <= jj; kk++)
    {
      if (jj % 2 == 0 && kk % 2 == 0 ) // J and K even
      {
        for (int mm = -1*jj; mm <= jj; mm++)
        {
          if ( eKeJBasisSets.find(mm) == eKeJBasisSets.end() )
            eKeJBasisSets[mm] = std::make_shared<basisSubset>();
          eKeJBasisSets[mm]->push_back(basis(jj,kk,mm));
        }
      }
      else if (jj % 2 == 0 && kk % 2 == 1) // J even, K odd
      {
        for (int mm = -1*jj; mm <= jj; mm++)
        {
          if ( oKeJBasisSets.find(mm) == oKeJBasisSets.end() )
            oKeJBasisSets[mm] = std::make_shared<basisSubset>();
          oKeJBasisSets[mm]->push_back(basis(jj,kk,mm));
        }
      }
      else if (jj % 2 == 1 && kk % 2 == 1 ) // J and K odd
      {
        for (int mm = -1*jj; mm <= jj; mm++)
        {
          if ( oKoJBasisSets.find(mm) == oKoJBasisSets.end() )
            oKoJBasisSets[mm] = std::make_shared<basisSubset>();
          oKoJBasisSets[mm]->push_back(basis(jj,kk,mm));
        }
      }
      else // J odd, K even
      {
        for (int mm = -1*jj; mm <= jj; mm++)
        {
          if ( eKoJBasisSets.find(mm) == eKoJBasisSets.end() )
            eKoJBasisSets[mm] = std::make_shared<basisSubset>();
          eKoJBasisSets[mm]->push_back(basis(jj,kk,mm));
        }
      }
    }
  }

  // Append all sets to a single list
  for (auto &set : oKoJBasisSets)
    basis_set->push_back(set.second);
  for (auto &set : oKeJBasisSets)
    basis_set->push_back(set.second);
  for (auto &set : eKoJBasisSets)
    basis_set->push_back(set.second);
  for (auto &set : eKeJBasisSets)
    basis_set->push_back(set.second);

  return basis_set;
}

std::shared_ptr<matrices> symmetricTopMolecule::createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets)
{
  auto ffHams = std::make_shared<matrices>();
  double A,C;
  if (symmetry_ == MOLSYM::SYMMETRIC_OBLATE)
  {
    A = rot_.Ae_;
    C = rot_.Ce_ - rot_.Ae_;
  }
  else
  {
    A = rot_.Ce_;
    C = rot_.Ae_ - rot_.Ce_;
  }

  for (auto &set : *sets)
  {
    int N = set->size();
    ffHams->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      int J = set->at(ii).J;
      int K = set->at(ii).K;
      ffHams->back()->element(ii,ii) = A*J*(J+1) + C*K*K;
    }
  }
  return ffHams;
}

std::shared_ptr<arrays> symmetricTopMolecule::initializePopulations(std::shared_ptr<basisSubsets> sets, std::shared_ptr<matrices> ffHam, double temperature)
{
  double temp;
  temperature == 0.0 ? temp = 1.0e-30 : temp = temperature; ///< Avoids divide by zero error
  auto pops = std::make_shared<arrays>();
  partition_function_ = 0.0;
  for (int ii = 0; ii < sets->size(); ii++)
  {
    auto Hamiltonian = ffHam->at(ii);
    int N = Hamiltonian->nr();
    pops->push_back(std::make_shared<std::vector<double>>(N,0.0));

    for (int jj = 0; jj < N; jj++)
    {
      double popTemp = exp(-1.0*Hamiltonian->element(jj,jj).real()/(temp*CONSTANTS::BOLTZ));
      /// Spin degeneracy statistics can be included here if need be
      pops->back()->at(jj) = popTemp;
      partition_function_ += popTemp;
    }
  }
  // Scale everything by the partition function
  for (auto &p : *pops)
    std::transform(p->begin(), p->end(), p->begin(), [&](double a){return a/partition_function_;});

  return pops;
}

std::shared_ptr<matrices> symmetricTopMolecule::initializeDensities(std::shared_ptr<arrays> pops)
{
  auto DMs = std::make_shared<matrices>();
  for (auto &p : *pops)
  {
    int N = p->size();
    DMs->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      DMs->back()->element(ii,ii) = p->at(ii);
  }
  return DMs;
}

std::shared_ptr<matrices> symmetricTopMolecule::createInteractionHamiltonians(std::shared_ptr<basisSubsets> sets)
{
  auto intHams = std::make_shared<matrices>();
  double coeff = (-1.0/4.0)*(pol_.aZZ_ - pol_.aXX_);
  if (symmetry_ == MOLSYM::SYMMETRIC_PROLATE) coeff *= -1.0;

  for (auto &set : *sets)
  {
    int N = set->size();
    intHams->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = 0; jj < N; jj++)
      {
        if (abs(set->at(ii).J - set->at(jj).J) > 2) // Skips unnecessary zero terms
          continue;
        double coupling = FMIME( set->at(ii).J, set->at(ii).K, set->at(ii).M,0,0,set->at(jj).J, set->at(jj).K ,set->at(jj).M);
        intHams->back()->element(ii,jj) += (2.0/3.0)*coeff*coupling;
      }
    }
  }
  return intHams;
}


/**************************
  Asymmetric Top Molecules
***************************/

asymmetricTopMolecule::asymmetricTopMolecule(inputParameters &IP) :
  moleculeBase(IP),
  Us_(nullptr),
  invUs_(nullptr)
{
  rot_.Ae_ = IP.rotational_constants_[2];
  pol_.aXX_ = IP.polarizabilities_[2];
  rot_.Be_ = IP.rotational_constants_[1];
  pol_.aYY_ = IP.polarizabilities_[1];
  rot_.Ce_ = IP.rotational_constants_[0];
  pol_.aZZ_ = IP.polarizabilities_[0];

  // Figure out which coordinate system to use. See Zare, pg 268-269
  if ( std::abs(IP.rotational_constants_[0] - IP.rotational_constants_[1])/IP.rotational_constants_[0] < 0.01 ) //prolate top case
  {
    std::cout << "Selecting Prolate-like coordinate system" << std::endl;
    Xe_ = IP.rotational_constants_[0]; pol_.aXX_ = IP.polarizabilities_[0];
    Ye_ = IP.rotational_constants_[2]; pol_.aYY_ = IP.polarizabilities_[2];
    Ze_ = IP.rotational_constants_[1]; pol_.aZZ_ = IP.polarizabilities_[1];
  }
  else if ( std::abs(IP.rotational_constants_[1] - IP.rotational_constants_[2])/IP.rotational_constants_[1] < 0.01 ) //oblate top case
  {
    std::cout << "Selecting Oblate-like coordinate system" << std::endl;
    Xe_ = IP.rotational_constants_[0]; pol_.aXX_ = IP.polarizabilities_[0];
    Ye_ = IP.rotational_constants_[1]; pol_.aYY_ = IP.polarizabilities_[1];
    Ze_ = IP.rotational_constants_[2]; pol_.aZZ_ = IP.polarizabilities_[2];
  }
  else
  {
    std::cout << "Selecting General Asymmetric coordinate system" << std::endl;
    Xe_ = IP.rotational_constants_[2]; pol_.aXX_ = IP.polarizabilities_[2];
    Ye_ = IP.rotational_constants_[1]; pol_.aYY_ = IP.polarizabilities_[1];
    Ze_ = IP.rotational_constants_[0]; pol_.aZZ_ = IP.polarizabilities_[0];
  }
}

std::shared_ptr<basisSubsets> asymmetricTopMolecule::createBasisSets(int JMAX)
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

std::shared_ptr<matrices> asymmetricTopMolecule::createFieldFreeHamiltonians(std::shared_ptr<basisSubsets> sets)
{
  auto ffHams = std::make_shared<matrices>();
  for (auto &set : *sets)
  {
    int N = set->size();
    ffHams->push_back(std::make_shared<matrixComp>(N,N));

    //Diagonal, symmetric top terms
    for (int ii = 0; ii < N; ii++)
    {
      int J = set->at(ii).J;
      int K = set->at(ii).K;
      ffHams->back()->element(ii,ii) = (0.5*(Xe_+Ye_)*(J*(J+1)-K*K) + Ze_*K*K);
    }
    //Off-diagonal, asymmetric couping terms
    for (int ii = 0; ii < N; ii++)
    {
      int J = set->at(ii).J;
      int K = set->at(ii).K;
      for (int jj = 0; jj < N; jj++)
      {
        int j = set->at(jj).J;
        int k = set->at(jj).K;
        if (k == K+2 && J == j)
          ffHams->back()->element(ii,jj) += ( (Xe_-Ye_)*0.25*sqrt(double(J*(J+1)-K*(K+1)))*sqrt(double(J*(J+1)-(K+1)*(K+2))) );
        else if (k == K-2 && J == j)
          ffHams->back()->element(ii,jj) += ( (Xe_-Ye_)*0.25*sqrt(double(J*(J+1)-K*(K-1)))*sqrt(double(J*(J+1)-(K-1)*(K-2))) );
      }
    }
  }
  constructTransformationMatrices(ffHams);

  return ffHams;
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

void asymmetricTopMolecule::constructTransformationMatrices(std::shared_ptr<matrices> offDiagHamiltonians)
{
  Us_    = std::make_shared<matrices>();
  invUs_ = std::make_shared<matrices>();

  for (auto &set : *offDiagHamiltonians)
  {
    // Allocate space for transformation matrices
    int N = set->nr();
    Us_   ->push_back(std::make_shared<matrixComp>(N,N));
    invUs_->push_back(std::make_shared<matrixComp>(N,N));

    // diagonalize field-free Hamiltonian, store transformation matrix
    std::vector<double> tempVec(N,0.0);
    set->diagonalize(tempVec.data());
    std::copy_n(set->data(), set->size(), Us_   ->back()->data());
    std::copy_n(set->data(), set->size(), invUs_->back()->data());
    invUs_->back()->invert();

    set->zero();
    for (int ii = 0; ii < N; ii++)
      set->element(ii,ii) = tempVec[ii];
    std::cout << *set;
  }
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
