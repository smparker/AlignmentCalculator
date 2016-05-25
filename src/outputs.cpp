#include "outputs.hpp"
#include "utilities.hpp"

/**************
  Base Class
***************/

observable::observable(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{}

double observable::density_evaluate_(std::shared_ptr<matrices> densities_)
{
  double val = 0.0;
  for (int ii = 0; ii < densities_->size(); ii++)
  {
    auto temp = (*operator_matrix_->at(ii))*(*densities_->at(ii));
    val += temp.trace().real();
  }
  return val;
}

double observable::wvfxn_evaluate_(std::shared_ptr<matrices> densities_, std::shared_ptr<arrays> populations_)
{
  if (densities_->size() != populations_->size())
    throw std::runtime_error("Matrix array and population array size mismatch: " + std::to_string(densities_->size()) + " and " + std::to_string(populations_->size()));
  cplx val = 0.0;
  for (int ii = 0; ii < densities_->size(); ii++)
  {
    matrixComp operatorTemp((*operator_matrix_->at(ii))*(*densities_->at(ii)));
    for (int jj = 0; jj < operator_matrix_->at(ii)->nc(); jj++)
      val += (populations_->at(ii)->at(jj))*zdotc_(operatorTemp.nr(),&densities_->at(ii)->element(0,jj),1,&operatorTemp(0,jj),1);
  }
  return val.real();
}


/********************
  Derived Observables
*********************/

obsCosTheta3D::obsCosTheta3D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<cos^2 theta>_3D";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsCosTheta3D::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = 0; jj < N; jj++)
      {
        double interaction = FMIME(set->at(ii).J,set->at(ii).K,set->at(ii).M,0,0,set->at(jj).J,set->at(jj).K,set->at(jj).M);
        operator_matrix_->back()->element(ii,jj) = (2.0/3.0)*interaction;
        if (ii == jj)
          operator_matrix_->back()->element(ii,jj) += (1.0/3.0);
      }
    }
  }
}

obsEnergy::obsEnergy(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<H>(field-free)";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsEnergy::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = fieldFreeHamiltonians;
}

obsCosThetaAlt::obsCosThetaAlt(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<cos^2 chi>_3D";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsCosThetaAlt::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
    {
      for (int jj = 0; jj < N; jj++)
      {
        double interaction = (-1.0/3.0)     *   FMIME(set->at(ii).J,set->at(ii).K,set->at(ii).M,0,0,set->at(jj).J,set->at(jj).K,set->at(jj).M)
                            -(1.0/sqrt(6.0))*(  FMIME(set->at(ii).J,set->at(ii).K,set->at(ii).M,0,2,set->at(jj).J,set->at(jj).K,set->at(jj).M)
                                              + FMIME(set->at(ii).J,set->at(ii).K,set->at(ii).M,0,-2,set->at(jj).J,set->at(jj).K,set->at(jj).M));
        operator_matrix_->back()->element(ii,jj) = interaction;
        if (ii == jj)
          operator_matrix_->back()->element(ii,jj) += (1.0/3.0);
      }
    }
  }
}



obsJ::obsJ(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<J^2>";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsJ::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
        operator_matrix_->back()->element(ii,ii) = (set->at(ii).J+1)*(set->at(ii).J);
  }
}

obsK::obsK(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<K^2>";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsK::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
        operator_matrix_->back()->element(ii,ii) = pow(set->at(ii).K,2);
  }
}

obsM::obsM(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<M>";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsM::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
        operator_matrix_->back()->element(ii,ii) = set->at(ii).M;
  }
}


obsCosTheta2D::obsCosTheta2D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<cos^2 theta>_2D";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsCosTheta2D::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{
  operator_matrix_ = std::make_shared<matrices>();
  for (auto &set : *basisSets)
  {
    int N = set->size();
    operator_matrix_->push_back(std::make_shared<matrixComp>(N,N));
    for (int ii = 0; ii < N; ii++)
      for (int jj = 0; jj < N; jj++)
          operator_matrix_->back()->element(ii,jj) = cos2D(set->at(ii).J, set->at(ii).M, set->at(jj).J,set->at(jj).M);
  }
}



 // Prints out the basis set from a basisSubsets object
 // for (auto &a : *basis_set)
 //  {
 //    for (int ii = 0; ii < a->size(); ii++)
 //      std::cout << a->at(ii).J << " " << a->at(ii).K << " " << a->at(ii).M << std::endl;
 //    std::cout << std::endl;
 //  }



// Prints the field free Hamiltonians
  // for (auto &iter : *fieldFreeHamiltonians_)
  //   std::cout << *iter;
