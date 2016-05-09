#include "outputs.hpp"

observable::observable(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{

}



obsCosTheta3D::obsCosTheta3D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) :
  observable(basisSets,fieldFreeHamiltonians)
{
  id_tag_ = "<cos^2 theta>_3D";
  initialize_(basisSets,fieldFreeHamiltonians);
}

void obsCosTheta3D::initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians)
{

}

double obsCosTheta3D::evaluate_(std::shared_ptr<matrices> densities_)
{
  return 0;
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
