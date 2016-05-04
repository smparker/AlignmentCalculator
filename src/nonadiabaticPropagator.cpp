#include "nonadiabaticPropagator.hpp"

nonadiabaticPropagator::nonadiabaticPropagator()
{

}

nonadiabaticPropagator::nonadiabaticPropagator(inputParameters &IP)
{
  symmetry_ = determineSymmetry(IP);
  if (symmetry_ == MOLSYM::LINEAR)
    molecule_ = std::make_shared<linearMolecule>(IP);
  else if ( (symmetry_ == MOLSYM::SYMMETRIC_PROLATE) || (symmetry_ == MOLSYM::SYMMETRIC_OBLATE))
    molecule_ = std::make_shared<symmetricTopMolecule>(IP);
  else
    molecule_ = std::make_shared<asymmetricTopMolecule>(IP);
}

MOLSYM nonadiabaticPropagator::determineSymmetry(inputParameters &IP)
{
  // If only one rotational constant provided, assume linear
  if (IP.rotational_constants_.size() == 1)
  {
    std::sort(IP.polarizabilities_.begin(),IP.polarizabilities_.end());
    std::cout << "SYMMETRY: Linear Molecule" << std::endl;
    return MOLSYM::LINEAR;
  }
  // Sort based on rotational constant
  else
  {
    std::vector<std::pair<double,double>> symVec;
    for (int ii = 0; ii < 3; ii++)
      symVec.push_back(std::make_pair(IP.rotational_constants_[ii],IP.polarizabilities_[ii]));

    std::sort(symVec.begin(),symVec.end(),[](std::pair<double,double> a, std::pair<double,double> b){return a.first > b.first;});

    IP.rotational_constants_ = {symVec[0].first,symVec[1].first,symVec[2].first};
    IP.polarizabilities_ = {symVec[0].second,symVec[1].second,symVec[2].second};

    if (symVec[1].first == symVec[2].first)
    {
      std::cout << "SYMMETRY: Prolate Symmetric Top" << std::endl;
      return MOLSYM::SYMMETRIC_PROLATE;
    }
    else if (symVec[0].first == symVec[1].first)
    {
      std::cout << "SYMMETRY: Oblate Symmetric Top" << std::endl;
      return MOLSYM::SYMMETRIC_OBLATE;
    }
    else
    {
      std::cout << "SYMMETRY: Asymmetric Top" << std::endl;
      return MOLSYM::ASYMMETRIC;
    }
  }
  throw std::runtime_error("Cannot determine molecule symmetry. Not sure how this line of code was ever reached.");
}




// class nonadiabaticPropagator
// {
//   double partitionFxn_;
//   double t0_;
//   double tFinal_;
//   double noutputs_;
//   double temperature_;
//   std::shared_ptr<moleculeBase> molecule_;
//   std::shared_ptr<basisSubset> basisSets_;
//   std::shared_ptr<matrices> fieldFreeHamiltonians_;
//   std::shared_ptr<matrices> intHamiltonians_;
//   std::shared_ptr<matrices> densities_;
//   std::shared_ptr<arrays> populations_;
//   std::vector<observable> observables_;
//   std::vector<pulse> pulses_;

//   void setupOutputs();
//   void initializeCVODE();
//   void evalRHS();
//   void step();
//   void run();
// };




//   polarizability pol_;
//   rotationalConstants rot_;

//   moleculeBase();
//   moleculeBase(inputParameters &);
//   virtual std::shared_ptr<basisSubset> createBasisSets(int JMAX) = 0;
//   virtual std::shared_ptr<matrices> createFieldFreeHamiltonians(std::shared_ptr<basisSubset> sets) = 0;
//   virtual std::shared_ptr<arrays> initializePopulations() = 0;
//   virtual std::shared_ptr<matrices> initializeDensities(std::shared_ptr<arrays>) = 0;
//   virtual std::shared_ptr<matrices> createInteractionHamiltonians(std::shared_ptr<basisSubset> sets) = 0;
