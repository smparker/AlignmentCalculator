#include "nonadiabaticPropagator.hpp"

nonadiabaticPropagator::nonadiabaticPropagator()
{

}

nonadiabaticPropagator::nonadiabaticPropagator(inputParameters &IP)
{
  molecule_ = std::make_shared<moleculeBase>(IP);

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
