#include "propagatorBase.hpp"

propagatorBase::propagatorBase(inputParameters &IP) :
  temperature_(IP.rotational_temp_)
{
  symmetry_ = determineSymmetry(IP);
  if (symmetry_ == MOLSYM::LINEAR)
    molecule_ = std::make_shared<linearMolecule>(IP);
  else if ( (symmetry_ == MOLSYM::SYMMETRIC_PROLATE) || (symmetry_ == MOLSYM::SYMMETRIC_OBLATE))
    molecule_ = std::make_shared<symmetricTopMolecule>(IP);
  else
    molecule_ = std::make_shared<asymmetricTopMolecule>(IP);
  basisSets_             = molecule_->createBasisSets(IP.max_j);
  initialize_();
}

void propagatorBase::initialize_()
{
  fieldFreeHamiltonians_ = molecule_->createFieldFreeHamiltonians(basisSets_);
  intHamiltonians_       = molecule_->createInteractionHamiltonians(basisSets_);
  populations_           = molecule_->initializePopulations(basisSets_,fieldFreeHamiltonians_,temperature_);
  partition_function_    = molecule_->partition_function_;
  std::cout << "Partition function at " << temperature_ << "K is " << partition_function_ << std::endl;
  densities_             = molecule_->initializeDensities(populations_);
  removeSmallPopulations();
}

void propagatorBase::removeSmallPopulations()
{
  // if population for a basis subset is below a threshold, remove from calculation
  std::vector<bool> isPopulated(populations_->size(),false);
  for (int ii = 0; ii < populations_->size(); ii++)
    if (std::accumulate(populations_->at(ii)->begin(), populations_->at(ii)->end(), 0.0) > 1.0e-6)
      isPopulated[ii] = true;

  // Removes starting from the end of the list to avoid index shifting
  for (int ii = isPopulated.size()-1; ii >= 0; ii--)
  {
    if (!isPopulated[ii])
    {
      basisSets_->erase(basisSets_->begin()+ii);
      fieldFreeHamiltonians_->erase(fieldFreeHamiltonians_->begin()+ii);
      intHamiltonians_->erase(intHamiltonians_->begin()+ii);
      populations_->erase(populations_->begin()+ii);
      densities_->erase(densities_->begin()+ii);
      if (molecule_->Us_ != nullptr)
      {
        molecule_->Us_->erase(molecule_->Us_->begin()+ii);
        molecule_->invUs_->erase(molecule_->invUs_->begin()+ii);
      }
    }
  }
}

MOLSYM propagatorBase::determineSymmetry(inputParameters &IP)
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

    std::sort(symVec.begin(),symVec.end(),[](std::pair<double,double> a, std::pair<double,double> b){return a.first < b.first;});

    IP.rotational_constants_ = {symVec[0].first,symVec[1].first,symVec[2].first};
    IP.polarizabilities_ = {symVec[0].second,symVec[1].second,symVec[2].second};

    if (symVec[0].first == symVec[1].first)
    {
      std::cout << "SYMMETRY: Prolate Symmetric Top" << std::endl;
      return MOLSYM::SYMMETRIC_PROLATE;
    }
    else if (symVec[1].first == symVec[2].first)
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

void propagatorBase::outputBasisStats()
{
  std::ofstream outFile;
  outFile.open("basisStatistics.txt");
  for (int ii = 0; ii < basisSets_->size(); ii++)
  {
    for (int jj = 0; jj < basisSets_->at(ii)->size(); jj++)
      outFile << basisSets_->at(ii)->at(jj).J << " "
              << basisSets_->at(ii)->at(jj).K << " "
              << basisSets_->at(ii)->at(jj).M << " "
              << fieldFreeHamiltonians_->at(ii)->element(jj,jj).real() << " "
              << densities_->at(ii)->element(jj,jj).real() << "\n";
  }
  outFile.close();
}

void propagatorBase::transformObservables()
{
  for (auto &obs : observables_)
    for (int ii = 0; ii < obs->operator_matrix_->size(); ii++)
      obs->operator_matrix_->at(ii) = std::make_shared<matrixComp>( *(molecule_->invUs_->at(ii)) * *(obs->operator_matrix_->at(ii)) * *(molecule_->Us_->at(ii)) );
}
