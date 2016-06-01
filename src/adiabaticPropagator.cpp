#include "adiabaticPropagator.hpp"
#include "utilities.hpp"
#include "constants.hpp"
#include <sstream>

adiabaticPropagator::adiabaticPropagator(inputParameters &IP) :
  propagatorBase(IP),
  intensity_(IP.initial_intensity_),
  dItn_(IP.intensity_increment_),
  itn_final_(IP.final_intensity_)
{
  print_n_vecs_ = IP.output_eigenvectors_;
  print_n_vals_ = IP.output_eigenvalues_;
  if (IP.add_increment_)
    intensity_stepper_ = [&](double a){return a+dItn_;};
  else
    intensity_stepper_ = [&](double a){return a*dItn_;};
  initializeOutputs(IP);
  if ( (molecule_->Us_ != nullptr) && (molecule_->invUs_ != nullptr) )
    transformObservables();
}

void adiabaticPropagator::initializeOutputs(inputParameters &IP)
{
  // Setup Output stream
  output_file_name_ = IP.molecule_name_ + "_ad.txt";
  out_file_.open(output_file_name_);

  // Initialize observable objects
  if (IP.output_cos3D_) observables_.push_back(std::make_shared<obsCosTheta3D>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_energy_) observables_.push_back(std::make_shared<obsEnergy>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_cos3DAlt_) observables_.push_back(std::make_shared<obsCosThetaAlt>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_J_) observables_.push_back(std::make_shared<obsJ>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_K_) observables_.push_back(std::make_shared<obsK>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_M_) observables_.push_back(std::make_shared<obsM>(basisSets_,fieldFreeHamiltonians_));
  if (IP.output_cos2D_ && (symmetry_ == MOLSYM::LINEAR))
    observables_.push_back(std::make_shared<obsM>(basisSets_,fieldFreeHamiltonians_));
  else if (IP.output_cos2D_ && (symmetry_ != MOLSYM::LINEAR))
    std::cout << "Output for the cos^2 projection in 2D not supported for molecules of nonlinear symmetry. Please set the output to false to supress this message." << std::endl;

  // Print data identifiers
  out_file_ << "# Intensity(W/cm^2)";
  for (auto obs : observables_)
    out_file_ << "\t" << obs->id_tag_;
  out_file_ << std::endl;
}


void adiabaticPropagator::printOutputs()
{
  out_file_ << intensity_*CONSTANTS::LASERINTEN << "\t";
  for (auto obs : observables_)
    out_file_ << "\t" << obs->wvfxn_evaluate_(densities_,populations_);
  if (print_n_vals_ > 0)
    out_file_ << print_lowest_eigenvals();
  out_file_ << std::endl;
}

void adiabaticPropagator::run()
{
  while (intensity_ < itn_final_ )
  {
    intensity_ = intensity_stepper_(intensity_);
    step();
    printOutputs();
  }
}

void adiabaticPropagator::step()
{
  eigenenergies_ = std::make_shared<arrays>();
  densities_ = std::make_shared<matrices>(); // Clear this array to use as storage space
  for (int ii = 0; ii < basisSets_->size(); ii++)
  {
    auto totalH = std::make_shared<matrixComp>(*(fieldFreeHamiltonians_->at(ii)) + *(intHamiltonians_->at(ii))*intensity_);

    auto energies = std::make_shared<std::vector<double>>(populations_->at(ii)->size(),0.0);
    totalH->diagonalize(energies->data());
    densities_->push_back(totalH);
    eigenenergies_->push_back(energies);
  }
}

std::string adiabaticPropagator::print_lowest_eigenvals()
{
  std::vector<double> vals;
  for (auto &ens : *eigenenergies_)
  {
// Grabs energies
    if (ens->size() < print_n_vals_)
      vals.insert(vals.end(), ens->begin(), ens->end());
    else
      vals.insert(vals.end(), ens->begin(), ens->begin()+print_n_vals_);
// Sort new global list
    std::sort(vals.begin(),vals.end());
// if too big, remove end
    if (vals.size() > print_n_vals_)
      vals.erase(vals.begin() + print_n_vals_, vals.end());
  }
  std::stringstream out_string;
  for (auto &v : vals)
    out_string << std::scientific << std::setprecision(6) << "\t" << v;
  return out_string.str();
}

