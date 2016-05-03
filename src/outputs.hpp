/**
 * \file "outputs.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_OUTPUTS
#define ALIGNMENTCALCULATOR_OUTPUTS

#include "molecules.hpp"
#include <fstream>

class observable
{
  observable();
  observable(std::shared_ptr<moleculeBase> molecule_);
  std::ofstream outFile_;
  std::string getName_();
  void initialize_();
  void evaluate_();
};


class obsCosTheta3D : public observable
{

};

class obsCosTheta2D : public observable
{

};

class obsEnergy : public observable
{

};

class obsIntensity : public observable
{

};

class obsCosChi : public observable
{

};

class obsJ : public observable
{

};

class obsK : public observable
{

};

class obsM : public observable
{

};
#endif