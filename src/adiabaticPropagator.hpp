/**
 * \file "adiabaticPropagator.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_ADIABATICPROP
#define ALIGNMENTCALCULATOR_ADIABATICPROP

#include "propagatorBase.hpp"

class adiabaticPropagator : public propagatorBase
{
public:
  double intensity_, dItn_, itn_final_;
  std::string output_file_name_;
  std::ofstream out_file_;
  std::function<double(double)> intensity_stepper_;
  std::shared_ptr<arrays> eigenenergies_;
  adiabaticPropagator(inputParameters &IP);
  void initializeOutputs(inputParameters &IP);
  void step();
  void run();
  void printOutputs();
};

#endif