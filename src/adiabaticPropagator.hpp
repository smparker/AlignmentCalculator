/**
 * \file "adiabaticPropagator.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_ADIABATICPROP
#define ALIGNMENTCALCULATOR_ADIABATICPROP

#include "propagatorBase.hpp"

/**
 * @brief Adiabatic Calculation Propagator
 * @details Class for managing data and outputs for aligning molecules in an adiabatic field (constant intensity)
 *
 */
class adiabaticPropagator : public propagatorBase
{
public:
  double intensity_; ///< Initial intensity and variable for holding the current field strength
  double dItn_; ///< intensity increment
  double itn_final_; ///< Final intensity
  std::string output_file_name_; ///< Output filename for observables
  std::ofstream out_file_; ///< Output file stream
  std::function<double(double)> intensity_stepper_; ///< Function to step intensity (addition and multiplication are currently supported)
  std::shared_ptr<arrays> eigenenergies_; ///< Storage for all eigenvalues at a particular intensity

  /**
   * @brief Constructor
   * @details Constructor function requiring input parameters to initialize molecule data and outputs
   *
   * @param IP input parameters
   */
  adiabaticPropagator(inputParameters &IP);

  /**
   * @brief output initializer
   * @details Creates output streams and matrix representations of all observables
   *
   * @param IP input parameters
   */
  void initializeOutputs(inputParameters &IP);

  /**
   * @brief calculation stepper
   * @details Evaluate the alignment for the current value of the intensity
   */
  void step();

  /**
   * @brief run calculation
   * @details Run all intensities
   */
  void run();

  /**
   * @brief print outputs
   * @details Calculates all observables and prints them to the output file
   */
  void printOutputs();
};

#endif