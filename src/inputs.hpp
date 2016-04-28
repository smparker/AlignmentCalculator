/**
 * \file "inputs.hpp"
 * \author  J. Szekely
 */

#ifndef ALIGNMENTCALCULATOR_INPUTS
#define ALIGNMENTCALCULATOR_INPUTS

#include <string>
#include <iostream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "pulses.hpp"


enum class JOBTYPE {ADIABATIC,NONADIABATIC};
enum class MOLSYM {LINEAR,SYMMETRIC,ASYMMETRIC};

class inputParameters
{

public:
/** @name General Inputs
 *  Input Parameters needed for all calculations.
 */
///@{
  std::string filename_;
  std::string molecule_name_;
  std::string library_file_;
  std::string library_molecule_;
  double rotational_temp_;
  std::vector<double> rotational_constants_;
  std::vector<double> polariabilities_;
  double odd_j_degeneracy_;
  double even_j_degeneracy_;
  int max_j;
  bool output_basis_list_;
  bool output_coupling_matrix_;
  bool oupput_cos2D_;
  bool output_cos3D_;
  bool output_energy_;
  bool output_intensity_;

  /**
   * @brief Parse inputs
   * @details Gets all calculation inputs from a json file
   *
   * @param infile filename
   */
  inputParameters(std::string infile);

  /**
   * @brief Strips all comments
   * @details Removes all comments in input file that begin with '//'
   */
  void stripComments_();
  void parseAllInputs_();
  void parseJobType(boost::property_tree::ptree &);
  void parseMoleculeInfo(boost::property_tree::ptree &);
  void parseFieldInfo(boost::property_tree::ptree &);
  void parseNumericalParams(boost::property_tree::ptree &);
  void parseOutputsInfo(boost::property_tree::ptree &);
///@}

/** @name Adiabatic Inputs
 *  Input Parameters needed for adiabatic calculations.
 */
///@{
  double intial_intensity_;
  double final_intensity_;
  double intensity_increment_;
  bool add_increment_;
  bool output_density_;
  bool output_eigenvectors_;

///@}


/** @name Nonadiabatic Inputs
 *  Input Parameters needed for nonadiabatic calculations.
 */
///@{
  std::vector<pulse> pulses_;
  int n_outputs_;

///@}

};

#endif