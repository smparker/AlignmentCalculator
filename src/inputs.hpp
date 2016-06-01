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

/**
 * @brief Jobtypes
 * @details Specifies adiabatic or nonadiabatic calculation
 *
 */
enum class JOBTYPE {ADIABATIC,NONADIABATIC};

/**
 * @brief Symmetry Class
 * @details Determines the symmetry influencing the basis and elements of the Hamiltonian
 *
 */
enum class MOLSYM {LINEAR,SYMMETRIC_OBLATE,SYMMETRIC_PROLATE,ASYMMETRIC};

class inputParameters
{
public:

/** @name General Inputs
 *  Input Parameters needed for all calculations.
 */
///@{
  JOBTYPE jobtype_; ///< Specifies nonadiabatic or adiabatic calculation
  std::string filename_; ///< JSON file containing all input parameters
  std::string molecule_name_; ///< Tag for output files, default is "Molecule"
  std::string library_file_; ///< Optional file containing polarizability and rotational constants
  std::string library_molecule_; ///< Name of molecule in the library file
  double rotational_temp_; ///< Initial rotational temperature
  std::vector<double> rotational_constants_; ///< Array of one (linear) or three (asymmetric or symmetric) rotational constants
  std::vector<double> polarizabilities_; ///< Three polarizability elements
  double odd_j_degeneracy_; ///< Term to modify thermal distribution for bosonic or fermionic nuclei
  double even_j_degeneracy_;///< Same as odd_j_degeneracy
  int max_j; ///< Maximum J state to include in calculation
  // bool is_full_calc_; ///< Forces full matrix calculation
  bool output_basis_list_; ///< Outputs list of basis functions
  bool output_coupling_matrix_; ///< Outputs couplings
  bool output_cos2D_; ///< Outputs projection of cosine squared onto a unit disk
  bool output_cos3D_; ///< Outputs cosine squared
  bool output_cos3DAlt_; ///< Outputs chi squared expectation value
  bool output_energy_; ///< Outputs total energy
  bool output_J_; ///< Output expectation value for J quantum number
  bool output_K_; ///< Output expectation value for K quantum number
  bool output_M_; ///< Output expectation value for M quantum number

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

  void parseAllInputs_(); ///< Read input file into property tree
  void parseJobType(boost::property_tree::ptree &); ///< Gets jobtype
  void parseMoleculeInfo(boost::property_tree::ptree &); ///< Parses molecule information
  void parseFieldInfo(boost::property_tree::ptree &); ///< Parses field information
  void parseNumericalParams(boost::property_tree::ptree &); ///< Parses numerical parameters
  void parseOutputsInfo(boost::property_tree::ptree &); ///< Parses output requests
///@}

/** @name Adiabatic Inputs
 *  Input Parameters needed for adiabatic calculations.
 */
///@{
  double initial_intensity_; ///< Initial intensity
  double final_intensity_; ///< Final intensity
  double intensity_increment_; ///< Intensity step value
  bool add_increment_; ///< Adds step rather than multiplies if true
  bool output_density_; ///< Outputs probability density of ground state in theta and phi (chi = 0 by default)
  double chi_; ///< change chi to number other than zero
  int output_eigenvectors_; ///< print eigenvector information (first 10)
  int output_eigenvalues_; ///< number of eigenvalues to output
///@}


/** @name Nonadiabatic Inputs
 *  Input Parameters needed for nonadiabatic calculations.
 */
///@{
  std::vector<pulse> pulses_; ///< List of laser pulses, all assumed to have same polarization
  int n_outputs_; ///< Total number of data points to collect
  double max_time_; ///< Time of nonadiabatic simulation
  double atol_; ///< CVODE Absolute tolerance
  double rtol_; ///< CVODE Relative tolerance

///@}

};


/**
 * @brief boost json to vector
 * @details helper function to convert boost json object into std::vector
 *
 * @param pt property tree
 * @param key property tree key
 *
 * @return vector
 */
template <typename T>
std::vector<T> as_vector(boost::property_tree::ptree const &pt, boost::property_tree::ptree::key_type const &key)
{
    std::vector<T> r;
    for (auto& item : pt.get_child(key))
        r.push_back(item.second.get_value<T>());
    return r;
}

#endif