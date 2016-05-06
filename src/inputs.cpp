#include "inputs.hpp"
#include <boost/filesystem.hpp>
#include <algorithm>
#include <cmath>
#include "constants.hpp"

inputParameters::inputParameters(std::string fn) : filename_(fn)
{
  try
  {
    stripComments_();
  }
  catch(...)
  {
    throw std::runtime_error("Make sure directory \"output_data\" exists in working directory. This is an probably from an error in the boost library that will be fixed in the future.");
  }
  parseAllInputs_();
}

void inputParameters::stripComments_()
{
//  Open input and output file
  // FIXME
  #if 0
  std::cout << boost::filesystem::current_path() << std::endl;
  boost::filesystem::path p ("output_data/");
  if (!boost::filesystem::exists(p))
    boost::filesystem::create_directory("output_data/");
  #endif

  std::string newfn = "output_data/" + filename_;
  std::fstream inputfile;
  std::ofstream inputcopy;
  inputfile.open(filename_);
  inputcopy.open(newfn);

//  search for '//', delete everything following, print remainder to new file
  std::string line;
  int found, found2;
  while (getline(inputfile,line))
  {
    found  = line.find('/');
    found2 = line.find('/', found+1);
    if (found != line.npos && found2 == found+1)
      inputcopy << line.erase(found, line.length()) << std::endl;
    else
      inputcopy << line << std::endl;
  }
  inputcopy.close();
  inputfile.close();

//  update filename;
  filename_ = newfn;
}


void inputParameters::parseAllInputs_()
{
  boost::property_tree::ptree IP;
  boost::property_tree::json_parser::read_json(filename_,IP);
  parseJobType(IP);
  parseMoleculeInfo(IP);
  parseFieldInfo(IP);
  parseNumericalParams(IP);
  parseOutputsInfo(IP);
}


void inputParameters::parseJobType(boost::property_tree::ptree &IP)
{
  std::string job = IP.get<std::string>("Jobtype","none");
  std::transform(job.begin(), job.end(), job.begin(), tolower);
  if (job == "adiabatic")
    jobtype_ = JOBTYPE::ADIABATIC;
  else if (job == "nonadiabatic")
    jobtype_ = JOBTYPE::NONADIABATIC;
  else
    throw std::runtime_error("Incorrect jobtype, only keywords 'nonadiabatic' and 'adiabatic' are supported");
}

void inputParameters::parseMoleculeInfo(boost::property_tree::ptree &IP)
{
// Get tag which will be used to name output files
  molecule_name_ = IP.get<std::string>("Molecule.molecule_name","Molecule");

// Check if library file and molecule name have been provided
// Proceed to import molecule data or read from input file
  library_file_ = IP.get<std::string>("Molecule.library_file","NONE");
  library_molecule_ = IP.get<std::string>("Molecule.library_molecule","NONE");
  if ((library_file_ == "NONE") || (library_molecule_ == "NONE"))
  {
    // Read Polarizability
    polarizabilities_ = as_vector<double>(IP,"Molecule.polarizability_components");
    if (polarizabilities_.size() != 3)
      throw std::runtime_error("Error reading the correct number of polarizability components");

    // Read Rotational Constants
    rotational_constants_ = as_vector<double>(IP,"Molecule.rotational_constants");
    if ((rotational_constants_.size() != 3) && (rotational_constants_.size() != 1) )
      throw std::runtime_error("Error: the number of rotational constants and polarizability components do not match: " + std::to_string(rotational_constants_.size()));
  }
  else
  {
    // Open molecular information json file
    boost::property_tree::ptree molIP;
    boost::property_tree::json_parser::read_json(library_file_,molIP);

    // Get polarization in either list format or listing "axx, ayy, azz"
    try
    {
      polarizabilities_ = as_vector<double>(molIP,library_molecule_+".polarizability_components");
    }
    catch(...)
    {
      if (polarizabilities_.size() == 0)
        polarizabilities_ = {
          molIP.get<double>(library_molecule_ + ".Polarizability.XX",-1.0),
          molIP.get<double>(library_molecule_ + ".Polarizability.YY",-1.0),
          molIP.get<double>(library_molecule_ + ".Polarizability.ZZ",-1.0) };
      if ((polarizabilities_[0] == -1.0) || (polarizabilities_[1] == -1.0) || (polarizabilities_[2] == -1.0))
        throw std::runtime_error("Error: Cannot import polarizability components from requested molecule library file");
    }
    // Get rotational constants in either list format or listing "axx, ayy, azz"
    try
    {
      rotational_constants_ = as_vector<double>(molIP,library_molecule_+".rotational_constants");
    }
    catch(...)
    {
    if (rotational_constants_.size() == 0)
      rotational_constants_ = {
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Ae",-1.0),
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Be",-1.0),
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Ce",-1.0) };
    if ((rotational_constants_[0] == -1.0) || (rotational_constants_[1] == -1.0) || (rotational_constants_[2] == -1.0))
      throw std::runtime_error("Error: Cannot import rotational constants from requested molecule library file");
    }
  }
// Get initial temperature of the rotational ensemble
  rotational_temp_ = IP.get<double>("Molecule.rotational_temperature",-1.0);
  if (rotational_temp_ < 0)
    throw std::runtime_error("Either your forgot to specify the \"rotational_temperature\" variable or you gave a negative number. Please fix it.");
// Get optional degeneracy terms
  odd_j_degeneracy_  = IP.get<double>("Molecule.odd_j_degeneracy",1.0);
  even_j_degeneracy_ = IP.get<double>("Molecule.even_j_degeneracy",1.0);

}

void inputParameters::parseFieldInfo(boost::property_tree::ptree &IP)
{
  if (jobtype_ == JOBTYPE::ADIABATIC)
  {
  // Import parameters
    initial_intensity_   = IP.get<double>("Field.initial_intensity",1.0e8);
    final_intensity_     = IP.get<double>("Field.final_intensity",1.0e12);
    intensity_increment_ = IP.get<double>("Field.intensity_increment",10.0);
    add_increment_       = IP.get<bool>("Field.add_increment",false);
  // Change to atomic units
    initial_intensity_ /= CONSTANTS::LASERINTEN;
    final_intensity_   /= CONSTANTS::LASERINTEN;
    if (add_increment_)
      intensity_increment_ /= CONSTANTS::LASERINTEN;
  }
  else if (jobtype_ == JOBTYPE::NONADIABATIC)
  {
    double FWHMFactor = 2.0*sqrt(2.0*log(2.0));
    for (auto &p : IP.get_child("Field.pulses"))
    {
      double sig = p.second.get<double>("pulse_fwhm",100.0)*CONSTANTS::AUperFS/ FWHMFactor;
      pulses_.push_back( pulse(
        p.second.get<double>("pulse_max_intensity",1.0e12)/CONSTANTS::LASERINTEN,
        sig,
        p.second.get<double>("pulse_center",10.0)*CONSTANTS::AUperFS));
    }
  }
}

void inputParameters::parseNumericalParams(boost::property_tree::ptree &IP)
{
  max_j = IP.get<int>("Numerical.maximum_J",50);
  if (jobtype_ == JOBTYPE::NONADIABATIC)
  {
    n_outputs_ = IP.get<int>("Numerical.n_outputs",1000);
    max_time_  = IP.get<double>("Numerical.maximum_time",100.0)*CONSTANTS::AUperFS*1000.0;
    atol_ = IP.get<double>("Numerical.atol", 1.0e-8);
    rtol_ = IP.get<double>("Numerical.rtol", 1.0e-6);
  }
}

void inputParameters::parseOutputsInfo(boost::property_tree::ptree &IP)
{

}