#include "inputs.hpp"
#include <boost/filesystem.hpp>
#include <algorithm>
#include "constants.hpp"

inputParameters::inputParameters(std::string fn) : filename_(fn)
{
  try
  {
    stripComments_();
  }
  catch(...)
  {
    throw std::runtime_error("Make sure directory \"output_data\" exists in workign directory. This is an probably from an error in the boost library that will be fixed in the future.");
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
    polarizabilities_ = as_vector<double>(molIP,library_molecule_+".polarizability_components");
    if (polarizabilities_.size() == 0)
      polarizabilities_ = {
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Ae",-1.0),
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Be",-1.0),
        molIP.get<double>(library_molecule_ + ".RotationalConstants.Ce",-1.0) };
    if ((polarizabilities_[0] == -1.0) || (polarizabilities_[1] == -1.0) || (polarizabilities_[2] == -1.0))
      throw std::runtime_error("Error: Cannot import polarizability components from requested molecule library file");
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
    constexpr double FWHMFactor = 2.0*sqrt(2.0*log(2.0));
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

}

void inputParameters::parseOutputsInfo(boost::property_tree::ptree &IP)
{

}



// inputParameters::inputParameters(std::string fn) : filename_(fn)
// {
//   stripComments_();
  // boost::property_tree::ptree IP;
  // boost::property_tree::json_parser::read_json(filename_,IP);
//   minPower_ = IP.get<double>("Field.min_laser_power",1.0e6);
//   maxPower_ = IP.get<double>("Field.max_laser_power",1.0e16);
//   factor_   = IP.get<double>("Field.power_factor",10.0);
//   if (minPower_ > maxPower_) maxPower_ = minPower_;
//   linearMoleculeFile_ = IP.get<string>("MoleculeParameters.LinearParamFile","Molecules.json");
//   symTopMoleculeFile_ = IP.get<string>("MoleculeParameters.SymTopParamFile","SymTopMolecules.json");
//   molecule_           = IP.get<string>("MoleculeParameters.Molecule","Cl2");
//   moleculeType_       = IP.get<int>("MoleculeParameters.RotationalType",0);
//   nMolecules_         = IP.get<int>("MoleculeParameters.nMolecules",10);
//   temperature_        = IP.get<double>("MoleculeParameters.RotTemp",1.0e-8);
//   transTemp_          = IP.get<double>("MoleculeParameters.TransTemp",1.0);
//   trajRepetition_     = IP.get<int>("MoleculeParameters.trajRepetition",1);
//   transField_         = IP.get<double>("MoleculeParameters.transField",1.0e10);
//   rampIntensity_      = IP.get<bool>("MoleculeParameters.rampIntensity",false);
//   rampMax_            = IP.get<double>("MoleculeParameters.rampMax",2.0);

//   jStates_ = IP.get<int>("CalculationParameters.JStates",30);
//   useM_    = IP.get<int>("CalculationParameters.UseM",0);
//   useOdd_  = IP.get<int>("CalculationParameters.UseOdd",0);

//   outFile_         = IP.get<string>("Outputs.outfile","output_data/AlignmentData.txt");
//   outJBasis_       = IP.get<int>("Outputs.JBasis",0);
//   outCouplings_    = IP.get<bool>("Outputs.Couplings",false);
//   outEigenvectors_ = IP.get<bool>("Outputs.EigenVectors",false);
//   outAngularRep_   = IP.get<bool>("Outputs.AngularRep",false);
//   cos2D_           = IP.get<bool>("Outputs.cos2D",false);
//   n_out_           = IP.get<int>("Outputs.n_out",100);

//   pulseDuration_ = IP.get<double>("PulseParameters.PulseDuration",100.0);
//   peakIntensity_ = IP.get<double>("PulseParameters.PeakIntensity",1.0);
//   t0_            = IP.get<double>("PulseParameters.InitialTime",10.0);
//   ramp_width_    = IP.get<double>("PulseParameters.ramp_width",41.0*1000.0);
//   ramp_start_    = IP.get<double>("PulseParameters.ramp_start",41.0*10000.0);


// //  Copies the json data to a file to check for debugging
//   ofstream ss;
//   ss.open("output_data/inputs_check.json");
//   boost::property_tree::write_json(ss,IP);
//   ss.close();

//   boost::property_tree::ptree mol;
//   if (moleculeType_ == 1 || moleculeType_ == 2)
//   {
// //  Defaults to HCCCCl3
//     boost::property_tree::json_parser::read_json(symTopMoleculeFile_,mol);
//     boost::property_tree::ptree pol;
//     pol  = mol.get_child(molecule_ + ".Polarizability");
//     aaa_ = pol.get<double>("XX", 65.37282978);
//     abb_ = pol.get<double>("YY", 52.12184222);
//     acc_ = pol.get<double>("ZZ", 52.09178890);

//     boost::property_tree::ptree rotConsts;
//     rotConsts = mol.get_child(molecule_ + ".RotationalConstants");
//     Ae_       = rotConsts.get<double>("Ae", 2.36291695e-07);
//     Be_       = rotConsts.get<double>("Be", 2.32510239e-07);
//     Ce_       = rotConsts.get<double>("Ce", 2.32605387e-07);
//   }
//   else
//   {
//     boost::property_tree::json_parser::read_json(linearMoleculeFile_,mol);
//     amass1_ = 1836.0*mol.get<double>(molecule_ + ".Mass1",14.0);
//     amass2_ = 1836.0*mol.get<double>(molecule_ + ".Mass2",14.0);
//     bond_   = (1.0/0.52918)*mol.get<double>(molecule_ + ".BondLength",1.0);

//     boost::property_tree::ptree pol;
//     pol    = mol.get_child(molecule_ + ".Polarizability");
//     apara_ = pol.get<double>("ZZ",12.49002644);
//     aperp_ = pol.get<double>("XX",7.06531958);

// // Try to import asymptotic coefficients if they are provided
//     try
//     {
//       zion_ = mol.get<double>(molecule_ + ".Zion",1.0);
//       ionizationEnergy_ = mol.get<double>(molecule_ + ".IonizationEnergy",0.0);
//       asympCoeffs_ = vector<double>();
//       for (auto& iter : mol.get_child(molecule_ + ".AsympCoeffs"))
//       {
//         assert(iter.first.empty());
//         asympCoeffs_.push_back(iter.second.get_value<double>());
//       }
//     }
//     catch (exception& e) {}

//   }
// }