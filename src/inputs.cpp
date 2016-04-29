#include "inputs.hpp"
#include <boost/filesystem.hpp>
#include <algorithm>

inputParameters::inputParameters(std::string fn) : filename_(fn)
{
  stripComments_();
  parseAllInputs_();
}

void inputParameters::stripComments_()
{
//  Open input and output file
  if (!boost::filesystem::exists("output_data/"))
    boost::filesystem::create_directory("output_data");

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

}

void inputParameters::parseFieldInfo(boost::property_tree::ptree &IP)
{

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