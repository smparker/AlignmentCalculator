#include "nonadiabaticPropagator.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
  fstream inputfile;
  if (argc < 2)
  {
    cout << "Please provide the name of the json inputfile" << endl;
    return(1);
  }
  inputParameters inputs(argv[1]);
  if (inputs.jobtype_ == JOBTYPE::NONADIABATIC)
  {
    nonadiabaticPropagator calculation(inputs);
    calculation.run();
  }
  // else if (inputs.jobtype_ == JOBTYPE::ADIABATIC)
  //   adiabaticPropagator calculation(inputs);

  return(0);
}