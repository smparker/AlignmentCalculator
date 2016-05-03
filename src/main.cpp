#include <iostream>
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
  nonadiabaticPropagator calculation(inputs);


  return(0);
}