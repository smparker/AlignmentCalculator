#include "nonadiabaticPropagator.hpp"

nonadiabaticPropagator::nonadiabaticPropagator(inputParameters &IP) :
  propagatorBase(IP),
  t0_(0.0),
  tFinal_(IP.max_time_),
  noutputs_(IP.n_outputs_),
  pulses_(IP.pulses_)
{}

void nonadiabaticPropagator::run()
{
  // if population is small, skip basis Subset
}

void initializeCVODE()
{

}

void evalRHS()
{

}

void step()
{

}

void run()
{

}