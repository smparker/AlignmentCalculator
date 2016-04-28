#include <cmath>
#include "pulses.hpp"

pulse::pulse(double I, double s, double t) : peakIntensity_(I), sigma_(s), t0_(t)
{}

double pulse::evaluate(double t)
{
  return peakIntensity_*exp(-1.0*pow(t - t0_,2)/(2.0*pow(sigma_,2)));
}