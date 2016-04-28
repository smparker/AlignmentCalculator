/**
 * \file "pulses.hpp"
 * \author  J. Szekely
 */

#ifndef ALIGNMENTCALCULATOR_PULSES
#define ALIGNMENTCALCULATOR_PULSES

class pulse {
	double peakIntensity_;
	double sigma_;
	double t0_;
	pulse (double, double, double);
	double evaluate(double t);
};

#endif
