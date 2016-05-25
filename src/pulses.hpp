/**
 * \file "pulses.hpp"
 * \author  J. Szekely
 */

#ifndef ALIGNMENTCALCULATOR_PULSES
#define ALIGNMENTCALCULATOR_PULSES

class pulse
{
public:
	double peakIntensity_; ///< Maximum intensity of Gaussian pulse
	double sigma_; ///< Variance of Gaussian function, determined from FWHM
	double t0_; ///< Time of peak maximum
	pulse (double, double, double); ///< Constructor
	double evaluate(double t); ///< Evaluate the field intensity at time t
};

#endif
