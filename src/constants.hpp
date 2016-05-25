/**
 * \file "constants.hpp"
 * \author  J. Szekely
 */

#ifndef ALIGNMENT_CONSTANTS
#define ALIGNMENT_CONSTANTS

namespace CONSTANTS {
	const double HBAR       = 1.0; ///< hbar in atomic units
	const double EMASS      = 1.0; ///< mass of an electron in atomic units
	const double ECHARGE    = 1.0; ///< charge of an electron in atomic units
	const double VACPERM    = (1.0/(4.0*M_PI)); ///< Vacuum permitivity in atomic units
	const double LEN        = 0.0529177; ///< nanometers in an atomic unit of length
	const double VEL        = 2.18e8; ///< 1 atomic unit of velocity in cm/s
	const double EN         = 27.21; ///< 1 atomic unit of energy in eV
	const double TIME       = 2.42e-17; ///< 1 atomic unit of time in s
	const double AUperFS    = 41.34137333656137; ///< atomic units of time in 1 fs
	const double FREQ       = 4.13e16; ///< 1 atomic unit of frequency in Hz
	const double ELECFIELD  = 5.14e9; ///< 1 atomic unit of electric field amplitude in V/cm^2
	const double LASERINTEN = 3.51e16; ///< 1 atomic unit of intensity in W/cm^2 (0.5*vacuum_permitivity*speed_of_light*EField^2)
	const double C          = 2.998e10/VEL; ///< speed of light in atomic units
	const double BOLTZ      = 3.1669e-6; ///< Boltzmann constant in atomic units (energy) per Kelvin
}

#endif