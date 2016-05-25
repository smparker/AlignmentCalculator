/**
 * \file "outputs.hpp"
 * \author J. Szekely
 */
#ifndef ALIGNMENTCALCULATOR_OUTPUTS
#define ALIGNMENTCALCULATOR_OUTPUTS

#include "molecules.hpp"
#include <fstream>
#include <gsl/gsl_sf_legendre.h> // Spherical harmonic functions

/**
 * @brief      Base class for observables that are obtained from a density matrix as Tr(O*rho). Other outputs are obtained from other members of the propagator class, such as the wavefunction density or the list of basis states used in the calculation.
 */
class observable
{
public:
  std::string id_tag_; ///< Unique identifier tag for printing to first line of output file
  std::shared_ptr<matrices> operator_matrix_; ///< Matrix representation of observable
  observable(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians); ///< Constructor
  virtual void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians) = 0; ///< Initialize the observable matrix
  virtual double density_evaluate_(std::shared_ptr<matrices> densities_); ///< Calculate expectation value given a density matrix
  virtual double wvfxn_evaluate_(std::shared_ptr<matrices> densities_, std::shared_ptr<arrays> populations_); ///< Calculate expectation value given a set of wavefunctions stored in a matrix
};


/// \f$  \langle \cos^2 \theta_{3D} \rangle = \langle \cos^2 (\hat z \cdot \hat Z ) \rangle \f$
class obsCosTheta3D : public observable
{
public:
  obsCosTheta3D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle \cos^2 \theta_{2D} \rangle \f$
class obsCosTheta2D : public observable
{
public:
  obsCosTheta2D(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle \hat H \rangle \f$
class obsEnergy : public observable
{
public:
  obsEnergy(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle \cos^2 \theta^\prime_{3D} \rangle = \langle \cos^2 (\hat x \cdot \hat Z ) \rangle \f$
class obsCosThetaAlt : public observable
{
public:
  obsCosThetaAlt(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle J^2 \rangle \f$
class obsJ : public observable
{
public:
  obsJ(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle K^2 \rangle \f$
class obsK : public observable
{
public:
  obsK(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};

/// \f$  \langle M \rangle \f$
class obsM : public observable
{
public:
  obsM(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
  void initialize_(std::shared_ptr<basisSubsets> basisSets,std::shared_ptr<matrices> fieldFreeHamiltonians);
};


/**
 * @brief cos^2 theta (2D) matrix elements
 * @details Computes the cos^2 2D matrix elements between two spherical harmonic functions
 *
 * @param J angular momentum of state 1
 * @param M angular momentum projection of state 1
 * @param j angular momentum of state 2
 * @param m angular momentum projection of state 2
 * @return overlap
 */
inline double cos2D(int J, int M, int j, int m)
{
  int Nx = 300;
  int Ny = 300;
  cplx sum = 0.0;
  for (double theta = 0; theta < M_PI; theta += M_PI/Nx)
    for (double phi = 0.0; phi < 2.0*M_PI; phi += 2.0*M_PI/Ny)
    {
      double c = cos(theta);
      sum += sin(theta)*std::conj(gsl_sf_legendre_sphPlm(J,abs(M),c))*gsl_sf_legendre_sphPlm(j,abs(m),c)/(1.0+pow(tan(theta)*cos(phi),2.0));
    }
    return real(sum*(M_PI/Nx)*(2.0*M_PI/Ny));
}

#endif