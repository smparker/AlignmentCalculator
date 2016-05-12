/**
 * \file "matrix.hpp"
 * \author J. Szekely
 */

#ifndef QLALIB_MATRIX
#define QLALIB_MATRIX

#include <memory>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <assert.h>
#include <random>
#include <complex>
#include <vector>
typedef std::complex<double> cplx;

/**
 * @brief Base matrix storage class
 * @details Matrix base class containing routines applicable to all matrices, such as
 * submatrix extraction, filling the matrix, and calculating the trace. Templated to
 * accept many data types
 *
 */

class matrixComp; ///< Need forward declaration for matrixReal*MatrixComp operation

template <typename T> class matrixBase
{
protected:
    static unsigned int memSize; ///< Total memory allocated to all matrix objects
    size_t nrows, ncols; ///< Number of rows and columns in matrix
    std::unique_ptr<T[]> vals; ///< Data containment array

public:
  /// Constructor
  matrixBase(const int nr, const int nc) : nrows(nr), ncols(nc), vals(std::unique_ptr<T[]>(new T[nc*nr]))
  {
    zero();
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }

  /// Copy Constructor
  matrixBase(const matrixBase& o) : nrows(o.nrows), ncols(o.ncols), vals(std::unique_ptr<T[]>(new T[nrows*ncols]))
  {
    std::copy_n(o.vals.get(), nrows*ncols, vals.get());
    memSize += sizeof(T)*size(); //number of bytes allocated to instance of class
  }

  /// Move Constructor
  matrixBase(matrixBase&& o) : nrows(o.nrows), ncols(o.ncols), vals(std::move(o.vals)) { o.nrows = 0; o.ncols = 0; };

  /// Destructor
  ~matrixBase(){memSize -= sizeof(T)*size();} //number of bytes allocated to instance of class

  ///  Returns the size of the matrix
  size_t size() const { return nrows * ncols; }

  /// Returns pointer to the beginning of the data array
  T* data() { return vals.get(); }
  const T* data() const { return vals.get(); }

  /// Fill entire matrix with zeroes
  void zero()
  {
    std::fill_n(vals.get(), nrows*ncols, T(0.0));
  }

  size_t nr() const {return nrows;} ///< Return number of rows
  size_t nc() const {return ncols;} ///< Return number of columns

  /// Fill matrix with randomly generated numbers
 void random()
 {
   std::random_device rd;
   std::mt19937 gen(rd());
   std::uniform_real_distribution<double> dis(0, 1);
   std::generate_n(vals.get(), nrows*ncols, [&dis, &gen](){return T(dis(gen));});
 }

  /// Get or set array element at (row,col) position
  T& element(const int row, const int col)
  {
    return vals[col*nrows + row];
  }

  /// Const version of previous function
  const T& element(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

  /// Operator definition calling on the element function
  T& operator()(const int row, const int col)
  {
    return element(row,col);
  }

  /// Const version of previous function
  const T& operator()(const int row, const int col) const
  {
    return vals[col*nrows + row];
  }

  /// Set Diagonal elements to unity
  void makeIdentity()
  {
    zero();
    for (int ii = 0; ii < std::min(ncols,nrows); ii++)
      vals[ii+ii*nrows] = T(1.0);
  }

  /// Compute the trace
  T trace()
  {
    T sum = T(0.0);
    for (int ii = 0; ii < std::min(int(ncols),int(nrows)); ii++)
      sum += element(ii,ii);
    return sum;
  }

  // Apply a scaling factor to all elements
  void scale(const T a)
  {
    std::for_each(data(), data()+size(), [&a](T& p){p*=a;});
  }

  /// Print memory usage for all matrices
  void printMem() const
  {
    std::cout << "Current memory allocated to this matrix: " << size()*sizeof(T) << " bytes." << std::endl;
    std::cout << "Total memory allocated for matrix storage: " << memSize << " bytes." << std::endl;
    return;
  }

  /// Extract a portion of the matrix starting at element(r,c), get (nr x nc matrix)
  template <class U>
  std::shared_ptr<U> getSub_impl(int r, int c, int nr, int nc) const
  {
    assert(r + nr <= nrows && c + nc <= ncols);
    auto out = std::make_shared<U>(nr,nc);
    for (int jj = 0; jj < nc; jj++)
      std::copy_n(&element(r,c+jj),nr,&out->element(0,jj));
    return out;
 }

  /// Place matrix o at position (r,c)
  template <typename U>
  void setSub(int r, int c, U o)
  {
    assert(r + o.nrows <= nrows && c + o.ncols <= ncols);
    for (int jj = 0; jj < o.ncols; jj++)
      std::copy_n(&o(0,jj), o.nrows, &element(r,c+jj));
    return;
  }

  /// Overload of the << operator to print out a portion of the matrix
  template <typename U> friend std::ostream &operator<<(std::ostream &out, const matrixBase <U> &o);
};

/**
 * @brief Matrix of real numbers
 * @details Matrix base class using doubles, added mathematical funcitonality
 *
 */
class matrixReal : public matrixBase<double>
{
public:

  matrixReal(const int nr, const int nc); ///< Default constructor
  matrixReal(const matrixReal&);          ///< Copy constructor
  matrixReal(matrixReal&&);               ///< Move constructor

/** @name Matrix-Matrix Operations
 *  Binary Operations accepting two real matricies
 */

///@{
  matrixReal& operator=(const matrixReal&);       ///< \f$ A = B \f$
  matrixReal operator*(const matrixReal&) const;  ///< \f$ A * B \f$
  matrixReal& operator*=(const matrixReal&);      ///< \f$ A = (A*B) \f$
  matrixReal operator+(const matrixReal&) const;  ///< \f$ A + B \f$
  matrixReal& operator+=(const matrixReal&);      ///< \f$ A = (A+B) \f$
  matrixReal operator-(const matrixReal&) const;  ///< \f$ A - B \f$
  matrixReal& operator-=(const matrixReal&);      ///< \f$ A = (A-B) \f$
  matrixReal operator|(const matrixReal&) const;  ///< \f$ A^T * B \f$
  matrixReal operator^(const matrixReal&) const;  ///< \f$ A * B^T \f$

///@}

/** @name Matrix(Real)-Matrix(Complex) Operations
 *  Binary Operations accepting two real matricies
 */
///@{
  matrixComp operator*(const matrixComp&) const;  ///< \f$ A * B \f$
///@}

/** @name Scalar-Matrix Operations
 *  Binary Operations accepting a matrix and a constant, only rhs operators at the moment
 */
///@{

  matrixReal operator*(const double&) const; ///< \f$ cA \f$
  matrixReal operator/(const double&) const; ///< \f$ \frac{1}{c}A \f$
  matrixReal& operator*=(const double&);     ///< \f$ A = cA \f$
  matrixReal& operator/=(const double&);     ///< \f$ A = \frac{1}{c}A \f$

///@}

/** @name BLAS and LAPACK
 *  Other routines that require BLAS and LAPACK libraries to function, assumes symmetric matrices
 */
///@{
  void diagonalize(double* eigVals); ///< Full Diagonalization with dsyev
  void diagonalize_alt(double* eigVals); ///< Alternate diagonalization wuth dsyevd
  void diagonalize(double* eigVals, bool getLowEigVal, int keepNum, double abstol); ///< Partial diagonaliation returning lowest several eigenvectors, uses dsyevr
  std::shared_ptr<matrixReal> transpose() const; ///< \f$ A^T \f$
  std::tuple<std::shared_ptr<matrixReal>, std::shared_ptr<matrixReal>>svd(std::vector<double>&); ///< Single Value Decomposition for non square matrices

///@}

/** @name Other Operations
 *
 */
///@{
  double dot_product(const matrixReal& o) const; ///< \f$ c = A \cdot B \f$
  double norm() const; ///< \f$ |A| = \sqrt{A \cdot A} \f$
  double rms() const; ///< \f$ \frac{|A|}{\sqrt{N}} \f$
  double variance() const; ///<\f$ \frac{|A|^2}{N} \f$
  double operator%(const matrixReal& o) const; ///< Dot product for vectors with one column

  matrixReal kron(matrixReal &o) const; ///< \f$ A \otimes A\f$, very inefficient

  /// Returns pointer to submatrix
  std::shared_ptr<matrixReal> getSub(int ii, int jj, int kk, int ll) const
  {
    return getSub_impl<matrixReal>(ii,jj,kk,ll);
  }

  void ax_plus_y(const double a, matrixReal &o); ///< \f$ cA+B \f$
  void invert();

///@}
};


/// Overload the << operator to print a matrix
template <typename T>
std::ostream &operator<<(std::ostream &out, const matrixBase <T> &o)
{
  for (int row = 0; row < std::min(10,int(o.nrows)); row++)
  {
    for (int col = 0; col < std::min(10,int(o.ncols)); col++)
    {
      out << std::setprecision(3) << o(row,col) << "\t";
    }
    out << "\n";
  }
  out << std::endl;
  return out;
};

template <typename T> unsigned int matrixBase<T>::memSize = 0;


/**
 * @brief Matrix of complex numbers
 * @details Matrix base class using std::complex<double>, added mathematical funcitonality
 *
 */
class matrixComp : public matrixBase<cplx>
{
public:
  matrixComp(const int nr, const int nc);
  matrixComp(const matrixComp&);
  matrixComp(matrixComp&&);

//Matrix-Matrix operations
  matrixComp& operator=(const matrixComp&);
  matrixComp operator*(const matrixComp&) const;
  matrixComp& operator*=(const matrixComp&);
  matrixComp operator+(const matrixComp&) const;
  matrixComp& operator+=(const matrixComp&);
  matrixComp operator-(const matrixComp&) const;
  matrixComp& operator-=(const matrixComp&);
  // matrixComp operator|(const matrixComp&) const;
  // matrixComp operator^(const matrixComp&) const;
  void getEigvals(double* eigVals);

  friend matrixComp matrixReal::operator*(const matrixComp& o) const;

//Scalar-Matrix Operations
  matrixComp operator*(const cplx&) const;
  matrixComp operator/(const cplx&) const;
  matrixComp& operator*=(const cplx&);
  matrixComp& operator/=(const cplx&);

  std::shared_ptr<matrixComp> getSub(int ii, int jj, int kk, int ll) const
  {
    return getSub_impl<matrixComp>(ii,jj,kk,ll);
  }
  void diagonalize(double* eigVals); ///< Full Diagonalization with dsyev
  void invert();
};

void printMatrix(matrixComp &o, std::string filename, double *x = nullptr, double *y = nullptr);
void printMatrix(matrixReal &o, std::string filename, double *x = nullptr, double *y = nullptr);


#endif
