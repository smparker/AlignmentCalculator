/**
 * \file "utilities.hpp"
 * \author  J. Szekely
 */

#ifndef QLALIB_UTILITIES
#define QLALIB_UTILITIES

#include <complex>
#include <memory>

//BLAS
extern "C"
{
  void dzgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
    const std::complex<double>* alpha, const double* a, const int* lda, const std::complex<double>* b, const int* ldb,
    const std::complex<double>* beta, std::complex<double>* c, const int* ldc);

  void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
    const double* alpha, const double* a, const int* lda, const double* b, const int* ldb,
    const double* beta, double* c, const int* ldc);

  void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);

  void dsyevd_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda, double *w, double *work, int *lwork, int *iwork, int *liwork, int* info);

  void zheev_(const char*, const char*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*, const int*, double*, int*);

  double ddot_(const int*, const double*, const int*, const double*, const int*);

#ifndef ZDOT_RETURN
 void zdotc_(std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*);
#else
 std::complex<double> zdotc_(const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*);
#endif

  void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);

  void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*);

  void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
               const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
               const std::complex<double>* beta, std::complex<double>* c, const int* ldc);

  void zhemm_(const char* side, const char* uplo, const int *m, const int *n, const std::complex<double> *alpha, const std::complex<double> *a, const int *lda, const std::complex<double> *b, const int *ldb, const std::complex<double> *beta, std::complex<double> *c, const int *ldc);

  void zgemv_(const char*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*,
             const std::complex<double>*, std::complex<double>*, const int*);
  int izamax_(const int*, const std::complex<double>*, const int*);

  int izamin_(const int*, const std::complex<double>*, const int*);

  int idamax_(const int*, const double *, const int*);

  int idamin_(const int*, const double *, const int*);

 void mkl_ddnscsr_(const int*, const int*, const int*, const double *, const int*, const double *, const int *, const int *, int *);

 void mkl_domatcopy_(const char*, const char*, const int *, const int *, const double *, const double* , const int *, double* , const int *);

 void mkl_zomatcopy_(const char*, const char*, const int *, const int *, const std::complex<double> *, const std::complex<double>* , const int *, std::complex<double>* , const int *);

 void mkl_dcsrgemv_(const char*, const int *, const double *, const int *, const int *, const double *, const double *);

 void dgetrf_(const int*, const int *, double* , int* , int* , int*); // LU decomoposition of a general matrix

 void dgetri_(const int*, double*, int*, int*, double*, int*, int*); // generate inverse of a matrix given its LU decomposition

 void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info);

 void dswap_(const int* n, double* x, const int* incx, double* y, const int* incy);

 void zswap_(const int* n, std::complex<double>* x, const int* incx, std::complex<double>* y, const int* incy);

 void zgetrf_(const int*, const int *, std::complex<double>* , int* , int* , int*);

 void zgetri_(const int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);

}

//LAPACK
extern "C"
{
  void dgesvd_(const char*, const char*, const int*, const int*, double*, const int*, double*, double*, const int*, double*, const int*,  double*, const int*, int*);

  void dsyevr_(const char*, const char*, const char*, const int*, double*, const int*, const double*, const double*, const int*, const int*, const double*,  int*, double*, double*, const int*, int*, double*, int*, int*, int*, int*);
}

//AlignmentTool interface
namespace
{
  void dgesv_(const int n, const int nrhs, double* a, const int lda, int* ipiv, double* b, const int ldb, int& info)
             { ::dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info); }
  void dgesv_(const int n, const int nrhs, std::unique_ptr<double[]>& a, const int lda, std::unique_ptr<int[]>& ipiv,
             std::unique_ptr<double[]>& b, const int ldb, int& info) { ::dgesv_(&n, &nrhs, a.get(), &lda, ipiv.get(), b.get(), &ldb, &info); }

  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const double* a, const int lda, const double* b, const int ldb,
              const double beta, double* c, const int ldc)
              { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

  void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb,
              const double beta, std::unique_ptr<double []>& c, const int ldc)
             { ::dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }


  void dzgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const std::complex<double> alpha, const double* a, const int lda, const std::complex<double>* b, const int ldb,
              const std::complex<double> beta, std::complex<double>* c, const int ldc)
              { ::dzgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

  void dzgemm_(const char* transa, const char* transb, const int m, const int n, const int k,
              const std::complex<double> alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<std::complex<double> []>& b, const int ldb,
              const std::complex<double> beta, std::unique_ptr<std::complex<double> []>& c, const int ldc)
             { ::dzgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

  void dsyev_(const char* a, const char* b, const int c, double* d, const int e, double* f, double* g, const int h, int& i)
             { ::dsyev_(a,b,&c,d,&e,f,g,&h,&i);}

  void dsyev_(const char* a, const char* b, const int c, std::unique_ptr<double []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<double []>& g, const int h, int& i)
             { ::dsyev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,&i);}

  void dsyevd_(const char* jobz, const char* uplo, const int n, double *a, const int lda, double *w, double *work, int lwork, int *iwork, int liwork, int info)
            { ::dsyevd_(jobz, uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info);}


  void dgesvd_(const char* a, const char* b, const int c, const int d, double* e, const int f, double* g, double* h, const int i, double* j, const int k,
              double* l, const int m, int& n) { ::dgesvd_(a,b,&c,&d,e,&f,g,h,&i,j,&k,l,&m,&n); }


  void mkl_domatcopy_(const char* ordering, const char* trans, const int r, const int c, const double alpha,
                    const double* A, const int nr, double* B, const int nc)
                    {::mkl_domatcopy_(ordering,trans,&r,&c,&alpha,A,&nr,B,&nc);}

  void mkl_zomatcopy_(const char* ordering, const char* trans, const int r, const int c, const std::complex<double> alpha,
                    const std::complex<double>* A, const int nr, std::complex<double>* B, const int nc)
                    {::mkl_zomatcopy_(ordering,trans,&r,&c,&alpha,A,&nr,B,&nc);}

  double ddot_(const int a, const double* b, const int c, const double* d, const int e) { return ::ddot_(&a,b,&c,d,&e);}

  double ddot_(const int a, const std::unique_ptr<double []>& b, const int c, const std::unique_ptr<double []>& d, const int e)
             { return ::ddot_(&a,b.get(),&c,d.get(),&e); }

#ifndef ZDOT_RETURN
 std::complex<double> zdotc_(const int b, const std::complex<double>* c, const int d, const std::complex<double>* e, const int f) {
   std::complex<double> a;
   ::zdotc_(&a,&b,c,&d,e,&f);
   return a;
 }
 std::complex<double> zdotc_(const int b, const std::unique_ptr<std::complex<double> []>& c, const int d, const std::unique_ptr<std::complex<double> []>& e, const int f) {
   std::complex<double> a;
   ::zdotc_(&a,&b,c.get(),&d,e.get(),&f);
   return a;
 }
#else
 std::complex<double> zdotc_(const int a, const std::complex<double>* b, const int c, const std::complex<double>* d, const int e) { return ::zdotc_(&a,b,&c,d,&e); }
 std::complex<double> zdotc_(const int a, const std::unique_ptr<std::complex<double> []>& b, const int c, const std::unique_ptr<std::complex<double> []>& d, const int e)
                             { return ::zdotc_(&a,b.get(),&c,d.get(),&e); }
#endif

  void dsyevr_(const char* jobz, const char* range, const char* uplo, const int n, double* a, const int lda, double vl,
                 double vu, const int il, const int iu, const double abstol, int m, double* w, double* z, const int ldz,
                int* isuppz, double* work, int lwork, int* iwork, int liwork, int& info)
                { ::dsyevr_(jobz, range, uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);}

  void daxpy_(const int a, const double b, const double* c, const int d, double* e, const int f) { ::daxpy_(&a,&b,c,&d,e,&f); }

  void zaxpy_(const int a, const std::complex<double> b, const std::complex<double>* c, const int d, std::complex<double>* e, const int f) { ::zaxpy_(&a,&b,c,&d,e,&f); }

  void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
                 const std::complex<double> alpha, const std::complex<double>* a, const int lda, const std::complex<double>* b, const int ldb,
                 const std::complex<double> beta, std::complex<double>* c, const int ldc) { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }

  void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
                 const std::complex<double> alpha, const std::unique_ptr<std::complex<double>[]>& a, const int lda, const std::unique_ptr<std::complex<double>[]>& b, const int ldb,
                 const std::complex<double> beta, std::unique_ptr<std::complex<double>[]>& c, const int ldc)
                 { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }

  void zhemm_(const char* side, const char* uplo, const int m, const int n, const std::complex<double> alpha, const std::complex<double> *a, const int lda, const std::complex<double> *b, const int ldb, const std::complex<double> beta, std::complex<double> *c, const int ldc)
                { ::zhemm_(side, uplo, &m, &n, &alpha, a, &lda, b, &ldb, &beta, c, &ldc); };

  void zgemv_(const char* a, const int b, const int c, const std::complex<double> d, const std::complex<double>* e, const int f, const std::complex<double>* g, const int h,
             const std::complex<double> i, std::complex<double>* j, const int k) { ::zgemv_(a,&b,&c,&d,e,&f,g,&h,&i,j,&k); }

  void zgemv_(const char* a, const int b, const int c, const std::complex<double> d, const std::unique_ptr<std::complex<double> []>& e, const int f,
             const std::unique_ptr<std::complex<double> []>& g, const int h, const std::complex<double> i, std::unique_ptr<std::complex<double> []>& j, const int k)
             { ::zgemv_(a,&b,&c,&d,e.get(),&f,g.get(),&h,&i,j.get(),&k); }

  int izamax_(const int n, const std::complex<double> *x, const int inc)
              { return ::izamax_(&n,x,&inc);}

  int izamin_(const int n, const std::complex<double> *x, const int inc)
              { return ::izamin_(&n,x,&inc);}

  int idamax_(const int n, const double *x, const int inc)
              { return ::idamax_(&n,x,&inc);}

  int idamin_(const int n, const double *x, const int inc)
              { return ::idamin_(&n,x,&inc);}

  void zheev_(const char* a, const char* b, const int c, std::complex<double>* d, const int e, double* f, std::complex<double>* g, const int h, double* i, int& j)
             { ::zheev_(a,b,&c,d,&e,f,g,&h,i,&j); }

  void zheev_(const char* a, const char* b, const int c, std::unique_ptr<std::complex<double> []>& d, const int e,
             std::unique_ptr<double []>& f, std::unique_ptr<std::complex<double> []>& g, const int h, std::unique_ptr<double[]>& i, int& j)
             { ::zheev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,i.get(),&j); }

 void mkl_ddnscsr_(const int* job, const int nr, const int nc, const double *data, const int nz, const double *csrdata, const int *cCdata, const int *cRdata, int info)
            { ::mkl_ddnscsr_(job, &nr, &nc, data, &nz, csrdata, cCdata, cRdata, &info); };

 void mkl_dcsrgemv_(const char* transa, const int nr, const double *data, const int *csrRow, const int *csrCol, const double *o, const double *out)
            { ::mkl_dcsrgemv_(transa,&nr,data,csrRow,csrCol,o,out);}

  void dgetrf_(const int a, const int b, double *c, int d, int *e, int f)
            {::dgetrf_(&a, &b, c, &d, e, &f);}

  void dgetri_(const int a, double *b, int c, int* d, double *e, int f, int g) // generate inverse of a matrix given its LU decomposition
            {::dgetri_(&a, b, &c, d, e, &f, &g);}

  void zgetrf_(const int a, const int b, std::complex<double>* c, int d, int *e, int f)
            {::zgetrf_(&a, &b, c, &d, e, &f);}

  void zgetri_(const int a, std::complex<double> *b, int c, int* d, std::complex<double> *e, int f, int g)
            {::zgetri_(&a, b, &c, d, e, &f, &g);}

  void dswap_(const int n, double* x, const int incx, double* y, const int incy)
            {::dswap_(&n,x,&incx,y,&incy);}

  void zswap_(const int n, std::complex<double>* x, const int incx, std::complex<double>* y, const int incy)
            {::zswap_(&n,x,&incx,y,&incy);}


}

#endif