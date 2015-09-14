// CIFlow is a very flexible configuration interaction program
// Copyright (C) Ghent University 2014-2015
//
// This file is part of CIFlow.
//
// CIFlow is developed by Mario Van Raemdonck <mario.vanraemdonck@ugent.be>
// a member of the Ghent Quantum Chemistry Group (Ghent University).
// See also : http://www.quantum.ugent.be
//
// At this moment CIFlow is not yet distributed.
// However this might change in the future in the hope that
// it will be useful to someone.
//
// For now you have to ask the main author for permission.
//
//--
#ifndef MATRIX_H
#define MATRIX_H

#include <memory>
#include <complex>
#include <vector>

// all kind of BLAS/LAPACK functions
#ifndef __INTEL_COMPILER__
extern "C" {
    void dcopy_(int *n,double *x,int *incx,double *y,int *incy);
    void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* iwork, int* liwork, int* info);
    double ddot_(int *n,double *x,int *incx,double *y,int *incy);
    void dscal_(int *n,double *alpha,double *x,int *incx);
    void dsymv_(char *uplo, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy);
    void daxpy_(int *n,const double *alpha,double *x,int *incx,double *y,int *incy);
    void dstev_( const char* jobz, const int* n, double* d, double* e, double* z, const int* ldz, double* work, int* info );
    void dgesvd_( char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );
    void dgemm_(char *transA,char *transB,const int *m,const int *n,const int *k,double *alpha,double *A,const int *lda,double *B,const int *ldb,double *beta,double *C,const int *ldc);

    void zheevd_( const char* jobz, const char* uplo, const int* n, std::complex<double>* a, const int* lda, double* w, std::complex<double>* work, const int* lwork, double* rwork,
            const int* lrwork, int* iwork, const int* liwork, int* info );

    void dsaupd_(int *ido, char *bmat, int *n, char *which,
            int *nev, double *tol, double *resid, int *ncv,
            double *v, int *ldv, int *iparam, int *ipntr,
            double *workd, double *workl, int *lworkl, int *info);

    void dseupd_(int *rvec, char *All, int *select, double *d,
            double *z, int *ldz, double *sigma,
            char *bmat, int *n, char *which, int *nev,
            double *tol, double *resid, int *ncv, double *v,
            int *ldv, int *iparam, int *ipntr, double *workd,
            double *workl, int *lworkl, int *info);
   void dsyev_(char *jobz,char *uplo,int *n,double *A,int *lda,double *W,double *work,int *lwork,int *info);
   void dpotrf_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dpotri_(char *uplo,int *n,double *A,int *lda,int *INFO);
   void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* IWORK, int* INFO);
   void dgeev_(char*jobvl,char*jobvr,int*n,double*A,int*lda,double* wr,double* wi,double * vl,int *ldvl,double *vr,int*ldvr,double*work,int*lwork,int*info);
   double dlamch_(char * ch);
   void dlasrt_(char* id, int* n, double* vec, int* info);
   double dlansy_(char * norm, char * uplo, int * dimR, double * mx, int * lda, double * work);
   double dlange_(char * norm, int * m, int * n, double * mx, int * lda, double * work);
   void dgesv_( int*  colmaj, int* num , double * A, int* nn, int* ipiv, double * x, int* nnn, int* d);
   void dgetrf_( int * m, int * n, double* a, int* lda, int * ipiv, int * info );
   void dgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );
}
#else
#include <mkl.h>
extern "C" { 
    void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info); 
    void dseupd_(int *rvec, char *All, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info); 
    void dgemm_(char *transA,char *transB,const int *m,const int *n,const int *k,double *alpha,double *A,const int *lda,double *B,const int *ldb,double *beta,double *C,const int *ldc);
    void dcopy_(int *n,double *x,int *incx,double *y,int *incy);
    void dsymv_(char *uplo, const int *n, const double *alpha, const double *a, const int *lda, const double *x, const int *incx, const double *beta, double *y, const int *incy); 
    double ddot_(int *n,double *x,int *incx,double *y,int *incy);
    //void dgesv_( int*  num, int* rhs , double * A, int* nn, int* ipiv, double * x, int* nnn, int* d);
    }
#endif


/**
 * Helper class, wrapper around a double array. Has methods to get the number of rows
 * and columns. A simple matrix class.
 */
class matrix
{
    public:
        matrix();

        matrix(int n_, int m_);

        matrix(const matrix &orig);

        matrix(matrix &&orig);

        //Make sure the vector vec contains the columns of the matrix in column major format.
        matrix(std::vector<double> vec, int _n , int _m);

        //Creates symmetric matrix from vector (which contains a diis system.)
        matrix(std::vector<double> vec, int _n );

        //Creates total matrix from blockdiagonal elements (from symmetry)
        matrix(std::vector<matrix> symvec);

        virtual ~matrix() { }

        matrix& operator=(const matrix &orig);

        matrix& operator=(double val);

        int getn() const;

        int getm() const;

        double operator()(int x,int y) const;

        double& operator()(int x,int y);

        void set(int i,int j,double val);

        double get(int i,int j);

        double& operator[](int x);

        double operator[](int x) const;

        double* getpointer() const;

        void add(matrix & a, double mult = 1);

        void subtract(matrix & a);

        void scale_column(int col, double val);

        double inproduct(int which);// returns inproduct of a particular(which) column of this matrix.

        double rms();

        // Transforms matrix a to space defined by u (the new basis should be saved in the columns of u) 
        void transform(matrix const & a , matrix const & u);

        /*The blas and lapack routines are in column major memory order, this matrix class
         * uses those routines frequently therefore it saves the memory also in column major memory order in a contiguous array. As opposed to the standard C order which is row major, therefore we don't have to transpose our memory layout before providing it to the blas routines. So we put the trans keywords standard to 'N'. But we still keep it as a keyword in case you do want to transpose the matrix.*/
        matrix& prod(matrix const &A, matrix const &B, char transa = 'N', char transb = 'N');

        void mvprod(matrix const &x, matrix &y, double alpha = 0) const;

        void mvprod(double const *x, double *y, double alpha = 0) const;

        std::unique_ptr<double []> svd();

        void Print() const;

        double trace() const;

        double matrix_trace(const matrix & a, bool tran) const;

        void diagonalize(matrix &eigs , matrix& eigvec);

        void transpose(matrix& transpose);

        void set_diagonal(matrix& diag);

        int linsolve(double * B, double *x);

        //Sets this matrix to zero.
        void reset();

        matrix invert(int overwrite);

        matrix matrix_square_root();
        matrix matrix_inv_square_root();

    private:
        //!n by m array of double
        std::unique_ptr<double []> mat;
        //! number of rows
        int n;
        //! number of columns
        int m;
};

/**
 * Helper class, wrapper around a complex double array. Has methods to get the number of rows
 * and columns. A simple matrix class.
 */
class cmatrix
{
    public:
        cmatrix();

        cmatrix(int n_, int m_);

        cmatrix(const cmatrix &orig);

        cmatrix(cmatrix &&orig);

        virtual ~cmatrix() { }

        cmatrix& operator=(const cmatrix &orig);

        cmatrix& operator=(std::complex<double> val);

        int getn() const;

        int getm() const;

        std::complex<double> operator()(int x,int y) const;

        std::complex<double>& operator()(int x,int y);

        std::complex<double>& operator[](int x);

        std::complex<double> operator[](int x) const;

        std::complex<double>* getpointer() const;

        matrix& prod(cmatrix const &A, cmatrix const &B);

        void Print() const;

    private:
        //!n by m array of complex<double>
        std::unique_ptr<std::complex<double> []> mat;
        //! number of rows
        int n;
        //! number of columns
        int m;
};

#endif /* MATRIX_H */

/* vim: set ts=8 sw=4 tw=0 expandtab :*/
