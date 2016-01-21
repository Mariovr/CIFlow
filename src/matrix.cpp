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
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cmath>
#include <vector>
#include <assert.h>
#include "scpp_assert.h"
#include "matrix.h"

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) {                 \
    if(status < 0)                                  \
    std::cerr << __FILE__ << ":" << __LINE__ <<     \
    ": Problem with writing to file. Status code="  \
    << status << std::endl;                         \
}

void vector_copy(int n, double* X, double* Y) {
    /*
        Copies a vector X of lenght n to vector Y
    */
    int i;
    for (i=0; i<n; i++) Y[i]=X[i];
}


/**
 * Empty matrix constructor. Don't use it
 * unless you know what you're doing
 */
matrix::matrix()
{
    this->n = -1;
    this->m = 0;
}

/**
 * @param n_ number of rows
 * @param m_ number of columns
 */
matrix::matrix(int n_, int m_)
{
    assert(n_ && m_);
    this->n = n_;
    this->m = m_;
    mat.reset(new double [n*m]());
}

/**
 * @param orig matrix to copy
 */
matrix::matrix(const matrix &orig)
{
    n = orig.n;
    m = orig.m;
    mat.reset(new double [n*m]());
    std::memcpy(mat.get(), orig.getpointer(), n*m*sizeof(double));
}

/**
 * move constructor
 * @param orig matrix to copy (destructive)
 */
matrix::matrix(matrix &&orig)
{
    n = orig.n;
    m = orig.m;
    mat = std::move(orig.mat);
}

matrix::matrix(std::vector<double>  vec , int _n , int _m)
{
   mat.reset(new double [vec.size()]);
   n = _n;
   m = _m ;
   for(int i = 0 ; i < n*m ; i++)
   {
       mat[i] = vec[i];
   }
}

//Creates DIIS system.
matrix::matrix(std::vector<double>  vec , int _n )
{
   mat.reset(new double [(_n+1) * (_n+1)]() );
   n = _n+1;
   m = _n+1 ;
   for(int i = 0 ; i < _n ; i++)
   {
       for(int j = i ; j < _n ; j++)
       {
           (*this)(i,j) = vec[(int) j*(j+1)/2. + i];
           if( i != j)
               (*this)(j,i) =(*this)(i,j);
       }
   }
   for( int i = 0 ; i< _n;i ++ )
   {
        (*this)(_n , i ) = -1.;
        (*this)(i , _n ) = -1.;
   }
   (*this)(_n , _n ) = 0; //Actually not necessary because full matrix is already initialized to zero.
}

matrix& matrix::operator=(const matrix &orig)
{
    n = orig.n;
    m = orig.m;
    mat.reset(new double [n*m]());
    std::memcpy(mat.get(), orig.getpointer(), n*m*sizeof(double));
    return *this;
}

/**
 * Set all matrix elements equal to a value
 * @param val the value to use
 */
matrix& matrix::operator=(double val)
{
    for(int i=0;i<n*m;i++)
        mat[i] = val;

    return *this;
}

/**
 * @return number of rows
 */
int matrix::getn() const
{
    return n;
}

/**
 * @return number of columns
 */
int matrix::getm() const
{
    return m;
}

double matrix::operator()(int x,int y) const
{
    assert(x<n && y<m);
    return mat[x+y*n];
}

double& matrix::operator()(int x,int y)
{
    //SCPP_ASSERT(x< n && y < m , "Boundary check error: x " << x << " should be smaller than " << "n " << n << " and y " << y << " should be smaller than m" << m );
    SCPP_TEST_ASSERT(x< n && y < m , "Boundary check error: x " << x << " should be smaller than n " << n<< " and y " << y << " should be smaller than m " << m );
    return mat[x+y*n];
}

double& matrix::operator[](int x)
{
    assert(x<n*m);
    return mat[x];
}

double matrix::operator[](int x) const
{
    assert(x<n*m);
    return mat[x];
}

void matrix::set(int i,int j,double val)
{
    mat[i+j*n] = val;
}

double matrix::get(int i,int j)
{
    return mat[i+j*n];
}

std::vector<double> matrix::get_column( int index)
{
    std::vector<double> col(getn());
    for(int i = 0; i < getn() ; i++)
        col[i] = (*this)(i,index);

    return col;
}

std::vector<double> matrix::get_row( int index)
{
    std::vector<double> row(getm()) ;
    for(int i = 0; i < getm() ; i++)
        row[i] = (*this)(index,i);

    return row;
}

double* matrix::getpointer() const
{
    return mat.get();
}

void matrix::add(matrix & a, double mult)
{
   for(int i = 0 ; i < n * m ; i ++) 
       mat[i] += mult * a[i];
}

void matrix::subtract(matrix & a)
{
   for(int i = 0 ; i < n * m ; i ++) 
       mat[i] -= a[i];
}

void matrix::multiply(double mul)
{
   for(int i = 0 ; i < n * m ; i ++) 
       mat[i] *= mul;
}

void matrix::scale_column(int col , double val)
{
   for(int i = 0 ; i < n  ; i ++) 
       (*this)(i,col) *= val;
}

double matrix::rms()
{
    double som = 0;
    for(int i = 0 ; i < n * m ; i ++) 
        som +=  mat[i] * mat[i];
    som /= n*m;
    return sqrt(som);
}

/*
 * Transforms matrix a to space defined by u (the new basis should be saved in the columns of u) 
 * */
void matrix::transform(const matrix & a ,const  matrix & u)
{
   matrix temp = matrix(a.getn() , a.getm()); //copy of a and we need extra space.
   temp.prod(a , u, 'N' , 'N');
   prod(u , temp, 'T' , 'N');
}


/*
 * Returns the trace of the product if tran == true: trace(This * a^T), else: trace(This * a)
 * REMARK: matrices should have same dimensions.
 */
double matrix::matrix_trace(const matrix & a , bool tran) const
{
    matrix temp = matrix(this->getm() , a.getn());
    if( tran)
        temp.prod(*this , a , 'N' , 'T');
    else
        temp.prod(*this , a , 'N' , 'N');
    return temp.trace();
}


/**
 * Matrix-Matrix product of A and B. Store result in this
 * @param A first matrix
 * @param B second matrix
 */
matrix& matrix::prod(matrix const &A, matrix const &B, char transa , char transb)
{
    //char trans = 'N';

    double alpha = 1.0;
    double beta = 0.0;

    assert(A.n == n && B.m == m);

    dgemm_(&transa,&transb,&A.n,&B.m,&A.m,&alpha,A.mat.get(),&A.n,B.mat.get(),&B.n,&beta,mat.get(),&A.n);

    return *this;
}

/**
 * Matrix-Vector product of this matrix and x. Store result in y. =>
 * y = A*x + alpha*y
 * We store the vector also in the matrix class (number of columns = 1)
 * @param x first vector
 * @param y second vector
 * @param alpha the constant factor to add y to y (optional, defaults to zero)
 */
void matrix::mvprod(matrix const &x, matrix &y, double alpha) const
{
    assert(m == x.n && n == y.n && "Dimension of vectors does not match!");

//    double beta = 1;
//    int incx = 1;
//    char uplo = 'U';

//    dsymv_(&uplo,&n,&beta,mat.get(),&n,x.getpointer(),&incx,&alpha,y.getpointer(),&incx);
    mvprod(x.getpointer(), y.getpointer(), alpha);
}

/**
 * Matrix-Vector product of this matrix and x. Store result in y. =>
 * y = A*x + alpha*y
 * Pass the vectors as double array's.
 * @param x first vector
 * @param y second vector
 * @param alpha the constant factor to add y to y (optional, defaults to zero)
 */
void matrix::mvprod(double const *x, double *y, double alpha) const
{
    assert(n == m && "Currently only works for square matrices");

    double beta = 1;
    int incx = 1;
    char uplo = 'U';

    dsymv_(&uplo,&n,&beta,mat.get(),&n,x,&incx,&alpha,y,&incx);
}

/**
 * Do a SVD on this matrix and store left singular values in
 * this. Changes the size of the matrix!
 * @return list of singular values
 */
std::unique_ptr<double []> matrix::svd()
{
    char jobu = 'A';
    char jobvt = 'N';

    int count_sing = std::min(n,m);

    std::unique_ptr<double []> sing_vals(new double[count_sing]);

    // MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
    int lwork = std::max( 3*count_sing + std::max(n,m), 5*count_sing);
    std::unique_ptr<double []> work(new double[lwork]);

    std::unique_ptr<double []> vt(new double[n*n]);

    int info;

    dgesvd_(&jobu,&jobvt,&n,&m,mat.get(),&n,sing_vals.get(),vt.get(),&n,0,&m,work.get(),&lwork,&info);

    if(info)
        std::cerr << "svd failed. info = " << info << std::endl;

    // overwrite the matrix with the right singular vectors
    m = n;
    mat = std::move(vt);

    return sing_vals;
}

void matrix::Print() const
{
        std::cout <<"#matrix: with dim( " <<getn() <<" , " << getm() << " )"  << std::endl;
        for(int i=0;i<n;i++){
            for(int j=0;j<m;j++)
                std::cout <<"    " <<std::setprecision(10) << (*this)(i,j) << "    "; 
            std::cout <<std::endl;
        }
}

double matrix::trace() const
{
    double result = 0;

    for(int i=0;i<std::min(m,n);i++)
        result += (*this)(i,i);

    return result;
}



void matrix::diagonalize(matrix &eigs , matrix &eigvec)
{
    
    assert(n==m && "Only square matrices for the moment");
    assert(eigs.n == n && "eigs has incorrect size, needs to be: (n,1)");

    char jobz;
    jobz = 'V'; //'N' eigenvalues only, 'V' eigenvectors also.

    char uplo = 'U'; //save in up or down triangular part.
    int lwork, liwork;

    lwork = 1 + 6*n + 2*pow(n,2);//minimal values of work arrays check lapack info.
    liwork = 3 + 5*n;

    double *work = new double[lwork];

    int *iwork = new int[liwork];

    // initialized to avoid uninitialized value errors
    int info = 1;

    int numelements = getn() * getm();
    int Bstep =1;
    dcopy_(&numelements,mat.get(),&Bstep,eigvec.getpointer(),&Bstep); //copy hamiltonian to eigenvectors because dsyev works in place and we don't want to loose the hamiltonian data
    //std::cout << eigvec(3,0)<< std::endl ;
    //std::cout << (*this)(0,0)<< std::endl ;
    dsyevd_(&jobz, &uplo, &n, eigvec.getpointer(), &n, eigs.getpointer(), work, &lwork, iwork, &liwork, &info);

    if(info != 0)
        std::cerr << "Calculating eigenvalues failed..." << std::endl;

    delete [] work;
    delete [] iwork;
}


void matrix::set_diagonal(matrix& diag)
{
    assert(diag.n == n);
    for(int i = 0 ; i < n ; i++)
    {
        (*this)(i,i) = diag(i,0);
    }
}

void matrix::transpose(matrix& trans)
{
    for(int i = 0 ; i < n ; i++)
    {
        for(int j = 0 ; j < m ; j++)
        {
            trans(j,i) = mat[i+j*n]; //mat is column major.
        }
    }
}

//Takes inproduct of a column of this matrix
double matrix::inproduct(int which)
{
    double som = 0;
    for(int i = 0 ; i < n ; i++)
    {
        som += (*this)(i,which) * (*this)(i,which);
    }
    return som;
}

void matrix::reset()
{
    std::fill(mat.get(), mat.get()+n*m, 0.);
}

int matrix::linsolve(double *B, double *x) {
    /*
        Solves the matrix equation A*x = B for x
        this is an n*n matrix in Fortran vector notation, B and x are vectors of lenght n
        A and B are not changed. x will contain the solution.
        x can equal B if B may be overwritten.
    */
    
    int info=0;
    
    /*  When x!=B we want to keep B unchanged. Copy to x and use this as working variable  */
    if (x!=B) vector_copy(n, B, x);
    
    /*  We need space for the pivot matrix  */
    int *ipiv;
    ipiv= (int *) malloc(n*sizeof(int));
    if (ipiv==NULL) {
        printf("malloc failed in linsolve\n"); 
        return 2;
    }
    
    int option = 1;//the number of right hand sides.
    /*  Now we can call the Lapack function  */
    dgesv_(&n, &option, mat.get() , &n, ipiv, x, &n, &info);
    
    /*  Cleanup and exit  */
    free(ipiv);
    return info;
}

matrix matrix::invert(int overwrite){
    /*  
        Calculates the inverse of the n*n matrix in this: Y = this^-1
        Does not change the value of this, unless overwrite = 1.
    */
    
    int info=0;
    double * result;
    
    matrix resultmat; 
    /*  When overwrite != 1 we want to keep this unchanged, so create new matrix object  */
    if (overwrite != 1)
    {
        resultmat = matrix(n,m);
        result = resultmat.getpointer();
    }
    else
    {
        result = getpointer();
    }
    
    /*  We need to store the pivot matrix obtained by the LU factorisation  */
    int *ipiv;
    ipiv=(int *) malloc(n*sizeof(int));
    if (ipiv==NULL) {
        printf("malloc failed in matrix_invert\n"); 
        //exit(EXIT_FAILURE);
    }
    
    /*  Turn Y into its LU form, store pivot matrix  */
    dgetrf_(&n, &n, result, &n, ipiv, &info);
    
    /*  Don't bother continuing when illegal argument (info<0) or singularity (info>0) occurs  */
    if (info!=0) printf("illegal element in invert\n"); 
        
    int workl = n*2;
    double work[workl];
    /*  Feed this to the lapack inversion routine.  */
    dgetri_(&n, result, &n, ipiv, work, &workl, &info);
    
    /*  Cleanup and exit  */
    free(ipiv);
    if (overwrite != 1)
    {
        return resultmat;
    }
    else
    {
        return *this;
    }
}

matrix::matrix(std::vector<matrix> symvec)
{
    int dim = 0;
    for(int i = 0 ; i < symvec.size() ; i++) dim += symvec[i].getn();
    this->n = dim;
    this->m = dim;
    mat.reset(new double [n*m]());

    int totn = 0;
    int totm = 0;
    for(int i = 0 ; i < symvec.size() ; i++) 
    {
        for(int l = 0 ; l < symvec[i].getm() ; l++) 
        {
            //We run first over columns because this is comp. more efficient, because stored in columnmajor.
            for(int y = 0 ; y < symvec[i].getn() ; y++) 
            {
                (*this)(y+totn,l + totm) = symvec[i](y,l);
            }
        }
        totm += symvec[i].getm();
        totn += symvec[i].getn();
    }
}

matrix matrix::matrix_square_root() {
    /*  
        This function calculates one of the square roots of the matrix in this and returns it:
        this needs to be a symmetric positive definite matrix of dimension n*n in Fortran vector
        format.         
        The variable I is a vector of length n containing +1 and -1 elements. It can be used 
        to select one of the 2^n different square roots of X
        
        This function first calculates the eigenvalue decomposition of X: X = U*D*U^T
        A new matrix F is then calculated with on the diagonal the square roots of D, with
        signs taken from I.
        The square root is then obtained by calculating U*F*U^T
    */
    matrix prod = matrix(n, m);
    matrix eigvec= matrix(n, m);
    matrix eigval= matrix(n, 1);
    (*this).diagonalize(eigval , eigvec);
    for(int i = 0 ; i < n ; i++)
    {
        if(eigval(i,0) > 1e-20)
        {
            prod(i,i) = sqrt(eigval(i,0));
        }
        else
        {
            std::cout << "Warning: the eigenvalues are to low in sqrt of matrix. "<<std::endl;
        }
    }
    matrix prod2 = matrix(n, m);
    prod2.prod( eigvec, prod );
    prod.prod( prod2, eigvec , 'N' , 'T');
    return prod;
}

matrix matrix::matrix_inv_square_root() {
    /*  
        This function calculates one of the square roots of the matrix in this and returns it:
        this needs to be a symmetric positive definite matrix of dimension n*n in Fortran vector
        format.         
        The variable I is a vector of length n containing +1 and -1 elements. It can be used 
        to select one of the 2^n different square roots of X
        
        This function first calculates the eigenvalue decomposition of X: X = U*D*U^T
        A new matrix F is then calculated with on the diagonal the square roots of D, with
        signs taken from I.
        The square root is then obtained by calculating U*F*U^T
    */
    matrix prod = matrix(n, m);
    matrix eigvec= matrix(n, m);
    matrix eigval= matrix(n, 1);
    (*this).diagonalize(eigval , eigvec);
    for(int i = 0 ; i < n ; i++)
    {
        if(eigval(i,0) > 1e-20)
        {
            prod(i,i) = 1./sqrt(eigval(i,0));
        }
        else
        {
            std::cout << "Warning: the eigenvalues are to low in sqrt of matrix. "<<std::endl;
        }
    }
    matrix prod2 = matrix(n, m);
    prod2.prod( eigvec, prod );
    prod.prod( prod2, eigvec , 'N' , 'T');
    return prod;
}

Vector2d::Vector2d(size_type num_rows, size_type num_cols) : _rows(num_rows), _cols(num_cols) ,_data(num_elem())
{
    SCPP_TEST_ASSERT(num_rows > 0, "Number of rows in a matrix must be positive: " << num_rows);
    SCPP_TEST_ASSERT(num_cols > 0, "Number of columns in a matrix must be positive: " << num_cols);
}

Vector2d::Vector2d(size_type num_rows, size_type num_cols, const double & init_value) : _rows(num_rows), _cols(num_cols), _data(num_elem())
{   
    SCPP_TEST_ASSERT(num_rows > 0, "Number of rows in a matrix must be positive: " << num_rows);
    SCPP_TEST_ASSERT(num_cols > 0, "Number of columns in a matrix must be positive: " << num_cols); 
}

Vector2d::Vector2d():_rows(0 ), _cols(0) , _data(0){}

Vector2d::Vector2d(const Vector2d & vec):_rows(vec.num_rows() ), _cols(vec.num_cols() )  ,_data(num_elem() )
{
        for(int j = 0 ; j < _cols ; j ++)
            for(int i = 0 ; i < _rows ; i ++)
                (*this)(i,j) = vec(i,j);
}


/*  
unsigned Vector2d::index(size_type row, size_type col) const 
{
	SCPP_TEST_ASSERT(row < _rows, "Row " << row  << " must be less than " << _rows);
 	SCPP_TEST_ASSERT(col < _cols, "Column " << col  << " must be less than " << _cols);
	return _rows * col + row; //Column major format.
}

unsigned Vector2d::num_elem() const
{
    return num_rows() * num_cols();
}*/

unsigned Vector2d::index(size_type row, size_type col) const 
{
    SCPP_TEST_ASSERT(row < this->num_rows(), "Row " << row  << " must be less than " << this->num_rows());
    SCPP_TEST_ASSERT(col < this->num_cols(), "Column " << col  << " must be less than " << this->num_cols());
    if( row > col)
    {
        double temp = col;
        col = row;
        row = temp;
    }
    return this->num_cols() * row + col - row*(row+1)/2; //Column major format.
}

unsigned Vector2d::num_elem() const
{
    return this->num_cols()* (this->num_cols()+1)/2;
}

inline std::ostream& operator << (std::ostream& os, const Vector2d & m) 
{
    for( unsigned r =0; r<m.num_rows(); ++r ) {
	for( unsigned c=0; c<m.num_cols(); ++c ) {
            os << m(r,c);
	    if( c + 1 < m.num_cols() )
	        os << "\t";
        }
        os << "\n";
    }
    return os;
} 


/* vim: set ts=8 sw=4 tw=0 expandtab :*/
