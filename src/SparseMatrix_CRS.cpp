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
#include <cmath>
#include <cassert>
#ifdef __OMP__
#include <omp.h>
#endif
#ifdef __INTEL_COMPILER__
#include <mkl.h>
#endif
#include "MyHDF5.h"
#include "SparseMatrix_CRS.h"


/**
 * Construct SparseMatrix_CRS object for n x m matrix
 * @param n the number of rows
 * @param m the number of columns
 */
SparseMatrix_CRS::SparseMatrix_CRS(unsigned int n, unsigned int m)
{
    this->n = n;
    this->m = m;
    data.resize(0);
    col.resize(0);
    row.resize(0);
    data.reserve(n*2);
    col.reserve(n*2);
    row.reserve(n);
    _zero = 1e-14;
}

SparseMatrix_CRS::SparseMatrix_CRS()
{
    this->n = 0;
    this->m = 0;
    data.resize(0);
    col.resize(0);
    row.resize(0);
    _zero = 1e-14;
}

/**
 * Copy constructor
 * @param origin the object to copy
 */
SparseMatrix_CRS::SparseMatrix_CRS(const SparseMatrix_CRS &origin)
{
    this->n = origin.n;
    this->m = origin.m;
    this->data = origin.data;
    this->col = origin.col;
    this->row = origin.row;
}

/**
 * constructor out of CCS object
 * @param ccs the CCS object to use
 */
/*  
SparseMatrix_CRS::SparseMatrix_CRS(const SparseMatrix_CCS &ccs)
{
    this->n = ccs.getn();
    this->m = ccs.getm();

    for(unsigned int i=0;i<n;i++)
    {
       NewRow();

       for(unsigned int k=0;k<ccs.col.size()-1;k++)
       {
          for(unsigned int l=ccs.col[k];l<ccs.col[k+1];l++)
             if( ccs.row[l] == i )
                PushToRow(k,ccs.data[l]);
             else if ( ccs.row[l] > i )
                break;
       }
    }

    NewRow();
}
*/

SparseMatrix_CRS::~SparseMatrix_CRS()
{

}

/**
 * Read only access operator
 * @param i the row number
 * @param j the column number
 * @return the matrix element
 */
double SparseMatrix_CRS::operator()(unsigned int i,unsigned int j) const
{
    for(unsigned int k=row[i];k<row[i+1];k++)
       if( col[k] == j )
          return data[k];
    return 0;
}

/**
 * @return the number of rows
 */
int SparseMatrix_CRS::getn() const
{
    return n;
}

/**
 * @return the number of columns
 */
int SparseMatrix_CRS::getm() const
{
    return m;
}

/**
 * Copy operator
 * @param origin the object to copy
 * @return this object
 */
SparseMatrix_CRS &SparseMatrix_CRS::operator=(const SparseMatrix_CRS &origin)
{
    this->n = origin.n;
    this->m = origin.m;
    this->data = origin.data;
    this->col = origin.col;
    this->row = origin.row;

    return *this;
}

/**
 * Convert a dense matrix to CRS format
 * @param dense the matrix to convert
 */
void SparseMatrix_CRS::ConvertFromMatrix(const matrix &dense)
{
   this->n = dense.getn();
   this->m = dense.getm();
   row.resize(n+1);

   data.clear();
   col.clear();

   row[0] = 0;

   for(unsigned int i=0;i<n;i++)
   {
      for(unsigned int j=0;j<m;j++)
         if( fabs(dense(i,j)) > _zero )
         {
            data.push_back(dense(i,j));
            col.push_back(j);
         }

      row[i+1] = col.size();
   }

   row.back() = col.size();
}

/**
 * Convert this CRS matrix to a dense matrix (only works for square matrices)
 * @param dense the matrix to fill
 */
void SparseMatrix_CRS::ConvertToMatrix(matrix &dense) const
{
   dense = 0;
   unsigned int dn = dense.getn();

   if(dn != n || n != m)
      return;

   for(unsigned int i=0;i<row.size() -1 ;i++)
      for(unsigned int k=row[i];k<row[i+1];k++)
         dense(i,col[k]) = data[k];
}

/**
 * Print the raw CRS data to stdout
 */
void SparseMatrix_CRS::PrintRaw() const
{
    std::cout << "Data(" << data.size() << "):" << std::endl;
    for(unsigned int i=0;i<data.size();i++)
        std::cout << data[i] << " ";
    std::cout << std::endl;

    std::cout << "Col indices:" << std::endl;
    for(unsigned int i=0;i<col.size();i++)
        std::cout << col[i] << " ";
    std::cout << std::endl;

    std::cout << "Row indices:" << std::endl;
    for(unsigned int i=0;i<row.size();i++)
        std::cout << row[i] << " ";
    std::cout << std::endl;
}

/**
 * Print sparse matrix to output
 * @param output the ostream to print to
 * @param matrix_p the matrix to print
 * @return the filled ostream (with the matrix)
 */
ostream &operator<<(ostream &output,SparseMatrix_CRS &matrix_p)
{
   for(unsigned int i=0;i<matrix_p.row.size()-1;i++)
      for(unsigned int k=matrix_p.row[i];k<matrix_p.row[i+1];k++)
         output << i << "\t" << matrix_p.col[k] << "\t" << matrix_p.data[k] << std::endl;

   return output;
}

/**
 * Adds a new column element to the current row.
 * To use this, first call NewRow() to start a row and then
 * use PushToRow() to add elements to that row. Always end
 * with calling NewRow() again.
 * @param j column
 * @param value the matrix element value
 */
void SparseMatrix_CRS::PushToRow(unsigned int j, double value)
{
    if( fabs(value ) > _zero)
    {
        if(col.empty() || row.back() == col.size() || col.back() < j)
        {
            data.push_back(value);
            col.push_back(j);
        }
        else
        {
            unsigned int begin = row.back();
            for(unsigned int i=begin;i<col.size();i++)
            {
               if( col[i] > j )
               {
                   col.insert(col.begin() + i,j);
                   data.insert(data.begin() + i,value);
                   break;
               } else if (col[i] == j)
               {
                   data[i] += value;
                   if(fabs(data[i])< _zero)//Check if we didn't cancel the previous value so this, has become zero and needs to be removed.
                   {
                      data.erase(data.begin() + i);
                      col.erase(col.begin() + i);
                   }
                   break;
               }
            }
        }
    }
}

void SparseMatrix_CRS::PushToRownoadd(unsigned int j, double value)
{
    if( fabs(value ) > _zero)
    {
        if(col.empty() || row.back() == col.size() || col.back() < j)
        {
           data.push_back(value);
           col.push_back(j);
        }
        else
        {
           unsigned int begin = row.back();
           for(unsigned int i=begin;i<col.size();i++)
           {
              if( col[i] > j )
              {
                  col.insert(col.begin() + i,j);
                  data.insert(data.begin() + i,value);
                  break;
              } else if (col[i] == j)
              {
                 data[i] = value;

                 break;
              }
           }
        }
    }
}

/**
 * Resets the data, row and col vectors.
 */
void SparseMatrix_CRS::reset_vecs(){
    data.resize(0);
    col.resize(0);
    row.resize(0);
}

/**
 * Adds the next row to the sparsematrix
 */
void SparseMatrix_CRS::NewRow()
{
   if(row.size() == (n+1))
      return;

   row.push_back(data.size());
}

/**
 * Do the matrix vector product y = A * x
 * @param xmat a m component vector
 * @param ymat a n component vector
 */
void SparseMatrix_CRS::mvprod(const matrix &xmat, matrix &ymat) const
{
   double *x = xmat.getpointer();
   double *y = ymat.getpointer();

   mvprod(x,y);
}

/**
 * Do the matrix vector product y = A * x
 * @param xmat a m component vector
 * @param ymat a n component vector
 */
void SparseMatrix_CRS::mvprod(const double *x, double *y) const
{
    #ifdef __INTEL_COMPILER__
    char uplo = 'U';
    mkl_cspblas_dcsrsymv(&uplo, ( int *) &n, (double *) data.data(), (int *) row.data(), (int *) col.data(), (double *) x, y);
    #else

        #pragma omp parallel for 
        for(unsigned int i=0;i<n;i++)
        {
            y[i] = 0;
            for(unsigned int k=row[i];k<row[i+1];k++)
                y[i] += data[k] * x[col[k]];
        }
    #endif    
}


void SparseMatrix_CRS::arpackDiagonalize(int nve , matrix & eigval , matrix &  eigvec) const
{
    // dimension of the matrix
    int n = this->n;

    // number of eigenvalues to calculate
    int nev = nve;

    // reverse communication parameter, must be zero on first iteration
    int ido = 0;
    // standard eigenvalue problem A*x=lambda*x
    char bmat = 'I';
    // calculate the smallest algebraic eigenvalue
    char which[] = {'S','A'};
    // calculate until machine precision
    double tol = 0;

    // the residual vector
    double *resid = new double[n];

    // the number of columns in v: the number of lanczos vector
    // generated at each iteration, ncv <= n
    // We use the answer to life, the universe and everything, if possible
    int ncv = 42;

    if( n < ncv )
        ncv = n;

    // v containts the lanczos basis vectors
    int ldv = n;
    double *v = new double[ldv*ncv];

    int *iparam = new int[11];
    iparam[0] = 1; // Specifies the shift strategy (1->exact)
    iparam[2] = 3*n; // Maximum number of iterations
    iparam[6] = 1; /* Sets the mode of dsaupd.
                      1 is exact shifting,
                      2 is user-supplied shifts,
                      3 is shift-invert mode,
                      4 is buckling mode,
                      5 is Cayley mode. */

    int *ipntr = new int[11]; /* Indicates the locations in the work array workd
                                 where the input and output vectors in the
                                 callback routine are located. */

    // array used for reverse communication
    double *workd = new double[3*n];
    for(int i=0;i<3*n;i++)
        workd[i] = 0;

    int lworkl = ncv*(ncv+8); /* Length of the workl array */
    double *workl = new double[lworkl];

    // info = 0: random start vector is used
    int info = 0; /* Passes convergence information out of the iteration
                     routine. */

    // rvec == 0 : calculate only eigenvalue
    // rvec > 0 : calculate eigenvalue and eigenvector
    int rvec = 1;

    // how many eigenvectors to calculate: 'A' => nev eigenvectors
    char howmny = 'A';

    int *select;
    // when howmny == 'A', this is used as workspace to reorder the eigenvectors
    if( howmny == 'A' )
        select = new int[ncv];

    // This vector will return the eigenvalues from the second routine, dseupd.
// double *d = new double[nev];

 //   double *z = 0;

  //  if( rvec )
   //     z = new double[n*nev];
    // not used if iparam[6] == 1
    double sigma;

    // first iteration
    dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    while( ido != 99 )
    {
        // matrix-vector multiplication
        mvprod(workd+ipntr[0]-1, workd+ipntr[1]-1);

        dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    }

    if( info < 0 )
        std::cerr << "Error with dsaupd, info = " << info << std::endl;
    else if ( info == 1 )
        std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
    else if ( info == 3 )
        std::cerr << "No shifts could be applied during implicit Arnoldi update, try increasing NCV." << std::endl;

    dseupd_(&rvec, &howmny, select, eigval.getpointer(), eigvec.getpointer(), &ldv, &sigma, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    if ( info != 0 )
        std::cerr << "Error with dseupd, info = " << info << std::endl;

    // use something to store the result before deleting...
    //sigma = d[0];

    delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    //delete [] d;

    //if( rvec )
        //delete [] z;

    if( howmny == 'A' )
        delete [] select;

    //return sigma; // lowest eigenvalue
}


int SparseMatrix_CRS::write_to_file(const char *filename, const char *name, bool append) const
{
    hid_t file_id, group_id, dataset_id, attribute_id, dataspace_id, scalar_id;
    herr_t status;

    if(append)
       file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    else
       file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    group_id = H5Gcreate(file_id, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dimblock = data.size();

    scalar_id = H5Screate(H5S_SCALAR);

    dataspace_id = H5Screate_simple(1, &dimblock, NULL);

    dataset_id = H5Dcreate(group_id, "data", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data() );
    HDF5_STATUS_CHECK(status);

    unsigned int size = data.size();
    attribute_id = H5Acreate (dataset_id, "size", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_STATUS_CHECK(attribute_id);
    status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &size );
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    dataset_id = H5Dcreate(group_id, "col", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, col.data() );
    HDF5_STATUS_CHECK(status);

    size = col.size();
    attribute_id = H5Acreate (dataset_id, "size", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite (attribute_id, H5T_NATIVE_UINT, &size );
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    dimblock = row.size();

    dataspace_id = H5Screate_simple(1, &dimblock, NULL);

    dataset_id = H5Dcreate(group_id, "row", H5T_STD_U64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, row.data() );
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(dataspace_id);
    HDF5_STATUS_CHECK(status);

    dataset_id = H5Dcreate(group_id, "n", H5T_STD_U64LE, scalar_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n );
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Sclose(scalar_id);
    HDF5_STATUS_CHECK(status);

    status = H5Gclose(group_id);
    HDF5_STATUS_CHECK(status);

    status = H5Fclose(file_id);
    HDF5_STATUS_CHECK(status);

    return 0;
}

/**
* Read a SparseMatrix_CRS from a HDF5 file
* @param filename the name of the file to read from
* @param name the name of the group in the HDF5 file
*/
int SparseMatrix_CRS::read_from_file(const char *filename, const char *name)
{
    hid_t file_id, group_id, dataset_id, attribute_id;
    herr_t status;

    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    HDF5_STATUS_CHECK(file_id);

    group_id = H5Gopen(file_id, name, H5P_DEFAULT);
    HDF5_STATUS_CHECK(group_id);

    dataset_id = H5Dopen(group_id, "n", H5P_DEFAULT);
    HDF5_STATUS_CHECK(dataset_id);

    status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &n);
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    row.resize(n+1);

    unsigned int size;

    dataset_id = H5Dopen(group_id, "data", H5P_DEFAULT);
    HDF5_STATUS_CHECK(dataset_id);

    attribute_id = H5Aopen(dataset_id, "size", H5P_DEFAULT);
    HDF5_STATUS_CHECK(attribute_id);

    status = H5Aread(attribute_id, H5T_NATIVE_UINT, &size);
    HDF5_STATUS_CHECK(status);


    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    data.resize(size);

    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);


    dataset_id = H5Dopen(group_id, "col", H5P_DEFAULT);
    HDF5_STATUS_CHECK(dataset_id);

    attribute_id = H5Aopen(dataset_id, "size", H5P_DEFAULT);
    HDF5_STATUS_CHECK(attribute_id);

    status = H5Aread(attribute_id, H5T_NATIVE_UINT, &size);
    HDF5_STATUS_CHECK(status);

    status = H5Aclose(attribute_id);
    HDF5_STATUS_CHECK(status);

    col.resize(size);

    status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, col.data());
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);


    dataset_id = H5Dopen(group_id, "row", H5P_DEFAULT);
    HDF5_STATUS_CHECK(dataset_id);

    status = H5Dread(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, row.data());
    HDF5_STATUS_CHECK(status);

    status = H5Dclose(dataset_id);
    HDF5_STATUS_CHECK(status);

    status = H5Gclose(group_id);
    HDF5_STATUS_CHECK(status);

    status = H5Fclose(file_id);
    HDF5_STATUS_CHECK(status);

    return 0;
}

void SparseMatrix_CRS::add_from_list(std::vector<std::unique_ptr<SparseMatrix_CRS>> &list){
    unsigned int total_size = 0;

    for(auto& smat : list)
        total_size += smat->data.size();
 
    data.reserve(total_size);
    col.reserve(total_size);
 
    data.clear();
    col.clear();
 
    row.reserve(n+1);
    row.clear();
 
    for(auto &smat: list)
    {
        auto prev_start = data.size();
 
        data.insert(data.end(), smat->data.begin(), smat->data.end());
        col.insert(col.end(), smat->col.begin(), smat->col.end());
 
        for(unsigned int i=0;i<smat->row.size();++i)
           row.push_back(prev_start + smat->row[i]);
 
        smat.reset();
    }
 
    assert(data.size() == total_size);
    assert(col.size() == total_size);
    assert(row.size() == n) ;
 
    row.push_back(data.size());
};

//----------------------------class SparseMatrix_CRS_Sym -------------------------------------------------------------------------------

void SparseMatrix_CRS_Sym::mvprod(const double *x, double *y) const{
    #pragma omp parallel for 
        for(unsigned int i=0;i<n;i++){
            y[i] = 0;
        }
    #ifdef __INTEL_COMPILER__
        char uplo = 'U';
        mkl_cspblas_dcsrsymv(&uplo, (int * ) &n, (double *) data.data(), (int *) row.data(), (int *) col.data(), (double * )x, y);
    #else
            //MPI::Init();
            //int size = MPI::COMM_WORLD.Get_size();
            //int rank = MPI::COMM_WORLD.Get_rank();
            for(unsigned int i = 0 ; i<n ; i++){
                for(unsigned int k=row[i];k<row[i+1];k++){
                    y[i] += data[k] * x[col[k]];
                    //std::cout << "i" << 0 << " " << y[0] << std::endl;
                    //std::cout << data[0] ;
                    if( i != col[k])
                            y[col[k]] += data[k] * x[i];
                }
            }
            //MPI::Finalize();
    #endif    
}

double SparseMatrix_CRS_Sym::operator()(unsigned int i,unsigned int j) const
{
    if( i > j){
        unsigned int k = i;
        i = j;
        j = k;
    }
    for(unsigned int k=row[i];k<row[i+1];k++)
        if( col[k] == j )
            return data[k];
    return 0;
}

/**
 * Convert this CRS_Sym matrix to a dense matrix (only works for square matrices)
 * @param dense the matrix to fill
 */
void SparseMatrix_CRS_Sym::ConvertToMatrix(matrix &dense) const
{
   dense = 0;
   unsigned int dn = dense.getn();

   if(dn != n || n != m)
      return;

   for(unsigned int i=0;i<row.size()-1;i++)
      for(unsigned int k=row[i];k<row[i+1];k++){
         dense(i,col[k]) = data[k];
         if( i != col[k]){
             dense(col[k],i) = data[k];
            }
      }
}

/**
 * Convert a dense matrix to CRS format
 * @param dense the matrix to convert
 */
void SparseMatrix_CRS_Sym::ConvertFromMatrix(const matrix &dense)
{
   this->n = dense.getn();
   this->m = dense.getm();
   row.resize(n+1);

   data.clear();
   col.clear();

   row[0] = 0;

   for(unsigned int i=0;i<n;i++)
   {
      for(unsigned int j=i;j<m;j++)
         if( fabs(dense(i,j)) >  _zero )
         {
            data.push_back(dense(i,j));
            col.push_back(j);
         }

      row[i+1] = col.size();
   }

   row.back() = col.size();
}
/* vim: set ts=4 sw=4 expandtab :*/
