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
#ifndef SPARSEMATRIX_CRS_H
#define SPARSEMATRIX_CRS_H

#include <iostream>
#include <vector>

#include "matrix.h"

using std::ostream;

/**
 * @author Ward Poelmans, Mario Van Raemdonck
 * This is a class written for sparse n x m matrices. It uses the CRS format to store
 * a matrix. This type of sparse matrix format is very suitable for an efficient matrix vector product.
 */

class SparseMatrix_CRS
{
   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << matrix << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << matrix << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param matrix_p the SparseMatrix_CRS you want to print
    */
   friend ostream &operator<<(ostream &output,SparseMatrix_CRS &matrix_p);

   //friend class SparseMatrix_CCS;

   public:

      //constructor
      SparseMatrix_CRS(unsigned int n, unsigned int m);
      SparseMatrix_CRS();

      //copy constructor
      SparseMatrix_CRS(const SparseMatrix_CRS &);

      //SparseMatrix_CRS(const SparseMatrix_CCS &);

      //destructor
      virtual ~SparseMatrix_CRS();

      //overload equality operator
      SparseMatrix_CRS &operator=(const SparseMatrix_CRS &);

      //easy to access the numbers
      double operator()(unsigned int i,unsigned int j) const;

      int getn() const;

      int getm() const;

      void ConvertFromMatrix(const matrix &dense);

      void ConvertToMatrix(matrix &dense) const;

      void PrintRaw() const;

      void PushToRow(unsigned int j, double value);

      int datasize(){return data.size();}

      void NewRow();

      void mvprod(const matrix &, matrix &) const;

      virtual void mvprod(const double *x, double *y) const;

      void arpackDiagonalize(int nve , matrix & eigval , matrix & eigvec) const;

      void reset_vecs();

      int write_to_file(const char *filename , const char *name , bool append) const;

      int read_from_file(const char *filename , const char *name);

      void add_from_list(std::vector<std::unique_ptr <SparseMatrix_CRS> > & list);

   protected:

      //! Array that holds the non zero values
      std::vector<double> data;
      //! Array that holds the column indexes
      std::vector<unsigned int> col;
      //! Array that holds the row index of data
      std::vector<unsigned int> row;

      //!dimension of the matrix (number of rows)
      unsigned int n;
      //!dimension of the matrix (number of columns)
      unsigned int m;
};

class SparseMatrix_CRS_Sym:public SparseMatrix_CRS{
   public:
      //constructor
      SparseMatrix_CRS_Sym(unsigned int n , unsigned int m ):SparseMatrix_CRS(n, m){}
      SparseMatrix_CRS_Sym();

      void mvprod(const double *x, double *y) const;
      double operator()(unsigned int i,unsigned int j) const;
      void ConvertToMatrix(matrix &dense) const;
      void ConvertFromMatrix(const matrix &dense);
};

#endif /* SPARSEMATRIX_CRS_H */
/* vim: set ts=3 sw=3 expandtab :*/
