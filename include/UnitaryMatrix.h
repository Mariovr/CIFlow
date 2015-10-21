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
#ifndef _UnitaryMatrix_H
#define _UnitaryMatrix_H

#include <memory>
#include <vector>
#include <iostream>

class matrix;
class OptIndex;

/** unitary class.
    The UnitaryMatrix class is a storage and manipulation class for the unitary matrix. This matrix is blockdiagonal in the irreducible representations. */
class UnitaryMatrix
{
    public:
      
        //! Constructor
        /** \param _hamindexIn */
        UnitaryMatrix(OptIndex& index);
        //! creates a UnitaryMatrix object from a matrix with dimension (L,L), where L is the number of orbitals. REMARK: if the unitary matrix behaves according to a particular point group symmetry. mat should be block diagonal with the nonzero irreps ranked in ascending order.
        UnitaryMatrix(matrix & mat, OptIndex & opt);

     	//!Copy constructor
        UnitaryMatrix(const UnitaryMatrix & unit);

        UnitaryMatrix(UnitaryMatrix && unit);

        UnitaryMatrix(OptIndex & index , std::istream & file);
        UnitaryMatrix(OptIndex & index , const std::string & file);

        //! Destructor
        virtual ~UnitaryMatrix() = default;

        UnitaryMatrix& operator=(const UnitaryMatrix &unit);

        UnitaryMatrix& operator=(UnitaryMatrix &&unit);

        //! Returns the inverse of this unitary matrix.
        UnitaryMatrix get_inverse() const;
        
        //! Get the number of variables in the x-parametrization of the unitary update
        /** \return The number of unique variables in the x-matrix */
        unsigned int getNumVariablesX() const;
        
        //! Get the unitary rotation for block irrep
        /** \param irrep The irreducible representation
            \return Pointer to the desired unitary block */
        double * getBlock(const int irrep) const;

        double get_element(int i , int j );

        //returns a vector of double that contains the unitary transformation, without irrep information, and in column major format (and new basis is in rows (same as unitary matrix))
        std::vector<double> get_full_transformation();

        //! Copy the x solution back (NR, augmented NR, ...)
        /** \param vector The x-solution */
        void copyXsolutionBack(double * vector);
        
        //! Update the unitary transformation based on the new vector and the previous unitary
        /** \param workmem1 Work memory
            \param workmem2 Work memory */
        void updateUnitary(double * workmem1, double * workmem2, double * vector , const bool multiply);

        //! Rotate the unitary matrix to a given basis (eg. NO) eigenbasis
        /** \param eigenvecs The eigenbasis
            \param work Work memory */
        void rotate_active_space_vectors(double * eigenvecs, double * work);
        
        //! Calculate the two-norm of U^T*U - I
        /** \param work Work memory */
        void CheckDeviationFromUnitary() const;
        
        //! Save the unitary to disk
        void saveU(const std::string savename) const;
        
        //! Load the unitary from disk
        void loadU(const std::string loadname);

        //! Load the unitary from a textfile
        void load_unitary(const std::string filename);

        void load_unitary(std::istream & file );
        
        //! Delete the stored unitary (on disk)
        void deleteStoredUnitary(std::string name) const;

        //! Implements a simple jacobi rotation of the i, and j th orbital.
        void jacobi_rotation(int irrep, int i, int j, double angle);

        //! Resets the unitary matrix
        void reset_unitary(bool allzero=0);

        void print_unitary(std::ostream & out) const;

        void build_skew_symm_x(const int irrep, double * xblock , const double * Xelem) const;

        int get_norb(int irrep);

        int get_nstart(int irrep);

        int get_Nirrep()const;

        //returns irrep from orbital
        int get_orbital_irrep(int orbital);

        //! U = U_2 * U
        void multiply_left_with(UnitaryMatrix & unit2) const;

        //! sets this unitary to a random unitary.
        void set_random_unitary();

        //! make skew symmetric matrix.
        void make_skew_symmetric();

        //!Checks difference with unit2
        double check_difference(UnitaryMatrix & unit2);

        //!add UnitaryMatrix unit2 (multiplied with mult) to this UnitaryMatrix
        void add_unit(UnitaryMatrix & unit2 , double mult);
        

        //!Returns _unitary *_unitary^T, each irrep in a separate matrix contained in the vector
        std::vector<matrix> calc_overlap();


    private:
      
        //Externally created and destroyed index handler
        std::unique_ptr<OptIndex> _index;
      
        //Number of variables in the x-matrix
        unsigned int x_linearlength;
        
        //The unitary matrix (e^x * previous unitary): unitary[irrep][row + size_irrep * col]
        std::vector< std::unique_ptr<double []> > unitary;
        
   };
#endif

/* vim: set ts=4 sw=4 expandtab :*/
