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
#ifndef TWOINDEX_H
#define TWOINDEX_H

#include <string>
#include "Irreps.h"

/** TwoIndex class.
 \author Sebastian Wouters <sebastianwouters@gmail.com>
 \date February 8, 2013
 
 Container class for symmetric two-index tensors with Abelian point group symmetry (real character table; see Irreps.h): 1DMs and 1-particle matrix elements. */
class TwoIndex{

   public:
   
      //! Constructor
      /** \param nGroup The symmetry group number (see Irreps.h)
          \param IrrepSizes Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep */
      TwoIndex(const int nGroup, const int * IrrepSizes);

      //!Copy Constructor
      TwoIndex(const TwoIndex & tindex);
      
      //! Destructor
      virtual ~TwoIndex();
      
      //! Set an element
      /** \param irrep The irrep number (see Irreps.h)
          \param i The first index (within the symmetry block)
          \param j The second index (within the symmetry block)
          \param val The value to which the element of the matrix should be set */
      void set(const int irrep, const int i, const int j, const double val);

      //! Get an element
      /** \param irrep The irrep number (see Irreps.h)
          \param i The first index (within the symmetry block)
          \param j The second index (within the symmetry block)
          \return The value of the matrix element */
      double get(const int irrep, const int i, const int j) const;
      
      //! Save the TwoIndex object
      /** \param name filename */
      void save(const std::string name) const;
      
      //! Load the TwoIndex object
      /** \param name filename */
      void read(const std::string name);

      int get_group()const {return SymmInfo.getGroupNumber();}

      int * get_irrepsizes()const {return Isizes;}

      double * get_storage_irrep(int cnt )const { return storage[cnt ];}
   
   private:
   
      //Contains the group number, the number of irreps, and the multiplication table
      Irreps SymmInfo;
      
      //Array with length the number of irreps of the specified group, containing the number of orbitals of that irrep
      int * Isizes;
      
      //storage[I_i][i+j*(j+1)/2] = mx_ij = mx_ji
      double ** storage;

};
#endif
