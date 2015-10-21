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
#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include<vector>
#include <memory>
#include <string>
#include <iostream>

//Includes
#include "Irreps.h"

//forward declarations
class OptIndex;
class FourIndex;
class TwoIndex;
class UnitaryMatrix;

/** Hamiltonian class.
    Container class for the Hamiltonian matrix elements.
    
    \section ham_info Specific Hamiltonian information
    
    Class containing all Hamiltonian information:\n
    - L: the number of orbitals
    - groupNumber (in SymmInfo): the number of the Abelian point group symmetry with real-valued character table (see Irreps.h)
    - orb2irrep: array with the irrep number for each orbital
    - Econst: nuclear repulsion energy; or any constant part of the energy not contained in the 1- or 2-particle matrix elements
    - Tmat: 1-particle matrix elements; Tmat\f$_{a,b}\f$ = 0 if \f$I_a\f$ is different from \f$I_b\f$
    - Vmat: 2-particle matrix elements; Vmat\f$_{a,b,c,d}\f$ = 0 if \f$I_a \otimes I_b\f$ is not equal to \f$I_c \otimes I_d\f$; the matrix elements are not antisymmetrized and are stored with the convention that both (a & c) and (b & d) have the same spatial variable for the nuclear repulsion integral (physics notation).
    
    The targeted spin, particle number and point group symmetry are not defined here. For convenience, the second quantized formulation of the Hamiltonian is given here: \n
    \f$ \hat{H} = E_{const} + \sum\limits_{ij\sigma} T_{ij} \delta_{I_i,I_j} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \sigma} + \frac{1}{2} \sum\limits_{ijkl\sigma\tau} V_{ijkl} \delta_{I_i \otimes I_j \otimes I_k \otimes I_l, I_{trivial}} \hat{a}_{i \sigma}^{\dagger} \hat{a}_{j \tau}^{\dagger} \hat{a}_{l \tau} \hat{a}_{k \sigma} \f$\n
    where the latin letters denote site-indices and the greek letters spin projections. This Hamiltonian preserves spin, spin projection, particle number, and Abelian point group symmetry (if its character table is real at least).
 */
class Hamiltonian{

    public:
      
        //! Constructor
        /** \param Norbitals The number of orbitals (L)
             \param nGroup The group number
             \param OrbIrreps Pointer to array containing the orbital irreps */
         Hamiltonian(const int Norbitals, const int nGroup, const int * OrbIrreps, int nup , int ndown, double Econstant);
         Hamiltonian(const int Norbitals, const int nGroup, int nup , int ndown, double Econstant);
         
         //! load an output file created by mointegrals/mointegrals.cc; which can be used as a plugin in psi4 beta5,
         /** \param filename filename of the psi4 mointegral output file 
	  * If this file has .h5 extension we load in the hdf5 file otherwise we assume it is in plain text format and we read data. */
         Hamiltonian(const string filename);

	 //!Copy constructor used in the orbital transformator
         Hamiltonian(const Hamiltonian & hamin);
         //Hamiltonian(Hamiltonian &&) = default;

         Hamiltonian& operator=(const Hamiltonian &);

         //Hamiltonian& operator=(Hamiltonian &&) = default;
         
         //! Destructor
         virtual ~Hamiltonian();
         
         //! Get the number of orbitals
         /** \return The number of orbitals */
         int getL() const;

         //! Get the group number
         /** \return The group number */
         int getNGroup() const;

         //! Get an orbital irrep number
         /** \param nOrb The orbital number
             \return The irrep of orbital nOrb */
         int getOrbitalIrrep(const int nOrb) const;

	     //!returns information about the Hamiltonian
         std::string get_ex_info();
         
         //! Set the constant energy
         /** \param val The new constant energy */
         void setEconst(const double val);
         
         //! Set a Tmat element
         /** \param index1 The first index
             \param index2 The second index
             \param val The new Tmat element */
         void setTmat(const int index1, const int index2, const double val);
         
         //! Set a Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \param val The new Vmat element */
         void setVmat(const int index1, const int index2, const int index3, const int index4, const double val);
         
         //! Add to Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \param val The value which should be added */
         void addToVmat(const int index1, const int index2, const int index3, const int index4, const double val);
         
         //! Get the constant energy
         /** \return The constant part of the Hamiltonian (nuclear repulsion & condensed orbitals) */
        double getEconst() const;
         
         //! Get a Tmat element
         /** \param index1 The first index
             \param index2 The second index
             \return \f$T_{index1,index2}\f$ */
         double getTmat(const int index1, const int index2) const;
         
         //! Get a Vmat element
         /** \param index1 The first index
             \param index2 The second index
             \param index3 The third index
             \param index4 The fourth index
             \return \f$V_{index1,index2,index3,index4}\f$ */
         double getVmat(const int index1, const int index2, const int index3, const int index4) const;
         
         //! Save the Hamiltonian
         void save(const string filename = "") const;
         
         //! Load the Hamiltonian
         void read(const string filename, bool newoverlap = false);
         
         //! Load the Hamiltonian
         void read_file(const string filename);

         //! Save the Hamiltonian to a textfile
	     virtual void save_file(const std::string & filename = "");

         //! Debug check certain elements and sums
         void debugcheck() const;

         //! Get a specific interaction matrix element
         /** \param alpha The first index (0 <= alpha < L)
             \param beta The second index
             \param gamma The third index
             \param delta The fourth index
             \return \f$ h_{\alpha \beta ; \gamma \delta} = \left(\alpha \beta \mid V \mid \gamma \delta \right) + \frac{1}{N-1} \left( \left( \alpha \mid T \mid \gamma \right) \delta_{\beta \delta} + \delta_{\alpha \gamma} \left( \beta \mid T \mid \delta \right) \right) \f$ */
         double gMxElement(const int alpha, const int beta, const int gamma, const int delta) const;

	     int getNORB(int irrep) const { return irrep2num_orb[irrep]; }
	     int getNirreps() const{ return  SymmInfo.getNumberOfIrreps(); }
	     int getNstart(int irrep)const { return _norbcum[irrep]; }
	     int * get_orb2irrep()const { return orb2irrep; }
	     void set_cum();

         int getnup() const;
         int getndown() const;
         void setnup(int nup) { _nup = nup;} 
         void setndown(int ndown) {_ndown = ndown;} 
	     double get_hf_energy(){return _hf_energy;}
	     //std::vector<double> get_moealpha() {return moealpha;}
	     //std::vector<int> get_perm_order();

	     double get_difference(Hamiltonian * ham_in);
	     void set_zero();

	     OptIndex get_index_object();

	     TwoIndex get_tmat_pointer()const; 
	     FourIndex get_vmat_pointer()const; 

	     std::vector<double> get_overlap() const;
         double get_overlap(int irrep , int i , int j )const; //in irrep ordering
         double get_overlap( int i , int j )const; //in orbital index

         void load_overlap( const std::string & filename);
	     void load_overlap(std::istream & file);
	     void print_overlap( std::ostream & file);
	     void set_overlap(int irrep, int i, int j , double val);
	     void set_overlap(int i, int j , double val);
         UnitaryMatrix * get_unitary() const {return _unit.get();};
         void load_unitary(const std::string & filename);
     	 std::string get_filename() const{return _filename;}
     	 void set_filename(std::string file ){_filename = file;}
     	 std::string get_short_filename() const;
     	 virtual std::string get_info() const;

         bool get_modham()const {return _modham;}
         void set_modham(bool mod){_modham = mod;}

    protected:
         //keyword that indicates if this instance is a model hamiltonian.
         bool _modham;
     
    private:

        //number of orbitals
        int L;

         //Overlap matrix (relevant for atomic orbitals), is public because it is not an important variable. Experimental
         std::vector<double>  _overlap;

         //Transformation that defines the transformation of the atomic orbitals defined by the basisset to the orbitals from which the current matrix elements are formed. Experimental
         std::unique_ptr<UnitaryMatrix>  _unit;
         //UnitaryMatrix *  _unit;

        //symmetry info
        Irreps SymmInfo;
        
        //irrep of each orbital
        int * orb2irrep;

        //number of orbitals before the index of present irrep starts
        int * _norbcum;
        
        //number of orbitals per irrep
        int * irrep2num_orb;
        
        //index of an orbital within irrep block
        int * orb2indexSy;
        
        //1-particle matrix elements
        TwoIndex * Tmat;
        
        //2-particle matrix elements
        FourIndex * Vmat;
        
        //Constant part of the Hamiltonian
        double Econst;

	    //1/(N-1)
	    double _oneOverNMinusOne;

	    double _hf_energy;
	    int _nup;
	    int _ndown;
	    //std::vector<double> moealpha;
	    std::string _filename;
};


class OptIndex
{
    public:
        OptIndex(const int L, const int Group, const int * NORBin);
        OptIndex(const OptIndex & index);
        int getL() const; 
        int getNirreps() const;
        int getNORB(const int irrep) const;
        int getNstart(const int irrep) const;
        int * get_irrep_each_orbital();
        void Print() const;
    private:
        int Nirreps;
        int L;
        std::vector<int> NORB;
        std::vector<int> NORBcumulative;
        std::vector<int> irrep_each_orbital;
};



#endif
