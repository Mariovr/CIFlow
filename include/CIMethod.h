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
#ifndef CIMethod_H
#define CIMethod_H

#include <string>
#include <utility>

//Includes (Only do this in header files for parent classes or classes contained as an object in the present class)
#include "binomial.h"
#include "matrix.h"
#include "SparseMatrix_CRS.h"
#include "Options.h" //definition of the integer TYPE of the bit representation based on the number of orbitals.

//Forward declarations
class Hamiltonian;
class Permutator;
class Permutator_File;
class Permutator_Bit;
class CIDens;
class OutputSingleFile;

class CIMethod{
	public:
		CIMethod(Hamiltonian * ham, const int nup , const int ndown);
		CIMethod(Hamiltonian * ham );
		virtual ~CIMethod(); 
		virtual void build_ci() = 0;
		virtual void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl) = 0;
		virtual int get_dim() = 0;
		virtual std::string get_name() = 0;
		virtual std::unique_ptr<CIDens> get_density() = 0;
        void construct_density(unsigned state = 0, bool trdm = true);
        void reset_density();


		//generate eigenvalues and eigenvectors
		virtual void solve(int neigval = NEIGVAL);
		double get_ci_energy(int neigval = 0);
		double get_hf_energy();

		//get and set functions
		Hamiltonian * get_ham();//{ return _ham ;}
		Permutator * get_perm();//{return _perm ;}
		double get_eigvec(int numstate , int numvec);
        std::vector<double> get_eigvec(int numvec);
		matrix get_eigval(){return _eigval;}
		matrix get_mat();
		double get_mat(int i , int j);
		int gUpDim()const;// { return _updim; }
		int gDownDim()const;// { return _downdim; }
		/** \return the number of particles with up spin.*/
		int gNup()const;// {return _nup; }
		/** \return The number of particles with down spin. */
		int gNdown()const;// { return _ndown; }
		int get_l() const ;//{return _ham->getL();}
		double get_econst() const;
		void set_ham(Hamiltonian * ham, bool build = true);   
		void reset_mat(){_mat.reset_vecs(); }
		bool is_solved(){return _solved;}
		void build_parallel();
		//loads ham from file hamname
		void load_ham(std::string hamname , bool build = true, bool reset = false);
		void delete_ham(); //Needs to be called when load_ham is used.
		void new_ham_pointer(); //Needs to be called when load_ham is used.
        double get_sz(){return (_nup - _ndown) / 2.;}
        double get_spin_squared(unsigned state = 0);
        double get_spin(unsigned state = 0);
        double get_mulliken(std::vector<int> orbs, unsigned state = 0);
        void set_sparsity(double zero);

		//helper functions to construct CI matrices and density matrices.
		int * get_orbs(TYPE, int nset);
		double get_offdiag_elem(TYPE up1);
		double get_diag_cont_one_bitstring(int * orbs, int npar);
		double get_diag_cont_jmix(int * orbsup, int *  orbsdown);
		void arrange_order(int * orbs , int num , TYPE up1 , TYPE up2);
		int get_nstacked(int * orbs, TYPE up);

		//functions necessary for output
		std::string get_ex_info();
		void print_output(std::vector<std::string> vec = {} , int num = 0 );
		void print_ham();
		void print_rdm(unsigned state = 0, bool trdm = false);
        void reset_output(std::string sort, bool partial = true);

		//test and check functions
		void check_hermiticity();

		matrix _eigval;
		matrix _eigvec;

		std::unique_ptr<Permutator> _perm;
        std::string get_perm_info();

	protected:
		bool _solved;
		int _nup;
		int _ndown;
		int _updim;
		int _downdim;
        std::unique_ptr<CIDens> _cid; //The density of an eigenstate .
		double _oneOverNMinusOne;
		Hamiltonian * _ham;
		SparseMatrix_CRS_Sym _mat;
        std::unique_ptr<OutputSingleFile> _output;
		//expectation values of the nonrelativistic quantumchemical Hamiltonian on different types of slaterdeterminants.
		double diagonal(TYPE up1 , TYPE down1);
		double one_diff_orbital(TYPE up1 , TYPE up2, TYPE down1);
		double two_diff_orbitals(TYPE up1, TYPE up2);
		double two_diff_orbitals(TYPE up1 , TYPE up2, TYPE down1, TYPE down2);
		//analyses the relationships between two different slaterdeterminants and returns the projection on the nonrel. quantumchem. Ham.
		double get_ham_element(TYPE up1 , TYPE down1 , TYPE up2 , TYPE down2);
};

class DOCI: public CIMethod{
	public:
		DOCI(Hamiltonian* ham);
		void build_ci(){build_parallel() ; }
		void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl);
		int get_dim(){ return gUpDim(); } 
		std::string get_name() {return "DOCI";}
        std::unique_ptr<CIDens> get_density();

    private:
        unsigned int determine_weight(TYPE string , const std::vector<std::vector<int> > & vw);
        void setup_vertex_weights(std::vector<std::vector<int>> & vw);
};

class FCI: public CIMethod{
	public :
		FCI(Hamiltonian*  ham); //Set print False if you are not interested to print FCI data.
        //!REMARK is not used at the moment. It is here for compatibility with the baseclass.
		void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl);	
        //!REMARK: this function doesn't use a sparse matrix, to fill the Hamiltonian, to maximally benefit from the extra speedup by using spinup and down symmetries for which it is easier to use a normal matrix. So it is very memory constrained use it only to test, debug and for reference data of small systems.
		void build_ci(){build_parallel() ; }
		int get_dim(){ return gUpDim() * gDownDim(); } 
		std::string get_name() {return "FCI";}
		std::unique_ptr<CIDens> get_density();

    private:
        unsigned int determine_weight(TYPE string , const std::vector<std::vector<int> > & vw);
        void setup_vertex_weights(std::vector<std::vector<int>> & vw);
};

class FCI_File: public CIMethod{
	public :
		FCI_File(Hamiltonian*  ham, std::string permfile);
		void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl); 
		void build_ci(){build_parallel() ; }
		int get_dim();
		std::string get_name() ;
        std::unique_ptr<CIDens> get_density();
        void print_dets();
};


class CI_Big: public CIMethod{
	public :
		CI_Big(Hamiltonian*  ham, std::string permfile);
		void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl);
		void build_ci(){build_parallel() ; }
		void solve(int neigval = 1);
		void arpackDiagonalize(int nve , matrix & eigval , matrix & eigvec);
		int get_dim();
		void mvprod(const matrix &xmat, matrix &ymat);
		void mvprod(const double *x, double *y); 
		std::string get_name() ;
        std::unique_ptr<CIDens> get_density();
        void print_dets();
};


/*  
class CI_BigDOCI: public CI_Big{
	public :
		CI_Big(Hamiltonian*  ham, std::string permfile);
		void construct_CI_matrix(SparseMatrix_CRS & mat , int startl , int endl);
		void build_ci(){build_parallel() ; }
		void construct_dict();
		void solve(int neigval = 1);
		void arpackDiagonalize(int nve , matrix & eigval , matrix & eigvec);
		int get_dim();
		void mvprod(const matrix &xmat, matrix &ymat);
		void mvprod(const double *x, double *y); 
		std::string get_name() ;
};*/
#endif
