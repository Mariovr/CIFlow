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
#ifndef __CIDENS_H
#define __CIDENS_H

#include <valarray>
#include <vector>

#include "matrix.h"
#include "Options.h"

//forward declaration to break the loop
class CIMethod;

class CIDens
{
	/* \brief This class constructs the 1-RDM and 2-RDM, starting from a CI-coefficient vector.
	 * The RDMs can be seen as the expectation values over excitation operators.
         * Four sorts of relevant 2rdm: 0 -> up, up, up ,up ,1 -> up, down, up, down, 2 -> down, up, down, up, 3-> down, down, down, down (where the first two are for creation and the last two for annihilation.) 
	 * ALL other combinations are irrelevant because when contracted with the matrixelements they give zero because they integrate down with up spinorbitals or they dont keep the number of ups and downs constant
	 * REMARK a_i+ a_j+ a_k a_l -> 2rdm(ij , lk) because 2rdm corresponds to physicist notation of inproduct.
    // In memory we keep only the 0,1,3 cases , because the 2 is completely equal to the 1. 
    // Because other notation of 2rdm = (pq,rs) with p < q , r < s (this for spin orbitals and we assume up < down for corresponding spatial orbitals (first L ups than L downs) )
    // So for case 1 we keep all elements because p < q and r < s always (u,d,u,d) , and for 0, and 3 we keep only the elements with p<q , r<s and the other ones are generated by changing fermion sign.  */
	public:
		CIDens(CIMethod * cim, unsigned state = 0); //state determines from which eigenvector of the CImatrix the density matrix is going to be constructed.
		virtual ~CIDens(){};
		//1rdm and 2rdm routines.
        void construct_density(bool twordm = true);
        void test_invariants(bool twordm = true);
		double get_one_electron_energy(int maxorb = -1);
		double get_two_electron_energy(int maxorb = -1);
		double get_dens_energy(int maxorb = -1);
		double get_one_pare(std::vector<int> vw); //energy decomposition of orbitals defined by vw
		double get_two_pare(std::vector<int> vw);
		double get_dens_pare(std::vector<int> vw);
        CIMethod * get_cim()const{return _cim;};
        std::vector<Vector2d >  get_two_dm(){return _two_dens;} 
        std::valarray<std::valarray<double> > get_one_dm(){return _one_dens;}

		void print_one_dens(std::ostream & os, const std::valarray<std::valarray<double>> & onerdm ) const;
		void print_one_dens(std::ostream & os) const {print_one_dens(os , _one_dens); }
	    void print_two_dens(std::ostream & os, const std::vector<Vector2d > & twordm ) const;
	    void print_two_dens(std::ostream & os) const{print_two_dens(os, _two_dens); }
		void compare_one_dens(CIDens *cid);
		void compare_two_dens(CIDens *cid);
		double get_one_rdm(int spin , int orb1 , int orb2, const std::valarray<std::valarray<double>> & onerdm) const;
		double get_one_rdm(int spin , int orb1 , int orb2) const {return get_one_rdm(spin , orb1,orb2 , _one_dens);}
        void add_one_rdm(int spin , int orb1 , int orb2 , double val, std::valarray<std::valarray<double>> & onerdm);
		double get_two_rdm(int spin , int orb1 , int orb2, int orb3 , int orb4,const std::vector<Vector2d > & twordm) const ; 
		double get_two_rdm(int spin , int orb1 , int orb2, int orb3 , int orb4) const {return get_two_rdm(spin ,orb1 , orb2 , orb3 , orb4 , _two_dens); } 
        void add_two_rdm(int spin, int orb1 , int orb2 , int orb3 , int orb4, double value, std::vector<Vector2d > & twordm);
        double operator()(int spin, int i , int j , int k ,int l) const {return get_two_rdm(spin, i , j , k , l);}
        double operator()(int spin, int i , int j ) const {return get_one_rdm(spin, i , j );}

        void reset_1rdm(); //resizes 1rdm to zero.
		void reset_2rdm(); //resizes 2rdm to zero.
        void reset_density();
		double get_seniority()const;
        double get_spin_squared()const;
        double get_spin()const;
        double get_mulliken(std::vector<int> orbs);//orbs defines the indices of the orbitals over which the partial trace runs.
        unsigned dim() const;

        void set_state(int state){if (state != _state ){ _state = state; reset_density(); }  }
        unsigned get_state() const {return _state;}

        std::pair<bool, bool> is_constructed(){return _solved; } //Check if we already have a two dens.
        void set_constructed(std::pair<bool,bool> solved){_solved = solved ;}

        //Puts in occupations the occupations of the orbitals, and in no the natural orbitals.
		//Get the NO in decreasing order of occupations.
		//New basis is saved in n columns.
		void get_no(matrix & occupations , matrix & no) const;
		matrix spin_summed_1rdm()const ;

        virtual std::pair<double,bool> find_min_angle(int k, int l, double start_angle, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const {return std::pair<double, bool> {-1000. , true} ; }
        virtual double calc_rotate(int k, int l, double theta, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const {return -1000. ;}

 
        void set_one_rdm(int spin ,int orb1 , int orb2, double value);
        void set_two_rdm(int spin ,int orb1 , int orb2, int orb3, int orb4, double value);
        void transform_ordm(const matrix & unitary);
        void transform_trdm(const matrix & unitary);
        //transforms the rdms to the AO basis, REMARK if this is used implicitely as is the case for the mulliken function, the density matriceces contained are not anymore in the MO, but in the AO so pay attention if one wants to use the density matrices later (so keep track of the basis) after a call to for example mulliken. To revert this transformation call this function with the revert bool = true.
        void transform_to_ao(bool trdm = true, bool revert = false);

        matrix get_oblock(int spin); //returns the square density matrix associated to a particular spin block of the ordm.
        void set_oblock(int spin, const matrix & dens);
        //matrix get_tblock(int spin);//returns the square density matrix associated to a particular spin block of the trdm (i,j , k,l). 
        //void set_tblock(int spin);//returns the square density matrix associated to a particular spin block of the trdm (i,j , k,l). 


	protected:
        unsigned _norb;
        unsigned _state;
        std::pair<bool, bool>  _solved;
		CIMethod * _cim;
        std::valarray<std::valarray<double >> _one_dens;
        std::vector< Vector2d >  _two_dens;
        void allocate_one_memory(std::valarray<std::valarray<double>> & onerdm);
        void allocate_two_memory(std::vector<Vector2d > & twordm);
		virtual void construct_CI_one_dens(unsigned start, unsigned end, std::valarray<std::valarray<double>> & onerdm) = 0;
		virtual void construct_CI_two_dens(unsigned start, unsigned end, std::vector<Vector2d > & twordm) = 0;
		//up_or_down < 2 fill only one_dense , > 2 fill 2rdm (3 = ups, 4 = downs).
        void fill_one_dense(TYPE bitstring1 , TYPE bitstring2,TYPE same ,  double contr, int up_or_down, std::valarray<std::valarray<double>> & onerdm, std::vector<Vector2d > & twordm);
	    //up_and_down -> 0 if up1 and up2 are upstring, 3 if up1 and up2 are downstrings.
        void fill_two_dense(TYPE up1 , TYPE up2, double contr, int up_or_down, std::vector<Vector2d > & twordm);
        void fill_two_dense_upanddown(TYPE up1 , TYPE up2, TYPE down1, TYPE down2, double contr, std::vector<Vector2d > & twordm );
        void fill_diagonal(TYPE up1 , TYPE down1, double contr, bool twordm, std::valarray<std::valarray<double>> & onerdm, std::vector<Vector2d > & two_rdm);
		void fill_density_matrix(TYPE up1 , TYPE down1 , TYPE up2, TYPE down2, double contr, std::vector<Vector2d > & twordm);
		void fill_density_matrix_one(TYPE up1 , TYPE down1 , TYPE up2 , TYPE down2, double contr,std::valarray<std::valarray<double>> & onerdm );
        //Functions for parallel calculation
        void build_parallel(bool twordm);
        void add_from_vector_one(std::vector< std::valarray<std::valarray<double>> > & parts);
        void add_from_vector_two(std::vector<std::vector<Vector2d >  > & parts);
};

class DensDOCI:public CIDens
{
	public:
		DensDOCI(CIMethod * cim):CIDens(cim){};

        //Functions necessary for the local optimizer.
        std::pair<double,bool> find_min_angle(int k, int l, double start_angle, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const;
        double calc_rotate(int k, int l, double theta, std::function<double(int,int)> &T, std::function<double(int,int,int,int)> &V) const;

    private:
		void construct_CI_one_dens(unsigned start, unsigned end, std::valarray<std::valarray<double>> & onerdm);
		void construct_CI_two_dens(unsigned start, unsigned end, std::vector<Vector2d > & twordm);
        //helper functions for the density matrices.
        void setup_vertex_weights(std::vector<std::vector<int> > & vw);
        /** This function determines the weight of a given string.
         *
         * This weight can be used to index contiguous arrays.
         * The vertex_weights are set up by a helper function.
         */
        unsigned int determine_weight(TYPE string, const std::vector<std::vector<int> > & vw);

};


class DensFILE:public CIDens
{
	public:
		DensFILE(CIMethod * cim ): CIDens(cim){};

    private:
		void construct_CI_one_dens(unsigned start, unsigned end, std::valarray<std::valarray<double>> & onerdm);
		void construct_CI_two_dens(unsigned start, unsigned end, std::vector<Vector2d > & twordm);
};

class DensFCI:public CIDens
{
	public:
		DensFCI(CIMethod * cim): CIDens(cim){};
    private:
		void construct_CI_one_dens(unsigned start, unsigned end, std::valarray<std::valarray<double>> & onerdm);
		void construct_CI_two_dens(unsigned start, unsigned end, std::vector<Vector2d > & twordm);
};

#endif
