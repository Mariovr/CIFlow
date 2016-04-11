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
#ifndef __OrbitalOptimization__
#define __OrbitalOptimization__

#include <memory>
#include <random>
#include <vector>
#include <string>

#include "matrix.h"

class Hamiltonian;
class UnitaryMatrix;
class OrbitalTransform;
class CIMethod;
class CIDens;
class DIIS;

class OrbitalOptimization
{
    public:
        OrbitalOptimization(CIMethod * cim);
        virtual ~OrbitalOptimization(); //Orbital optimizations contains an orbital transform and when it is deleted. We multiply the resulting unitary matrix of the orbital optimization in the unitary matrix contained in _optham of the original hamiltonian, which defines the original hamiltonian. So at the end the unitary matrix contained in _optham represents the transformation from the atomic orbitals to _optham. This can be used for visualisation in a molden file.
        virtual double optimize() = 0 ;
        void get_back_original_ham();
        double get_opt_ci_energy() { return _opt_energy; }
        void set_cim_optimal();
        void update_cim();
        //!Transforms _optham to the natural orbitals of the 1rdm in cid.
        double transform_to_no(CIDens & cid);
        //!Transforms to the NO's of previous solution.
        double get_cim_no();
        //!Transforms to the FCI no's.
        double get_fci_no();
        //!Sets the ham in cim to _optham and calculates the cim energy in this basis.
        double set_basis(CIMethod & cim);

        //!get hamiltonian pointer
        Hamiltonian * get_ham(){return _optham ;}
        UnitaryMatrix * get_opt_unitary();

        void save_ham(std::string name);
        //Saves the unitary matrix that rotates the start Hamiltonian, to the Hamiltonian with 
        //the desired properties (if hdf5 is true in hdf5 format else in textfile format)
        void save_unitary(std::string name = "" , bool hdf5 = 0);
        //virtual void reset_values() = 0;

    protected:
        void  transform_orbitals(int irrep, int i, int j, double angle);
        void transform_withasx(double * step, double temp, double max_angle);
        double _opt_energy;
        std::unique_ptr<OrbitalTransform> _orbtrans;
        CIMethod * _cim;
        //This is commented out because you have to make sure the Hamiltonian pointer is not used by other instances anymore.
        Hamiltonian *  _optham;
        std::unique_ptr<UnitaryMatrix> _opt_unitary;
};

class SimulatedAnnealing: public OrbitalOptimization
{
    public:
        SimulatedAnnealing(CIMethod * cim, double temp, double change_temp, double angle, double change_angle);
        double optimize();
        void run_multiple(int num, bool saveopt = true);
        bool accept_function(double e_new, double e_old, double temp, std::mt19937 &mt, std::uniform_real_distribution<double> &dist);
        double cost_function();
        double get_cooling_steps() const { return _c_steps; }
        void set_cooling_steps(int c_step) { _c_steps = c_step; }
        void set_good_start();

    private:
        double _temp;
        double _change_temp;
        double _max_angle;
        double _change_angle;
        unsigned int _c_steps;
};

/* This class uses the iterative procedure of Subotnik to rotate the orbitals to get an extremal value of some quantity dependend on the same orbitals.
 * See: Subotnik, J. E. and Shao, Y. and Liang, W. Z. and Head-Gordon, M. : An efficient method for calculating maxima of homogeneous functions of orthogonal matrices: Applications to localized occupied orbitals.
 * This method can be used to localize orbitals to the Edminston Rudenberg orbitals, to obtain the orbitals which minimize the seniority number of the wavefunction, ...
 * In current implementation we minize the seniority number by maximizing the spinsummed 2rdm D^ii_ii
 * REMARK: this paper saves unitary transformations in columns.
 */
class Iterative_Subotnik: public OrbitalOptimization
{
    public:
        Iterative_Subotnik(CIMethod * cim, double crit , int type);//type = 0 (seniority minimization), type = 1 (maximizing self repulsion)
        void set_good_start();
        virtual double optimize(); 
        void fill_R_ER(int all= 1); //only occupied hf orbitals type = 0 , or all orbitals type =1
        void fill_R_mmin();
        void rotate_unitary_from_R();
        double check_converged();
        //! Returns the current value you want to maximize (Trace(R) ).
        double calc_value();
        int get_type(){return _type;}
        void set_type(int type){ _type = type;}

     protected:
        matrix _R;
        int _type;
        double _crit;
};

class Iterative_Subotnik_DIIS: public Iterative_Subotnik
{
    public: 
        Iterative_Subotnik_DIIS(CIMethod * cim, double crit, int type);
        ~Iterative_Subotnik_DIIS();
        void rotate_unitary_from_extrapolated_R();
        void interpolate_unitary();
        void interpolate_R();
        void transform_R();
        double optimize(); 
        matrix calc_error_matrix();
        void add_unitlist(UnitaryMatrix & val);

    private:
        std::vector<UnitaryMatrix *> _unitlist;
        DIIS * _diis;
};


class LineSearch: public OrbitalOptimization{
    public:
        LineSearch(CIMethod *cim): OrbitalOptimization(cim){};
        //work function that executes the actual optimization. At this moment Powell method.
        double optimize();

        double func(double angle , int i , int j);
        /*   
        Finds the brackets (a,b) of a minimum of the energy of cim
        The search starts downhill from the position x1 with step length h.
        Change x1 if you think the minimum is located near another point
        //Fills bracket with an interval in which a minimum is located.*/
        void find_bracket(double x1, double step, double * bracket, int i , int j);
        /*  
        std::vector<double> (xmin , fmin) = search(i,j,a,b,tol=1.0e-10) returns angle and optimized energy of rotating orbital i and j.
        With a Golden section method for determining the angle that minimizes
        the CI energy. The minimum must be bracketed in (a,b).*/
        std::vector<double> search(int i , int j , double a , double b , double tol = 1e-10);


    private:

};

#endif

/* vim: set ts=4 sw=4 expandtab :*/
