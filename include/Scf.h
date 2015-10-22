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
#ifndef __SCF_H
#define __SCF_H

#include <matrix.h>

class Hamiltonian;
class UnitaryMatrix;
class DIIS;

class RHF
{
    public:
        RHF(Hamiltonian * ham, int max_iter , double e_conv , double d_conv , int print, int diis); 
        ~RHF();
	    UnitaryMatrix load_overlap(std::string overlap_file);//Remark unitary matrix is not unitary :), because of nonorthogonal atomic orbitals.
        UnitaryMatrix get_transformation();

	    void symmetric_diag();
	    void construct_density();
	    double get_hf_energy(); //Depends on density (which depends on c the unitary transformation from ao to mo) and fock matrix.
	    double scf(); //returns scf energy

            matrix get_orbital_energies(int print = 0);
	    void diag_symmetric_matrix( matrix Fockt, matrix & val_part , matrix & part);
	    void sort_energies(matrix & orbe);

        int get_max_iter(){return _max_iter ; }
        double get_e_conv(){return _e_conv ; }
        double get_d_conv(){return _d_conv ; }
        void set_max_iter(int maxi){_max_iter  = maxi; }
        void set_e_conv(double ec){ _e_conv = ec; }
        void set_d_conv(double dc){ _d_conv = dc; }

    protected:

    private:
    	std::unique_ptr<Hamiltonian> _ham;
      	matrix _X; //transformation that makes the non-orthogonal atomic orbitals orthogonal.
      	matrix _P;
      	matrix _Fock;
      	matrix _Fockt ; //Contains transformed fock matrix to orthogonal basis
      	matrix _C;
      	matrix _overlap;
      	matrix _orbitale;

	    int _max_iter;
	    double _e_conv;
	    double _d_conv;
	    int _print;
	    int _dodiis;

	    DIIS * _diis;

	    std::vector<int> _orbarray; //Contains the map to the lowest (orbital energy) orbitals, this is necessary
	    //because we keep spatial symmetry and only order by energy in one irrep.
};

#endif
