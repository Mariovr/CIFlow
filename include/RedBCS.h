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
#ifndef RedBCS_H
#define RedBCS_H

#include <vector>

class Hamiltonian;
class CIMethod;
class DOCI;
class CIDens;

class RedBCS
{
    public:
        RedBCS(Hamiltonian & Ham);
        RedBCS(Hamiltonian & Ham, double g , std::vector<double> epsilon, std::vector<int> permarray);
        ~RedBCS();

      	double calculate_energy();
      	//Returns DOCI energy of Red BCS Hamiltonian.
      	double get_energy();
      	//Returns the energy of the current wavefunction multiplied with the matrix elements of ham and changed indices with _permarray;
      	double get_energy(Hamiltonian & ham);
      
        void set_ham(double g , std::vector<double> epsilon,  bool reset);
      	void set_permarray(std::vector<int> permarray){_permarray = permarray;}
      	void set_g(double g);
      	void set_epsilon(std::vector<double> epsilon);
      	void set_epsilon(double epsilon, int i);
      	void add_epsilon(double difepsilon, int i);
      	//Generates epsilons from ham.
      	std::vector<double> get_epsilon(Hamiltonian & ham);
      
      	CIDens * get_density();
      	//returns pair coupling constant of current Red BCS Hamiltonian.
      	double get_g(){return _hamrg->getVmat(0,0,0,0);}
      	double get_epsilon(int i ){return _hamrg->getTmat(i,i); }
      	//returns the single particle energies of current Red BCS Hamiltonian.
      	std::vector<double > get_epsilon();
      	//returns the permutation array.
      	std::vector<int> get_permarray(){return _permarray;}
      
        std::string get_information();

    private:
        std::unique_ptr<DOCI>  _doci;
        Hamiltonian * _hamrg;
        std::unique_ptr<CIDens> _cid;
        std::vector<int> _permarray;
};


#endif
