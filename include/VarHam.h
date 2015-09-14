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
#ifndef VarHam_H
#define VarHam_H

#include <vector>

class Hamiltonian;
class CIDens;
class RedBCS;
//class OrbitalOptimization;
class SimulatedAnnealing;

class VarHam //: public CIMethod
{
    public:
    	VarHam(Hamiltonian * Ham, int epsilontype);
	~VarHam(){};
	std::pair<double, double > optimize_g(double crit);
	void print_gscan(std::ostream & os);
	double energy_from_g(double g);
	double energy_from_epsilon(double epsilon , int num);
	double solve();

        //For simulated annealing to work.
	Hamiltonian * get_ham(){return _chemham.get() ;}
	void set_ham(Hamiltonian * ham){_chemham.reset(new Hamiltonian(*ham) ) ;}
	int get_l(){return _chemham->getL();}
	double get_ci_energy();
	bool is_solved();

        //Returns the energy of the current model hamiltonian posessed by the VarHam object.
	double get_mod_ham_energy();
	//Returns info about the current VarHamm object.
	std::string get_info();

	void set_epsilon(std::vector<double> epsilon){_redbcs->set_epsilon(epsilon) ;}
	void set_g(double g){_redbcs->set_g(g) ;}
	double get_g(){return _redbcs->get_g();}
	std::vector<double > get_epsilon(){return _redbcs->get_epsilon();};


    private:
    	//std::unique_ptr<SimulatedAnnealing>  _orbopt;
        std::unique_ptr<RedBCS> _redbcs;
        std::unique_ptr<Hamiltonian> _chemham;

};


#endif
