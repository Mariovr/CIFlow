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
#include "vector"

#include "Hamiltonian.h"
#include "CIMethod.h"
#include "CIDens.h"
#include "RedBCS.h"

template <typename T> 
std::vector<int> sort_indexes(const std::vector<T> &v) 
{
  // initialize original index locations
  std::vector<int> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2){return v[i1] < v[i2];});
  return idx;
}

RedBCS::RedBCS(Hamiltonian & Ham)
{
    _hamrg =  new Hamiltonian(Ham.getL() , Ham.getNGroup() ,Ham.get_orb2irrep() , Ham.getnup() , Ham.getndown(), 0.  ); //We set econst to zero, because reduced BCS ham doesn't have a econst to calculate the energy.
    std::vector<double> epsilon = get_epsilon(Ham);
    _doci.reset(new DOCI(_hamrg )); 
    set_ham(0.0001 , epsilon , 0);
}

RedBCS::RedBCS(Hamiltonian & Ham, double g , std::vector<double> epsilon, std::vector<int> permarray)
{
    _hamrg = new Hamiltonian(Ham.getL() , Ham.getNGroup() ,Ham.get_orb2irrep() , Ham.getnup() , Ham.getndown()  , 0.);
    _permarray = permarray;
    _doci.reset(new DOCI(_hamrg )); 
    set_ham(g , epsilon , 0);
    //_hamrg->debugcheck();
}

RedBCS::~RedBCS()
{
    delete _hamrg;
}

//we assume the epsilons are already sorted with permarray the sorted array of indices.
void RedBCS::set_ham(double g , std::vector<double> epsilon  , bool reset)
{
    if (reset)
    	_hamrg->set_zero();
    for(int i = 0 ; i < _hamrg->getL() ; i++)
    {
    	_hamrg->setTmat(i,i, epsilon[i]);
    }
    for(int i = 0 ; i < _hamrg->getL() ; i++)
	    for(int j = 0 ; j < _hamrg->getL() ; j++)
	    {
	         _hamrg->setVmat(i,i,j,j,g);	
	    }
    //Because we changed the Hamiltonian	    
    _doci->set_ham(_hamrg);
}

void RedBCS::set_g(double g)
{
    for(int i = 0 ; i < _hamrg->getL() ; i++)
	    for(int j = i ; j < _hamrg->getL() ; j++)
	    {
	         _hamrg->setVmat(i,i,j,j,g);	
	    }
    _doci->set_ham(_hamrg); //Because _hamrg is changed we need to update the CIMatrix.
}

void RedBCS::set_epsilon( std::vector<double> epsilon)
{
    for(int i = 0 ; i < _hamrg->getL() ; i++)
	 _hamrg->setTmat(i,i,epsilon[i]);	
    _doci->set_ham(_hamrg); //Because _hamrg is changed we need to update the CIMatrix.
}

void RedBCS::set_epsilon( double  epsilon, int i)
{
    _hamrg->setTmat(i,i,epsilon);	
    _doci->set_ham(_hamrg); //Because _hamrg is changed we need to update the CIMatrix.
}

void RedBCS::add_epsilon( double  epsilon, int i)
{
    _hamrg->setTmat(i,i,epsilon+_hamrg->getTmat(i,i) );	
    _doci->set_ham(_hamrg); //Because _hamrg is changed we need to update the CIMatrix.
}

double RedBCS::calculate_energy()
{
    _doci->solve();
    _cid = _doci->get_density() ; //We automatically calculate the density from the ground-state solution.
    _cid->construct_density();
    return _doci->get_ci_energy(); //The Hamiltonian presumes that the nuclear repulsion is in the epsilon.
}

double RedBCS::get_energy()
{
    return _doci->get_ci_energy();
}

//Returns a vector containing epsilon guesses for the reduced BCS Hamiltonian.
std::vector<double> RedBCS::get_epsilon(Hamiltonian & ham)
{
    std::vector<double> epsilon(_hamrg->getL() ); //Saves speed, because vector doesn't have to reallocate regularly.
    for(int i = 0 ; i < _hamrg->getL() ; i++)
    {
    	double som = ham.getTmat(i,i);
	//for(int j = 0 ; j < _hamrg->getnup() ; j++)
	//{
          //som += 1./2. * ( ham.getVmat(i,j,i,j)*2 - ham.getVmat(i,j,j,i));
          //som += ham.getVmat(i,j,i,j)*2 - ham.getVmat(i,j,j,i);
	//}
	epsilon[i] = som;
     }
    _permarray = (std::vector<int>) sort_indexes(epsilon);
    std::sort(epsilon.begin(), epsilon.end());
    return epsilon;
}


CIDens * RedBCS::get_density()
{
    //return _doci->get_density(); //Creates always a new pointer.
    return _cid.get(); // Uses the same pointer, so we don't need to calculate the 2rdm again, and if it changes all instances have acces to the same changed variant.
}

double RedBCS::get_energy(Hamiltonian & ham)
{
    double energy = 0.0;
    //one body energy terms.
    for (int i = 0; i < _doci->get_l(); i++) {
        energy += 2.*ham.getTmat(_permarray[i],_permarray[i]) * _cid->get_one_rdm(0, i,i  );
        //std::cout << ham.getTmat(_permarray[i],_permarray[i]) << " " << energy << std::endl;
    }
    //two body energy terms.
    for (int i = 0; i < _doci->get_l(); i++) {
        for (int j = 0; j < _doci->get_l(); j++) {
            double vmat = ham.getVmat(_permarray[i],_permarray[i],_permarray[j],_permarray[j]);
            double vmat2 = ham.getVmat(_permarray[i],_permarray[j],_permarray[i],_permarray[j]);
            double vmat3 = ham.getVmat(_permarray[i],_permarray[j],_permarray[j],_permarray[i]);
            //std::cout << vmat << " " << vmat2 << " " << vmat3 << std::endl;
            if( i != j )
            {
        	    energy +=  (2*vmat2-vmat3) *_cid->get_two_rdm(1,i,j,i,j);
        	    energy +=  vmat * _cid->get_two_rdm(1,i,i,j,j);
            }
        }
        double vmat = ham.getVmat(_permarray[i],_permarray[i],_permarray[i],_permarray[i]);
        energy += vmat * _cid->get_two_rdm(1,i,i,i,i);
    }
    energy += ham.getEconst(); //nuclear repulsion constant.
    return energy;
}

std::vector<double> RedBCS::get_epsilon()
{
    std::vector<double> epsilons(_doci->get_l());
    for (int i = 0; i < _doci->get_l(); i++) 
        epsilons[i] = _hamrg->getTmat(i,i); 
    return epsilons;
}

std::string RedBCS::get_information()
{
    std::ostringstream info;
    info << std::setprecision(16);
    info << "#g = " << get_g() << "     epsilon = ";
    for(int i = 0 ; i < _doci->get_l(); i++)
    {
        info << get_epsilon(i) << " ";
    }
    info << "   npair = " << _doci->gNup();
    info << "    E = " << get_energy() <<std::endl; 
    return info.str();
}

/*  
double RedBCS::get_energy(Hamiltonian & ham)
{
    double energy = 0.0;
    //one body energy terms.
    for (int i = 0; i < _doci->get_l(); i++) {
        for (int j = 0; j < _doci->get_l(); j++) {
            energy += ham.getTmat(_permarray[i],_permarray[j]) * _cid->get_one_rdm(0, i,j );
            energy += ham.getTmat(_permarray[i],_permarray[j]) * _cid->get_one_rdm(1, i , j );
	    std::cout << " " << energy << std::endl;
        }
    }
    //two body energy terms.
    double energytwo = 0;
    for (int i = 0; i < _doci->get_l(); i++) {
        for (int j = 0; j < _doci->get_l(); j++) {
            for (int k = 0; k < _doci->get_l(); k++) {
                for (int l = 0; l < _doci->get_l(); l++) {
                    double vmat = ham.getVmat(_permarray[i],_permarray[j],_permarray[k],_permarray[l]);
                    energytwo +=  vmat * _cid->get_two_rdm(0,i,j,k,l);
                    energytwo +=  vmat * _cid->get_two_rdm(1,i,j,k,l);
                    energytwo +=  vmat * _cid->get_two_rdm(2,i,j,k,l);
                    energytwo +=  vmat * _cid->get_two_rdm(3,i,j,k,l);
                }
            }
        }
    }
    energytwo *= 0.5;
    energy += energytwo;
    energy += ham.getEconst(); //nuclear repulsion constant.
    return energy;
}*/

