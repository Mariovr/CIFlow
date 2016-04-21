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
#include <iomanip>
#include <iostream>
#include <cmath>
#include <utility>

#include "Hamiltonian.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"
#include "RedBCS.h"
#include "VarHam.h"


VarHam::VarHam(Hamiltonian * Ham , int epsilontype)
{
    _redbcs.reset(new RedBCS(*Ham));
    _chemham.reset(new Hamiltonian(*Ham));
}


void VarHam::print_gscan(std::ostream & os)
{
    double old_g = _redbcs->get_g();
    os <<"#g\tE_MOD_HAM\tE_VarRG" <<std::endl;
    os.precision(10);
    for( double g = 0.0001 ; g < 1. ; g += 0.005)
    {
        _redbcs->set_g(g);
	double e_redbcs = _redbcs->calculate_energy();
	double e_var = _redbcs->get_energy(*_chemham);

        os <<g << "\t"<<  e_redbcs<< "\t"<< e_var <<std::endl;
    }
    _redbcs->set_g(old_g); //We set our reduced BCS Hamiltonian back to the previous g.
    _redbcs->calculate_energy();
}

//Simple implementation of the golden search method.
std::pair<double,double> VarHam::optimize_g(double crit)
{
    double gstart = 0.0001;
    double gend = 1.+ _redbcs->get_g() ;
    
    double ref = _redbcs->calculate_energy();
    int it = 0;
    while(energy_from_g(gend) < ref && it < 10)
    {
    	gend += 0.5;
	    it +=1 ;
    }

    int numit = std::ceil(std::log(crit/(gend-gstart))/std::log(0.618));
    for(int i = 0 ; i < numit ; i++)
    {
    	double x1 = gend - 0.618*(gend -gstart);
    	double x2 = gstart + 0.618*(gend -gstart);
	if(energy_from_g(x1) < energy_from_g(x2) )
	{
	    gend = x2;  
	}
	else
	{
	   gstart = x1;
	}

    }
    double val = energy_from_g(gstart) ;
    return std::pair<double,double>(gstart, val);

}

double VarHam::energy_from_g(double gval)
{
    _redbcs->set_g(gval);
    _redbcs->calculate_energy();
    return  _redbcs->get_energy(*_chemham);
}

double VarHam::energy_from_epsilon(double epsilon , int which)
{
    _redbcs->set_epsilon(epsilon , which);
    _redbcs->calculate_energy();
    return  _redbcs->get_energy(*_chemham);
}

double VarHam::get_ci_energy()
{
    return _redbcs->get_energy(*_chemham);
}

double VarHam::solve()
{
    int max_iter = 120;
    int it = 0;
    double conv = 0;

    while(it < max_iter && conv <  get_l())
    {
        conv = 0;
        std::pair<double,double> g_opt = optimize_g(1e-6);
	std::cout << "Energy after g opt. in solve : " << g_opt.second << " g is : " << g_opt.first << std::endl;
        for(int epsnum = 0 ; epsnum < get_l()-1 ; epsnum++)
        {
            double check = 0.005;
            double epsdif = 0.1;
            //epsdif /= ((it+1) );
	    double gend = 0;
	    double gstart = 0;
            double cur_epsilon = _redbcs->get_epsilon(epsnum);
	    double cur_e =0;
            if( energy_from_epsilon(cur_epsilon - check, epsnum) <  energy_from_epsilon(cur_epsilon, epsnum) -1e-7 )
            {   
	        double cur_e = energy_from_epsilon(cur_epsilon - check, epsnum);
		double prev_e = energy_from_epsilon(cur_epsilon - epsdif, epsnum);
                while(prev_e <  cur_e)
		{
                	epsdif += 0.01;
		        std::cout << "Epsdif: " <<  epsdif <<std::endl;
		        std::cout << "E: " <<  energy_from_epsilon(cur_epsilon - epsdif, epsnum) <<std::endl;
			if( epsdif > 1. )
			{
			   if( energy_from_epsilon(cur_epsilon - epsdif, epsnum) > prev_e -1e-4)
				    break;	
			}
			prev_e = energy_from_epsilon(cur_epsilon - epsdif, epsnum);
	        }
                
                gstart = cur_epsilon  - epsdif;
                gend = cur_epsilon ;
            }
            else if( energy_from_epsilon(cur_epsilon + check, epsnum) <  energy_from_epsilon(cur_epsilon, epsnum) -1e-7)
            {
	        double cur_e = energy_from_epsilon(cur_epsilon + check, epsnum);
		double prev_e = energy_from_epsilon(cur_epsilon + epsdif, epsnum);
                while(prev_e <  cur_e)
		{
                	epsdif += 0.01;
		        std::cout << "Epsdif: " <<  epsdif <<std::endl;
		        std::cout << "E: " <<  energy_from_epsilon(cur_epsilon + epsdif, epsnum) <<std::endl;
			if( epsdif > 1.)
			{
			   if( energy_from_epsilon(cur_epsilon + epsdif, epsnum) > prev_e -1e-4 )
				    break;	
			}
			prev_e = energy_from_epsilon(cur_epsilon + epsdif, epsnum);
	        }
		
                gstart = cur_epsilon;
                gend = cur_epsilon + epsdif;
            }
            else
            {
                conv += 1; //One epsilon converged, looking left and right doesn't reduce the energy.
		std::cout << "Augmented conv: " <<  conv <<std::endl;
                continue;
            }
            if(epsdif >= 1. &&  get_ci_energy() < g_opt.second - 1e-3) 
	    {
	    	break;
	    }

            int numit = std::ceil(std::log(1e-6/(epsdif))/std::log(0.618));
	    std::cout << "started linesearch" <<  std::endl;
            for(int i = 0 ; i < numit ; i++)
            {
                double x1 = gend  - 0.618*epsdif;
                double x2 = gstart + 0.618*epsdif;

                if(energy_from_epsilon(x1, epsnum) < energy_from_epsilon(x2, epsnum) )
                {
                    gend = x2;  
                }
                else
                {
                   gstart = x1;
                }
            }
	    std::cout << "ended linesearch" <<  std::endl;
            std::pair<double,double> g_opt = optimize_g(1e-5);
        }
	//Energy doesn't get any lower.
	std::cout << get_ci_energy() << std::endl;
	if( fabs(g_opt.second - get_ci_energy() )< 1e-6)
            break;
	std::cout << "iteration: " << it << " finished. " << std::endl;
	it +=1; //One iteration finished.
    }
    std::pair<double,double> g_opt = optimize_g(1e-6);
    return g_opt.second;
}

bool VarHam::is_solved()
{
    try 
    {
        _redbcs->get_energy();
	return true;
    }
    catch(...)
    {
    	std::cout << "VarHam is associated to a model Hamiltonian that is not solved yet." << std::endl;
	return false;
    }

}

double VarHam::get_mod_ham_energy()
{
    return _redbcs->get_energy();
}

std::string VarHam::get_info()
{
    return _chemham->get_info() + _redbcs->get_information();
}






