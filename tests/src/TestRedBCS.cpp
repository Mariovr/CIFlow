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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE TestRedBCS
#include <iostream>
#include <sys/stat.h>
#include <assert.h>
#include <fstream>
#include <boost/test/unit_test.hpp>

#include "Hamiltonian.h"
#include "CIDens.h"
#include "RedBCS.h"
//compile with icpc -std=c++11 -I../include -I/opt/local/include -o redbcs.x redbcs.cpp -L../lib -lciflow -L/usr/lib -lhdf5 -lmpi -lmpi_cxx -lboost_unit_test_framework -lm -larpack -mkl -openmp

BOOST_AUTO_TEST_CASE( test_RedBCS_VarRG_energy)
{
    string matrixelements = "../data/psi0_6-31g0.60.mout";
    //Create objects to do the work, and execute the work.
    Hamiltonian ham = Hamiltonian(matrixelements);
    //std::vector<double > epsilon=  {-1.34976654, -0.51770481, -0.26592237,  0.33503513};
    //std::vector<int> permarray = { 0 , 2, 1, 3};
    //RedBCS red = RedBCS(ham, 0.119034351059 , epsilon, permarray);
    RedBCS red { ham } ;
    red.set_g(0.119034351059);
    double e_red = red.calculate_energy();
    double e_var = red.get_energy(ham);
    std::cout <<  "Red BCS energy: " << e_red << " var RG energy: " << e_var <<std::endl;

    BOOST_CHECK_CLOSE(-2.5978142225  , e_red , pow(10,-7));
    BOOST_CHECK_CLOSE(-1.1247132522  , e_var , pow(10,-7));
    
    //red.get_density()->print_one_dens(std::cout);
    //red.get_density()->print_two_dens(std::cout);
    //std::vector<double> epsilon = red.get_epsilon();
    //std::vector<int> permarray = red.get_permarray();
    //std::cout << "g " << red.get_g() <<std::endl;
    //for(int i = 0; i < ham.getL() ; i++)
    //{
    //    std::cout << epsilon[i] << " " << permarray[i]  <<std::endl;
    //}
}
