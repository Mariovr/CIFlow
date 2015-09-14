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
#define BOOST_TEST_MODULE TestTransformations
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cstdlib>
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"
#include "CIMethod.h"

//using namespace std;


/*  
void transform_unitary(double max_angle, OrbitalTransform * _orbtrans, OptIndex index)
{
    std::random_device rd();
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist_real(0, 1);
    for(unsigned int irrep = 0 ; irrep < index.getNirreps() ; irrep++)
    {	    
        for(int i = 0 ; i < index.getNORB(irrep)-1 ; i++)
	{	
            for(int j = i+1 ; j < index.getNORB(irrep) ; j++)
            {
                double cur_angle =  max_angle * (dist_real(mt)-dist_real(mt));
	        _orbtrans->get_unitary().jacobi_rotation(irrep, i + index.getNstart(irrep), j+index.getNstart(irrep), cur_angle);
            }        
	}

    }	
    //_orbtrans->get_unitary().print_unitary();
    _orbtrans->CheckDeviationFromUnitary();
}*/

void transform(OrbitalTransform * _orbtrans){
    _orbtrans->get_unitary().jacobi_rotation(0, 0, 1, 1.2);
    _orbtrans->get_unitary().jacobi_rotation(0, 0, 2, -0.4);
    _orbtrans->get_unitary().jacobi_rotation(5, 3, 4, 3);
    //_orbtrans->get_unitary().print_unitary();
    //_orbtrans->CheckDeviationFromUnitary();
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  main
 *  Description:  Test main for unitary transformations.
 * =====================================================================================
 */
BOOST_AUTO_TEST_CASE( test_transformations)
{
    //This test transforms a Hamiltonian, then transforms it back with the inverse transformation.
    //Checks if the difference of the original with the (transformed and backtransformed) is zero.
    //Initialisation
    Hamiltonian * ham = new Hamiltonian("../data/beh2_sto3g_1_5.dat");
    OptIndex index = ham->get_index_object();
    DOCI * cim = new DOCI(ham);
    //Start solving for the initial values.
    cim->solve();
    std::cout << "E_before transformation: " << cim->get_ci_energy() <<std::endl; 
    //Generate the orbital transformed and calculate new energies.
    OrbitalTransform * orbtrans = new OrbitalTransform(ham);
    transform(orbtrans );
    Hamiltonian *hamorig = new Hamiltonian(*ham);
    orbtrans->fillHamCI(*hamorig);
    cim->set_ham(hamorig);
    cim->solve();
    std::cout << "E_after transformation: " << cim->get_ci_energy() <<std::endl; 
    std::cout << "after transformation the Hamiltonian differs: " << ham->get_difference(hamorig) << std::endl;
    UnitaryMatrix old =orbtrans->get_unitary();

    orbtrans->do_jacobi_rotation(0,1,1., hamorig);
    cim->set_ham(hamorig);
    cim->solve();
    std::cout << "E_after do_jacobi_rotation: " << cim->get_ci_energy() <<std::endl; 

    old.jacobi_rotation(0,0,1,1.);
    orbtrans->set_unitary(old);
    Hamiltonian * ham2 = new Hamiltonian("../data/beh2_sto3g_1_5.dat");
    orbtrans->fillHamCI(*ham2);
    cim->set_ham(ham2);
    cim->solve();

    std::cout << "E_after fillHamCI should be the same as the previous for consistency: " << cim->get_ci_energy() <<std::endl; 
    std::cout << "Hamdifference: " <<hamorig->get_difference(ham2) <<   std::endl; 
    delete ham2;


    //check reset_function and difference between ham and original should be zero.
    UnitaryMatrix inverse = old.get_inverse();
    //cout << "U^T = " <<endl;
    //inverse.print_unitary();
    old.multiply_left_with(inverse);
    //cout << "U^T * U = " <<endl;
    //old.print_unitary();
    orbtrans->set_unitary(old);
    orbtrans->fillHamCI(*hamorig);
    cim->set_ham(hamorig);
    cim->solve();
    std::cout << "E_after inverse transformation: " << cim->get_ci_energy() <<std::endl; 
    std::cout << "Should give back the first energy: " <<std::endl; 

    BOOST_CHECK_CLOSE( 0.0, ham->get_difference(hamorig), pow(10,-10));
    delete ham;
    delete hamorig;
    delete cim;
    delete orbtrans;
}				/* ----------  end of function main  ---------- */
