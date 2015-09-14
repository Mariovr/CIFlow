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
#define BOOST_TEST_MODULE TestER
#include <boost/test/unit_test.hpp>
#include "iostream"
#include "Hamiltonian.h"
#include "Permutator.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"

BOOST_AUTO_TEST_CASE(test_edmiston_ruedenberg){
    Hamiltonian * ham = new Hamiltonian("../data/but_2_5.dat");
    DOCI * cim = new DOCI(ham);
    cim->solve();
    std::cout << "DOCI E in c1 MO: " << cim->get_ci_energy() <<std::endl;
    //Iterative_Subotnik * er_calc = new Iterative_Subotnik_DIIS(cim  , 1e-5, 1);
    Iterative_Subotnik * er_calc = new Iterative_Subotnik(cim  , 1e-5, 1);
    er_calc->optimize();
    std::cout << "DOCI E in Edmiston Ruedenberg orbitals: " << er_calc->get_opt_ci_energy() <<std::endl;
    //SimulatedAnnealing * sim_anneal = new SimulatedAnnealing(cim  , 0.1 , 0.99 , 1.3 , 0.999 );
    //sim_anneal->run_multiple(3, false);
    //std::cout << "DOCI E after sim an starting from ER orbitals: " << sim_anneal->get_opt_ci_energy() <<std::endl;
    //BOOST_CHECK_CLOSE(-15.307318, sim_anneal->get_opt_ci_energy() , pow(10,-4) );
    delete er_calc;
    //delete sim_anneal;
    delete cim;
    delete ham;
}
