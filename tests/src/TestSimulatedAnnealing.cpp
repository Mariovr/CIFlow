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
#define BOOST_TEST_MODULE TestSimulatedAnnealing
#include <boost/test/unit_test.hpp>
#include "Hamiltonian.h"
#include "Permutator.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"

BOOST_AUTO_TEST_CASE(test_SimulatedAnnealing_energy){
    Hamiltonian * ham = new Hamiltonian("../data/beh2_sto3g_1_5.dat");
    DOCI * cim = new DOCI(ham);
    SimulatedAnnealing * sim_anneal = new SimulatedAnnealing(cim  , 0.1 , 0.99 , 1.3 , 0.999 );
    sim_anneal->run_multiple(3, false);
    BOOST_CHECK_GE(-15.328, sim_anneal->get_opt_ci_energy() );
    delete sim_anneal;
    delete ham;
    delete cim;
}

