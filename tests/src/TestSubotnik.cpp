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
#define BOOST_TEST_MODULE TestSubotnik
#include <boost/test/unit_test.hpp>
#include "Hamiltonian.h"
#include "Permutator.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"

BOOST_AUTO_TEST_CASE(test_iterative_subotnik){
    Hamiltonian * ham = new Hamiltonian("../data/beh2_sto3g_1_5.dat");
    FCI * cim = new FCI(ham);
    Iterative_Subotnik * min_sen = new Iterative_Subotnik(cim  , 1e-5, 0);
    min_sen->optimize();
    DOCI* cim2 = new DOCI(min_sen->get_ham()); //Calculates DOCI energy in Mmin basis
    cim2->solve();
    BOOST_CHECK_CLOSE(-15.307318, cim2->get_ci_energy(), pow(10,-4) );
    std:: cout << "DOCI energy in FCI mmin (subotnik ): " << cim2->get_ci_energy() <<std::endl;
    delete min_sen;
    delete cim;
    delete cim2;
    delete ham;
}

