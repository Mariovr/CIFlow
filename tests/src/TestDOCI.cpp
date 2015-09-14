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
#define BOOST_TEST_MODULE TestDOCI
#include <boost/test/unit_test.hpp>
#include "Hamiltonian.h"
#include "Permutator.h"
#include "CIMethod.h"

BOOST_AUTO_TEST_CASE( test_DOCI_energy)
{
  
    Hamiltonian * ham = new Hamiltonian("../data/beh2_sto3g_1_5.dat");
    DOCI * cim = new DOCI(ham);
    cim->solve();

    BOOST_CHECK_CLOSE(-15.225362047972723 , cim->get_ci_energy() , pow(10,-10));
    delete ham;
    delete cim;
}

