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
#define BOOST_TEST_MODULE Testmultipleprecision
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cstdlib>
#include "Hamiltonian.h"
#include "Permutator.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"
#include "Options.h"


BOOST_AUTO_TEST_CASE(test_multipleprecision91orbitals_energy){
    #ifdef TL
	std::cout<< "Recompile the binaries with TLL defined in include/Options.h"  << std::endl;
	std::cout<< "This because we need at least a 128 bit integer type to perform this test, because we use 70 orbitals (lowest 70 of the 91 HF cc-pv5z)."  << std::endl;
	exit(1);
    #endif
    Hamiltonian * ham = new Hamiltonian("../data/mpbe5zmatrixelements.dat");
    FCI_File * cim = new FCI_File(ham,  "../data/mpbe5zcut70determinants.dat");
    //cim->print_dets();
    cim->solve();
    BOOST_CHECK_CLOSE( -14.346369472944886, cim->get_ci_energy(), pow(10,-10) );
    delete ham;
    delete cim;
}

