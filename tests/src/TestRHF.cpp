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
#define BOOST_TEST_MODULE TestSCF
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "Scf.h"

BOOST_AUTO_TEST_CASE( test_scf_energy)
{
    //Hamiltonian * ham = new Hamiltonian("../data/h2otestpsi4site.dat"); //uses symmetric integrals
    //ham->load_overlap("../data/overlaph2o.dat");
    Hamiltonian * ham = new Hamiltonian("../data/atomicintegralsbeh2_sto3g.h5"); //uses symmetric integrals
    ham->load_overlap("../data/overlap.txt");
    //Hamiltonian * ham = new Hamiltonian("../data/beh2atomic.dat"); //uses c1 integrals
    //ham->load_overlap("../data/beh2atomicoverlap.txt");

    //Test overlap functionality.
    //std::vector<double> vec = ham->get_overlap();
    //matrix mat = matrix(vec, ham->getL() , ham->getL());
    //mat.Print();

    std::cout << "We do RHF without DIIS " << std::endl;
    RHF rhf(ham ,1000 , 1e-10 , 1e-12 , 1 ,0);
    double ehf = rhf.scf();

    std::cout << "We do RHF with DIIS " << std::endl;
    RHF rhf2(ham ,1000 , 1e-10 , 1e-12 , 1 ,1);
    ehf = rhf2.scf();

    BOOST_CHECK_CLOSE(-15.12681016931100863 , ehf , pow(10,-7));//For the beh2 tests

    //Remark the HF determinants are not dependend on unitary rotations of the occupied with occupied and nonoccupied with nonoccupied.
    OptIndex opt =  ham->get_index_object();
    UnitaryMatrix unit(opt);
    unit.load_unitary("../data/unitary-mounit_ao_to_mo.txt");
    UnitaryMatrix unit2 = rhf.get_transformation();
    std::cout << "Mo from psi4." << std::endl;
    unit.print_unitary(std::cout );
    std::cout << "Mo from ciflow." << std::endl;
    unit2.print_unitary(std::cout );
    //BOOST_CHECK_CLOSE(unit.check_difference(unit2) , 1e-10 , pow(10,-6));
    delete ham;
}

