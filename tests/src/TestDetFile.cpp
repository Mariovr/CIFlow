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
#define BOOST_TEST_MODULE TestDetFile
#include <boost/test/unit_test.hpp>
#include "Hamiltonian.h"
#include "CIMethod.h"

BOOST_AUTO_TEST_CASE( test_det_file)
{
 
  // Be_ccpVDZ test case with determinant list.
  Hamiltonian * ham = new Hamiltonian("../data/Be_ccpVDZ.dat");
  FCI_File * cim = new FCI_File(ham, "../data/GKCI_Be_cntSp_cc-pVDZ_0100.txt");

  cim->solve();

  BOOST_CHECK_CLOSE(-14.5951673370889, cim->get_ci_energy() , pow(10,-10));
  delete cim;
  delete ham;
}

