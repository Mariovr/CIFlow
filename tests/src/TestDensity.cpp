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
#define BOOST_TEST_MODULE TestDensity
#include <boost/test/unit_test.hpp>
#include <iostream>
#include "Hamiltonian.h"
#include "CIMethod.h"
#include "CIDens.h"
#include "OrbitalTransform.h"
#include "UnitaryMatrix.h"


BOOST_AUTO_TEST_CASE( test_density)
{
 
  Hamiltonian * ham = new Hamiltonian("../data/beh2_sto3g_1_5.dat");

  FCI_File * cim = new FCI_File(ham, "../data/fci_beh2_sto-3gdeterminants.dat");
  cim->solve();
  CIDens * densmatfile = new DensFILE(cim);
  densmatfile->construct_density();
  //densmatfile->print_one_dens(std::cout);
  //densmatfile->print_two_dens(std::cout);
  //densmatfile->get_NO();
  std::cout << "FILE-DENS one-electron energy: " << densmatfile->get_one_electron_energy() << std::endl;
  std::cout << "FILE-DENS two-electron energy: " << densmatfile->get_two_electron_energy() << std::endl;
  std::cout << "FILE-DENS total energy: " << densmatfile->get_dens_energy() << std::endl;

  FCI * fci = new FCI(ham );
  fci->solve();
  CIDens * densmatfci = new DensFCI(fci);
  densmatfci->construct_density();
  std::cout << "FCI-DENS one-electron energy: " << densmatfci->get_one_electron_energy() << std::endl;
  std::cout << "FCI-DENS two-electron energy: " << densmatfci->get_two_electron_energy() << std::endl;
  std::cout << "FCI-DENS total energy: " << densmatfci->get_dens_energy() << std::endl;
  //densmatfci->print_one_dens(std::cout);
  //densmatfci->print_two_dens(std::cout);
  densmatfci->compare_two_dens(densmatfile);
  densmatfci->compare_one_dens(densmatfile);


  BOOST_CHECK_CLOSE(densmatfci->get_dens_energy(), fci->get_ci_energy() , pow(10,-10));

  delete densmatfile;
  delete densmatfci;
  delete cim;
  delete fci;
  delete ham;
}
