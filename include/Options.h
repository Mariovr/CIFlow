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
#ifndef OPTIONS
#define OPTIONS

#include <string>
#include <bitset>

//#define _HPC_ //Define this if the compiler has c++11 support and the random headers. Otherwise you need boost. (as the hpc)

using std::string;

const string TMPpath                       = "/tmp/";

const bool Orbopt_debugPrint  = false;
const bool   HAMILTONIAN_debugPrint        = false;
const string HAMILTONIAN_TmatStorageName   = "Ham_Tmat.h5";
const string HAMILTONIAN_VmatStorageName   = "Ham_Vmat.h5"; 
const string HAMILTONIAN_ParentStorageName = "Ham_parent.h5";

const bool CIMethod_debug = false;

#define NEIGVAL 1 //determines the number of eigenvectors to calculate.

#include <boost/multiprecision/cpp_int.hpp>
namespace mp = boost::multiprecision;     // Reduce the typing a bit later...

#define CISIZE 14 //Maximum space for the CImatrix in gigabyte (14 for my ubuntu machine, 80 on delcatty)

#define TL

#if defined TL
	typedef unsigned long TYPE;
//#if defined TL
	//typedef std::bitset<64> TYPE;
#elif defined TLL
	typedef mp::uint128_t TYPE;
#elif defined TLLL
	typedef mp::uint256_t TYPE;
#elif defined TLLLL
	typedef mp::uint512_t TYPE;
#else
	cout << "ERROR: make sure number of orbitals is smaller than 512 _norb is now and define the corresponding type: L ,LL, LLL , LLLL" << endl;

#endif

#endif
