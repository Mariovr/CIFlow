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

#include <sys/stat.h>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>


#include "Properties.h"
#include "Permutator.h"
#include "CIMethod.h"
#include "scpp_assert.h"


using namespace std;

//Properties::Properties(vector<double> eigvec, int norb , Permutator * perm)
//{
//    //When a full copy of the cimethod is not necessary to extract the properties, and one only needs the eigenvector (shannon entropy, mpd's, ...)
//    _perm = perm;
//    _wf = eigvec;
//    _cim = nullptr;
//}

Properties::Properties(CIMethod * cim)
{
    _cim = cim;
    _cim->construct_density();
}

double Properties::shannon_ic(int num)
{
    double som = 0;
    std::vector<double> wf = _cim->get_eigvec(num);
    for(auto & x : wf )
    {
        if (fabs(x) > 1e-15)
            som += x*x*log(x * x)/log(2.); // log_2 (x) = log(x) / log(2) 
    }
    return -1.* som;
}

std::pair<std::string , double >  Properties::get_max_det(int num)
{
    std::vector<double> wf = _cim->get_eigvec(num);
    auto val = std::max_element(wf.begin() , wf.end() ); //returns iterator
    int index = std::distance(wf.begin() , val);

    std::pair<std::string , double> values { "d", *val};
    return values;
}

std::string Properties::get_property(std::string prop, int num)
{
    if(prop ==  "shannon")
        return std::to_string(shannon_ic(num) );
    else if(prop == "spin_squared")
        return std::to_string(_cim->get_spin_squared(num));
    else if(prop == "spin")
        return std::to_string(_cim->get_spin(num ));
    else if(prop == "sz")
        return std::to_string(_cim->get_sz());
    else if(prop == "mulliken")
        return std::to_string(_cim->get_mulliken({0,1,2,3,4}, num) ) ;
    else if(prop == "maxdet")
    {
        return "d";
    }
    else
        SCPP_ASSERT(1 , "Property " << prop << " is not one of shannon, spin, spin_squared, sz, mulliken, ..." << std::endl);
    return "";
} 


