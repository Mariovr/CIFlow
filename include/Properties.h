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

#ifndef __Properties__
#define __Properties__

#include <vector>
#include <string>
#include <memory>

class Permutator;
class Permutator_Bit;
class Permutator_File;
class CIMethod;


//This class calculates properties of CI wavefunctions.
class Properties
{
    public:
        //When a full copy of the cimethod is not necessary to extract the properties, and one only needs the eigenvector (shannon entropy, mpd's, ...)
        Properties(std::vector<double> eigvec, int norb , Permutator * perm);
        //Extract properties from the solutions this CIMethod contains.
        Properties(CIMethod * cim);
        double shannon_ic();
        std::pair<std::string , double >  get_max_det();
        std::string get_property(std::string prop); 

    private:
        std::vector<double> _wf;
        Permutator *  _perm;
        CIMethod * _cim;
};







#endif 
