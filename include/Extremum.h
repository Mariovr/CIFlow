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
#ifndef EXTREMUM_H
#define EXTREMUM_H

#include <string>
#include <memory>

class Hamiltonian;

class Extremum
{
    public:
        Extremum(std::string filename);

    private:
        Hamiltonian * _ham;

};





#endif //EXTREMUM_H
