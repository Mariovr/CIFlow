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

#include <string>
#include <memory>

#ifndef HAMCONSTRUCT_H
#define HAMCONSTRUCT_H


class Hamiltonian;

class HamConstruct
{
    public:
        std::shared_ptr<Hamiltonian> generate_ham(std::string filename);
        std::string get_type(std::string filename);

    private:
        std::string _filename;

};











#endif
