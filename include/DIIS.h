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

#ifndef DIIS_H
#define DIIS_H

#include <vector>
#include "matrix.h"

class DIIS
{
    public:
        DIIS(int max);
     	~DIIS();
     	void add_val_matrix(const matrix & val);
        void add_error_matrix(const matrix & R);
     	void add_matrices(const matrix & val, const matrix & R);
     	void construct_diis_matrix(matrix & diismat);
        void do_diis(const matrix & val , const matrix & er);
        int solve_diis(double * coef);
     	int get_current_size(){return _errorlist.size();}
     	std::vector<double> get_coef(); 
     	matrix get_extrapolated(std::vector<matrix *> & matrixlist);
     	matrix get_extrapolated();
     	int get_max_size(){return _max_size;}
     	void reset();

    private:
        std::vector<matrix *> _errorlist;
	    std::vector<matrix * > _vallist;
	    std::unique_ptr<double []> _coef;
	    int _max_size; 
	    int _print;
};

#endif
