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

#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include "DIIS.h"

DIIS::DIIS(int max)
{
    _errorlist.resize(0);
    _vallist.resize(0);
    _coef.reset(new double[100]);
    _max_size = max;
    _print = 1;
    std::cout << "We created DIIS, with " << max << " errormatrices kept in memory." << std::endl;
}

DIIS::~DIIS()
{
    for(int i = 0 ; i < _errorlist.size() ; i++)
    {
    	delete _errorlist[i];
	    delete _vallist[i];
    }
}

void DIIS::add_val_matrix(const matrix & val)
{
    _vallist.push_back(new matrix(val));
    if(_vallist.size() > _max_size)
    {
        delete _vallist[0];
        _vallist.erase(_vallist.begin());
    }
}

void DIIS::add_error_matrix(const matrix & R)
{
    _errorlist.push_back(new matrix(R));
    if(_errorlist.size() > _max_size)
    {
        delete _errorlist[0];
        _errorlist.erase(_errorlist.begin());
    }
}

void DIIS::add_matrices(const matrix & val , const matrix & error)
{
    add_val_matrix(val);
    add_error_matrix(error);
}

void DIIS::construct_diis_matrix(matrix & diismat)
{
    //Remark the DIIS matrix is by definition symmetric, when each error matrix is added we only have to calculate the extra columns.
    //
    for(int i = 0 ; i < get_current_size() ; i++)
    {
        for(int j = i ; j < get_current_size() ; j++)
        {
            diismat(i,j) = _errorlist[i]->matrix_trace(*_errorlist[j], true);
	    if(i != j)
                diismat(j,i) = diismat(i,j);
	}
    }	
    for(int i = 0 ; i < get_current_size() ; i++)
    {
        //_diisb->set(get_current_size() ,i, -1.);  
        //_diisb->set(i,get_current_size() , -1.);  
        diismat(get_current_size() ,i) = -1.;
        diismat(i,get_current_size() ) = -1.;
    }
    diismat(get_current_size() , get_current_size() ) = 0.; 
}

int DIIS::solve_diis(double * coef)
{
    matrix diissystem = matrix(get_current_size() +1, get_current_size() +1);
    construct_diis_matrix(diissystem);
    //diissystem.Print();
    std::fill(coef, coef+get_current_size() +1, 0.);
    coef[get_current_size() ] = -1.;
    int info = diissystem.linsolve(coef , coef); //Solves the linear system defined by diissystem * x = coef, and saves the solution x in coef (Remark the coef rhs gets overwritten).
    return info;
}

void DIIS::do_diis(const matrix & val , const matrix & er)
{
    add_matrices(val , er); //This is done separately in the routine that contains the diis.
    _coef.reset( new double[get_current_size() +1]() ) ; 
    int info =  solve_diis(_coef.get());

    if(_print > 2)
    {
        for(int i = 0 ; i < get_current_size()   ; i++)
        {
	    std::cout << "coef " << _coef[i] << " " << std::endl       ;
        }
    }	    
    if (info != 0)
    {
        std::cout << "Something went wrong in DIIS, do_diis: info = " << info << std::endl;
	exit(EXIT_FAILURE);
    }
}

matrix DIIS::get_extrapolated( std::vector<matrix *> & matrixlist)
{
    //We only extrapolate with the get_step last matrices in matrixlist.
    matrix mat(matrixlist[0]->getn() , matrixlist[0]->getm()); 
    for(int i = 0 ; i < get_current_size()   ; i++)
    {
        mat.add(*matrixlist[matrixlist.size()-get_current_size()  +i], _coef[i]);
    }
    return mat;
}

matrix DIIS::get_extrapolated()
{
    //We only extrapolate with the get_step last matrices in matrixlist.
    matrix mat(_vallist[0]->getn() ,_vallist[0]->getm()); 
    assert(_vallist.size() == get_current_size() );
    for(int i = 0 ; i < get_current_size()   ; i++)
    {
        mat.add(*_vallist[i], _coef[i]);
    }
    return mat;
}

std::vector<double> DIIS::get_coef()
{

    std::vector<double> coefs(get_current_size() +1 );
    for(int i = 0 ; i < get_current_size() +1 ; i++)
    {
    	coefs[i] = _coef[i];
    }
    return coefs;
}

void DIIS::reset()
{
    _errorlist.resize(0);
    _vallist.resize(0);
    _coef.reset();
}




