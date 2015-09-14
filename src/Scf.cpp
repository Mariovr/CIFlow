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

#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdexcept>
#include <exception>

#include "Options.h"
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "matrix.h"
#include "Scf.h"
#include "DIIS.h"

using namespace std;

RHF::RHF(Hamiltonian * ham, int max_iter , double e_conv , double d_conv , int print , int diis)
{
    _ham.reset(new Hamiltonian(*ham) );
    _max_iter = max_iter; 
    _e_conv = e_conv;
    _d_conv = d_conv ;
    _print = print;

    //Density matrix. (C * C^T)
    _P = matrix(ham->getL() , ham->getL() );
    //Initialize th-> symmetrizat->on transformation of the atomic orbitals.
    _X = matrix(ham->getL() , ham->getL() );
    _C = matrix(ham->getL() , ham->getL() );
    //Initialize matrix that saves the orbital energies.
    _orbitale = matrix(_ham->getL() , 1 );

    _orbarray = vector<int>(_ham->getL(), 0.);

    //Initialize the Fock matrix from Core Hamiltonian guess
    _Fock = matrix(ham->getL() , ham->getL() );
    _Fockt = matrix(ham->getL() , ham->getL() );
    for (int h=0; h< ham->getL(); ++h)
        for (int l= 0 ; l < ham->getL(); ++l)
	    _Fock(h,l) = ham->getTmat(h,l);
    
   if(print > 1)
   {
       cout << "Fock matrix from core hamiltonian guess." << endl;
       _Fock.Print();
   }

   vector<double> overl = ham->get_overlap();
   _overlap = matrix(overl, ham->getL() , ham->getL());
   symmetric_diag(); //Sets the _X matrix.
   if(_print > 1)
   {
       cout << "_x : symmetric diagonalization matrix. " << endl;
       _X.Print();
   }

    _Fockt.transform(_Fock, _X);

    matrix evecs = matrix(ham->getL() , ham->getL() );
    diag_symmetric_matrix(_Fockt, _orbitale, evecs ); //The first guess for the MO comes from diagonalizing the transformed Fock matrix, that comes from the core hamiltonian.
    //_Fockt.diagonalize(_orbitale, evecs ); //The first guess for the MO comes from diagonalizing the transformed Fock matrix, that comes from the core hamiltonian.

    _C.prod( _X, evecs, 'N' , 'N'); //The C is the transformation from the nonorthognonal atomic orbitals, to the first HF mo guess so first X (to orthogonalize), and then rightmultiply with Evecs (right multiply because new basis is saved in columns.)
    construct_density();

    _dodiis = diis;
    if (diis)
        _diis = new DIIS(8);

    if(_print > 1){
        printf("Transformed Fock matrix.:\n");
        _Fockt.Print();
        printf("MO Coefficients and density from Core Hamiltonian guess:\n");
        _C.Print();
        _P.Print();
    }
}

void RHF::symmetric_diag()
{
    matrix eigvec = matrix(_ham->getL(), _ham->getL() );
    matrix eigval = matrix(_ham->getL(), 1 );
    _overlap.diagonalize(eigval, eigvec);//vectors are in column format.
    // Convert the eigenvales to 1/sqrt(eigenvalues)
    double min_S = fabs(eigval(0,0));
    for (int h=0; h< _ham->getL(); ++h)
    {
            if (min_S > eigval(h,0))//eigvals are always positive but watch out for to small eigenvalues which point to linear dependencies.
                min_S = eigval(h,0);
            double scale = 1.0 / sqrt(eigval(h, 0));
            eigval(h, 0) = scale;
     }

    cout << "Lowest eigenvalue of overlap S = " <<  min_S << endl;

    //If symmetric diagonalization of the overlap matrix didn't work, start with canonical diagonalization.
    //REMARK canonical orthonormalization is not working yet, TODO: debug this.
    if(min_S < 1e-10 )
    {
        cout << "WARNING: Min value of overlap below treshold!!!!" << endl;
        cout << "So we use canonical orthogonalization to create the orthonormal single-particle orbitals.!!!!" << endl;

        _overlap.diagonalize(eigval, eigvec);//vectors are in column format. reset eigenvectors and overlaps.
        int delta_mos = 0;
        for (int h=0; h < _ham->getL(); ++h) 
	{
            //in each irrep, scale significant cols i  by 1.0/sqrt(s_i)
            if ( eigval(h,0) > 1e-9) 
	    {
                double scale = 1.0 / sqrt(eigval(h, 0));
                eigvec.scale_column(h,  scale);
            } 
	    else 
	    {
                delta_mos++;
            }
        }
        if (_print>2)
	    cout << delta_mos << " of " <<  _ham->getL() << " possible MOs eliminated.\n";

        //Initialize th-> symmetrizat->on transformation of the atomic orbitals.
        _X = matrix(_ham->getL(), _ham->getL()-delta_mos );
	int xvec = 0;
        for (int h=0; h< _ham->getNirreps(); ++h) {
            //Copy significant columns of eigvec into X in
            //descending order, remarkt lapack diagonalization ordered the eigenvectors according to increasing eigenvalue (lowest first).
            int start_index = 0;
            for (int i=0; i< _ham->getNORB(h); ++i) {
                if ( 1e-9< eigval(h,0)) 
		{
                } 
		else 
		{
                    start_index++;
                }
            }
            for (int i=0; i< _ham->getNORB(h) -start_index; ++i) 
	    {
                for (int m = 0; m < _ham->getL(); m++)
                    _X(m ,xvec ) = eigvec(m,_ham->getNstart(h+1)-i-1);
		xvec ++;
            }
        }

        //Density matrix. (C * C^T)
        _P = matrix(_ham->getL() -delta_mos, _ham->getL() -delta_mos );
        _C = matrix(_ham->getL() , _ham->getL() -delta_mos);
        //Initialize matrix that saves the orbital energies.
        _orbitale = matrix(_ham->getL()-delta_mos , 1 );
        _orbarray = vector<int>(_ham->getL()-delta_mos, 0.);

        //Initialize the Fock matrix from Core Hamiltonian guess
        _Fock = matrix(_ham->getL()-delta_mos , _ham->getL()-delta_mos );
        _Fockt = matrix(_ham->getL()-delta_mos , _ham->getL()-delta_mos );
    
    }

    matrix temp2 = matrix(eigvec.getn(), eigvec.getm() );//initialize all elements to zero.
    // Create a vector matrix from the converted eigenvalues
    temp2.set_diagonal(eigval);

    _X.prod(temp2, eigvec, 'N' , 'T');
    temp2.prod(eigvec, _X , 'N' , 'N');
    _X = matrix(temp2);
}

UnitaryMatrix RHF::get_transformation()
{
    OptIndex opt =  _ham->get_index_object();
    UnitaryMatrix unit = UnitaryMatrix(_C, opt);
    return unit;
}

void RHF::diag_symmetric_matrix(matrix Fockt, matrix & eigval , matrix & eigvec)
{
    for(int irrep = 0 ; irrep < _ham->getNirreps() ; irrep++){
        int nstart =_ham->getNstart(irrep);
        int norb = _ham->getNORB(irrep);
        //handles the case if there are no orbitals associated to a given irrep.
        if( norb > 0){
            matrix part = matrix(norb ,norb );
            matrix val_part = matrix(norb ,1 );
            for(int i = 0 ; i < norb; i++){
                for(int j = i ; j < norb; j++){
                    part(i,j) = Fockt(i+nstart,j+nstart);
                    if(j!= i)
                        part(j,i) = Fockt(j+nstart,i+nstart);
                }
            }
            part.diagonalize(val_part,part);
            if(_print > 3 )
	    {
                part.Print();
		val_part.Print();
            }
            for(int i = nstart ; i < norb+nstart; i++){
                eigval(i,0) = val_part(i-nstart , 0);
                for(int j = nstart ; j < norb+nstart ; j++){
                    eigvec(j,i) = part(j -nstart , i - nstart);
                }
		//eigval.Print();
		//eigvec.Print();
            }

        } //end if norb > 0.
    }//end for irrep

    sort_energies(eigval); //To know the indices of the nup orbitals with the lowest orbital energies.
    //This is necessary when you diagonalize in irreps, because the orbitals are only ordered according to lowest fock matrix energy in an irrep.
}


//Only for one dimensional matrices.
void RHF::sort_energies(matrix & orbe){
    // initialize original index locations
    for (int i = 0; i != _orbarray.size(); ++i) _orbarray[i] = i;
    
    // sort indexes based on comparing values in v
    sort(_orbarray.begin(), _orbarray.end(),
                 [&orbe](size_t i1, size_t i2) {return orbe[i1] < orbe[i2];});
    
}



void RHF::construct_density()
{
    for(int p = 0; p < _ham->getL(); ++p){
        for(int q = 0; q < _ham->getL(); ++q){
            double val = 0.0;
            for(int i = 0; i < _ham->getnup(); ++i){
                val += _C(p, _orbarray[i]) * _C(q, _orbarray[i]);
            }
            _P(p, q ) =  val;
        }
    }   
}

double RHF::get_hf_energy()
{
    double som = 0;
    for (int h=0; h < _ham->getL(); ++h)
        for (int l= 0 ; l < _ham->getL(); ++l)
            som += _P(h,l) * ( _ham->getTmat(l,h) + _Fock(l, h) ) ;
    
    return som;
}

matrix RHF::get_orbital_energies(int print )
{
    matrix evecs = matrix(_ham->getL() , _ham->getL() );
    diag_symmetric_matrix(_Fockt, _orbitale, evecs ); 
    if( print)
        _orbitale.Print();
    return _orbitale;
}

double RHF::scf()
{
   int iter = 1;
   bool converged = false;
   double e_old;
   double e_new = _ham->getEconst() + get_hf_energy();

   printf("\tEnergy from core Hamiltonian guess: %20.16f\n\n", e_new);

   printf("\t*=======================================================*\n");
   printf("\t* Iter       Energy            delta E    ||gradient||  *\n");
   printf("\t*-------------------------------------------------------*\n");

   matrix evals = matrix(_ham->getL() , 1 );
   matrix evecs = matrix(_ham->getL() , _ham->getL() );

   matrix Temp1 = matrix(_ham->getL() , _ham->getL() );
   matrix SPF = matrix(_ham->getL() , _ham->getL() );
   matrix FPS = matrix(_ham->getL() , _ham->getL() );

   while(!converged && iter < _max_iter){
       e_old = e_new;

       // Add the core Hamiltonian term to the Fock operator
       for (int h=0; h< _ham->getL(); ++h)
           for (int l= 0 ; l < _ham->getL(); ++l)
	           _Fock(h,l) = _ham->getTmat(h,l);

       // Add the two electron terms to the Fock operator
       for(int p = 0; p < _ham->getL(); ++p){
           for(int q = 0; q < _ham->getL() ; ++q){
               double J = 0.0;
               double K = 0.0;
               for(int r = 0; r < _ham->getL(); ++r){
                   for(int s = 0; s < _ham->getL() ; ++s){
                       J += _ham->getVmat(p, s, q, r ) * _P(r, s); //Hamiltonian is in physicist notation.
                       K += _ham->getVmat(p,   q , r ,s  ) * _P(r, s);

                   }
               }
	       //cout << " i " << p << " j " << q << " J" << J <<"K"  << K << endl;
               _Fock(p, q) += (2.0 * J - K);
           }   
       }

       // Transform the Fock operator and diagonalize it
       _Fockt.transform(_Fock, _X);
       //_Fockt.diagonalize(_orbitale, evecs ); //The first guess for the MO comes from diagonalizing the transformed Fock matrix, that comes from the core hamiltonian.
       diag_symmetric_matrix(_Fockt, _orbitale, evecs ); 

       _C.prod( _X, evecs, 'N' , 'N'); 


       // Compute the energy
       e_new = _ham->getEconst() + get_hf_energy();
       double dE = e_new - e_old;
 
       // Compute the orbital gradient, FDS-SDF
       Temp1.prod(_P, _overlap, 'n' , 'n');
       FPS.prod( _Fock, Temp1, 'n' , 'n');
       Temp1.prod(_P, _Fock, 'n' , 'n');
       SPF.prod(_overlap, Temp1, 'n' , 'n' );
       //FPS.Print();
       //SPF.Print();
       Temp1= matrix(FPS);
       Temp1.subtract(SPF);
       //SPF.prod(Temp1, _X , 'n' ,'t');
       //Temp1.prod(_X , SPF , 'n' ,'n');
       //Temp1.Print();

       // Rebuild the density using the new orbitals. After construction of orbital gradient, because it presumes the density of the orbitals of wich the fock matrix is constructed.
       construct_density();

       if(_dodiis)
       {
           if(iter <= 2)
	       _diis->add_matrices(_Fock,Temp1);
           if(iter > 2) 
	   {
               _diis->do_diis(_Fock, Temp1);
    	       _Fock = _diis->get_extrapolated();
    
               _Fockt.transform(_Fock, _X);
               diag_symmetric_matrix(_Fockt, _orbitale, evecs ); 
               _C.prod( _X, evecs, 'N' , 'N'); 
               construct_density();
               //e_new = _ham->getEconst() + get_hf_energy();
               //dE = e_new - e_old;
	   }
       }
       double dRMS = Temp1.rms();

       if(_print > 1){
           _Fockt.Print();
           evecs.Print();
           _orbitale.Print();
           _C.Print();
           _P.Print();
           FPS.Print();
           SPF.Print();
       }
     
       converged = (fabs(dE) < _e_conv) && (dRMS < _d_conv);

       printf("\t* %3d %20.14f    %9.2e    %9.2e    *\n", iter, e_new, dE, dRMS);
       iter++;
   }
   printf("\t*=======================================================*\n");

   if(!converged)
       throw std::runtime_error("The SCF iterations did not converge.");

   cout << "Orbital Energies" << endl;
   _orbitale.Print();
   cout << "HF energy: " << e_new <<  endl;
   return e_new;
}

