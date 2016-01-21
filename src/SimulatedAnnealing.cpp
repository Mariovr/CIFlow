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
#include <iomanip>
#include <cassert>
#include <stdexcept>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#ifdef __OMP__
#include <omp.h>
#endif

#include "Options.h"
#include "Hamiltonian.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"
#include "CIMethod.h"
#include "SimulatedAnnealing.h"
#include "matrix.h"
#include "CIDens.h"
#include "DIIS.h"
#include "Properties.h"

using namespace std;

OrbitalOptimization::OrbitalOptimization(CIMethod * cim)
{
    _cim = cim;
    // REMARK because optham points to the same adress, as the ham in cim
    // any change of _optham will also affect the ham in cim.
    _optham = new Hamiltonian (*cim->get_ham());

    // we save a copy of the original ham in the orbital transformation object.
    _orbtrans.reset(new OrbitalTransform(cim->get_ham() ) );
    //Check if _cim is already solved.
    if(not _cim->is_solved())
        _cim->solve();
    _opt_energy = _cim->get_ci_energy();
    OptIndex index = _optham->get_index_object();
    _opt_unitary.reset(new UnitaryMatrix(index));

    if(Orbopt_debugPrint)
        assert(fabs(_optham->get_difference(_orbtrans->get_ham())) < 1e-12);
}

//This is commented out because you have to make sure the Hamiltonian object is different.
OrbitalOptimization::~OrbitalOptimization()
{
    if(_optham->get_unitary() )//Check if _optham contains already a unitary transformation. This can be the case if the original Hamiltonian is generated from an old psi4 outputfile.
        _optham->get_unitary()->multiply_left_with(*_opt_unitary);
}

void OrbitalOptimization::set_cim_optimal()
{
    _orbtrans->set_unitary(*_opt_unitary);
    _orbtrans->fillHamCI(*_optham);
    update_cim();
}

void OrbitalOptimization::get_back_original_ham()
{
    _orbtrans->get_unitary().reset_unitary();
    _orbtrans->fillHamCI(*_optham);
    update_cim();
}

void OrbitalOptimization::update_cim()
{
    _cim->set_ham(_optham);
    _cim->solve();
}

void OrbitalOptimization::transform_orbitals(int irrep, int i, int j, double cur_angle)
{
    /*  
    _orbtrans->get_unitary().jacobi_rotation(irrep, i, j, cur_angle);

    if(Orbopt_debugPrint)
    {
        _orbtrans->get_unitary().print_unitary(cout);
        _orbtrans->CheckDeviationFromUnitary();
    }*/
    _orbtrans->do_jacobi_rotation(i,j , cur_angle, _optham);
}

void OrbitalOptimization::transform_withasx(double * step, double temp, double max_angle)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist_real(0, 1);
    for(unsigned int index = 0 ; index < _orbtrans->get_unitary().getNumVariablesX() ;index ++)
    {
        step[index] = max_angle * (dist_real(mt)-dist_real(mt));
    }    
    _orbtrans->update_unitary(step);

    if(Orbopt_debugPrint)
    {
        _orbtrans->get_unitary().print_unitary(cout);
        _orbtrans->CheckDeviationFromUnitary();
    }
}

double OrbitalOptimization::transform_to_no(CIDens & cid)
{
    //This method the Hamiltonian to the no from cim or fci.
    matrix NO = matrix(_cim->get_l() , _cim->get_l());
    matrix occupations = matrix(_cim->get_l() ,1);
    cid.construct_density(false); //False makes sure we don't construct the 2rdm.
    cid.get_no( occupations, NO);
    if(CIMethod_debug){
        NO.Print();
        occupations.Print();
    }	    
    _orbtrans->rotate_active_space_vectors(NO.getpointer());
    _orbtrans->fillHamCI(*_optham);
    _cim->set_ham(_optham);//funny, but necessary because this function also resets the eigenvector and values, and reconstructs the ci matrix.
    _cim->solve();
    //std::cout  << "E in natural orbitals: " << _cim->get_ci_energy() << std::endl;
    return _cim->get_ci_energy();
}

double OrbitalOptimization::get_cim_no()
{
    std::unique_ptr<CIDens> cid = _cim->get_density();
    double noe = transform_to_no(*cid);
    return noe;
}

double OrbitalOptimization::get_fci_no()
{
    FCI * fci = new FCI(_cim->get_ham());
    fci->solve();
    std::cout << "FCI total E: " << fci->get_ci_energy() << std::endl;
    std::unique_ptr<CIDens> densmatfci {new DensFCI(fci)} ;
    double cim_fcino_energy = transform_to_no(*densmatfci);
    std::cout  << "E in FCI natural orbitals: " << cim_fcino_energy << std::endl;
    delete fci; //OPMERKING delete fci always after the density object generated from it.
    return cim_fcino_energy;
}

double OrbitalOptimization::set_basis(CIMethod &cim)
{
    
    cim.set_ham(_optham);
    cim.solve();
    return cim.get_ci_energy();
}

//Saves the Hamiltonian with the matrixelements of the basis that maximizes the quantity you wanted to maximize.
void OrbitalOptimization::save_ham(std::string name)
{
    set_cim_optimal(); //To make sure the optimal hamiltonian is saved.
    std::cout << "we save ham with this E:" << _cim->get_ci_energy() << endl;
    //_optham->save(name);
    _optham->save_file(name);
}

void OrbitalOptimization::save_unitary(std::string name, bool hdf5)
{
    string filename = "unitary_" + _optham->get_short_filename() + _cim->get_perm_info() +_cim->get_name() + name ;
    if( hdf5)
    {
        _opt_unitary->saveU(filename + ".h5");
    }
    else
    {
        ofstream file(filename + ".dat");
        _opt_unitary->print_unitary(file );
        file.close();
    }
}

UnitaryMatrix * OrbitalOptimization::get_opt_unitary(){return new UnitaryMatrix(*_opt_unitary);}


//-------------------------------------------SimulatedAnnealing-----------------------------------

SimulatedAnnealing::SimulatedAnnealing(CIMethod * cim, double temp, double change_temp, double angle, double change_angle):
    OrbitalOptimization(cim)
{
    _temp = temp;
    _change_temp = change_temp;
    _max_angle = angle;
    _change_angle = change_angle;
    _c_steps = 20000;
    std::cout << "We constructed a SA optimizer." << std::endl;
}

void SimulatedAnnealing::set_good_start()
{
    std::cout << "We search for a good starting guess for sa." << std::endl;
    double stemp = _temp;
    double sma = _max_angle;
    double sct = _change_temp;
    double sca = _change_angle;
    double e_new;
    try{
        get_fci_no(); //If possible we start from the FCI natural orbitals which is a very good starting guess.
    }
    catch(...){
        std::cerr << "It was not possible to generate the FCI natural orbitals of this system, so we generate the natural orbitals of this CIMethod (as starting guess)." << std::endl;
        get_cim_no(); //If not possible we start from the natural orbitals from this method.
    }   
    _temp = 0.01;
    _max_angle = 0.1;
    _change_temp = 0.9;
    _change_angle = 0.99;
    optimize(); //restart calculation.
    _temp = stemp;
    _max_angle = sma;
    _change_temp = sct;
    _change_angle = sca;
    cout << "Start E from guess is: " << _opt_energy  << endl;
    set_cim_optimal(); //We set our Hamiltonian to the lowest minima found.
}

void SimulatedAnnealing::run_multiple(int num, bool saveopt)
{
    //Play in this loop with different starting simulated annealing parameters.
    set_good_start();
    //First we do num simmulated annealing runs. To try to escape local minima.
    for(int i = 0 ; i < num; i++)
    {
        double e_new = optimize(); //restart calculation.
        cout << "#run : " << i << " E found: " << e_new << endl;
        _temp *= 1. ;//0.8;
        //_change_temp -= (1-_change_temp);//*(1 + 0.001);
        //_change_angle -= (1-_change_angle);//*(1+0.0001) ;
        set_cim_optimal(); //We set our Hamiltonian to the lowest minima found.
        //_orbtrans->fillHamCI(*_optham);
        //get_back_original_ham();
        //This does an extra scan around the optimal point.
        if(i == num-2)
        {
            //set_cim_optimal(); //We set our Hamiltonian to the lowest minima found.
            _temp = 0.01;
            _max_angle = 0.1;
            _change_temp = 0.9;
            _change_angle = 0.99;
            //cout << "around minima ->" << endl;
        }
    }

    if(Orbopt_debugPrint)
    {
        cout << "THE MINIMUM ENERGY VALUE WE FOUND IS: " << _opt_energy << endl;
        cout << "With the following unitary transformation of the starting orbitals." << endl;
        _opt_unitary->print_unitary(cout);
    }    

    if (saveopt)
        save_unitary();

    set_cim_optimal(); //We set our Hamiltonian to the lowest minima found.
}

bool SimulatedAnnealing::accept_function(double e_new, double e_old, double temp, std::mt19937 &mt, std::uniform_real_distribution<double> &dist_real)
{
    if(e_new < e_old)
    {
        /*  
            double cost = exp((e_new - e_old) / temp);
            if ( dist_real(mt) *(1+cost) > cost ){
            return true;
            }// biggest chance for accepting
            else 
            return false;*/
        return true ;
    }
    else
    {
        double cost = exp((e_old - e_new) / temp);
        if ( dist_real(mt) *(1+cost) > cost)
            //biggest chance for failing
            return false;
        else 
            return true;
    }
}

double SimulatedAnnealing::cost_function()
{
    return _cim->get_mulliken({0,1,2,3}, 0);
    //return _cim->get_ci_energy(); //For energy based orbital optimization.
    //Properties prop { _cim->get_eigvec(0)  , _cim->get_l() , _cim->get_perm() };
    //return prop.shannon_ic();
}

double SimulatedAnnealing::optimize()
{
    //int size = _orbtrans->_unitary->getNumVariablesX();
    //double step[size];
    std::cout << setprecision(15);
    double temp = _temp;
    double max_angle = _max_angle;
    double cur_angle;
    double e_old =  cost_function();
    if(Orbopt_debugPrint)
        cout << "start E is:" << e_old << endl;
    //Setup random number generators. (making use of c++11 standards.)
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<int> dist(0, _optham->getL() -1); //for benzene -> 27, 38 ); 
    std::uniform_real_distribution<double> dist_real(0, 1);
    //to get random numbers: dist(mt);
    int unaccept = 0;
    for(unsigned int ll=0;ll<_c_steps;ll++)
    {
        int rani = dist(mt);
        int start = _optham->getNstart(_optham->getOrbitalIrrep(rani));
        int end = _optham->getNstart(_optham->getOrbitalIrrep(rani)+1);
        std::uniform_int_distribution<int> dist2(start, end -1); //for benzene -> 27, 38 ); 
        int ranj = dist2(mt);
        if (rani != ranj)
        {
            //We simulate a more normal like distribution.
            cur_angle = max_angle * (dist_real(mt) - dist_real(mt));
            //cout << "Orbitals we mix: index " << rani << " and index " << ranj << "from 0 till " << _optham->getL() << "angle : " << cur_angle << endl;
            transform_orbitals(_optham->getOrbitalIrrep(rani) , rani, ranj , cur_angle); 
            //transform_withasx(step, temp, maxangle);
            update_cim();

            double e_new = cost_function();

            if (accept_function(e_new, e_old, temp, mt, dist_real))
            {
                e_old = e_new;
                cout << "new E: " << e_new << "old E: " << e_old << " at Temp: " << temp << "Angle : " << cur_angle <<endl;
                if(e_new < _opt_energy)
                {
                    _opt_energy = e_new;
                    _opt_unitary.reset(new UnitaryMatrix(_orbtrans->get_unitary()));
                }
                if(ll <  1000 && ll % 5 == 0)
                {
                    save_ham("sim"+ to_string(ll));

                }
                unaccept = 0;
                if(_c_steps % 1000 == 0)
                {
                    //get_cim_no(); //If not possible we start from the natural orbitals from this method.
                }
            }
            else
            {
                //Go back to the old solution.
                transform_orbitals(_optham->getOrbitalIrrep(rani) , rani, ranj , -1.*cur_angle); 
                unaccept ++;
                /*  
                    for(int ii = 0; ii < size ; ii++){
                    step[ii] *= -1.;
                    }        
                    _orbtrans->update_unitary(step);*/
                //cout << "Step not accepted.  E remains : " << e_old << "Old energy is : " << e_new;
            }
            temp *= _change_temp;
            max_angle *= _change_angle;
            if(unaccept > 2000)
            {
                //we reached a minimum from which we cant escape.
                cout << "Stopped through the unaccept criterium." << endl;
                break;
            }
            if(max_angle < 1e-7)
            {
                break;
            }
        }// endif test good orbitals
    }//endloop coolingsteps

    return e_old;
}

//------------------------------------------------Iterative Subotnik-------------------------------
Iterative_Subotnik::Iterative_Subotnik(CIMethod * cim, double crit, int type): OrbitalOptimization(cim)
{
   _R = matrix(cim->get_l(), cim->get_l());
   _crit = crit;
   _type = type;
}

//Fills R to minimize the seniority contribution of the wavefunction. The sum of 2rdm_iiii is maximized : \sigma_i 2rdm_iiii, which minimizes seniority number and orbital optimizes DOCI a bit.
void Iterative_Subotnik::fill_R_mmin()
{
    std::unique_ptr<CIDens> _cid = _cim->get_density(); //sets 2rdm back to zero. So we can reconstruct it again.
    _cid->construct_density();
    std::vector<int> vals {1,1}; //0->u,u,u,u ; 1-> u,d,u,d ; 2=1 -> d,u,d,u ; 2-> d,d,d,d
    for(int i = 0 ; i < _cim->get_l() ; i++)
    {
        for(int j = 0 ; j < _cim->get_l() ; j++)
        {
            _R(i,j) = 0;
            for(auto spin : vals)//spin-summed
            {
                _R(i,j) +=  _cid->get_two_rdm(spin,i,j,j,j);
            }
        }
    }
    if(CIMethod_debug){
        double sen = _cid->get_seniority();
        std::cout  << "Trace(R): " << calc_value() << " (should increase) " <<std::endl;
        std::cout  << "Minimum seniority is : " << sen << std::endl;
    }
}

//Fills R to calculate Edmiston-Ruedenberg orbitals. Which are a set of localized orbitals, because the self-repulsion of (the HF occupied orbitals or all orbitals) is maximized. \sigma_i(1_max(nalpha,nbeta) <i i | i i>. 
void Iterative_Subotnik::fill_R_ER(int all)
{
   int max_orb;
   if(all = 0)//HF det is invariant through unitary transformations of occupied orbitals so we can keep the HF energy, but localize the occupied orbitals.
   {
       max_orb =  _cim->gNup() > _cim->gNdown()?  _cim->gNup(): _cim->gNdown() ;
   }
   else
   {
       max_orb = _cim->get_l();
   }
   _R.reset();
   for(int i = 0 ; i < max_orb ; i++)
   {
       for(int j = 0 ; j < max_orb ; j++)
       {
           _R(i,j) =  _cim->get_ham()->getVmat(i,j,j,j);
       }
   }
   if(CIMethod_debug)
   {
       std::cout  << "Trace(R): " << calc_value() << " (should increase) " <<std::endl;
   }
}

//Watch out this doesn't take into account the second derivative and sometimes you can go to a minimum of the \sigma_i 2rdm_iiii instead of the maximum and you maximizes the seniority number of the state.
//This is based on the fact that: Tr(R^TU) has a unique global maximum at U = R(R^TR)^(-1/2)
//But this is only for eta which has the same points where the gradient is zero as \sigma_i 2rdm_iiii. This is when R is symmetric.
void Iterative_Subotnik::rotate_unitary_from_R()
{
   //Makes use of the polar form  of a normal matrix A: A = UP with U unitary and P Hermitian positiv-semidefinite. (P = sqrt(A * A^T) )
   matrix prod = matrix(_cim->get_l(), _cim->get_l());
   prod.prod( _R, _R, 'T' , 'N' );
   matrix eigvec= matrix(_cim->get_l(), _cim->get_l());
   matrix eigval= matrix(_cim->get_l(), 1);
   prod.diagonalize(eigval , eigvec);
   matrix tran = matrix(_cim->get_l(), _cim->get_l());
   for(int i = 0 ; i < _cim->get_l() ; i++)
   {
       if(eigval(i,0) > 1e-20)
       {
           tran(i,i) = 1./ sqrt(eigval(i,0));
       }
       else{
           tran(i,i) = 1.;
           _R(_cim->get_l()-i-1, _cim->get_l()-i-1) = 1.; //because lapack orders eigenvectors and eigenvalues in increasing order; virtual orbitals which correspond to eigenvalue zero have low indices (in eigenvector basis) but in determinant order they have high indices, and _R is in determinant order.
       }
   }
   prod.prod( eigvec,tran );
   tran.prod( prod, eigvec , 'N' , 'T');
   prod.prod(_R , tran); //Contains the unitary in columnformat
   if(CIMethod_debug){
       cout << "Should be unitary: " <<endl;
       prod.Print(); 
   }
   _orbtrans->rotate_active_space_vectors(prod.getpointer()); //Transposes the unitary so , it is saved in rowformat and multiplies in from the left with the existing unitaries.
}

double Iterative_Subotnik::check_converged()
{
    //We are at the desired maximum when _R is symmetric.
    double som = 0;
    for(int i = 0 ; i < _cim->get_l() ; i++)
    {
        for(int j = i+1 ; j < _cim->get_l() ; j++)
        {
            som += 2*fabs(_R(i,j) - _R(j,i));
        }
    }
    som /= (double) pow(_cim->get_l(),2);

    if(Orbopt_debugPrint)
        std::cout  << "R deviates this amount from a symmetric matrix: " << som << std::endl;
    return som;
}

double Iterative_Subotnik::calc_value()
{
    return _R.trace();
}

void Iterative_Subotnik::set_good_start()
{
    try{
        get_fci_no(); //If possible we start from the FCI natural orbitals which is a very good starting guess.
    }
    catch(...){
        get_cim_no(); //If not possible we start from the natural orbitals from this method.
    }   
}

   
double Iterative_Subotnik::optimize()
{   
    if( _type == 0 )
        fill_R_mmin();//Fills the R matrix with the current density matrix of cim (if cim changes the density matrix will also change automatically because it points to this cim)
    else if( _type == 1)
        fill_R_ER(1);
    bool notconverged = true;
    int i = 0;
    //set_good_start();
    //double prevtrace = -100;
    while(notconverged)
    {
        rotate_unitary_from_R();
        _orbtrans->fillHamCI(*_optham);
        update_cim();//sets ham in cim to _optham and solves.        
        if( _type == 0 )
            fill_R_mmin();//Fills the R matrix with the current density matrix of cim (if cim changes the density matrix will also change automatically because it points to this cim)
        else if( _type == 1)
            fill_R_ER(1);
        if( i > 10)
        {//Always 10 steps.
            double difsym =check_converged(); 
            notconverged = (difsym > _crit );
            //if(calc_value() < prevtrace)       
                //notconverged = false;
            //else
                //prevtrace = calc_value();
        }

        if(i % 200 == 0)
        {
            std::cout  << "Trace(R): " << calc_value() << " (should increase) " <<std::endl;
        }
        /*  
        if(!notconverged)
        {
            double e_bef = _cim->get_ci_energy();
            double enommin = get_cim_no();
            if( _type == 0 )
                fill_R_mmin();//Fills the R matrix with the current density matrix of cim (if cim changes the density matrix will also change automatically because it points to this cim)
            else if( _type == 1)
                fill_R_ER(1);
            std::cout  << "Trace(R)_NO: " << calc_value() << " (should increase) " <<std::endl;
            cout << " E after R convergence: " << e_bef << " E after no trans: "  << enommin <<std::endl;
            //notconverged = !check_converged(_crit);
        }*/
        i++;
    }
    if(_type == 0)
    {
        std::unique_ptr<CIDens> _cid = _cim->get_density(); //sets 2rdm back to zero. So we can reconstruct it again.
        _cid->construct_density();
        double sen = _cid->get_seniority();
        std::cout  << "Trace(R): " << calc_value() << " (should increase) " <<std::endl;
        std::cout  << "Minimum seniority is : " << sen << std::endl;
        //double e_bef = _cim->get_ci_energy();
        //double enommin = get_cim_no();
        //cout << " E after R: " << e_bef << " energy after no trans: "  << enommin <<std::endl;
    }
    else if(_type ==1)
    {
        std::cout  << "maximal self repulsion is : " << calc_value() <<std::endl;

    }
    _opt_unitary.reset(new UnitaryMatrix(_orbtrans->get_unitary()));
    _cim->solve();
    _opt_energy = _cim->get_ci_energy();
    return calc_value();
}

//---------------start of DIIS subclass ----------------------------------------------------

Iterative_Subotnik_DIIS::Iterative_Subotnik_DIIS(CIMethod * cim, double crit, int type): Iterative_Subotnik(cim, crit, type)
{
    _diis = new DIIS(8);
    _unitlist.resize(0);
}

Iterative_Subotnik_DIIS::~Iterative_Subotnik_DIIS()
{
    delete _diis;
    for(int i = 0 ; i < _unitlist.size() ; i++)
    {
        delete _unitlist[i];
    }
}

void Iterative_Subotnik_DIIS::interpolate_unitary()
{
    std::vector<double> coef = _diis->get_coef();
    _opt_unitary->reset_unitary(1); //we put allzero to 1, so all elements of unitarymatrix are zero. In contrast to the standard reset when on the diagonal are all 1.
    for(int i = 0 ; i < _unitlist.size() ; i++)
    {
        _opt_unitary->add_unit(*_unitlist[i], coef[i]);
        //cout << coef[i] << endl;
    }
    //_opt_unitary->print_unitary(std::cout);
}

void Iterative_Subotnik_DIIS::interpolate_R()
{
    _R = _diis->get_extrapolated();
}

void Iterative_Subotnik_DIIS::rotate_unitary_from_extrapolated_R()
{
    vector<matrix> overlap = _opt_unitary->calc_overlap();
    for(int i = 0 ; i < overlap.size() ; i++)
    {
        //overlap[i].Print();
        overlap[i].invert(1);
        //overlap[i].Print();
    }
    matrix inverseoverlap = matrix(overlap);
    matrix prod = matrix(_cim->get_l(), _cim->get_l());
    matrix prod2 = matrix(_cim->get_l(), _cim->get_l());

    prod.prod( inverseoverlap,_R, 'N' , 'N' );
    prod2.prod( _R, prod, 'T' , 'N' );
    prod2 = prod2.matrix_inv_square_root();
    inverseoverlap.prod( prod,prod2 ); //Contains V matrix of subotnik article.
    //inverseoverlap.Print();


    OptIndex index = _optham->get_index_object();
    UnitaryMatrix V = UnitaryMatrix(inverseoverlap, index);
    //V.CheckDeviationFromUnitary();

    //V.multiply_left_with(*_opt_unitary);
    _opt_unitary->multiply_left_with(V);
    //_opt_unitary->print_unitary(std::cout);
    //_opt_unitary->CheckDeviationFromUnitary();
    _orbtrans->set_unitary(*_opt_unitary); 
    //_orbtrans->get_unitary().multiply_left_with(*_opt_unitary); 
}

void Iterative_Subotnik_DIIS::add_unitlist(UnitaryMatrix & val)
{
    _unitlist.push_back(new UnitaryMatrix(val));
    if(_unitlist.size() > _diis->get_max_size() )
    {
        delete _unitlist[0];
        _unitlist.erase(_unitlist.begin());
    }
}

double Iterative_Subotnik_DIIS::optimize()
{   
    double difsym = 100;
    bool notconverged = true;
    int i = 0;
    double prev = -100.;
    //set_good_start();
    while(notconverged)
    {
        if( _type == 0 )
            fill_R_mmin();//Fills the R matrix with the current density matrix of cim (if cim changes the density matrix will also change automatically because it points to this cim)
        else if( _type == 1)
            fill_R_ER(1);
        if(i > _diis->get_max_size() )
        {
          difsym =check_converged(); 
          notconverged = (difsym > _crit );
        }
        if(notconverged)
        {
            if(_unitlist.size() < _diis->get_max_size() ) //Start collecting matrices for the interpolation, this between the 2nd and maximum number of error matrices kept, this to make sure that we start to interpolate when we are already close enough to our final point (because the R interpolation is approximative).
            {
                add_unitlist(_orbtrans->get_unitary()); //We are rotated after rotate_unitary_from_R or extrapolated R.
                _diis->add_matrices(_R , calc_error_matrix());
                rotate_unitary_from_R();
                _orbtrans->fillHamCI(*_optham);
                update_cim();//sets ham in cim to _optham and solves.        
            }
            else
            {
                _diis->do_diis(_R , calc_error_matrix());
                add_unitlist(_orbtrans->get_unitary()); //We are rotated after rotate_unitary_from_R or extrapolated R.
                interpolate_unitary();
                interpolate_R();
                rotate_unitary_from_extrapolated_R();
                _orbtrans->fillHamCI(*_optham);
                update_cim();//sets ham in cim to _optham and solves.        
                double calcval = calc_value();
                //Make sure we maximize the trace -> minimize the seniority.
                _diis->reset();
                _unitlist.resize(0);
                if (calcval > prev)
                {
                    std::cerr << "DIIS decreased the trace: prev =  " <<  prev << " new = " << calcval <<std::endl;
                    //_orbtrans->set_unitary(*_unitlist[_unitlist.size() -1]); 
                    //_orbtrans->fillHamCI(*_optham);
                    //update_cim();//sets ham in cim to _optham and solves.        
                    //_diis->reset();
                    //_unitlist.resize(0);
                }
                else
                    prev = calcval;
            }
            //If it takes to long to converge, we reset the DIIS and do a few normal steps, to stop the oscillation of the approximate DIIS procedure.
            if(i % 50 == 0)
            {
                //_diis->reset();
                //_unitlist.resize(0);
            }
            else if(i % 60 == 0. )
            {
                //_crit *= 10.;
            }
            std::cout  << "iteratie: "<<i << " Trace(R): " << calc_value() << " (didnt increase) prev: " << prev << "crit : " << _crit << "current difsym "<< difsym << std::endl;
            std::cout << "dif smaller" <<std::endl;
            //UnitaryMatrix unit = *_unitlist.back();
            //_orbtrans->set_unitary(unit.get_inverse() );
            //_orbtrans->fillHamCI(*_optham);
            //update_cim();//sets ham in cim to _optham and solves.        
            //_diis->reset();
            //_unitlist.resize(0);
        }/*      
        if(!notconverged) //Extra final check with stable subotnik solver.
        {
            int extranum = _diis->get_max_size();
            for(int j = 0 ; j < extranum ; j++)
            {
                if( _type == 0 )
                    fill_R_mmin();//Fills the R matrix with the current density matrix of cim (if cim changes the density matrix will also change automatically because it points to this cim)
                else if( _type == 1)
                    fill_R_ER(1);
                add_unitlist(_orbtrans->get_unitary()); //We are rotated after rotate_unitary_from_R or extrapolated R.
                _diis->add_matrices(_R , calc_error_matrix());
                rotate_unitary_from_R();
                _orbtrans->fillHamCI(*_optham);
                update_cim();//sets ham in cim to _optham and solves.        
                double difsym =check_converged(); 
                notconverged = (difsym > _crit );

            }     
            i+= extranum -1;
        }*/
        i++;
    }
    if(_type == 0)
    {
        std::unique_ptr<CIDens> _cid = _cim->get_density(); //sets 2rdm back to zero. So we can reconstruct it again.
        _cid->construct_density();
        double sen = _cid->get_seniority();
        std::cout  << "Trace(R): " << calc_value() << " (should increase) " <<std::endl;
        std::cout  << "Minimum seniority is : " << sen << std::endl;
        //double e_bef = _cim->get_ci_energy();
        //double enommin = get_cim_no();
        //cout << " E after R: " << e_bef << " energy after no trans: "  << enommin <<std::endl;
    }
    else if(_type ==1)
    {
        std::cout  << "maximal self repulsion is : " << calc_value() <<std::endl;

    }
    _opt_unitary.reset(new UnitaryMatrix(_orbtrans->get_unitary()));
    //_cim->solve();
    _opt_energy = _cim->get_ci_energy();
    return calc_value();
}

matrix Iterative_Subotnik_DIIS::calc_error_matrix()
{
    int n = _R.getn(); 
    matrix D = matrix(n, n);
    for(int i = 0 ; i < n ; i++)
    {
        for(int j = i+1 ; j < n ; j++)
        {
            D(i,j) = _R(i,j) - _R(j,i);
            D(j,i) = -1*D(i,j);
        }
    }	
    return D;
}

//-------------------------------------------LineSearch-----------------------------------------
void LineSearch::find_bracket(double x1, double step , double * bracket, int i, int j){
    double c = 1.61803398875;//magical line search  number. see -> Kiusalaas_J._Numerical_Methods_in_Engineering_with_Python_(2005)(en)(424s).pdf p. 383-385 based on the golden ratio.
    _orbtrans->do_jacobi_rotation(i,j , x1, _optham);
    update_cim();
    double f1 = _cim->get_ci_energy();
    double x2 = x1 + step;
    _orbtrans->do_jacobi_rotation(i,j , step, _optham);//REMARK _optham is already in x1 so if we add step we go to x2.
    update_cim();
    double f2 = _cim->get_ci_energy();
    //To make sure we advance in lowering direction.
    if(f2>f1){
        step *= -1;
        x2 = x1 + step;
        _orbtrans->do_jacobi_rotation(i,j , 2*step, _optham);//REMARK _optham is already in x1+step so if we want to go to x1-step, we have to change 2*step.
        update_cim();
        f2 = _cim->get_ci_energy();
        //Check if minimum between x1-h and x1+h
        if(f2> f1){
            bracket[0] = x2;
            bracket[1] = x1-step; //reminder step is negative now.
            return;
        }
    }
    for(int l = 0 ; l < 100 ; l++)
    {
        step *= c;
        double x3 = x2+step;
        _orbtrans->do_jacobi_rotation(i,j , step, _optham);//REMARK _optham is already in x1 so if we add step we go to x2.
        update_cim();
        double f3 = _cim->get_ci_energy();
        if(f3> f2)
        {
            bracket[0] = x1;
            bracket[1] = x3;
            return;
        }
        x1 = x2;
        x2 = x3;
        f1 = f2;
        f2 = f3;
    }
    throw std::logic_error( "Bracket did not find a minimum.");
}

double LineSearch::func(double angle, int i , int j)
{
    _orbtrans->do_jacobi_rotation(i,j , angle, _optham);
    update_cim();
    return _cim->get_ci_energy();
}

std::vector<double> LineSearch::search(int i , int j, double a , double b ,double tol)
{
    std::vector<double> results(2,0.);
    //a , b contains the interval of the angle in which the minimum is situated.
    int niter = -2.078087*std::log(tol/abs(b-a)); // Eq. (10.4) in book of Kiusalaas J.
    double R = 0.618033989;
    double C = 1.0 - R;
    //first telescoping.
    double x1 = R*a+ C*b;
    double x2 = C*a+ R*b;
    _orbtrans->do_jacobi_rotation(i,j , x1,_optham);
    update_cim();
    double f1 = _cim->get_ci_energy();//if you don't give a pointer to a Hamiltonian object it will copy the original ham in orbtrans and rotates from there.
    _orbtrans->do_jacobi_rotation(i,j , x2-x1,_optham); //_optham is in x1 -> to get x2 rotate with (x2-x1)
    update_cim();
    double f2 = _cim->get_ci_energy();
    int last = 1; //If last=1 ham is in current x2 , if last = 0 ham is in current x1
    for(int nloop = 0 ; nloop < niter ; nloop++)
    {
        if(f1 < f2)
        {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = C*a+R*b;
            if(last)
                _orbtrans->do_jacobi_rotation(i,j , x2-x1,_optham);//_optham is in x2_prev
            else
                _orbtrans->do_jacobi_rotation(i,j , x2-a,_optham);//_optham is in x1_prev
                
            update_cim();
            f2 = _cim->get_ci_energy();
            cout << "E after optimizing orbitals i " << i<< "and j " << j << " is :" <<f2 <<endl;
            last = 1;
        }
        else
        {
            b = x2;
            x2 =x1 ;
            f2 = f1;
            x1 = R*a+C*b;
            if(last)
                _orbtrans->do_jacobi_rotation(i,j , x1-b,_optham);//_optham is in x2_prev
            else
                _orbtrans->do_jacobi_rotation(i,j , x1-x2,_optham);//_optham is in x1_prev
            update_cim();
            f1 = _cim->get_ci_energy();
            cout << "E after optimizing orbitals i " << i<< "and j " << j << " is :" <<f1 <<endl;
            last = 0;
        }
    }
    if(last)
    {
        results[0] = x2;
        results[1] = f2;
        return results;
    }
    else
    {
        results[0] = x1;
        results[1] = f1;
        return results;
    }
}

double LineSearch::optimize()
{
    int numpair = 0;
    //Only rotate between orbitals in same irrep to keep the symmetry.
    double x1 = -1.5; //start_angle for interval search
    double step = 1.5;
    double bracket[2];//array that contains the interval that contains a minimum.
    double tol = 1e-10;
    double e_prev = 100;
    std::vector<double> results;
    results.resize(2);
    OptIndex index = _optham->get_index_object();
    while(true)
    {
        for(int n = 0; n < index.getNirreps() ; n++)
        {
            for(int i = 0 ; i < index.getNORB(n) ; i++)
            {
                for(int j =i+1 ; j < index.getNORB(n) ; j++)
                {
                    
                    find_bracket( x1, step , bracket, i+index.getNstart(n), j+index.getNstart(n));
                    results = search(i + index.getNstart(n),j+index.getNstart(n), bracket[0], bracket[1] ,1e-11);
                }
            }
            //cout << "E after optimizing irrep " << n << " is :" << results[1] <<endl;
        }
        cout << "E after full run: " << results[1] <<endl;
        if(abs(abs(e_prev) - abs(results[1])) < tol)
        {
            break;
        }
        e_prev = results[1];
    }    
    return 0;
                
}

/*
def powell(F,x,h=0.1,tol=1.0e-6):
def f(s): return F(x + s*v) # F in direction of v
n = len(x) # Number of design variables
df = zeros((n),type=Float64) # Decreases of F stored here
u = identity(n)*1.0 # Vectors v stored here by rows
for j in range(30): # Allow for 30 cycles:
xOld = x.copy()
# Save starting point
fOld = F(xOld)
# First n line searches record decreases of F
for i in range(n):
v = u[i]
a,b = bracket(f,0.0,h)
s,fMin = search(f,a,b)
df[i] = fOld - fMin
fOld = fMin
x = x + s*v
# Last line search in the cycle
v = x - xOld
a,b = bracket(f,0.0,h)
s,fLast = search(f,a,b)
x = x + s*v
394
Introduction to Optimization
# Check for convergence
if sqrt(dot(x-xOld,x-xOld)/n) < tol: return x,j+1
# Identify biggest decrease & update search directions
iMax = int(argmax(df))
for i in range(iMax,n-1):
u[i] = u[i+1]
u[n-1] = v
print ’’Powell did not converge’’
*/
/* vim: set ts=4 sw=4 expandtab :*/
