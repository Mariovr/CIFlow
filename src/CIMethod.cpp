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
#include <cmath>
#include <string>
#include <cstdlib>
#include<limits.h> // for CHAR_BIT
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <exception>
#include <vector>
#include <bitset>
#include <boost/multiprecision/cpp_int.hpp>

#include "Permutator.h"
#include "Output.h"
#include "CIDens.h"
#include "Hamiltonian.h"
#include "CIMethod.h"
#include "Options.h"
#include "scpp_assert.h"
#ifdef __OMP__
#include <omp.h>
#endif 

using namespace std;
namespace mp = boost::multiprecision;     // Reduce the typing a bit later...


//Show the progress of the calculation.
void printProgBar( int percent ){
	string bar;
	for(int i = 0; i < 50; i++) {
		if( i < (percent/2))	
			bar.replace(i,1,"=");
		else if( i == (percent/2))
			bar.replace(i,1,">");
		else
			bar.replace(i,1," "); 
	}
	cout<< "\r" "[" << bar << "] ";
	cout.width( 3 );
	cout<< percent << "%     " << std::flush;
}

//Constructor
CIMethod::CIMethod(Hamiltonian * hams, const int nup , const int ndown){
    _ham = hams;
    _nup = nup;
    _ndown = ndown;
    _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    _updim = bico(_ham->getL(), _nup);	
    _downdim = bico(_ham->getL(), _ndown);
    _solved = false;
}

//Constructor
CIMethod::CIMethod(Hamiltonian * hams){
    _ham = hams;
    _nup = hams->getnup();
    _ndown = hams->getndown();
    _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    _updim = bico(_ham->getL(), _nup);	
    _downdim = bico(_ham->getL(), _ndown);
    _solved = false;
}

//Destructor
CIMethod::~CIMethod()
{
}

int CIMethod::gUpDim()const { return _updim; }
int CIMethod::gDownDim()const { return _downdim; }
int CIMethod::gNup()const {return _nup; }
int CIMethod::gNdown()const { return _ndown; }
Hamiltonian * CIMethod::get_ham(){ return _ham ;}
Permutator * CIMethod::get_perm(){return _perm.get() ;}

int CIMethod::get_l() const {return _perm->get_norb();}

double CIMethod::get_eigvec(int numstate , int numvec){ return _eigvec(numstate,numvec); }
std::vector<double> CIMethod::get_eigvec(int numvec){ return _eigvec.get_column(numvec); }
double CIMethod::get_econst() const{ return _ham->getEconst(); }

matrix CIMethod::get_mat(){
    matrix d(_mat.getn(),_mat.getm());
     _mat.ConvertToMatrix(d);
    return d;
}

double CIMethod::get_mat(int i , int j ){
    return _mat(i,j);
}

void CIMethod::set_sparsity(double zero)
{
    return _mat.set_zero(zero);
}

void CIMethod::set_ham(Hamiltonian * ham, bool build){
    _ham = ham; //actually not necessary because ham and _ham point to the same adres in the orbital optimization routines, but we keep it because in the future this function could also be used for applications where they differ.
    reset_mat();
    if(build)
    {
        build_ci();
    }
    _solved = false;
    _cid->reset_density();
}   

void CIMethod::load_ham( std::string hamname , bool build, bool reset)
{
    if(reset)
        delete_ham();
    _ham = new Hamiltonian(hamname);
    reset_mat();
    if(build)
    {
        build_ci();
    }
    _solved = false;
	std::cout  << "Loaded Hamiltonian from : " << hamname << std::endl;
}

void CIMethod::new_ham_pointer()
{
    _ham = new Hamiltonian(*_ham);
}
void CIMethod::delete_ham()
{
    delete _ham;
}

void CIMethod::print_output(std::vector<std::string> vec , int num , bool all )
{
    _output->print_output(vec, num, all);
}

void CIMethod::print_ham()
{
    _output->print_ham();
}

void CIMethod::construct_density(unsigned state, bool trdm)
{
    //_cid->print_one_dens(std::cout);
    _cid->set_state(state);
    if(  ! _cid->is_constructed().first ||  (trdm  && !_cid->is_constructed().second ) )
    {
        _cid->construct_density(trdm);
        //cout << "Density is constructed from the solution in CIMethod because density was not yet constructed." << std::endl;
    }
    //_cid->print_one_dens(std::cout);
}

void CIMethod::reset_density()
{
    _cid->reset_density();
}

void CIMethod::print_rdm(unsigned state, bool trdm)
{
    construct_density(state, trdm);
    _output->print_rdm( trdm);
}

void CIMethod::print_one_rdm(std::ostream & os)
{
    os << "#The eigenstate from which the density matrix is constructed: " << _cid->get_state() << endl;
    os << "#THE 1RDM:"<<  endl;
    _cid->print_one_dens(os);
}

void CIMethod::print_two_rdm(std::ostream & os)
{
    os << "#THE 2RDM:" << endl;
    _cid->print_two_dens(os);
}
        
double CIMethod::get_spin_squared(unsigned state)
{
    construct_density(state, true);
    return _cid->get_spin_squared();
}

double CIMethod::get_seniority(unsigned state)
{
    construct_density(state, true);
    return _cid->get_seniority();
}

double CIMethod::get_pare(int maxorb, unsigned state)
{
    construct_density(state, true);
    return _cid->get_dens_energy(maxorb);
}

double CIMethod::get_mulliken(std::vector<int> orbs, unsigned state )
{
    construct_density(state, false);
    //SCPP_ASSERT(_cid->is_constructed().first , "density is not constructed but we extract variables.");
    double mul = _cid->get_mulliken(orbs);
    return mul; 
}


double CIMethod::get_spin(unsigned state)
{
    construct_density(state, true);
    return _cid->get_spin();
}

void CIMethod::reset_output(std::string sort, bool partial)
{
    _output->change_filename(sort, partial);
}

std::string CIMethod::get_perm_info()
{
    return _perm->get_clean_detfile();
}


void CIMethod::solve(int neigval){
    if(_solved != true){
        if(CIMethod_debug)
            cout << "We solve the Hamiltonian with dimension: " << get_dim() << endl;
        _eigval = matrix(neigval , 1);
        _eigvec = matrix(get_dim() , neigval);
        _mat.arpackDiagonalize(neigval , _eigval , _eigvec);
        //matrix mat = matrix(get_dim() , get_dim() );
        //_mat.ConvertToMatrix(mat);
        //_eigval = matrix(get_dim(), 1);
        //_eigvec = matrix(get_dim() , get_dim());
        ////mat.Print();
        //cout << "We solve the Hamiltonian with dimension: " << get_dim() << endl;
        //mat.diagonalize(_eigval, _eigvec);
        //_eigval.Print();
        _solved = true;
        SCPP_TEST_ASSERT( fabs(_eigvec.inproduct(0) - 1.) < 1e-9  , "Error in solve routine of CIMethod: the eigenvectors are not normalized: the inproduct is : " << _eigvec.inproduct(0)<< "with eigenvalue : " << _eigval[0]);
    }    
    else{
        cout << "The system is unchanged after the last call of solve." << endl;
    }
    _mat.reset_vecs();
    _cid->reset_density();
}

double CIMethod::get_ci_energy(int neigval){
    return _eigval[neigval] + get_econst();
}

double CIMethod::get_hf_energy(){
    /*  
    //REMARK we have assumed that the orbitals are ordered according to increasing eigenvalue of fock matrix.
    mp::uint512_t int _last_bit = 1;
    _last_bit <<= _nup;
    _last_bit --; //all up electrons occupy nup lowest sp orbitals
    mp::uint512_t int _last_bit2 = 1;
    _last_bit2 <<= _ndown;
    _last_bit2 --;// all down electrons occupy ndown lowest sp orbitals
    return diagonal(_last_bit, _last_bit2)+ get_econst();// doesn't work ofcourse if the HF determinant is not included 
    */
    return _ham->get_hf_energy();
}

void CIMethod::check_hermiticity(){
    /* 
     * Checks the hermiticity of the Hamiltonian matrix
     * */
    for(int i = 0 ; i < _mat.getn() ; i ++){
        for(int j = 0 ; j < _mat.getm() ; j ++){
            if(_mat(i,j) != _mat(j,i)){
                cout << "matrix is not hermitian: i,j: "<< i << "," << j <<"mat(i,j) " << _mat(i,j) << "mat(j,i)" << _mat(j,i) ;  
                exit(1);
            }
        }// end for i
    }//end for j
    cout << "Hamiltonian matrix is hermitian.\n" ; 
}// end check hermiticity

std::string CIMethod::get_ex_info(){
    return _ham->get_info();
}


double CIMethod::get_offdiag_elem(TYPE up ){
    //THIS function is deprecated remove in due time.
    long ret = sizeof(TYPE ) * CHAR_BIT - 1;
    //long ind1 = up ? ret - __builtin_clzl(up) : ret; // get highest set bit  (__builtin_clzl(up) for mp::uint512_t)
    //long ind2 = __builtin_ffsl(up) - 1; // get lowest set bit (__builtin_ffsl(up) for long int)
    //return _ham->gMxElement(ind1 , ind1 , ind2 , ind2);
    return 0;
}

int * CIMethod::get_orbs(TYPE up, int nset){
	// returns array with the occupied orbitals determined by the bit string
	int * array = new int [nset];
	int i = 0;
	int count = 0;
	while( up  != 0){
		if ((up &1) == 1 ){
		array[i] = count;
		i++;
		}
		up = up >> 1;
		count += 1;
	} 
	return  array;
}

double CIMethod::get_diag_cont_one_bitstring(int * _orbs , int npar){
    double somup = 0;
    for(int k = 0 ; k < npar ; k++){
        somup += _ham->getTmat(_orbs[k], _orbs[k]);
    }
    for(int k = 0 ; k < npar-1 ; k++){
        for(int l = k+1 ; l < npar ; l++){
            somup +=   _ham->getVmat(_orbs[k],_orbs[l] , _orbs[k] , _orbs[l]) -_ham->getVmat(_orbs[k], _orbs[l] , _orbs[l] , _orbs[k]) ;
        } // endl  
    }// end k
    return somup;

}

double CIMethod::get_diag_cont_jmix(int * orbsup , int * orbsdown){
    double som = 0;
    for(int k = 0 ; k < gNup() ; k++){
        for(int l = 0 ; l < gNdown() ; l++){
            som += _ham->getVmat(orbsup[k],orbsdown[l] , orbsup[k] , orbsdown[l]);
        }
    }
    return som;
}

void CIMethod::arrange_order(int * orbs , int num ,TYPE up1 , TYPE up2){
	int orb_reorder[num];
	int upin1,upin2;
	if( num == 4){
		upin1 = 0 ;
		upin2 = 2 ;
	}
	else if(num == 2) {
		upin1 = 0 ;
		upin2 = 1 ;
	}

    TYPE d = 1;
	for(int i = 0 ; i < num ; i ++){
		if(up1 & ( d << orbs[i])){
			orb_reorder[upin1] = orbs[i];
			upin1 ++;

		}
		else{
			orb_reorder[upin2] = orbs[i];
			upin2++ ;
		}
	}
	for( int i = 0 ; i < num ; i++){
		orbs[i] = orb_reorder[i];
	}
}

int CIMethod::get_nstacked(int * orbs , TYPE up)
{
	int sum = 0;
    int small = orbs[0] < orbs[1] ? orbs[0] : orbs[1];

    TYPE shiftbit = 1;
	for(int i = small ; i < fabs(orbs[1] - orbs[0]) + small ; i++)
    {
		if (up & (shiftbit << i)  )
        {
			sum += 1;
		}
	} 
	return sum;
}

double CIMethod::diagonal(TYPE up1 , TYPE down1){
	int *  orbsup = get_orbs(up1,gNup());
	double somtot = get_diag_cont_one_bitstring(orbsup ,gNup());
	int * orbsdown = get_orbs(down1,gNdown());
	somtot += get_diag_cont_one_bitstring(orbsdown,gNdown());
	somtot += get_diag_cont_jmix(orbsup , orbsdown);
	delete [] orbsup;
	delete [] orbsdown;
	return  somtot;
}

double CIMethod::one_diff_orbital(TYPE up1 , TYPE up2, TYPE down1){
    int * orbs = get_orbs(up1 ^ up2 , 2);
    double elemu = _ham->getTmat(orbs[0] , orbs[1]);
    int * otheruporbs = get_orbs(up1 & up2 , gNup() -1);
    for(int z = 0 ; z < gNup() -1 ; z++){
        elemu += _ham->getVmat(orbs[0],otheruporbs[z] , orbs[1] , otheruporbs[z]) - _ham->getVmat(orbs[0],otheruporbs[z] , otheruporbs[z] , orbs[1]) ;
    }
    int diff = get_nstacked(orbs, up1 & up2);
    int * downorbs = get_orbs(down1 ,gNdown());
    for(int z = 0 ; z <gNdown()  ; z++){
        elemu += _ham->getVmat(orbs[0],downorbs[z] , orbs[1] , downorbs[z]);
    }
    if (diff % 2 != 0 )
        elemu *= -1;
    delete [] downorbs;
    delete [] otheruporbs;
    delete [] orbs;
    return elemu;
}

double CIMethod::two_diff_orbitals(TYPE up1, TYPE up2){
	int * orbs = get_orbs(up1^ up2 , 4);
	arrange_order(orbs , 4 , up1 , up2);
	int checkorbs[2] = {orbs[0] , orbs[1]};
	int diff = get_nstacked(checkorbs, up1 & up2);
	checkorbs[0] = orbs[2]; checkorbs[1] = orbs[3];
	diff += get_nstacked(checkorbs,up1 & up2);
	double elem =_ham->getVmat(orbs[0],orbs[1] , orbs[2] , orbs[3]) - _ham->getVmat(orbs[0],orbs[1] , orbs[3] , orbs[2]) ;
	 if (diff % 2 != 0 ){
		 elem *= -1;
	 }
	delete [] orbs;
	return elem;
}

double CIMethod::two_diff_orbitals(TYPE up1 , TYPE up2, TYPE down1, TYPE down2){
    int * orbs = get_orbs(up1 ^ up2 , 2);
    int * downorbs = get_orbs(down1 ^ down2 , 2);
    arrange_order(orbs , 2 , up1 , up2);
    arrange_order(downorbs , 2 , down1 , down2);
    int diff = get_nstacked(orbs, up1 & up2);
    diff += get_nstacked(downorbs , down1 & down2);
    double elem =_ham->getVmat(orbs[0],downorbs[0] , orbs[1] , downorbs[1]);
    if (diff % 2 != 0 /*|| downorbs[0] > orbs[0]*/){
        elem *= -1;}
    delete [] downorbs;
    delete [] orbs;
    return elem;
}

double CIMethod::get_ham_element(TYPE up1 , TYPE down1 , TYPE up2 , TYPE down2){
    if(_perm->popcount( down1 ^ down2) > 4  || _perm->popcount(up1 ^ up2) > 4){ 
        return 0;
    }
	else if (_perm->popcount( up1 ^ up2) == 4 && _perm->popcount(down1 ^ down2) == 0){ // two different occupied orbitals
		return two_diff_orbitals(up1, up2);
	}
	else if(_perm->popcount(up1 ^ up2) == 2){ // one different occupied orbital-> two posibilities: 1) downs are the same , 2) downs have also one differen occupied orbital.
	    if(_perm->popcount(down1 ^ down2) == 0){
		//1)downs are the same:
	        return  one_diff_orbital(up1 , up2, down1);
	    }//end downs are the same

	    else if(_perm->popcount(down1 ^ down2) == 2){
	    //2) downs differ also by own orbital.
	     return  two_diff_orbitals(up1 , up2, down1, down2);
	    } // ends down and ups differ by one
	} // end ups differ by one.
	else if (_perm->popcount( down1 ^ down2) == 4 && _perm->popcount(up1 ^ up2) == 0){ // two different occupied orbitals of the downs
		return  two_diff_orbitals(down1, down2);
	} //end downs differ by two.
	else if(_perm->popcount(down1 ^ down2) == 2 &&_perm->popcount(up1 ^ up2) == 0 ){ // one different occupied orbital of downs->  ups are the same (downs and up differ by 1 is already taken into acount above)
		return one_diff_orbital(down1 , down2, up1);
	} // end down differ by one.
   else if(_perm->popcount( down1 ^ down2) == 0 && _perm->popcount(up1 ^ up2) == 0){ // diagonal contribution
	    return diagonal(up1 , down1);
	}
	return 0.; // when there are two orbitals different in up and one in down or the other way around.
}

void CIMethod::build_parallel()
{
    #ifdef __OMP__
    int num_t = omp_get_max_threads();
    #else
    int num_t = 1;
    #endif
    unsigned long long num_elems;
    if (get_name() != "FCI")
        num_elems = (get_dim()*1ull*(get_dim()+1ull))/2.;
    else 
        num_elems = (gUpDim()*1ull*(gUpDim()+1ull))/2.;
    unsigned long long size_part = num_elems/num_t + 1;

    // every thread should process the lines between i and i+1
    // with i the thread number
    std::vector<unsigned int> workload(num_t+1);
    workload.front() = 0;
    if (get_name() != "FCI")
        workload.back() = get_dim();
    else 
        workload.back() = gUpDim();

    for(int i=1;i<num_t;i++)
    {
        auto num_lines = workload[i-1];
        unsigned long long num_elems = 0;

        if (get_name() != "FCI")
        {
            while(num_elems < size_part)
               num_elems += get_dim() - num_lines++;

            if(num_lines > get_dim() )
               num_lines = get_dim();
        }
        else 
        {
            while(num_elems < size_part)
               num_elems += gUpDim() - num_lines++;

            if(num_lines > gUpDim() )
               num_lines = gUpDim();
        }

        workload[i] = num_lines;
    }
    vector< std::unique_ptr<SparseMatrix_CRS> > parts(num_t);
    if (CIMethod_debug){
        cout << "Running with " << num_t << " threads." << std::endl;
    }
    #pragma omp parallel
    {
        #ifdef __OMP__
        int me = omp_get_thread_num();
        #else
        int me = 0;
        #endif
        if (get_name() != "FCI")
            parts[me].reset(new SparseMatrix_CRS(workload[me+1] - workload[me], get_dim() ));
        else
            parts[me].reset(new SparseMatrix_CRS((workload[me+1] - workload[me]) * gDownDim() , get_dim() ));


        construct_CI_matrix((*parts[me]), workload[me], workload[me+1]);
    }// End parallel

    _mat.add_from_list(parts); //REMARK because we only filled the upper diagonal => _mat should be a SparseMatrix_CRS_Sym type.

    if (CIMethod_debug){
        cout << "Number of nonzero elements: " << _mat.datasize() << " dimension :  " << get_dim() << " number uppperdiagonal elements: " << get_dim() * (get_dim() +1) /2. << " sparsity : " << _mat.datasize() / ( get_dim() * (get_dim() +1) /2. )<<endl ;
    }
}


//---------------------------------SUBCLASSES-------------------------------------------------

DOCI::DOCI(Hamiltonian * ham ) : CIMethod(ham){
    _output.reset(new OutputSFDOCI(this, ham->get_short_filename() + "outputdoci.dat") );
    _perm.reset(new Permutator_Bit(ham->getL() ) );
    _mat = SparseMatrix_CRS_Sym(get_dim(),get_dim());
    _cid = get_density() ;
    build_parallel();
}

FCI::FCI(Hamiltonian * ham ):CIMethod( ham){
    _output.reset( new OutputSFFCI(this, ham->get_short_filename() + "outputfci.dat" ) );
    _perm.reset(new Permutator_Bit(ham->getL() ) ) ;
	_mat = SparseMatrix_CRS_Sym(get_dim(),get_dim());
    _cid = get_density() ;
    build_parallel(); //Because of the use of symmetry we have to use different parallelization.
}

FCI_File::FCI_File(Hamiltonian * ham ,  std::string permfile):CIMethod( ham){
    _perm.reset(new Permutator_File(permfile , ham->getL() ) ) ;
    _output.reset( new OutputSFFCI_File(this ,ham->get_short_filename() + _perm->get_clean_detfile() +"outputci_file.dat" ) );
	_mat = SparseMatrix_CRS_Sym(get_dim(),get_dim());
    _cid = get_density() ;
	cout << ham->get_short_filename() + _perm->get_clean_detfile() +"outputci_file.dat" << endl;
    build_parallel();
}

CI_Big::CI_Big(Hamiltonian * ham ,  std::string permfile):CIMethod( ham)
{ 
    _perm.reset(new Permutator_File(permfile, ham->getL() ) );
    _output.reset( new OutputSFFCI_File(this ,ham->get_short_filename() + _perm->get_clean_detfile() +"outputcibig_file.dat") );
    _cid = get_density() ;
    //CISIZE is defined in options and gives the size in the gigabytes of available storage for the CI matrix.
    long numelements = CISIZE*1024.*1024.*1024.*3/(64.*5)*10./2.; //We keep 4/5 of the available memory for the CI matrix, the extra factor 10 is because we assume a sparsity of less than 10 percent. So we have 10 times as much as possible spots of nonzero elements to obtain numelements nonzero elements.
    unsigned long dim = get_dim();
    unsigned long numlines = round((-1+sqrt(1+8*numelements))/2);
    cout << "the lines we keep in memory:"<< numlines << endl;
    if( numlines > get_dim()) {
        cout << "We can keep the total matrix in memory you should use FCI_File, but we'll give it a try." << endl;
        _mat = SparseMatrix_CRS_Sym(get_dim(),get_dim());
    }
    else{
        _mat = SparseMatrix_CRS_Sym(numlines,numlines);
    }
    build_parallel();
};

void DOCI::construct_CI_matrix(SparseMatrix_CRS & mat, int start_l , int end_l){
    //we make use of the hermicity of the matrix to add off diagonal elements. (we dont use get_ham_element because this implementation is a bit faster (and still short))
    vector< vector<int> > vw  ( get_l() +1,  vector<int> ( gNup() + 1 , 0));
    setup_vertex_weights(vw);
    TYPE shiftbit = 1;

    Permutator_Bit my_perm(_ham->getL());
    TYPE up1 = my_perm.get_start_int(gNup());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        up1 = my_perm.permutate_bit(up1);
    //TYPE up2 = up1;

    for( long i = start_l; i < end_l ; i++ ) { //add off diagonal elements.
        mat.NewRow();
        if (CIMethod_debug){
            if( i % 5000 == 0)
                cout << "line : " << i << " of lines: " << _perm->get_dim() << endl;
        }
        //First add diagonal contribution.
        mat.PushToRow(i , diagonal(up1 , up1)); //with DOCI, ups are always equal to downs

        // Then test all the single excitations you can think off.
        for (int j = 0; j < get_l() ; j++) //Loop over the L spatial orbitals.
            if(up1 & ( shiftbit << j))
            {
                TYPE up_interm = up1 ^ (shiftbit << j);
                for (int k = j + 1; k < get_l() ; k++) 
                {
                    if(up1 & ( shiftbit << k))
                    {
                        continue;
                    } 
                    else //Match is possible so we select all matches.
                    {
                        TYPE up2 = up_interm ^ (shiftbit << k);
                        unsigned int index = determine_weight(up2, vw);
                        SCPP_TEST_ASSERT(index > i, "We full only upper half of the DOCI Hamiltonian, -> row < col but : row = " << i << "col = " << index);
                        mat.PushToRow(index , two_diff_orbitals(up1 , up2, up1, up2));
                    } 
                }
            }//End jth orbital is occupied in current determinant.

        up1 = my_perm.permutate_bit(up1);
        //up2 = up1;
    } //end for rows
    //mat.NewRow(); //not necessary for parallel, we add the last newrow in add_from_list (in sparsematrix_crs)

    if (CIMethod_debug)
        cout << "We constructed the DOCI Hamiltonian matrix with "<< mat.datasize() << "nonzero elements. Between line:" << start_l << "and line: "<< end_l << endl;
} // end DOCI


void DOCI::setup_vertex_weights(vector<vector<int>> & vw) 
{
    /* Set up the vertex_weight vector */
    vw[0][0] = 1; /* Set the first element to one. */
    for (int k = 1; k <= ( get_l()  - gNup()) ; k++) 
    { /* Go down the first column; at most num_orbs-num_elecs can be set to one (must have N elecs in orbitals) */
        vw[k][0] = 1;
    }

    for (int m = 1; m <= gNup(); m++) 
    { /* Fill the remaining columns */
        for (int k = m; k <= ( get_l()  - (gNup() - m)); k++) 
        { /* At the most all previous orbitals are filled */
            vw[k][m] = vw[k-1][m] + vw[k-1][m-1];
        }
    }
    //for (int i = 0; i < _norb +1; i++) 
    //{
        //for (int j = 0; j < _cim->gNup() +1; j++) 
        //{
            //cout << vw[i][j] << " " ;
        //}
        //cout << std::endl;
    //}
}

unsigned int DOCI::determine_weight(TYPE string , const vector<vector<int> > & vw) 
{
    unsigned int weight = 1;
    int num_elecs_interm = 0;

    for (int k = 0; k < get_l(); k++) 
    {
        TYPE shiftbit = 1;
        if(string & ( shiftbit << k))
        {
            num_elecs_interm += 1;
            weight += vw[k][num_elecs_interm];
        }
    }
    return weight-1;
}

void FCI::setup_vertex_weights(vector<vector<int>> & vw) 
{
    /* Set up the vertex_weight vector */
    vw[0][0] = 1; /* Set the first element to one. */
    for (int k = 1; k <= ( get_l()  - gNup()) ; k++) 
    { /* Go down the first column; at most num_orbs-num_elecs can be set to one (must have N elecs in orbitals) */
        vw[k][0] = 1;
    }

    for (int m = 1; m <= gNup(); m++) 
    { /* Fill the remaining columns */
        for (int k = m; k <= ( get_l()  - (gNup() - m)); k++) 
        { /* At the most all previous orbitals are filled */
            vw[k][m] = vw[k-1][m] + vw[k-1][m-1];
        }
    }
    //for (int i = 0; i < _norb +1; i++) 
    //{
        //for (int j = 0; j < _cim->gNup() +1; j++) 
        //{
            //cout << vw[i][j] << " " ;
        //}
        //cout << std::endl;
    //}
}

unsigned int FCI::determine_weight(TYPE string , const vector<vector<int> > & vw) 
{
    unsigned int weight = 1;
    int num_elecs_interm = 0;

    for (int k = 0; k < get_l(); k++) 
    {
        TYPE shiftbit = 1;
        if(string & ( shiftbit << k))
        {
            num_elecs_interm += 1;
            weight += vw[k][num_elecs_interm];
        }
    }
    return weight-1;
}

void FCI::construct_CI_matrix(SparseMatrix_CRS & mat, int start_l , int end_l){
    //and make this function as general as possible (should also work for unrestricted orbitals)
    //the columns and rows of the Hamiltonian matrix are determined such that with a given up, first all the possible downs change and then we change an up again (the downs are fast changing and the ups slow)
    //REMARK: this is only intended for small molecules and debugging purposes, for FCI calculations on big systems use: Gaussian, psi4, games, ... 
    //REMARK:we only fill the upper diagonal elements.
    //we make use of the hermicity of the matrix to add off diagonal elements. (we dont use get_ham_element because this implementation is a bit faster (and still short))
    vector< vector<int> > vw  ( get_l() +1,  vector<int> ( gNup() + 1 , 0));
    vector< vector<int> > vwd  ( get_l() +1,  vector<int> ( gNdown() + 1 , 0));
    setup_vertex_weights(vw);
    setup_vertex_weights(vwd);
    TYPE shiftbit = 1;

    Permutator_Bit my_perm(_ham->getL());
    TYPE up1 = my_perm.get_start_int(gNup());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        up1 = my_perm.permutate_bit(up1);

    for( long i = start_l; i < end_l ; i++ ) { //ups between start_l and end_l.
        TYPE down1 = my_perm.get_start_int(gNdown());
        for(long k = 0; k < gDownDim(); k++){ //runs over rows.
            mat.NewRow();
            //First add diagonal contribution.
            mat.PushToRow(i*gDownDim()+ k , diagonal(up1 , down1)); //add diagonal elements.

            // Then test all the single excitations of the up you can think off.
            for (int j = 0; j < get_l() ; j++) //Loop over the L spatial orbitals.
                if(up1 & ( shiftbit << j))
                {
                    TYPE up_interm = up1 ^ (shiftbit << j);
                    for (int l = j + 1; l < get_l() ; l++) 
                    {
                        if(up1 & ( shiftbit << l))
                        {
                            continue;
                        } 
                        else //Match is possible so we select all matches.
                        {
                            TYPE up2 = up_interm ^ (shiftbit << l);
                            unsigned int index = determine_weight(up2, vw);
                            SCPP_TEST_ASSERT(index > i, "We full only upper half of the FCI Hamiltonian, -> row < col but : row = " << i << "col = " << index);
                            mat.PushToRow(index*gDownDim()+ k , one_diff_orbital(up1 , up2, down1)); //Ads single up different.
                            // Then test all the single excitations of the down you can think off.
                            for (int m = 0; m < get_l() ; m++) //Loop over the L spatial orbitals.
                                if(down1 & ( shiftbit << m))
                                {
                                    TYPE down_interm = down1 ^ (shiftbit << m);
                                    for (int w = 0; w < get_l() ; w++) 
                                    {
                                        if(down1 & ( shiftbit << w) )
                                        {
                                            continue;
                                        } 
                                        else //Match is possible so we select all matches.
                                        {
                                            TYPE down2 = down_interm ^ (shiftbit << w);
                                            unsigned int indexd = determine_weight(down2, vwd);
                                            mat.PushToRow(index*gDownDim()+ indexd , two_diff_orbitals(up1, up2 , down1, down2)); //Ads single down different and single up different.
                                        } 
                                    }
                                }//End jth orbital is occupied in current determinant.
                            // Then test all the single excitations of the single excited up, to get double excited up.
                            for (int m = 0; m < get_l() ; m++) //Loop over the L spatial orbitals.
                                if(up2 & ( shiftbit << m) && m != l)
                                {
                                    TYPE ups_interm = up2^ (shiftbit << m);
                                    for (int w = 0; w < get_l() ; w++) 
                                    {
                                        if(up2& ( shiftbit << w) || w == j)
                                        {
                                            continue;
                                        } 
                                        else //Match is possible so we select all matches.
                                        {
                                            TYPE updex = ups_interm ^ (shiftbit << w);
                                            unsigned int indexdex = determine_weight(updex, vw);
                                            if(indexdex > i)
                                            {
                                                SCPP_TEST_ASSERT(indexdex > i, "We full only upper half of the FCI Hamiltonian, -> row < col but : row = " << i << "col = " << index);
                                                mat.PushToRownoadd(indexdex*gDownDim()+ k , two_diff_orbitals(up1, updex )); //Ads single down different and single up different.
                                            }
                                        } 
                                    }
                                }//End mth orbital is occupied in current determinant.

                        } //Finished building on top of singly excited up1, namely singly excited downs, and doubly excited up.
                    }
                }//End jth orbital is occupied in current determinant.

            // Then test all the single excitations of the down you can think off.
            for (int j = 0; j < get_l() ; j++) //Loop over the L spatial orbitals.
                if(down1 & ( shiftbit << j))
                {
                    TYPE down_interm = down1 ^ (shiftbit << j);
                    for (int l = j + 1; l < get_l() ; l++) 
                    {
                        if(down1 & ( shiftbit << l))
                        {
                            continue;
                        } 
                        else //Match is possible so we select all matches.
                        {
                            TYPE down2 = down_interm ^ (shiftbit << l);
                            unsigned int index = determine_weight(down2, vwd);
                            //SCPP_TEST_ASSERT(index > i, "We full only upper half of the FCI Hamiltonian, -> row < col but : row = " << i << "col = " << index);
                            mat.PushToRow(i*gDownDim()+ index , one_diff_orbital(down1 , down2, up1)); //Ads single down different.
                            //Excite down2 for the second time.
                            for (int m = 0; m < get_l() ; m++) //Loop over the L spatial orbitals.
                                if(down2 & ( shiftbit << m) && m != l)
                                {
                                    TYPE downs_interm = down2^ (shiftbit << m);
                                    for (int w = 0; w < get_l() ; w++) 
                                    {
                                        if(down2& ( shiftbit << w) || w == j)
                                        {
                                            continue;
                                        } 
                                        else //Match is possible so we select all matches.
                                        {
                                            TYPE downdex = downs_interm ^ (shiftbit << w);
                                            unsigned int indexddex = determine_weight(downdex , vwd);
                                            if(indexddex> k)
                                            {
                                                SCPP_TEST_ASSERT(indexddex > k, "We full only upper half of the FCI Hamiltonian, -> row < col but : row = " << k << "col = " << indexddex);
                                                mat.PushToRownoadd(i*gDownDim()+ indexddex , two_diff_orbitals(down1, downdex )); //Ads single down different and single up different.
                                            }
                                        } 
                                    }
                                }//End mth orbital is occupied in current determinant.
                        } //End single excited down.
                    }
                }//End jth orbital is occupied in current determinant.

            down1 = my_perm.permutate_bit(down1);
        }
        up1 = my_perm.permutate_bit(up1);
    }

    if (CIMethod_debug)
        cout << "We constructed the FCI Hamiltonian matrix with "<< _mat.datasize() << "nonzero elements, of dimension: " << get_dim() << endl;
}

void FCI_File::construct_CI_matrix(SparseMatrix_CRS & mat , int start_l , int end_l){
	//for a file we don't expect any ordering information of the determinants-> much slower but more robust
    //cout << "Started construct_CI_matrix in FCI_File with determinantfile: " << _perm->get_filename() <<endl;
    Permutator_File my_perm(_perm->get_filename(),_ham->getL());
    if (CIMethod_debug)
        cout << "Started construct_CI_matrix in FCI_File with determinantfile: " << my_perm.get_filename() << "between lines: " << start_l << " and " << end_l << endl;
	TYPE dets_row[2] = {0};
	TYPE dets_col[2] = {0};
	for(unsigned long i = start_l ; i < end_l ; i ++){
        mat.NewRow();
        if (CIMethod_debug){
            if( i % 5000 == 0)
                cout << "line : " << i << " of lines: " << _perm->get_dim() << endl;
        }
		my_perm.permutate_bit(dets_row,i);
		TYPE up1 = dets_row[0];
		TYPE down1 = dets_row[1];
		for(unsigned long j = i ; j < get_dim() ; j ++){
			my_perm.permutate_bit(dets_col,j);
			TYPE up2 = dets_col[0];
		    TYPE down2 = dets_col[1];
			mat.PushToRow(j , get_ham_element(up1,down1 , up2, down2));
	} // end cols 
	} // end rows 
    //_mat.NewRow(); Last new row is added in sparsematrix class in function add_from_list, change because of parallel implementation.
    if (CIMethod_debug)
        cout << "We constructed the Hamiltonian matrix with "<< mat.datasize() << "nonzero elements. " << endl;
}

int FCI_File::get_dim(){ return _perm->get_dim(); } 

//----------class CI_Big-------------------------------------------------------------------

/**
 * Do the matrix vector product y = A * x
 * @param xmat a m component vector
 * @param ymat a n component vector
 */
void CI_Big::mvprod(const matrix &xmat, matrix &ymat){
    double *x = xmat.getpointer();
    double *y = ymat.getpointer();
 
    mvprod(x,y);
}

/**
 * Do the matrix vector product y = A * x
 * @param xmat a m component vector
 * @param ymat a n component vector
 * We calculate this on the fly so we dont keep the Hamiltonian in memory
 * WARNING: this is rather slow, for smaller systems where it is possible to
 * keep the Hamiltonian in memory use other CIMethods.
 */
void CI_Big::mvprod(const double *x, double *y){
    int msize = _mat.getn();
    #pragma omp parallel for firstprivate(msize)
    for(int i= msize;i< get_dim();i++)
        y[i] = 0;
 
    _mat.mvprod(x, y); //calculate part of the matrix vector product that fitted in memory.

    #pragma omp parallel firstprivate(msize)
    {
        int dim = get_dim();
        double * y_private = new double[dim](); //() initializes automatically to zero.

        Permutator_File my_perm(_perm->get_filename(),_ham->getL());
        TYPE dets_row[2] = {0};
        TYPE dets_col[2] = {0};
        #pragma omp for schedule(dynamic)
        for(int i=0;i< dim ;i++){
            my_perm.permutate_bit(dets_row,i);
            TYPE up1 = dets_row[0];
            TYPE down1 = dets_row[1];
            int startfly = (i < msize )? msize: i;
            for(int j= startfly ;j< dim;j++){
                my_perm.permutate_bit(dets_col,j);
                TYPE up2 = dets_col[0];
                TYPE down2 = dets_col[1];
                double ham_el = get_ham_element(up1,down1 , up2, down2);
                y_private[i] += ham_el * x[j];
                if(i != j ){
                    y_private[j] += ham_el * x[i];
                }
           }// end cols
        }// end rows
        
        for(int i=0; i< dim; i++){
            #pragma omp atomic
            y[i] += y_private[i];
        }
        delete [] y_private;
    }//end parallel

}

void CI_Big::solve(int neigval){
	_eigval = matrix(neigval , 1);
	_eigvec = matrix(get_dim() , neigval);
    cout << "We solve the Hamiltonian with dimension: " << get_dim() << endl;
    cout << "Started to diagonalize the CI Hamiltonian, while calculating the CI matrix elements on the fly." << endl;
    arpackDiagonalize(neigval , _eigval , _eigvec);
}


void CI_Big::arpackDiagonalize(int nve , matrix & eigval , matrix & eigvec)
{
    // dimension of the matrix
    int n = get_dim();

    // number of eigenvalues to calculate
    int nev = nve;

    // reverse communication parameter, must be zero on first iteration
    int ido = 0;
    // standard eigenvalue problem A*x=lambda*x
    char bmat = 'I';
    // calculate the smallest algebraic eigenvalue
    char which[] = {'S','A'};
    // calculate until machine precision
    double tol = 0;

    // the residual vector
    double *resid = new double[n];

    // the number of columns in v: the number of lanczos vector
    // generated at each iteration, ncv <= n
    // We use the answer to life, the universe and everything, if possible
    int ncv = 42;

    if( n < ncv )
        ncv = n;

    // v containts the lanczos basis vectors
    int ldv = n;
    double *v = new double[ldv*ncv];

    int *iparam = new int[11];
    iparam[0] = 1; // Specifies the shift strategy (1->exact)
    iparam[2] = 3*n; // Maximum number of iterations
    iparam[6] = 1; /* Sets the mode of dsaupd.
                      1 is exact shifting,
                      2 is user-supplied shifts,
                      3 is shift-invert mode,
                      4 is buckling mode,
                      5 is Cayley mode. */

    int *ipntr = new int[11]; /* Indicates the locations in the work array workd
                                 where the input and output vectors in the
                                 callback routine are located. */

    // array used for reverse communication
    double *workd = new double[3*n];
    for(int i=0;i<3*n;i++)
        workd[i] = 0;

    int lworkl = ncv*(ncv+8); /* Length of the workl array */
    double *workl = new double[lworkl];

    // info = 0: random start vector is used
    int info = 0; /* Passes convergence information out of the iteration
                     routine. */

    // rvec == 0 : calculate only eigenvalue
    // rvec > 0 : calculate eigenvalue and eigenvector
    int rvec = 1;

    // how many eigenvectors to calculate: 'A' => nev eigenvectors
    char howmny = 'A';

    int *select;
    // when howmny == 'A', this is used as workspace to reorder the eigenvectors
    if( howmny == 'A' )
        select = new int[ncv];

    // This vector will return the eigenvalues from the second routine, dseupd.
    //double *d = new double[nev];

    //double *z = 0;

    //if( rvec )
        //z = new double[n*nev];

    // not used if iparam[6] == 1
    double sigma;

    // first iteration
    dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    while( ido != 99 )
    {
        // matrix-vector multiplication
        mvprod(workd+ipntr[0]-1, workd+ipntr[1]-1);
        dsaupd_(&ido, &bmat, &n, &which[0], &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
    }

    if( info < 0 )
        std::cerr << "Error with dsaupd, info = " << info << std::endl;
    else if ( info == 1 )
        std::cerr << "Maximum number of Lanczos iterations reached." << std::endl;
    else if ( info == 3 )
        std::cerr << "No shifts could be applied during implicit Arnoldi update, try increasing NCV." << std::endl;

    dseupd_(&rvec, &howmny, select, eigval.getpointer(), eigvec.getpointer(), &ldv, &sigma, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    if ( info != 0 )
        std::cerr << "Error with dseupd, info = " << info << std::endl;

    // use something to store the result before deleting...
    //sigma = d[0];

    delete [] resid;
    delete [] v;
    delete [] iparam;
    delete [] ipntr;
    delete [] workd;
    delete [] workl;
    //delete [] d;

    //if( rvec )
        //delete [] z;

    if( howmny == 'A' )
        delete [] select;

    //return sigma; // lowest eigenvalue
}

void CI_Big::construct_CI_matrix(SparseMatrix_CRS & mat , int start_l , int end_l){
    //full matrix is to big so we construct only the part that fits in memory the rest is calculated on the fly.
    Permutator_File my_perm(_perm->get_filename(),_ham->getL());
    if (CIMethod_debug)
        cout << "Started construct_CI_matrix in FCI_File with determinantfile: " << my_perm.get_filename() << "between lines: " << start_l << " and " << end_l << endl;
	TYPE dets_row[2] = {0};
	TYPE dets_col[2] = {0};
	for(unsigned long i = start_l ; i < end_l ; i ++){
        mat.NewRow();
        if( i % 10000 == 0)
            cout << "line : " << i << " of lines: " << _perm->get_dim() << endl;
		my_perm.permutate_bit(dets_row,i);
		TYPE up1 = dets_row[0];
		TYPE down1 = dets_row[1];
		for(unsigned long j = i ; j < _mat.getn() ; j ++){
			my_perm.permutate_bit(dets_col,j);
			TYPE up2 = dets_col[0];
		    TYPE down2 = dets_col[1];
			mat.PushToRow(j , get_ham_element(up1,down1 , up2, down2));
	} // end cols 
	} // end rows 
    //_mat.NewRow(); Last new row is added in sparsematrix class in function add_from_list, change because of parallel implementation.
    if (CIMethod_debug)
        cout << "We constructed the Hamiltonian matrix with "<< mat.datasize() << "nonzero elements. " << endl;
}

int CI_Big::get_dim(){ return _perm->get_dim(); } 

std::string FCI_File::get_name() {return "FCI_File" + _perm->get_clean_detfile() ;}
std::string CI_Big::get_name() {return "CI_Big" + _perm->get_clean_detfile() ;}

std::unique_ptr<CIDens> DOCI::get_density(){return std::unique_ptr<CIDens> {new DensDOCI(this)}; }
std::unique_ptr<CIDens> FCI::get_density(){return std::unique_ptr<CIDens> {new DensFCI(this) };}
std::unique_ptr<CIDens> FCI_File::get_density(){return std::unique_ptr<CIDens> {new DensFILE(this) };}
std::unique_ptr<CIDens> CI_Big::get_density(){return std::unique_ptr<CIDens> {new DensFILE(this)  };}

void FCI_File::print_dets()
{
    dynamic_cast<Permutator_File *>( _perm.get())->print_dets();
}

void CI_Big::print_dets()
{
    dynamic_cast<Permutator_File *>( _perm.get())->print_dets();
}
/* vim: set ts=4 sw=4 expandtab :*/
