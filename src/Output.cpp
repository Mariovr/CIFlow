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
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <bitset>
#include <boost/multiprecision/cpp_int.hpp>
#include <regex>
#include <algorithm>

#include "Permutator.h"
#include "CIMethod.h"
#include "matrix.h"
#include "Output.h"
#include "CIDens.h"
#include "Properties.h"
#include "Hamiltonian.h"

namespace mp = boost::multiprecision;     // Reduce the typing a bit later...

using namespace std;

Output::Output(CIMethod * cim)
{
	_cim =cim;
    _num = 0;
}

Output::~Output()
{
}

//-------------------------------------------------------------------------------------
OutputSingleFile::OutputSingleFile(CIMethod *cim, string filename):Output(cim)
{
    _filename = filename;
}

void OutputSingleFile::change_filename( const std::string & filen, bool partial)
{
    if(_file.is_open() ) //Check if file is already open, then close it because we change the filename
        _file.close();	// Close file 
    if (partial)
    {
        unsigned found = _filename.find_last_of(".");
        std::string sub = _filename.substr(0,found);
        int pos = sub.find(filen);
        if ( pos != std::string::npos) 
        {
            sub  = _filename.substr(0,pos);
            _num += 1;
            _filename = sub + filen + to_string(_num)  + ".dat";
        }
        else
        {
            std::string name = _cim->get_name();
            transform(name.begin() , name.end(), name.begin() , ::tolower );
            _filename = _cim->get_ham()->get_short_filename() + "output" + name + filen + "0.dat";
        }
    }
    else
        _filename = filen;
}

void OutputSingleFile::prepare_for_writing()
{
    if (! _file.is_open() )
        _file.open(_filename.c_str(), ios::out);	// Opening of file.
	_file.precision(14);
 	assert(_file.is_open());
	print_psi4_input();
	_file << _cim->get_ex_info(); // Write the information of the problem.
}

void OutputSingleFile::print_energy(){
	_file << setprecision(16) << "#The groundstate energy = " << _cim->get_ci_energy() << endl ;
}

void OutputSingleFile::print_properties(vector<string> props ){
    Properties prop { _cim };
    for( int i = 0 ; i < _cim->_eigvec.getm() ; i++)
    {   
        for(std::string & which : props)
            _file << setprecision(16) << "#Property eigvec : " << i << " "  << which << " = " <<  prop.get_property(which, i) << endl ;
    }
}

void OutputSingleFile::print_output(std::vector<std::string> props, int num , bool all)
{
    prepare_for_writing() ; print_energy(); print_solutions(); print_properties(props);  
    if(all)//Print all eigenvectors.
        for( int i = 0 ; i < _cim->_eigvec.getm() ; i++)
        {   
            _file << "#Eigenvector : " << i << std::endl ;
            print_ci_vec(i ); 
        }
    else
    {
        print_ci_vec(num); 
    }
}

void OutputSingleFile::close_stream()
{
    _file.close();	// Close file 
}

void OutputSingleFile::print_ci_vec(int num , double nauw){ //If a CIMethod has no implementation of an output for a particular civec
	_file << "#start eigenvector " << num << endl;
	for(int i = 0 ; i < _cim->get_dim() ; i ++){
		if(nauw < fabs(_cim->get_eigvec(num,i)))
			_file << i << "\t" << _cim->get_eigvec(num,i) <<endl ;
	}
	_file << "#end eigenvector " << endl;
}

void OutputSingleFile::print_solutions()
{
	//print the energy spectrum
        _file << "#the obtained energy spectrum with the nuclear repulsion energy ("<<  _cim->get_econst() <<") added is: " << endl;
	_file.precision(16);
	_file << "#energies: " << endl;
	for(int i = 0 ; i < _cim->_eigval.getn() ;i ++){
		_file << _cim->_eigval(i,0) + _cim->get_econst()<< endl;
	}
}

void OutputSingleFile::print_rdm(bool twordm)
{
	_file.precision(16);
	_cim->print_one_rdm(_file);
    if(twordm)
    {
        _cim->print_two_rdm(_file);
    }
}

void OutputSingleFile::print_ham(){
	//Please do this only for moderate size Hamiltonians, because get_mat only works for small size Hamiltonians, and it makes no sense to write huge Hamiltonians to a file.
        _file << "#The Hamiltonian matrix is:" << endl;
	_file.precision(14);
        for(int i=0;i< _cim->get_dim();i++)
            for(int j=0;j<_cim->get_dim();j++)
                _file << "(" << i << "," << "\t" << j << ")\t" << _cim->get_mat(i,j) << std::endl;
	_file << "#end Hamiltonian matrix" << endl;
}

void OutputSingleFile::print_psi4_input(string psi_input){
	string line;
	_file <<"#The one and two electron integrals were generated by psi4 using the following inputfile: " << psi_input  <<endl;
	_file <<"--------------------------------------------------------------------------" <<endl;
	ifstream psifile (psi_input.c_str());
	psifile.precision(16);
	if (psifile.is_open()){
		while(getline(psifile , line)){
			if( line.find('#')  == string::npos){
			//line.insert(0,"#psi_input#\t");
			_file << line << endl;}
		}
		psifile.close();
	}
	else cerr << "Unable to open " << psi_input << endl;
	_file <<"--------------------------------------------------------------------------" <<endl;
	_file << endl;
}

OutputSingleFile::~OutputSingleFile()
{
    if( _file.is_open())
        _file.close();	// Close file when object is destroyed.
}
//------------------------------------------------------------------------------------
OutputSFDOCI::OutputSFDOCI(CIMethod * cim, string filename):OutputSingleFile(cim,filename){
}

void OutputSFDOCI::print_ci_vec(int numvec, double nauw){
	// Writing the header of the file.
	_file << "#CIMETHOD=  DOCI\n";
	_file << "#Start  DOCI vec " << numvec << " only significant determinants are kept (abs(coef) > " << nauw <<" )" << endl;
	_file << "#dimension DOCI " << " is: " << _cim->get_dim() << endl;
	_file << "#index\tupbitstring\tdoci_coef" << endl;
	TYPE up =  _cim->_perm->get_start_int( _cim->gNup());
	int count = 0;
	int index = 0;
        for( long i = 0; i < _cim->gUpDim() ; i++ ) { 
		if(fabs(_cim->get_eigvec(i,numvec)) >= nauw){
			_file << index << "\t" << setfill('0') << setw(_cim->get_perm()->get_norb()) <<_cim->_perm->binary(up) << "\t" <<_cim->get_eigvec(i,numvec) << endl;
			up =  _cim->_perm->permutate_bit(up);
			count += 1;
		}
		index ++;
	}
       _file << "#Significant determinants: " << count << endl;
       _file << "#End vecoutput" << endl;
}

//-------------------------------------------------------------------------------------
OutputSFFCI::OutputSFFCI(CIMethod * cim, string filename):OutputSingleFile(cim, filename){

}

//numvec determines the eigenvector (mostly 0-> groundstate , 1->first excited,...)
//nauw is the criterium for the size of the amplitude of the slaterdeterminant to be outputed.
void OutputSFFCI::print_ci_vec(int numvec, double nauw){
	_file <<"#CIMETHOD= FCI\n";
	_file << "#Start  FCI vec " << numvec << " only significant determinants are kept (abs(coef) > " << nauw <<" )" << endl;
	_file << "#dimension FCI" << " is: " << _cim->get_dim() << endl;
	_file << "#index\tupbitstring\tdownbitstring\tfci_coef" << endl;
	TYPE up =  _cim->_perm->get_start_int( _cim->gNup());
	TYPE down =  _cim->_perm->get_start_int( _cim->gNdown());
	int count = 0;
	int index = 0;
        for( long i = 0; i < _cim->gUpDim() ; i++ ) { 
		down =  _cim->_perm->get_start_int(_cim->gNdown());
		for( long j = 0; j < _cim->gDownDim() ; j++ ) { 
			if(fabs(_cim->get_eigvec(i*_cim->gDownDim()+j,numvec)) >= nauw){
				_file << index << "\t" << setfill('0') << setw(_cim->get_perm()->get_norb()) << _cim->_perm->binary(up) << "\t" << setfill('0') << setw(_cim->get_perm()->get_norb()) << _cim->_perm->binary(down) << "\t"<<  _cim->get_eigvec(i*_cim->gDownDim()+j,numvec) << endl;
				count += 1;
			}
		index ++;
		down =  _cim->_perm->permutate_bit(down);
		}
		up =  _cim->_perm->permutate_bit(up);
       }		
       _file << "#Significant determinants: " << count <<  " of " << _cim->get_dim() << endl;
       _file << "#End vecoutput" << endl;
}
//-------------------------------------------------------------------------------------
OutputSFFCI_File::OutputSFFCI_File(CIMethod * cim, string filename):OutputSingleFile(cim, filename){

}

void OutputSFFCI_File::print_ci_vec(int numvec, double nauw){
	_file <<"#CIMETHOD= configuration interaction of all the determinants in " << _cim->_perm->get_clean_detfile() << "\n";
	_file << "#Start  FCI_File vec " << numvec << " only significant determinants are kept (abs(coef) > " << nauw <<" )" << endl;
	_file << "#CI dimension in file" << " is: " << _cim->get_dim() << endl;
	_file << "#index\tupbitstring\tdownbitstring\tfci_coef" << endl;
	TYPE * determinant = new TYPE[2];
	int count = 0;
        for( int pos = 0; pos < _cim->get_dim() ; pos++ ) { 
		        _cim->_perm->permutate_bit(determinant, pos);
			if(fabs(_cim->get_eigvec(pos,numvec)) >= nauw){
				_file << pos << "\t" << setfill('0') << setw(_cim->get_perm()->get_norb()) << _cim->_perm->binary(determinant[0]) << "\t" << setfill('0') << setw(_cim->get_perm()->get_norb()) << _cim->_perm->binary(determinant[1]) << "\t"<<  _cim->get_eigvec(pos,numvec) << endl;
				count += 1;
			}
	}		
	delete [] determinant;
       _file << "#Significant determinants: " << count <<  " of " << _cim->get_dim() << endl;
       _file << "#End vecoutput" << endl;
}


OutputHDF5::OutputHDF5(CIMethod *cim, string filename):Output(cim){
	_file.open(filename.c_str());	// Open file
	assert(_file.is_open());
}

OutputHDF5::~OutputHDF5(){
	_file.close();
}

