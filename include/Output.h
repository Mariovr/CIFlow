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
#ifndef __OUTPUT_HPP__ 
#define __OUTPUT_HPP__ 

#include <fstream> 
#include <string>
#include <iostream> 

// forward declarations 
class CIMethod; 
class CIDens;

class Output {
public:
Output(CIMethod * cim);

virtual ~Output();


protected:
CIMethod * _cim;
    int _num ;
};

//______________________________________________________________________________
class OutputSingleFile : public Output
{
public:
	OutputSingleFile(CIMethod * cim, std::string filename);
	virtual ~OutputSingleFile();

	virtual void print_ci_vec(int num = 0 , double nauw = -1.);
	void print_energy();
    void print_properties(std::vector<std::string> props , int num = 0);
	void print_rdm(const CIDens & cid, bool twordm = true);
	void print_solutions();
	void print_ham();
	void print_output(std::vector<std::string> props = {} , int num = 0);
	void print_psi4_input(std::string psi_input = "input.dat");
    void prepare_for_writing();
    void change_filename(const std::string & filen, bool partial = true);
    void close_stream();
	//void printOrbitals();

protected:
	std::ofstream _file;
	std::string _filename;
};

class OutputSFDOCI: public OutputSingleFile
{
	public:
		OutputSFDOCI(CIMethod* cim, std::string filename);
		void print_ci_vec(int num = 0 , double nauw = -1.);
	
};

class OutputSFFCI : public OutputSingleFile
{
	public:
		OutputSFFCI(CIMethod* cim, std::string filename);
		void print_ci_vec(int num = 0 , double nauw = -1.);

};

class OutputSFFCI_File : public OutputSingleFile
{
	public:
		OutputSFFCI_File(CIMethod *cim , std::string filename );
		void print_ci_vec(int num = 0 , double nauw = -1.);
	private: 
		std::string _dets;
};

class OutputHDF5 : public Output
{
	public:
		OutputHDF5(CIMethod* cim, std::string filename);
		virtual ~OutputHDF5();

		void print_ci_vec(){std::cout << " d" ;}

	private:
		std::ofstream _file;
};

#endif // __OUTPUT_HPP__
