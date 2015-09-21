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
#ifndef PERMUTATOR_H
#define PERMUTATOR_H

#include <vector>
#include <string>

#include "Options.h" //definition of the integer TYPE of the bit representation based on the number of orbitals.

class Permutator{
	public:
		Permutator(int norb);
		virtual ~Permutator(){};
		virtual TYPE permutate_bit(TYPE) = 0;
		virtual TYPE get_start_int(int n) = 0;
		virtual void  permutate_bit(TYPE * det , int pos) = 0 ;
		virtual int get_dim() = 0;
		virtual std::string get_filename() = 0;
		std::string binary(TYPE x);
		unsigned popcount(TYPE num);
		std::string get_clean_detfile();
	        int get_norb(){return _norb;}
	private:
		unsigned int _norb;
};

class Permutator_Bit:public Permutator{
	//used for DOCI (1 permutator) and FCI (2 permutators)
	public:
		Permutator_Bit(int norb):Permutator(norb){}
		int bitcount(TYPE x);
		TYPE permutate_bit();
		TYPE permutate_bit(TYPE);
		TYPE get_start_int(int n); 
		void set_last_bit(TYPE newbit){_last_bit = newbit;}
		//void reset(){_last_bit = get_start_int(__builtin_popcountll(_last_bit));}
		TYPE get_bit(){return _last_bit;}
		// for compatibility (find a better solution)
		int get_dim(){return -1;}
		void  permutate_bit(TYPE * det, int pos );
		void  print_dets(int nup , int orbs);

		std::string get_filename(){return "permbit";}
		std::string binary_string(TYPE x);

	private:
		TYPE _last_bit;
};

class Permutator_File: public Permutator{
	public:
		Permutator_File(std::string file,int norb);
		void  permutate_bit(TYPE * det , int pos);
		void set_last_bit(TYPE up , TYPE down);
		void set_pos(int pos){ _pos = pos;}
		//void reset(){_pos = 0;}
		void read_file();
		TYPE get_up_bit(){return _up_dets[_pos];}
		TYPE get_down_bit(){return _down_dets[_pos];}
		void print_dets();
		int get_dim(){return _up_dets.size();}
		TYPE permutate_bit(TYPE d){return 0;}
		TYPE get_start_int(int n){return 0;} 
		std::string get_filename(){return _filename;}
		void set_filename(std::string filename){_filename = filename; read_file();}
	        std::string binary_string(int pos);

	
	private:
		std::vector<TYPE> _up_dets;
		std::vector<TYPE> _down_dets;
		unsigned int _pos; //length of ups and downs is the same, -> pos determines where we are in both vectors
		unsigned int _norb;
		std::string _filename;
};
#endif
