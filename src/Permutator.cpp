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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string>
#include <assert.h>
#include <bitset>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "Permutator.h"
using namespace std;

Permutator::Permutator(int norb):_norb(norb){
	/*
	if(norb <= 64){
		typedef unsigned long TYPE;
	}
	else if(norb <= 128){
		typedef mp::uint128_t TYPE;
	}
	else if(norb <= 256){
		typedef mp::uint256_t TYPE;
	}
	else{
		assert(_norb <= 512) , cout << "ERROR: make sure number of orbitals is smaller than 512 _norb is now :" << _norb << endl;
		typedef mp::uint512_t TYPE;
	}
	*/
}

std::string Permutator::binary(TYPE x){
    std::string s = "";
    do {
        if(x & 1){
	    s.push_back('1');
	}
	else{
            s.push_back('0');
	}
    } while (x >>= 1);
    std::reverse(s.begin(), s.end());
    return s;
}

unsigned Permutator::popcount(TYPE num){
	#if defined TL
	    return __builtin_popcountl(num);
	#else
	unsigned popcount = 0;
	unsigned i = 0;
	unsigned end = 0;
	try{
	    i = lsb(num);
	    end = msb(num);
	}
	catch(range_error){
	    i = 0;
	    end = 0;
	}    
	for (; i <= end; ++i) {
	    if (bit_test(num,i) != 0) {
		++popcount;
	    }
	}
	return popcount;
	#endif
}
    
std::string Permutator::get_clean_detfile()
{
    std::string detfile = get_filename();
    unsigned found = detfile.find_last_of("/\\");
    unsigned found2 = detfile.find_last_of(".");
    detfile = detfile.substr(found+1,found2-found-1);
    char chars[] = "_()-./\\";
    for (unsigned int i = 0; i < strlen(chars); ++i)
    {
	detfile.erase(std::remove(detfile.begin(), detfile.end(), chars[i]), detfile.end());
    }
    return detfile;
}

std::string Permutator_Bit::binary_string(TYPE x){
    //return binary(_last_bit);	
    return binary(x);	
}    

TYPE Permutator_Bit::permutate_bit(TYPE d ){
	//to keep compatibility with previous implementation
	_last_bit = d;
	return permutate_bit();
}

TYPE Permutator_Bit::permutate_bit(){
#if defined TL
        unsigned long t = _last_bit | (_last_bit - 1); // t gets v's least significant 0 bits set to 1
        // Next set to 1 the most significant bit to change,
        // set to 0 the least significant ones, and add the necessary 1 bits.
        _last_bit = (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctzll(_last_bit) + 1));

#else
        TYPE t = (_last_bit | (_last_bit - 1)) +1; // t gets v's least significant 0 bits set to 1
        // Next set to 1 the most significant bit to change,
        // set to 0 the least significant ones, and add the necessary 1 bits.
        unsigned int lowestsetb = 0;
        try{
            lowestsetb = lsb(t) ;
        }    
        catch(range_error){
            lowestsetb = 0;
        }    
        unsigned int lastb  = 0;
        try{
            lastb = lsb(_last_bit);
        }    
        catch(range_error){
            lastb = 0;
        }    
        TYPE one(1);
        TYPE val1 = one << lowestsetb;
        TYPE val2 = one << lastb;
        _last_bit = _last_bit==0 ? 0 : t | (((val1 / val2) >> 1) - 1);
#endif
    return _last_bit;
}

TYPE Permutator_Bit::get_start_int(int n){
    //n number of bits set to 1, returns the start integer for the permutator with the n first bits set to 1
    // save as long to be sure we have a long int
    _last_bit = 1;
    _last_bit = _last_bit << n;
    _last_bit -= 1;
    return _last_bit;
}

void  Permutator_Bit::permutate_bit(TYPE * det, int pos ){std::cout << "you shouldn't be here" ;} 

void Permutator_Bit::print_dets(int nup , int orbs){
	cout << "Determinants contained in the permutator Bit object:" <<endl;
    TYPE start = get_start_int(nup);
    TYPE shiftbit = 1;
    int it = 0;
	while ( !(start & (1<< orbs))  ){
		cout << binary_string(start)  << endl;
        start = permutate_bit(start);
        it += 1;
	}
	cout << "Number of determinants: " << it << endl;
}

//-------------------------class Permutator_File-----------------------------------------------


Permutator_File::Permutator_File(string file, int norb): Permutator(norb){
	_filename = file;
	read_file();
	_pos = 0;
}

void Permutator_File::read_file(){
	_up_dets.resize(0);
	_down_dets.resize(0);
	string line;
	ifstream detfile (_filename.c_str());
	if (detfile.is_open()){
		while(getline(detfile , line)){
			if( line.find('#')  == string::npos){
				int cut_pos = line.find("|");
				int num_convert = ceil(cut_pos/64.);//division of integers rounds automatically.
				int cur_pos = 0;
				//cout << line.substr(0,cut_pos).c_str() <<"   " <<line.substr(cut_pos+1,string::npos).c_str()<< endl;
				TYPE up = 0;
				TYPE down = 0;
				for(int i = num_convert-1 ; i > -1 ; i--){
					if( i== 0){
						//cout << " cup  " << cut_pos << " cp " << cur_pos << " nc " << num_convert << endl;
						up += TYPE( strtoull(line.substr(cur_pos,cut_pos-(num_convert-1)*64).c_str(),NULL,2) );
						down +=TYPE(strtoull(line.substr(cut_pos+1+(num_convert-1)*64,string::npos).c_str(),NULL,2) ) ;
						//cout << " up  " << up << " down " << down << endl;
					}
					else{
						up += TYPE (strtoull(line.substr(cur_pos,64).c_str(),NULL,2)) *TYPE (pow(2,(cut_pos - 64.*i)) );
						down += TYPE (strtoull(line.substr(cut_pos+1+64*(num_convert-1 - i),64).c_str(),NULL,2)) * TYPE(pow(2,(cut_pos - 64.*i)) );

					}
					cur_pos += 64;

				}
				_up_dets.push_back(up);
				_down_dets.push_back(down);

			}
		}
		detfile.close();
	}
	else cerr << "Unable to open " << _filename;
#ifdef CIMethod_debug
	cout << "We read in the determinants from file: "<< _filename << endl;
#endif
}

void  Permutator_File::permutate_bit(TYPE *  det, int pos){
	det[0] = _up_dets[pos]; 
	det[1] =  _down_dets[pos];
}

void Permutator_File::set_last_bit(TYPE up , TYPE down){
	for (unsigned int it = 0 ;it < _up_dets.size() ; ++it){
		if (_up_dets[it] == up && _down_dets[it] == down){
			_pos = it;
			break ;
		}

	}
}


std::string Permutator_File::binary_string(int pos)
{
    std::stringstream ups;  
    std::stringstream downs;
    ups << setfill('0') << setw(_norb) << binary(_up_dets[pos]);
    downs  << setfill('0') << setw(_norb) << binary(_down_dets[pos]);
    std::string total = ups.str()+ '|' + downs.str();
    return total;
}

void Permutator_File::print_dets(){
	cout << "Determinants contained in the permutator File object:" <<endl;
	for (int it = 0 ; it != _up_dets.size() ; it++){
		cout << binary_string(it)  << endl;
	}
	cout << "Number of determinants: " << _up_dets.size() << endl;
}
