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
#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm> //std::sort
#include <fstream>
#include <sstream>

#include "Irreps.h"
#include "TwoIndex.h"
#include "FourIndex.h"
#include "Hamiltonian.h"
#include "MyHDF5.h"
#include "Options.h"
#include "UnitaryMatrix.h"
#include "scpp_assert.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;

Hamiltonian::Hamiltonian(const int Norbitals, const int nGroup, const int * OrbIrreps, int nup , int ndown, double Econstant): _modham{false}{

    Econst = Econstant;
    L = Norbitals;
    if ((nGroup<0) || (nGroup>7)){
       cout << "Error at Hamiltonian::Hamiltonian : nGroup out of bound." << endl;
    }
    SymmInfo.setGroup(nGroup);
    
    orb2irrep = new int[L];
    orb2indexSy = new int[L];
    int nIrreps = SymmInfo.getNumberOfIrreps();
    irrep2num_orb = new int[nIrreps];
    for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
    for (int cnt=0; cnt<L; cnt++){
       if ((OrbIrreps[cnt]<0) || (OrbIrreps[cnt]>=nIrreps)){
          cout << "Error at Hamiltonian::Hamiltonian : OrbIrreps[cnt] out of bound." << endl;
       }
       orb2irrep[cnt] = OrbIrreps[cnt];
       orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
       irrep2num_orb[orb2irrep[cnt]]++;
    }
    set_cum(); //sets the cummulative int array to determine how many orb indices are before the start of an irrep.
    
    Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
    Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);

    _nup = nup;
    _ndown = ndown;
    //SCPP_ASSERT(_nup == _ndown , "Error: At this moment CIFlow only works when nup = ndown, make sure this is the case. nup = " << _nup << " ndown = " << _ndown << std::endl);
    _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    _filename = "";

    //Orthogonal atomic orbitals.(for example for model hamiltonians.)
    _overlap = std::vector<double>(  0. );
    _unit.reset(nullptr);
    _hf_energy = 0.;
}

Hamiltonian::Hamiltonian(const int Norbitals, const int nGroup, int nup , int ndown, double Econstant): _modham{false}{
    //Associateds all orbitals to 0 irrep.

    Econst = Econstant;
    L = Norbitals;
    if (nGroup != 0 ){
       cerr << "Error at Hamiltonian::Hamiltonian : nGroup out of bound." << endl;
       cerr << "Make sure the group with this constructor is 0 because we make all orbitals correspond to irrep 0." << endl;
       exit(1);
    }
    SymmInfo.setGroup(nGroup);
    
    orb2irrep = new int[L];
    orb2indexSy = new int[L];
    int nIrreps = SymmInfo.getNumberOfIrreps();
    irrep2num_orb = new int[nIrreps];
    for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
    for (int cnt=0; cnt<L; cnt++){
       orb2irrep[cnt] = 0;
       orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
       irrep2num_orb[orb2irrep[cnt]]++;
    }
    set_cum(); //sets the cummulative int array to determine how many orb indices are before the start of an irrep.
    
    Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
    Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);

    _nup = nup;
    _ndown = ndown;
    SCPP_ASSERT(_nup == _ndown , "Error: At this moment CIFlow only works when nup = ndown, make sure this is the case. nup = " << _nup << " ndown = " << _ndown << std::endl);
    _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    _filename = "";
    //Orthogonal atomic orbitals.(for example for model hamiltonians.)
    _overlap = std::vector<double>(  0. );
    _unit.reset(nullptr);
    _hf_energy = 0;
}

Hamiltonian::Hamiltonian(const Hamiltonian & HamIn){
	int numberOfIrreps = HamIn.getNirreps();
        int SyGroup = HamIn.getNGroup();
        SymmInfo.setGroup(SyGroup);
	L = HamIn.getL();
        Econst = HamIn.getEconst();
	_nup = HamIn.getnup();
	_ndown = HamIn.getndown();
    //SCPP_ASSERT(_nup == _ndown , "Error: At this moment CIFlow only works when nup = ndown, make sure this is the case. nup = " << _nup << " ndown = " << _ndown << std::endl);
	_oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    _hf_energy = HamIn.get_hf_energy();

        orb2irrep = new int[L];
	orb2indexSy = new int[L];
	irrep2num_orb = new int[numberOfIrreps];
	for (int cnt=0; cnt<numberOfIrreps; cnt++) irrep2num_orb[cnt] = 0;
        for (int cnt=0; cnt<L; cnt++){
		orb2irrep[cnt] = HamIn.getOrbitalIrrep(cnt); 
		if ((orb2irrep[cnt]<0) || (orb2irrep[cnt]>=numberOfIrreps)){
		    cout << "Error at Hamiltonian::Hamiltonian : OrbIrreps[cnt] out of bound." << endl;
		}
		orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
		irrep2num_orb[orb2irrep[cnt]]++;
        }	
        set_cum(); //sets the cummulative int array to determine how many orb indices are before the start of an irrep.
       Tmat = new TwoIndex(HamIn.get_tmat_pointer());
       Vmat = new FourIndex(HamIn.get_vmat_pointer());
       _filename = HamIn.get_filename();

       _modham = HamIn.get_modham(); 
       if(HamIn.get_overlap().size() == L *L)
       {
           _overlap = std::vector<double>( L*L, 0. );
           for(int i = 0 ; i < L*L ; i++)
               _overlap[i] = HamIn.get_overlap()[i];
       }
       if(HamIn.get_unitary())
           _unit.reset(new UnitaryMatrix(*HamIn.get_unitary()) );
       else
           _unit.reset(nullptr);
}

Hamiltonian::Hamiltonian(const string filename) : _modham{false}
{
    int last_index = filename.find_last_of('.');
    //Check if inputfile is in hdf5 format.
    if("h5" == filename.substr(last_index+1 , string::npos))
    {
        read(filename);
    }
    else
    {
    	read_file(filename);
    }
    _filename = filename; 
    SCPP_ASSERT(_nup == _ndown , "Error: At this moment CIFlow only works when nup = ndown, make sure this is the case. nup = " << _nup << " ndown = " << _ndown << std::endl);
}

Hamiltonian& Hamiltonian::operator=(const Hamiltonian &orig)  
{
    L = orig.getL();
    _nup = orig.getnup();
    _ndown = orig.getndown();
    _hf_energy = orig.get_hf_energy();
    SCPP_ASSERT(_nup == _ndown , "Error: At this moment CIFlow only works when nup = ndown, make sure this is the case. nup = " << _nup << " ndown = " << _ndown << std::endl);
    _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    int SyGroup = orig.getNGroup();
    SymmInfo.setGroup(SyGroup);
    Econst = orig.getEconst();

    orb2irrep = new int[L];
    orb2indexSy = new int[L];
    int nIrreps = SymmInfo.getNumberOfIrreps();
    irrep2num_orb = new int[nIrreps];

    for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
    for (int cnt=0; cnt<L; cnt++){
	orb2irrep[cnt] = orig.getOrbitalIrrep(cnt); 
	if ((orb2irrep[cnt]<0) || (orb2irrep[cnt]>=nIrreps)){
	   cout << "Error at Hamiltonian::Hamiltonian : OrbIrreps[cnt] out of bound." << endl;
	}
      orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
      irrep2num_orb[orb2irrep[cnt]]++;
    }

    set_cum();

    Tmat = new TwoIndex(orig.get_tmat_pointer() );
    Vmat = new FourIndex(orig.get_vmat_pointer() );
    _filename = orig.get_filename();
    _modham = orig.get_modham(); 
   if(orig.get_overlap().size() == L *L)
   {
        for(int i = 0 ; i < L*L ; i++)
            _overlap[i] = orig.get_overlap()[i];
   }
   if(orig.get_unitary())
       _unit.reset(new UnitaryMatrix(*orig.get_unitary()) );
   else
       _unit.reset(nullptr);
    return *this;
}

Hamiltonian:: ~Hamiltonian(){
   
   delete [] orb2irrep;
   delete [] orb2indexSy;
   delete [] irrep2num_orb;
   delete Tmat;
   delete Vmat;

   if(_norbcum)
       delete [] _norbcum;
}

std::string Hamiltonian::get_info() const
{
    std::ostringstream info;
    info << std::setprecision(16);
    info << "#nup = " << getnup() << "\n#ndown = " << getndown() << "\n#norbs = "<< getL() << "\n#point group symmetry = " << getNGroup() << "\n#constant energy (nuc. rep.) = " << Econst << "\n#HF energy (energy of 0..01..1 determinant): " <<  _hf_energy << "\n#filename = " << _filename << std::endl;
    return info.str();
}

std::string Hamiltonian::get_short_filename() const
{
    unsigned found = _filename.find_last_of("/\\");
    unsigned found2 = _filename.find_last_of(".");
    return _filename.substr(found+1,found2-found-1);
}

int Hamiltonian::getL() const{ return L; }

int Hamiltonian::getNGroup() const{ return SymmInfo.getGroupNumber(); }

int Hamiltonian::getOrbitalIrrep(const int nOrb) const{ return orb2irrep[nOrb]; }

void Hamiltonian::setEconst(const double val){ Econst = val; }

double Hamiltonian::getEconst() const{ return Econst; }

void Hamiltonian::setTmat(const int index1, const int index2, const double val){

   if (orb2irrep[index1]==orb2irrep[index2]){
      Tmat->set(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2], val);
      return;
   } else {
      if (val==0.0) return;
      cout << "Error at Hamiltonian::setTmat: sy1 != sy2 and val != 0.0" << endl;
   }

}

double Hamiltonian::getTmat(const int index1, const int index2) const{

   if (orb2irrep[index1]==orb2irrep[index2]){
      return Tmat->get(orb2irrep[index1], orb2indexSy[index1], orb2indexSy[index2]);
   }

   return 0.0;
   
}

void Hamiltonian::setVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   if ( SymmInfo.directProd(orb2irrep[index1],orb2irrep[index2]) == SymmInfo.directProd(orb2irrep[index3],orb2irrep[index4]) ){
      Vmat->set(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);
      return;
   } else {
      if (val==0.0) return;
      cout << "Error at Hamiltonian::setVmat: sy1 x sy2 != sy3 x sy4 and val != 0.0" << endl;
   }

}

void Hamiltonian::addToVmat(const int index1, const int index2, const int index3, const int index4, const double val){

   if ( SymmInfo.directProd(orb2irrep[index1],orb2irrep[index2]) == SymmInfo.directProd(orb2irrep[index3],orb2irrep[index4]) ){
      Vmat->add(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4], val);
      return;
   } else {
      if (val==0.0) return;
      cout << "Error at Hamiltonian::addToVmat: sy1 x sy2 != sy3 x sy4 and val != 0.0" << endl;
   }

}

double Hamiltonian::getVmat(const int index1, const int index2, const int index3, const int index4) const{

    if( !_modham)
    {
       if ( SymmInfo.directProd(orb2irrep[index1],orb2irrep[index2]) == SymmInfo.directProd(orb2irrep[index3],orb2irrep[index4]) ){
          return Vmat->get(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4]);
       }

       return 0.0;
    } //The pair Hamiltonians miss the <ij|ji> symmetry of the <ibia|jbja>, <jbja|ibia> , <iaib|jajb>, ... so we need an extra check, to avoid this symmetry, and put the corresponding matrix elements to zero.
    else
    {
        if ( SymmInfo.directProd(orb2irrep[index1],orb2irrep[index2]) == SymmInfo.directProd(orb2irrep[index3],orb2irrep[index4]) && index1 == index2){
           return Vmat->get(orb2irrep[index1], orb2irrep[index2], orb2irrep[index3], orb2irrep[index4], orb2indexSy[index1], orb2indexSy[index2], orb2indexSy[index3], orb2indexSy[index4]);
        }

        return 0.0;
    }
   
}

void Hamiltonian::save(const string filename) const{
    const string filenames = "ham" +  get_short_filename() + filename + ".h5";
   //The hdf5 file
   hid_t file_id = H5Fcreate(filenames.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
         //The chain length
         hsize_t dimarray       = 1;
         hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
         hid_t dataset_id       = H5Dcreate(group_id, "L", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);
         
         //The group number
         hsize_t dimarray2      = 1;
         hid_t dataspace_id2    = H5Screate_simple(1, &dimarray2, NULL);
         hid_t dataset_id2      = H5Dcreate(group_id, "nGroup", H5T_STD_I32LE, dataspace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         int nGroup = SymmInfo.getGroupNumber();
         H5Dwrite(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
         
         //orb2irrep
         hsize_t dimarray3      = L;
         hid_t dataspace_id3    = H5Screate_simple(1, &dimarray3, NULL);
         hid_t dataset_id3      = H5Dcreate(group_id, "orb2irrep", H5T_STD_I32LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);
         
         //Econst
         hsize_t dimarray4      = 1;
         hid_t dataspace_id4    = H5Screate_simple(1, &dimarray4, NULL);
         hid_t dataset_id4      = H5Dcreate(group_id, "Econst", H5T_IEEE_F64LE, dataspace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);

	 //nalpha
         hsize_t dimarray5      = 1;
         hid_t dataspace_id5    = H5Screate_simple(1, &dimarray5, NULL);
         hid_t dataset_id5      = H5Dcreate(group_id, "nup", H5T_STD_I32LE, dataspace_id5, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id5, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_nup);

         //nbeta
         hsize_t dimarray6      = 1;
         hid_t dataspace_id6    = H5Screate_simple(1, &dimarray6, NULL);
         hid_t dataset_id6      = H5Dcreate(group_id, "ndown", H5T_STD_I32LE, dataspace_id6, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id6, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_ndown);
    
         //overlap
         hsize_t dimarray7      = _overlap.size();
         hid_t dataspace_id7    = H5Screate_simple(1, &dimarray7, NULL);
         hid_t dataset_id7      = H5Dcreate(group_id, "overlap", H5T_IEEE_F64LE, dataspace_id7, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
         H5Dwrite(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_overlap[0]);

         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);
         H5Dclose(dataset_id2);
         H5Sclose(dataspace_id2);
         H5Dclose(dataset_id3);
         H5Sclose(dataspace_id3);
         H5Dclose(dataset_id4);
         H5Sclose(dataspace_id4);
         H5Dclose(dataset_id5);
         H5Sclose(dataspace_id5);
         H5Dclose(dataset_id6);
         H5Sclose(dataspace_id6);
         H5Dclose(dataset_id7);
         H5Sclose(dataspace_id7);

      H5Gclose(group_id);
      
   H5Fclose(file_id);

   Tmat->save(filenames);
   Vmat->save(filenames);

   if(_unit)
   {
       _unit->saveU(filenames+"unitary"); 
       //_unit->print_unitary(filenames+"unitary"); 
   }
}

void Hamiltonian::read(const string filename){

   //The hdf5 file
   hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      
      //The data
      hid_t group_id = H5Gopen(file_id, "/Data",H5P_DEFAULT);
       
         //The chain length
         hid_t dataset_id = H5Dopen(group_id, "L", H5P_DEFAULT);
         H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &L);

         //The group number
         hid_t dataset_id2 = H5Dopen(group_id, "nGroup", H5P_DEFAULT);
         int nGroup;
         H5Dread(dataset_id2, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nGroup);
	     SymmInfo.setGroup(nGroup);
         
         //orb2irrep
         hid_t dataset_id3 = H5Dopen(group_id, "orb2irrep", H5P_DEFAULT);
         orb2irrep = new int[L];
         H5Dread(dataset_id3, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, orb2irrep);//orb2irrep is already address
         orb2indexSy = new int[L];
         int nIrreps = SymmInfo.getNumberOfIrreps();
         irrep2num_orb = new int[nIrreps];
         for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
         for (int cnt=0; cnt<L; cnt++){
            orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
            irrep2num_orb[orb2irrep[cnt]]++;
         }
         
         //Econst
         hid_t dataset_id4 = H5Dopen(group_id, "Econst", H5P_DEFAULT);
         H5Dread(dataset_id4, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Econst);
         
	 //nalpha
         hid_t dataset_id5 = H5Dopen(group_id, "nup", H5P_DEFAULT);
         H5Dread(dataset_id5, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_nup);
	 //nbeta
         hid_t dataset_id6 = H5Dopen(group_id, "ndown", H5P_DEFAULT);
         H5Dread(dataset_id6, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_ndown);
         _oneOverNMinusOne = 1.0/((_ndown+_nup)-1);
    
        try
        {
             //overlap
             hid_t dataset_id7 = H5Dopen(group_id, "overlap", H5P_DEFAULT);
             _overlap = std::vector<double>( L*L, 0. );
             H5Dread(dataset_id7, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_overlap[0]);//orb2irrep is already address
             H5Dclose(dataset_id7);
         }
        catch(...)
         {
             _overlap = std::vector<double>( 0 );
         }



         H5Dclose(dataset_id);
         H5Dclose(dataset_id2);
         H5Dclose(dataset_id3);
         H5Dclose(dataset_id4);
         H5Dclose(dataset_id5);
         H5Dclose(dataset_id6);

      H5Gclose(group_id);
      
   H5Fclose(file_id);

   Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);

   Tmat->read(filename);
   Vmat->read(filename);

   set_cum(); //sets the cummulative int array to determine how many orb indices are before the start of an irrep.
   
   _filename = filename;

   try
   {
       OptIndex opt = OptIndex(L, nGroup ,irrep2num_orb);
       _unit.reset(new UnitaryMatrix(opt));
       _unit->loadU(filename+"unitary");
   }
   catch(...)
   {
       _unit.reset(nullptr);
   }


   if (HAMILTONIAN_debugPrint) debugcheck();

}

void Hamiltonian::load_overlap(const std::string & filename)
{
    ifstream file(filename.c_str() );
    load_overlap(file);
    file.close();
}

void Hamiltonian::load_unitary(const std::string & filename)
{
    int last_index = filename.find_last_of('.');
    //Check if inputfile is in hdf5 format.
    if("h5" == filename.substr(last_index+1 , string::npos))
    {
        _unit->loadU(filename);
    }
    else
    {
    	_unit->load_unitary(filename);
    }
}

void Hamiltonian::set_unitary(const std::vector<double> & unit)
{
    OptIndex opt = get_index_object();
    _unit.reset(new UnitaryMatrix(unit,opt ));
}

bool Hamiltonian::load_overlap(std::istream & file)
{
    //std::ifstream file(overlapname.c_str()); 
    bool roverlap = false;

    //Wind forward to relevant part
    std::string line; 
    int position = file.tellg();
    while( getline(file, line) )
    {
        if(   line.find("CIFlowOverlap") != std::string::npos)
        {
            roverlap = true;
            _overlap = std::vector<double>( L*L, 0. );
            int tot = 0;
            for( int irrep = 0 ; irrep < getNirreps() ; irrep ++)
            {
                int linsize = getNORB(irrep);
                if(linsize > 0)
                {
                    while(line.find("irrep") == std::string::npos && line.find("Irrep") == std::string::npos)
                    {
                        std::getline(file , line);
                    }
                    for( int j = 0 ; j < linsize ; j ++)
                    {
                        std::getline(file, line );
                        std::istringstream sin(line);
                        double coef;
                        for( int l = 0 ; l < linsize ; l ++)
                        {
                            sin >> coef;
                            _overlap[L* (j+tot) + (l+tot) ] =  coef;
                        }

                    }
                    tot += linsize;
                } //End linsize > 0
            }
            break;
        }
    }
    if(!roverlap)
    {
        file.clear();
        file.seekg(position);
    }

    //file.close();
    return roverlap;
}

void Hamiltonian::print_overlap(std::ostream & file)
{
    if(_overlap.size() > 0)
    {
        file << "#CIFlowOverlap:   "<< std::endl;
        int tot = 0;
        for( int irrep = 0 ; irrep < getNirreps() ; irrep ++)
        {
            int linsize = getNORB(irrep);
            if(linsize > 0)
            {
                file << "#irrep : " << irrep << std::endl;
                for( int j = 0 ; j < linsize ; j ++)
                {
                    for( int l = 0 ; l < linsize ; l ++)
                    {
                        file << get_overlap( (j+tot) ,  (l+tot) ) << " ";
                    }
                    file << std::endl;
                }
                tot += linsize;
            } //End linsize > 0
        }
    }
    else
    {
        file << "#### There is no overlap.   "<< std::endl;
    }
}

std::vector<double> Hamiltonian::get_overlap() const
{
    return _overlap;
}


UnitaryMatrix * Hamiltonian::get_unitary() const 
{ 
    return _unit.get();
}

void Hamiltonian::set_overlap(std::vector<double> overlap)
{
    _overlap = overlap;
}

void Hamiltonian::set_overlap(int irrep, int i , int j , double val)
{
    int rowindex = (getNstart(irrep) + i);
    int colindex = (getNstart(irrep) + j);
    _overlap[getL() * rowindex + colindex] = val; //c-style first rows than columns.
}

double Hamiltonian::get_overlap(int irrep , int i, int j )const
{
    int rowindex = (getNstart(irrep) + i);
    int colindex = (getNstart(irrep) + j);
    return _overlap[getL() * rowindex + colindex]; //c-style first rows than columns.
}

double Hamiltonian::get_overlap(int i, int j )const
{
    return _overlap[getL() * i + j]; //c-style first rows than columns.
}

void Hamiltonian::set_overlap(int i , int j , double val)
{
    _overlap[getL() * i+j] = val; //c-style first rows than columns.
}

void Hamiltonian::save_file(const string & filename)
{
   const string filenames = "ham" +  get_short_filename() + filename + ".dat";
   std::ofstream outputfile(filenames.c_str());
   //First go to the start of the integral dump.
   outputfile << "****  Molecular Integrals Start Here \n";
   outputfile << "Nalpha = ";
   outputfile << _nup << " \n";
   outputfile << "Nbeta = "<< _ndown << " \n";
   outputfile << "Symmetry Label = " << SymmInfo.getGroupName() << " \n"; 
   outputfile << "Nirreps = " <<  SymmInfo.getNumberOfIrreps() << " \n" ;
   outputfile.precision(15);
   outputfile << "Nuclear Repulsion Energy = " <<  Econst << " \n" ;
   outputfile << "Number Of Molecular Orbitals = " <<  L << " \n" ;
   outputfile << "Irreps Of Molecular Orbitals = \n" ;
   for(int i = 0 ; i < L; i ++)
   {
       outputfile << orb2irrep[i] << " ";
   }
   outputfile << "\n" ;
   outputfile << "DOCC = " << "" << " \n";
   outputfile << "SOCC = " << "" << " \n";
   print_overlap(outputfile);
   if(_unit)
   {
       _unit->print_unitary(outputfile);
   }
   outputfile << "****  MO OEI \n";
   for(int i = 0 ; i < L ; i++)
   {
       for(int j = 0 ; j < L ; j++)
       {
           if(getTmat(i,j) != 0.)
               outputfile << i << " " << j << " "<< getTmat(i,j) << " \n";
       }
   }

   //write two-electron integrals --> in file: chemical notation; in Vmat: physics notation
   outputfile << "****  MO TEI \n";
   for(int i = 0 ; i < L ; i++)
   {
       for(int j = 0 ; j < L ; j++)
       {
           for(int k = 0 ; k < L ; k++)
           {
               for(int l = 0 ; l < L ; l++)
               {
		   if(getVmat(i,j,k,l) != 0.)
		       outputfile << i << " " << k << " "<< j << " " << l << " " << getVmat(i,j,k ,l) << " \n";
               }
           }
       }
   }
   outputfile << "****  HF Energy = " << _hf_energy << " \n";
   outputfile << "****  Molecular Integrals End Here" << " \n" ; 
   outputfile <<"\n";
   
   outputfile.close();
}

//Works for the file mointegrals/mointegrals.cc_PRINT which can be used as a plugin in psi4 beta5
void Hamiltonian::read_file(const string filename){
   
   string line, part;
   int pos;
   
   ifstream inputfile(filename.c_str());
   _filename = filename;
   //First go to the start of the integral dump.
   string start = "****  Molecular Integrals";
   do{
      getline(inputfile,line);
   } while (line.find(start) == std::string::npos);
   /*  
   // Get the alpha molecular orbital energies. 
   getline(inputfile,line);
   pos = line.find("=");
   pos += 2;
	while(pos!= 0){
			  int pos2 = line.find(" ",pos);
			  double value = atof(line.substr(pos,pos2).c_str());
			  moealpha.push_back(value);
			  pos = pos2 +1;
	}
	//moealpha.pop_back();
	//std::sort(moealpha.begin(), moealpha.end());
	//for(std::vector<double>::iterator it =  moealpha.begin() ; it != moealpha.end() ; it++){
			  //cout << *it << endl;
	//}*/
   ///Get number of alpha electrons
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   _nup = atoi(part.c_str());

   ///Get number of beta electrons
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   _ndown = atoi(part.c_str());

   _oneOverNMinusOne = 1.0/(_nup+_ndown-1.0);

   //Get the group name and convert it to the group number
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   int nGroup = 0;
   bool stop = false;
   do {
      if (part.compare(SymmInfo.getGroupName(nGroup))==0) stop = true;
      else nGroup += 1;
   } while (!stop);
   SymmInfo.setGroup(nGroup);
   //cout << "The group was found to be " << SymmInfo.getGroupName() << " ." << endl;
   
   //This line says how many irreps there are: skip.
   getline(inputfile,line);
   
   //This line contains the nuclear energy part.
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   Econst = atof(part.c_str());

   //This line contains the number of MO's.
   getline(inputfile,line);
   pos = line.find("=");
   part = line.substr(pos+2,line.size()-pos-3);
   L = atoi(part.c_str());
   
   //This line contains only text
   getline(inputfile,line);
   
   //This line contains the irrep numbers --> allocate, read in & set
   getline(inputfile,line);
   
   orb2irrep = new int[L];
   orb2indexSy = new int[L];
   int nIrreps = SymmInfo.getNumberOfIrreps();
   irrep2num_orb = new int[nIrreps];
   
   pos = 0;
   do {
      orb2irrep[pos] = atoi(line.substr(2*pos,1).c_str());
      pos++;
   } while (2*pos < (int)line.size()-1);
   
   for (int cnt=0; cnt<nIrreps; cnt++) irrep2num_orb[cnt] = 0;
   for (int cnt=0; cnt<L; cnt++){
      orb2indexSy[cnt] = irrep2num_orb[orb2irrep[cnt]];
      irrep2num_orb[orb2irrep[cnt]]++;
   }
   Tmat = new TwoIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   Vmat = new FourIndex(SymmInfo.getGroupNumber(),irrep2num_orb);
   
   //print_overlap(std::cout );
   //Skip three lines --> number of double occupations, single occupations and test line
   //getline(inputfile,line);
   //getline(inputfile,line);
   //getline(inputfile,line);
   if( load_overlap( inputfile))
   {
       OptIndex opt = OptIndex(L, nGroup ,irrep2num_orb);
       _unit.reset(new UnitaryMatrix(opt, filename));
       //_unit->print_unitary(std::cout);
   }
   else
   {
       std::cout << "we didn't succeed in reading in an overlap and unitary matrix for this matrixelements." << endl;
       _overlap = std::vector<double> (0);
       _unit.reset(nullptr);
   }
   //Wind till integrals
   stop = false; 
   start = "**";
   do{
      getline(inputfile,line);
      //std::cout << line << endl;
      pos = line.find(start);
      if (pos==0) stop = true;
   } while (!stop);
   
   //Read in one-electron integrals
   getline(inputfile,line);
   int pos2, index1, index2;
   double value;
   while( (line.substr(0,1)).compare("*")!=0 ){
   
      pos = 0;
      pos2 = line.find(" ",pos);
      index1 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index2 = atoi(line.substr(pos,pos2-pos).c_str());
      
      value = atof(line.substr(pos2+1,line.size()-pos2-2).c_str());
      
      setTmat(index1,index2,value);
      
      getline(inputfile,line);
   }
   
   //Read in two-electron integrals --> in file: chemical notation; in Vmat: physics notation
   getline(inputfile,line);
   int index3, index4;
   while( (line.substr(0,1)).compare("*")!=0 ){
   
      pos = 0;
      pos2 = line.find(" ",pos);
      index1 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index2 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index3 = atoi(line.substr(pos,pos2-pos).c_str());
      
      pos = pos2+1;
      pos2 = line.find(" ",pos);
      index4 = atoi(line.substr(pos,pos2-pos).c_str());
      
      value = atof(line.substr(pos2+1,line.size()-pos2-2).c_str());
      
      setVmat(index1, index3, index2, index4, value);
      
      getline(inputfile,line);
   
   }
   pos = line.find("=");
   _hf_energy = atof(line.substr(pos+2,line.size()-pos -3).c_str());
   
   set_cum(); //sets the cummulative int array to determine how many orb indices are before the start of an irrep.
   if (HAMILTONIAN_debugPrint) debugcheck();
   
   inputfile.close();
}

void Hamiltonian::debugcheck() const{

   cout << "Econst = " << Econst << endl;
   
   double test = 0.0;
   double test2 = 0.0;
   double test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getTmat(i,i);
      for (int j=0; j<L; j++){
         test += getTmat(i,j);
         if (i<=j) test2 += getTmat(i,j);
      }
   }
   cout << "1-electron integrals: Trace                  : " << test3 << endl;
   cout << "1-electron integrals: Sum over all elements  : " << test << endl;
   cout << "1-electron integrals: Sum over Tij with i<=j : " << test2 << endl;
      
   test = 0.0;
   test2 = 0.0;
   test3 = 0.0;
   for (int i=0; i<L; i++){
      test3 += getVmat(i,i,i,i);
      for (int j=0; j<L; j++){
         for (int k=0; k<L; k++){
            for (int l=0; l<L; l++){
               test += getVmat(i,j,k,l);
               if ((i<=j) && (j<=k) && (k<=l)) test2 += getVmat(i,j,k,l);
            }
         }
      }
   }
   cout << "2-electron integrals: Trace                          : " << test3 << endl;
   cout << "2-electron integrals: Sum over all elements          : " << test << endl;
   cout << "2-electron integrals: Sum over Vijkl with i<=j<=k<=l : " << test2 << endl;

}

double Hamiltonian::gMxElement(const int alpha, const int beta, const int gamma, const int delta) const{

	return getVmat(alpha, beta, gamma, delta) + _oneOverNMinusOne*(((alpha==gamma)?getTmat(beta,delta):0) + ((beta==delta)?getTmat(alpha,gamma):0));

}


void Hamiltonian::set_cum(){
	 int *NORBcumulative  = new int[getNirreps()+1]; 
	 NORBcumulative[0]  = 0; 
	 for (int irrep=0; irrep< getNirreps(); irrep++){ 
	    if (getNORB(irrep) < 0){ cerr << "DMRGSCFindices::DMRGSCFindices : NVIRT[" << irrep << "] = " << getNORB(irrep) << std::endl; } 
	     
	    NORBcumulative[ irrep+1] = NORBcumulative[ irrep] + getNORB(irrep); 

	}
	 _norbcum = NORBcumulative;
}

int Hamiltonian::getnup() const
{
    return _nup;
}

int Hamiltonian::getndown() const
{
    return _ndown;
}

double Hamiltonian::get_difference(Hamiltonian * ham_in){
	double dif = 0;
	double count = 1;	
        dif += abs(this->getEconst() - ham_in->getEconst());
        for (int cnt=0; cnt< L; cnt++){
           int I1 = this->getOrbitalIrrep(cnt);
           for (int cnt2=cnt; cnt2< L; cnt2++){
              int I2 = this->getOrbitalIrrep(cnt2);
              if (I1==I2){
                 dif += abs(this->getTmat(cnt,cnt2) - ham_in->getTmat(cnt,cnt2));
		 count ++;
              }
              for (int cnt3=cnt; cnt3< L; cnt3++){
                 int I3 = this->getOrbitalIrrep(cnt3);
                 for (int cnt4=cnt2; cnt4< L; cnt4++){
                    int I4 = this->getOrbitalIrrep(cnt4);
                    if (SymmInfo.directProd(I1,I2) == SymmInfo.directProd(I3,I4)){
                       dif += abs( this->getVmat(cnt,cnt2,cnt3,cnt4) - ham_in->getVmat(cnt,cnt2,cnt3,cnt4));
		       count ++;
                    }
                 }
              }
           }
	}
	return dif /= count;
}	

void Hamiltonian::set_zero(){
        for (int cnt=0; cnt< L; cnt++){
           int I1 = this->getOrbitalIrrep(cnt);
           for (int cnt2=cnt; cnt2< L; cnt2++){
              int I2 = this->getOrbitalIrrep(cnt2);
              if (I1==I2){
                 this->setTmat(cnt,cnt2,0); 
              }
              for (int cnt3=cnt; cnt3< L; cnt3++){
                 int I3 = this->getOrbitalIrrep(cnt3);
                 for (int cnt4=cnt2; cnt4< L; cnt4++){
                    int I4 = this->getOrbitalIrrep(cnt4);
                    if (SymmInfo.directProd(I1,I2) == SymmInfo.directProd(I3,I4)){
                       this->setVmat(cnt,cnt2,cnt3,cnt4, 0); 
                    }
                 }
              }
           }
	}
}	

OptIndex Hamiltonian::get_index_object(){
	return OptIndex(L, SymmInfo.getGroupNumber() ,irrep2num_orb);
}

/*  
std::vector<int> Hamiltonian::get_perm_order(){
//This function is depcrecated because it happens automatically in the psi4 outputroutine
	std::vector<int> index(moealpha.size(), 0);
	for (int i = 0 ; i != index.size() ; i++) {
		    index[i] = i;
	}
	sort(index.begin(), index.end(),
			    [&](const int& a, const int& b) { return (moealpha[a] <moealpha[b]); } ); //make use of lambda functions std=c11 inherited from boost
	return index;
}
*/

TwoIndex Hamiltonian::get_tmat_pointer()const{return *Tmat;}
FourIndex Hamiltonian::get_vmat_pointer()const{return *Vmat;}

//-------------------------class OptIndex----------------------------------------------------

OptIndex::OptIndex(const int L, const int _Group, const int * NORBin)
{
    this->L = L;
    Irreps SymmInfo;
    SymmInfo.setGroup(_Group);
    this->Nirreps = SymmInfo.getNumberOfIrreps();

    NORB.resize(Nirreps);
    NORBcumulative.resize(Nirreps+1);

    int sum_check = 0;
    NORBcumulative[0] = 0;
    for (int irrep=0; irrep<Nirreps; irrep++)
    {
        NORB[irrep] = NORBin[irrep];

        sum_check += NORB[irrep];

        NORBcumulative[irrep+1] = NORBcumulative[irrep] + NORB[irrep];
    }

    if (sum_check != L)
        cerr << "OptIndex::OptIndex : Sum over all orbitals is not L." << endl;

    irrep_each_orbital.resize(NORBcumulative[Nirreps]);

    for (int irrep=0; irrep<Nirreps; irrep++)
        for (int cnt=0; cnt<NORB[irrep]; cnt++)
            irrep_each_orbital[ NORBcumulative[irrep] + cnt ] = irrep;
}

OptIndex::OptIndex(const OptIndex & index)
{
    this->L = index.getL();
    this->Nirreps = index.getNirreps();

    NORB.resize(Nirreps);
    NORBcumulative.resize(Nirreps+1);

    int sum_check = 0;
    NORBcumulative[0] = 0;
    for (int irrep=0; irrep<Nirreps; irrep++)
    {
        NORB[irrep] = index.getNORB(irrep);

        sum_check += NORB[irrep];

        NORBcumulative[irrep+1] = NORBcumulative[irrep] + NORB[irrep];
    }

    if (sum_check != L)
        cerr << "OptIndex::OptIndex : Sum over all orbitals is not L." << endl;

    irrep_each_orbital.resize(NORBcumulative[Nirreps]);

    for (int irrep=0; irrep<Nirreps; irrep++)
        for (int cnt=0; cnt<NORB[irrep]; cnt++)
            irrep_each_orbital[ NORBcumulative[irrep] + cnt ] = irrep;

}

int OptIndex::getL() const { return L; }

int OptIndex::getNirreps() const { return Nirreps; }

int OptIndex::getNORB(const int irrep) const { return NORB[irrep]; }

int OptIndex::getNstart(const int irrep) const { return NORBcumulative[irrep]; }

int * OptIndex::get_irrep_each_orbital() { return irrep_each_orbital.data(); }

void OptIndex::Print() const
{
    cout << "NORB  = [ ";
    for (int irrep=0; irrep<Nirreps-1; irrep++)
        cout << NORB[irrep] << " , ";
    cout << NORB[Nirreps-1] << " ]" << endl;
}
