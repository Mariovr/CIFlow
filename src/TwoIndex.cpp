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
#include <sstream>
#include <string>
#include <algorithm>
#include <assert.h>

#include "TwoIndex.h"
#include "MyHDF5.h"

using namespace std;

TwoIndex::TwoIndex(const int nGroup, const int * IrrepSizes){

   SymmInfo.setGroup(nGroup);
   
   Isizes = new int[SymmInfo.getNumberOfIrreps()];
   storage = new double*[SymmInfo.getNumberOfIrreps()];
   
   for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
      Isizes[cnt] = IrrepSizes[cnt];
      if (Isizes[cnt]>0) storage[cnt] = new double[Isizes[cnt]*(Isizes[cnt]+1)/2]();
   }

}

TwoIndex::TwoIndex(const TwoIndex & tindex){
   SymmInfo.setGroup(tindex.get_group());
   
   Isizes = new int[SymmInfo.getNumberOfIrreps()];
   storage = new double*[SymmInfo.getNumberOfIrreps()];
   
   copy(tindex.get_irrepsizes(), tindex.get_irrepsizes() + SymmInfo.getNumberOfIrreps() , Isizes ); //from algorithm
   for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
      int irrepsize = Isizes[cnt]*(Isizes[cnt]+1)/2;
      if (Isizes[cnt]>0){
	      storage[cnt] = new double[irrepsize];
	      copy(tindex.get_storage_irrep(cnt) , tindex.get_storage_irrep(cnt) + irrepsize, storage[cnt] ); 
      }
   }

}

TwoIndex::~TwoIndex(){
   
   for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++) if (Isizes[cnt]>0) delete [] storage[cnt];
   delete [] storage;
   delete [] Isizes;
   
}

void TwoIndex::set(const int irrep, const int i, const int j, const double val){

   if (i>j) storage[irrep][j + i*(i+1)/2] = val;
   else storage[irrep][i + j*(j+1)/2] = val;

}

double TwoIndex::get(const int irrep, const int i, const int j) const{

   if (i>j) return storage[irrep][j + i*(i+1)/2];
   return storage[irrep][i + j*(j+1)/2];

}

void TwoIndex::save(const std::string name) const{
    //The hdf5 file
    hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    hid_t group_id2 = H5Gcreate(file_id, "TwoIndex", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //The metadata
    hid_t group_id = H5Gcreate(group_id2, "MetaData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //The IrrepSizes
    hsize_t dimarray       = SymmInfo.getNumberOfIrreps();
    hid_t dataspace_id     = H5Screate_simple(1, &dimarray, NULL);
    hid_t dataset_id       = H5Dcreate(group_id, "IrrepSizes", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Isizes);

    //Attributes
    hid_t attribute_space_id1  = H5Screate(H5S_SCALAR);
    hid_t attribute_id1        = H5Acreate(dataset_id, "nGroup", H5T_STD_I32LE, attribute_space_id1, H5P_DEFAULT, H5P_DEFAULT);
    int nGroup                 = SymmInfo.getGroupNumber();
    H5Awrite(attribute_id1, H5T_NATIVE_INT, &nGroup); 

    hid_t attribute_space_id2  = H5Screate(H5S_SCALAR);
    hid_t attribute_id2        = H5Acreate(dataset_id, "nIrreps", H5T_STD_I32LE, attribute_space_id2, H5P_DEFAULT, H5P_DEFAULT);
    int nIrreps                = SymmInfo.getNumberOfIrreps();
    H5Awrite(attribute_id2, H5T_NATIVE_INT, &nIrreps); 

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);
    H5Sclose(attribute_space_id1);
    H5Sclose(attribute_space_id2);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    H5Gclose(group_id);

    //The object itself.
    for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){

        if (Isizes[cnt]>0) {

            std::stringstream sstream;
            sstream << "TwoIndex" << cnt ;
            hid_t group_id3 = H5Gcreate(group_id2, sstream.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hsize_t dimarray3       = Isizes[cnt]*(Isizes[cnt]+1)/2;
            hid_t dataspace_id3     = H5Screate_simple(1, &dimarray3, NULL);
            hid_t dataset_id3       = H5Dcreate(group_id3, "Matrix elements", H5T_IEEE_F64LE, dataspace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[cnt]);

            H5Dclose(dataset_id3);
            H5Sclose(dataspace_id3);

            H5Gclose(group_id3);

        }

    }

    H5Gclose(group_id2);

    H5Fclose(file_id);
 
}

void TwoIndex::read(const std::string name){
    //The hdf5 file
    hid_t file_id = H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    //The metadata
    hid_t group_id = H5Gopen(file_id, "/TwoIndex/MetaData",H5P_DEFAULT);

    //The IrrepSizes
    hid_t dataset_id = H5Dopen(group_id, "IrrepSizes", H5P_DEFAULT);

    //Attributes
    hid_t attribute_id1 = H5Aopen_by_name(group_id,"IrrepSizes", "nGroup", H5P_DEFAULT, H5P_DEFAULT);
    int nGroup;
    H5Aread(attribute_id1, H5T_NATIVE_INT, &nGroup);
    assert( nGroup==SymmInfo.getGroupNumber() );

    hid_t attribute_id2 = H5Aopen_by_name(group_id,"IrrepSizes", "nIrreps", H5P_DEFAULT, H5P_DEFAULT);
    int nIrreps;
    H5Aread(attribute_id2, H5T_NATIVE_INT, &nIrreps);
    assert( nIrreps==SymmInfo.getNumberOfIrreps() );

    H5Aclose(attribute_id1);
    H5Aclose(attribute_id2);

    int * IsizesAgain = new int[SymmInfo.getNumberOfIrreps()];
    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, IsizesAgain);
    for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){
        assert( IsizesAgain[cnt]==Isizes[cnt] );
    }
    delete [] IsizesAgain;
    H5Dclose(dataset_id);

    H5Gclose(group_id);

    //The object itself.
    for (int cnt=0; cnt<SymmInfo.getNumberOfIrreps(); cnt++){

        if (Isizes[cnt]>0) {

            std::stringstream sstream;
            sstream << "/TwoIndex/TwoIndex" << cnt ;
            hid_t group_id3 = H5Gopen(file_id, sstream.str().c_str(), H5P_DEFAULT);

            hid_t dataset_id3 = H5Dopen(group_id3, "Matrix elements", H5P_DEFAULT);
            H5Dread(dataset_id3, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, storage[cnt]);
            H5Dclose(dataset_id3);

            H5Gclose(group_id3);

        }

    }

    H5Fclose(file_id);
 
}



