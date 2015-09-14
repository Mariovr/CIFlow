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
#ifndef HDF5_DOCI_H
#define HDF5_DOCI_H

   //Force the use of the 1.8 API of HDF5
   #undef H5_USE_16_API
   #define H5_NO_DEPRECATED_SYMBOLS
   #define H5Acreate_vers 2
   #define H5Dcreate_vers 2
   #define H5Dopen_vers 2
   #define H5Gcreate_vers 2
   #define H5Gopen_vers 2

   #include <hdf5.h>

   // macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) \
    do { \
        if(status < 0) \
            std::cerr << __FILE__ << ":" << __LINE__ << \
            ": Problem with writing to file. Status code=" \
            << status << std::endl; \
    }while(false)

#endif /* HDF5_DOCI_H */
