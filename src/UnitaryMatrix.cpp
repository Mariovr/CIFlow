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
#include <random>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <fstream>

#include "MyHDF5.h"
//######include "Lapack.h"
#include "UnitaryMatrix.h"
#include "matrix.h"
#include "Options.h"
#include "Hamiltonian.h" //For optindex (should put this in a separate file)

using std::string;
using std::ifstream;
using std::cerr;
using std::cout;
using std::endl;

UnitaryMatrix::UnitaryMatrix(OptIndex& index)
{
    _index.reset(new OptIndex(index));

    //Allocate the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = 0;
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        x_linearlength += size;
        unitary[irrep].reset(new double[size]);
        memset(unitary[irrep].get(), 0, sizeof(double)*size);
        for (int cnt=0; cnt<linsize; cnt++)
            unitary[irrep][cnt*(1+linsize)] = 1.0;
    }
}

UnitaryMatrix::UnitaryMatrix(const UnitaryMatrix & unit)
{
    _index.reset(new OptIndex(*unit._index));
    //Copy the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = unit.x_linearlength;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        unitary[irrep].reset(new double[size]);
        for (int index=0; index<size; index++)
            unitary[irrep][index] = unit.getBlock(irrep)[index];
    }
}

UnitaryMatrix::UnitaryMatrix(OptIndex & index, std::istream & file):UnitaryMatrix(index)
{
  load_unitary(file);
}

UnitaryMatrix::UnitaryMatrix(OptIndex & index, const std::string & file):UnitaryMatrix(index)
{
  load_unitary(file);
}

UnitaryMatrix UnitaryMatrix::get_inverse() const{
    UnitaryMatrix unit = UnitaryMatrix(*_index);
    for (int irrep=0; irrep< _index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        for( int j = 0 ; j < linsize ; j ++)
        {
            for( int l = 0 ; l < linsize ; l ++)
                unit.getBlock(irrep)[j + linsize* l] =  unitary[irrep][l+linsize * j];
        }
    }
    return unit;
}

UnitaryMatrix::UnitaryMatrix(UnitaryMatrix && unit)
{
    _index = std::move(unit._index);
    unitary = std::move(unit.unitary);
}

UnitaryMatrix& UnitaryMatrix::operator=(const UnitaryMatrix &unit)
{
    _index.reset(new OptIndex(*unit._index));
    //Copy the unitary
    unitary.resize(_index->getNirreps());
    x_linearlength = unit.x_linearlength;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        unitary[irrep].reset(new double[size]);
        for (int index=0; index<size; index++)
            unitary[irrep][index] = unit.getBlock(irrep)[index];
    }

    return *this;
}

UnitaryMatrix& UnitaryMatrix::operator=(UnitaryMatrix &&unit)
{
    _index = std::move(unit._index);
    unitary = std::move(unit.unitary);

    return *this;
}

double UnitaryMatrix::get_element(int i , int j)
{
    int irrepi = get_orbital_irrep(i);
    int irrepj = get_orbital_irrep(j);
    if (irrepi != irrepj)
        return 0.;
    else
    {
        int start = get_nstart(irrepi);
        int norb = get_norb(irrepi);
        return unitary[irrepi][i-start + (j - start) * norb ];
    }

}


std::vector<double> UnitaryMatrix::get_full_transformation()
{
    int L = _index->getL();
    std::vector<double> transform(L*L , 0);
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        int nstart = get_nstart(irrep);
        for( int j = 0 ; j < linsize ; j ++)//Row
        {
            for( int l = 0 ; l < linsize ; l ++)//col
                transform[(j+nstart)+L*(l+nstart)]  =  unitary[irrep][j+linsize * l];
        }
    }
    return transform;
}


void UnitaryMatrix::jacobi_rotation(int irrep , int i , int j , double angle)
{
    //Simple implementation of jacobi rotation. 
    //REMARK new eigenvectors are saved in the rows.
    i -= _index->getNstart(irrep);
    j -= _index->getNstart(irrep);
    int linsize = _index->getNORB(irrep);
    double  work[2*linsize];

    //cout << i << j<< linsize << cos(angle) << sin(angle);
    for(int l = 0 ; l < linsize ; l++)
    {
        work[l] = unitary[irrep][i+l*linsize]*cos(angle) + unitary[irrep][j+l*linsize] * -1.*sin(angle);
        work[l+linsize] = unitary[irrep][j+l*linsize]*cos(angle) + unitary[irrep][i+l*linsize] * sin(angle);
    }

    for(int l = 0 ; l < _index->getNORB(irrep) ; l++)
    {
        unitary[irrep][i+l*linsize] = work[l];
        unitary[irrep][j+l*linsize] = work[l+linsize];
    }
}

void UnitaryMatrix::reset_unitary(bool allzero)
{
   for (int irrep=0; irrep<_index->getNirreps(); irrep++)
   {
      const int linsize = _index->getNORB(irrep);
      const int size = linsize * linsize;
      unitary[irrep].reset(new double[size]());
      //memmove(unitary[irrep].get(), 0, sizeof(double)*size);
      if(not allzero)
      {
          for (int cnt=0; cnt<linsize; cnt++)
              unitary[irrep][cnt*(1+linsize)] = 1.0;
      }
   }
}

unsigned int UnitaryMatrix::getNumVariablesX() const{ return x_linearlength; }

double * UnitaryMatrix::getBlock(const int irrep) const { return unitary[irrep].get(); }

void UnitaryMatrix::updateUnitary(double * temp1, double * temp2, double * vector, const bool multiply)
{
    //Per irrep
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        int size = linsize * linsize;

        //linsize is op z'n minst 2 dus temp1, temp1+size, temp1+2*size,temp1+3*size zijn zeker ok
        if (linsize>1)
        {         
            //Construct the anti-symmetric x-matrix
            double * xblock    = temp1;
            double * Bmat      = temp1 +   size;   //linsize*linsize
            double * work1     = temp1 + 2*size;   //linsize*linsize
            double * work2     = temp1 + 3*size;   //linsize*linsize

            double * workLARGE = temp2;  //4*size
            int     lworkLARGE = 4*size; //4*size = 4*linsize*linsize > 3*linsize-1

            //Construct the antisymmetric x-matrix
            build_skew_symm_x(irrep, xblock , vector);

            //Bmat <= xblock * xblock
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0; //SET !!!
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,xblock,&linsize,&beta,Bmat,&linsize); 

            //Calculate its eigenvalues and eigenvectors
            //Bmat * work1 * Bmat^T <= xblock * xblock
            char uplo = 'U';
            char jobz = 'V';
            int info;
            dsyev_(&jobz, &uplo, &linsize, Bmat, &linsize, work1, workLARGE, &lworkLARGE, &info); 

            if (info != 0 )
            {
                cerr << "A problem occured in the diagonalisation of the anti-symmetric matrix squared to generated the unitary -> exp(X).";
                //exit(EXIT_FAILURE);	 
            }

            //work2 <= Bmat^T * xblock * Bmat
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,xblock,&linsize,Bmat,&linsize,&beta,work1,&linsize);
            char trans = 'T';
            dgemm_(&trans,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work1,&linsize,&beta,work2,&linsize); //work2 = Bmat^T * xblock * Bmat

            if (Orbopt_debugPrint)
            {
                cout << "   UnitaryMatrix::updateUnitary : Lambdas of irrep block " << irrep << " : " << endl;
                for (int cnt=0; cnt<linsize/2; cnt++)
                {
                    cout << "      block = [ " << work2[2*cnt   + linsize*2*cnt] << " , " << work2[2*cnt   + linsize*(2*cnt+1)] << " ] " << endl;
                    cout << "              [ " << work2[2*cnt+1 + linsize*2*cnt] << " , " << work2[2*cnt+1 + linsize*(2*cnt+1)] << " ] " << endl;
                }
            }

            //work1 <= values of the antisymmetric 2x2 blocks
            for (int cnt=0; cnt<linsize/2; cnt++)
                work1[cnt] = 0.5*( work2[2*cnt + linsize*(2*cnt+1)] - work2[2*cnt+1 + linsize*(2*cnt)] );

            if (Orbopt_debugPrint)
            {
                for (int cnt=0; cnt<linsize/2; cnt++)
                {
                    work2[2*cnt + linsize*(2*cnt+1)] -= work1[cnt];
                    work2[2*cnt+1 + linsize*(2*cnt)] += work1[cnt];
                }

                int inc = 1;
                double RMSdeviation = ddot_(&size, work2, &inc, work2, &inc);
                RMSdeviation = sqrt(RMSdeviation);

                cout << "   UnitaryMatrix::updateUnitary : RMSdeviation of irrep block " << irrep << " (should be 0.0) = " << RMSdeviation << endl;
            }

            //Calculate exp(x)
            //work2 <= exp(Bmat^T * xblock * Bmat)
            memset(work2, 0, sizeof(double)*size);
            for (int cnt=0; cnt<linsize/2; cnt++)
            {
                double cosine = cos(work1[cnt]);
                double sine = sin(work1[cnt]);
                work2[2*cnt   + linsize*(2*cnt  )] = cosine;
                work2[2*cnt+1 + linsize*(2*cnt+1)] = cosine;
                work2[2*cnt   + linsize*(2*cnt+1)] = sine;
                work2[2*cnt+1 + linsize*(2*cnt  )] = - sine;
            }

            //Handles the case when there is an odd number of orbitals.
            for (int cnt=2*(linsize/2); cnt<linsize; cnt++)
                work2[cnt*(linsize + 1)] = 1.0;

            //work2 <= Bmat * exp(Bmat^T * xblock * Bmat) * Bmat^T = exp(xblock)
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,Bmat,&linsize,work2,&linsize,&beta,work1,&linsize);
            dgemm_(&notr,&trans,&linsize,&linsize,&linsize,&alpha,work1,&linsize,Bmat,&linsize,&beta,work2,&linsize); //work2 = exp(xblock)

            //U <-- exp(x) * U
            int inc = 1;
            if (multiply)
            {
                //U <-- exp(x) * U
                dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,work2,&linsize,unitary[irrep].get(),&linsize,&beta,work1,&linsize);
                dcopy_(&size, work1, &inc, unitary[irrep].get(), &inc);
            }
            else  //U <-- exp(x)
                dcopy_(&size, work2, &inc, unitary[irrep].get(), &inc);
        }
    }

    if (Orbopt_debugPrint){ CheckDeviationFromUnitary(); }
}

void UnitaryMatrix::build_skew_symm_x(const int irrep, double * result , const double * Xelem) const
{
    //should have size: linsize*linsize.	
    //Remark we build the antisymmetrix X matrix within an irrep.
    const int linsize = _index->getNORB(irrep);

    int jump = 0;
    for (int cnt=0; cnt<irrep; cnt++)
    {
        int linsizeCNT = _index->getNORB(cnt);
        jump += linsizeCNT * (linsizeCNT-1) / 2;
    }

    for (int row=0; row<linsize; row++)
    {
        result[ row + linsize * row ] = 0;

        for (int col=row+1; col<linsize; col++)
        {
            result[ row + linsize * col ] =   Xelem[ jump + row + col*(col-1)/2 ];
            result[ col + linsize * row ] = - Xelem[ jump + row + col*(col-1)/2 ];
        }
    }
}

void UnitaryMatrix::rotate_active_space_vectors(double * eigenvecs, double * work)
{
    int passed = 0;
    int norb = _index->getL();
    for (int irrep=0; irrep<_index->getNirreps(); irrep++){

        const int nirrep = _index->getNORB(irrep);
        if (nirrep> 1){

            int rotationlinsize =nirrep;

            double * temp1 = work;
            double * temp2 = work + rotationlinsize*rotationlinsize;
            double * BlockEigen = eigenvecs + passed * (norb+1); //after rotating the irrep before point to next eigenvector irrep block, take into account we have to pass a stride of norb although (distance between next columns in column major) the dimension of matrix is nirrep.

            for (int row = 0; row<rotationlinsize; row++){
                for (int col = 0; col<rotationlinsize; col++){
                    temp1[row + rotationlinsize*col] = unitary[irrep][ row + rotationlinsize * col ];
                }
            }

            char tran = 'T';
            char notr = 'N';
            double alpha = 1.0;
            double beta = 0.0;
            //eigvec_jk C_ji ph_i (new eigenbasis is saved in columns in column major format, but eigenbasis in unitary matrix is in rows in column major format.)
            //so we have to transform BlockEigen so it saves the new eigenvectors also in the rows and then multiply left with the unitary matrix.
            dgemm_(&tran,&notr,&rotationlinsize,&rotationlinsize,&rotationlinsize,&alpha,BlockEigen,&norb,temp1,&rotationlinsize,&beta,temp2,&rotationlinsize);

            for (int row = 0; row<rotationlinsize; row++){
                for (int col = 0; col<rotationlinsize; col++){
                    unitary[irrep][  row + rotationlinsize * col ] = temp2[row + rotationlinsize*col];
                }
            }

        }

        passed +=nirrep;

    }

    if (Orbopt_debugPrint){ CheckDeviationFromUnitary(); }
}

void UnitaryMatrix::CheckDeviationFromUnitary() const
{
    char tran = 'T';
    char notr = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        double work[linsize * linsize];
        // linsize > 0 (if the system has orbitals within the irrep)
        if(linsize > 0)
        {
            dgemm_(&tran,&notr,&linsize,&linsize,&linsize,&alpha,unitary[irrep].get(),&linsize,unitary[irrep].get(),&linsize,&beta,work,&linsize);

            double value = 0.0;
            for (int cnt=0; cnt<linsize; cnt++){
                value += (work[cnt*(1+linsize)]-1.0) * (work[cnt*(1+linsize)]-1.0); //diagonal elements should be one.
                for (int cnt2=cnt+1; cnt2<linsize; cnt2++){
                    value += work[cnt + cnt2*linsize] * work[cnt + cnt2*linsize] + work[cnt2 + cnt*linsize] * work[cnt2 + cnt*linsize];//off diagonal elements are zero, we use the triangle indices for the upper and down triangular matrix elements.
                }
            }
            value = sqrt(value);
            cout << "Two-norm of unitary[" << irrep << "]^(dagger) * unitary[" << irrep << "] = " << value << endl;
            if((value) > 1e-10)
            {
                cerr << "WARNING: we reseted the unitary because we lost unitarity." << endl;
                throw std::logic_error("The UnitaryMatrix is not unitary.");
            }
        }
    }
}

void UnitaryMatrix::multiply_left_with(UnitaryMatrix & unit2) const
{
    //puts in this unitary: unit2*unitary, corresponds to a basis transformation to new orbitals if the unitary matrix basisses are saved in the rows.
    char tran = 'T';
    char notr = 'N';
    double alpha = 1.0;
    double beta = 0.0;

    int inc =1;
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        int size = linsize*linsize;
        double * work= new double[size]();
        if(linsize > 1)
        {
            dgemm_(&notr,&notr,&linsize,&linsize,&linsize,&alpha,unit2.unitary[irrep].get(),&linsize,unitary[irrep].get(),&linsize,&beta,work,&linsize);
            dcopy_(&size, work, &inc, unitary[irrep].get(), &inc);

        }
        delete [] work;
    }
}

void UnitaryMatrix::saveU(const string savename) const
{
    hid_t file_id = H5Fcreate(savename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); //New format for unitarymatrices, compatible with Hamiltonian hdf5 format.

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int norb = _index->getNORB(irrep); 

        if(norb > 0)
        {
            std::stringstream irrepname;
            irrepname << "irrep_" << irrep;

            hsize_t dimarray      = norb * norb;
            hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
            hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep].get());

            H5Dclose(dataset_id);
            H5Sclose(dataspace_id);
        }
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
}

void UnitaryMatrix::loadU(const string unitname)
{
    hid_t file_id = H5Fopen(unitname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t group_id = H5Gopen(file_id, "/Data", H5P_DEFAULT);

    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        if(_index->getNORB(irrep) != 0)
        {
            std::stringstream irrepname;
            irrepname << "irrep_" << irrep;

            hid_t dataset_id = H5Dopen(group_id, irrepname.str().c_str(), H5P_DEFAULT);
            H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, unitary[irrep].get());

            H5Dclose(dataset_id);
        }
    }

    H5Gclose(group_id);
    H5Fclose(file_id);
}

void UnitaryMatrix::deleteStoredUnitary(std::string name) const
{
    std::stringstream temp;
   temp << "rm " << name;
   int info = system(temp.str().c_str());
   cout << "Info on CASSCF::Unitary rm call to system: " << info << endl;
}

void UnitaryMatrix::print_unitary(std::ostream & out ) const
{
    out << std::setprecision(12);
    out << "CIFlowTransformation" << std::endl; 
    for( int irrep = 0 ; irrep < _index->getNirreps() ; irrep ++)
    {
        int linsize = _index->getNORB(irrep);
        if(linsize > 0 ){
            out << "#irrep : " << irrep << "(new basis is in rows)" << std::endl;
            for( int j = 0 ; j < linsize ; j ++)
            {
                for( int l = 0 ; l < linsize ; l ++)
                    out << unitary[irrep][j + linsize* l] << " ";

                out << std::endl;
            }
            out << std::endl;
        }
    }
}

void UnitaryMatrix::print_unitary(const std::string & filename) const
{
    std::ofstream file(filename.c_str()); 
    print_unitary(file);
    file.close();
}

void UnitaryMatrix::load_unitary(const std::string filename)
{
    int last_index = filename.find_last_of('.');
    //Check if inputfile is in hdf5 format.
    if("h5" == filename.substr(last_index+1 , string::npos))
    {
        loadU(filename);
    }
    else
    {
        std::ifstream file(filename.c_str()); 
        load_unitary(file);
        file.close();
    }

}

void UnitaryMatrix::load_unitary(std::istream & file)
{
    std::string line;
    //Wind the file forward to the relevant part.
    int position = file.tellg();
    while(line.find("CIFlowTransformation") == std::string::npos )
    {
        std::getline(file, line);
    }
    for( int irrep = 0 ; irrep < _index->getNirreps() ; irrep ++)
    {
        int linsize = _index->getNORB(irrep);
        if(linsize > 0)
        {
            while(line.find("irrep") == std::string::npos && line.find("Irrep") == std::string::npos)
            {
                std::getline(file, line);
            }
            for( int j = 0 ; j < linsize ; j ++)
            {
                std::getline(file, line );
                std::istringstream sin(line);
                double coef;
                for( int l = 0 ; l < linsize ; l ++){
                    sin >> coef;
                    unitary[irrep][j + linsize* l] =  coef;
                }

            }
        }
    }
    file.seekg(position);
}

UnitaryMatrix::UnitaryMatrix(matrix & mat, OptIndex &opt){
    //const int norbin[1] = {mat.getn()};
    //_index.reset(new OptIndex(mat.getn() , 0, norbin)); //sets explicitly c1 symmetry.
    _index.reset(new OptIndex(opt)); 
    //Allocate the unitary
    unitary.resize(_index->getNirreps());
    for( int irrep = 0 ; irrep < _index->getNirreps() ; irrep ++)
    {
        int linsize = _index->getNORB(irrep);
        if (linsize > 0)
        {
            int nstart = _index->getNstart(irrep);
            unitary[irrep].reset(new double[linsize*linsize]);
            for( int j = 0 ; j < linsize ; j ++)
            {
                for( int l = 0 ; l < linsize ; l ++)
                    unitary[irrep][j + linsize* l] = mat(l+nstart,j+nstart);//mat contains new basis in columns but unitary contains them in the rows. (Both are in column major representation.)

            }
        }//end linsize > 1
    }
}    

UnitaryMatrix::UnitaryMatrix(const std::vector<double> & mat , OptIndex &opt){
    //const int norbin[1] = {mat.getn()};
    //_index.reset(new OptIndex(mat.getn() , 0, norbin)); //sets explicitly c1 symmetry.
    _index.reset(new OptIndex(opt)); 
    //Allocate the unitary
    unitary.resize(_index->getNirreps());
    for( int irrep = 0 ; irrep < _index->getNirreps() ; irrep ++)
    {
        int linsize = _index->getNORB(irrep);
        if (linsize > 0)
        {
            int nstart = _index->getNstart(irrep);
            unitary[irrep].reset(new double[linsize*linsize]);
            for( int j = 0 ; j < linsize ; j ++)
            {
                for( int l = 0 ; l < linsize ; l ++)
                    unitary[irrep][j + linsize* l] = mat[l+j*linsize];//mat contains new basis in columns but unitary contains them in the rows. (Both are in column major representation.)

            }
        }//end linsize > 1
    }
}    

int UnitaryMatrix::get_norb(int irrep){ return _index->getNORB(irrep); }

int UnitaryMatrix::get_nstart(int irrep){ return _index->getNstart(irrep); }

int UnitaryMatrix::get_orbital_irrep(int orbital){ return _index->get_irrep_each_orbital()[orbital]; }

int UnitaryMatrix::get_Nirrep() const { return _index->getNirreps(); }

/**
* Spread the UnitaryMatrix from rank orig to all
* other UnitaryMatrices in the world
*/
 /*
void UnitaryMatrix::sendreceive(int orig)
{
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        MPI_Bcast(unitary[irrep].get(), size, MPI_DOUBLE, orig, MPI_COMM_WORLD);
    }
}

*/  
void UnitaryMatrix::set_random_unitary(){
    std::random_device rd;
    auto mt = std::mt19937_64(rd());
    std::uniform_real_distribution<double> dist_angles(-1000, 1000);
    //std::uniform_real_distribution<double> dist_angles(-M_PI, M_PI);
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        if(linsize > 1)
        {
            int size = linsize * linsize;
            matrix eigvec(linsize , linsize);
            matrix eigval(linsize , 1);
            matrix sym(linsize,linsize);
            //fill symmetric matrix. Remember the eigenvectors of a real symmetric matrix are orthogonal.
            for(int i=0;i<linsize;i++)
                for(int j=i;j<linsize;j++)
                {
                    sym(i,j) = dist_angles(mt);
                    sym(j,i) = sym(i,j);
                }
            sym.diagonalize(eigval, eigvec);
            int inc = 1;
            dcopy_(&size, eigvec.getpointer(), &inc, unitary[irrep].get(), &inc);
        }
        else if(linsize ==1 )
        {
            unitary[irrep][0] = 1.;
        }
    }
}

void UnitaryMatrix::make_skew_symmetric()
{
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        for(int i=0;i<linsize;i++)
        {
            unitary[irrep][i+i*linsize] = 0;
            for(int j=i+1;j<linsize;j++)
                unitary[irrep][j+i*linsize] = -1 * unitary[irrep][i+j*linsize];
        }
    }
}

double UnitaryMatrix::check_difference(UnitaryMatrix & unit2)
{
    double som = 0 ;
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        for (int cnt=0; cnt< size; cnt++)
            som += pow( (fabs(unitary[irrep][cnt]) - fabs(unit2.getBlock(irrep)[cnt]) ) , 2);
    }
    som /=  x_linearlength;

    return sqrt(som);
}

void UnitaryMatrix::add_unit(UnitaryMatrix & unit2, double mult)
{
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        const int size = linsize * linsize;
        for (int cnt=0; cnt< size; cnt++)
            unitary[irrep][cnt] += mult * unit2.getBlock(irrep)[cnt];
    }
}

std::vector<matrix> UnitaryMatrix::calc_overlap()
{
    char tran = 'T';
    char notr = 'N';
    double alpha = 1.0;
    double beta = 0.0;
    std::vector<matrix> overlap;
    overlap.resize(0);
    for (int irrep=0; irrep<_index->getNirreps(); irrep++)
    {
        int linsize = _index->getNORB(irrep);
        // linsize > 0 (if the system has orbitals within the irrep)
        if(linsize > 0)
        {
            matrix work(linsize, linsize) ;
            dgemm_(&notr,&tran,&linsize,&linsize,&linsize,&alpha,unitary[irrep].get(),&linsize,unitary[irrep].get(),&linsize,&beta,work.getpointer(),&linsize);
            overlap.push_back(work);
        }
    }
    return overlap;
}

/* vim: set ts=4 sw=4 expandtab :*/
