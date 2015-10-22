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
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cstring>

#include "OrbitalTransform.h"
#include "Hamiltonian.h" 
#include "UnitaryMatrix.h"
#include "Lapack.h"

using std::min;
using std::max;

OrbitalTransform::OrbitalTransform(Hamiltonian *  HamIn)
{
    int L = HamIn->getL();
    SymmInfo.setGroup(HamIn->getNGroup());
    _hamorig = HamIn;

    _index.reset(new OptIndex(_hamorig->get_index_object()));
    numberOfIrreps = _index->getNirreps();
    _unitary.reset(new UnitaryMatrix(*_index));

    //Create the memory for the orbital transformations.
    unsigned long long maxlinsize = 0;
    for (int irrep=0; irrep< _index->getNirreps(); irrep++)
    {
        unsigned int linsize_irrep = _index->getNORB(irrep);
        if (linsize_irrep > maxlinsize)  
            maxlinsize  = linsize_irrep;
    }

    //Determine the blocksize for the 2-body transformation
    auto& maxBlockSize = maxlinsize;
    //Allocate 2-body rotation memory: One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
    auto maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!
    auto sizeWorkmem1 = max( max( maxBSpower4 , 3*maxlinsize*maxlinsize ) , 1uLL * L * L ); //For (2-body tfo , updateUnitary, calcNOON)
    auto sizeWorkmem2 = max( max( maxBSpower4 , 2*maxlinsize*maxlinsize ) , L*(L + 1uLL) ); //For (2-body tfo, updateUnitary and rotate_to_active_space, rotate2DMand1DM)
    mem1.reset(new double[sizeWorkmem1]);
    mem2.reset(new double[sizeWorkmem2]);
}

void OrbitalTransform::set_ham(Hamiltonian * ham){ _hamorig = ham ; }
UnitaryMatrix & OrbitalTransform::get_unitary(){ return (*_unitary); } //Watch out pass by reference.

void OrbitalTransform::fillHamCI(Hamiltonian& HamCI)
{
    assert(&HamCI != _hamorig);	
    buildOneBodyMatrixElements();	
    fillConstAndTmat(HamCI); //fill one body terms and constant part.	

    //Two-body terms --> use eightfold permutation symmetry in the irreps :-)
    for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++)
        for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++)
        {
            const int productSymm = SymmInfo.directProd(irrep1,irrep2);
            for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++)
            {
                const int irrep4 = SymmInfo.directProd(productSymm,irrep3);
                // Generated all possible combinations of allowed irreps
                if (irrep4>=irrep2)
                {
                    int linsize1 =_index->getNORB(irrep1);
                    int linsize2 =_index->getNORB(irrep2);
                    int linsize3 =_index->getNORB(irrep3);
                    int linsize4 =_index->getNORB(irrep4);

                    if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0))
                    {
                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ]
                                            = _hamorig->getVmat(_index->getNstart(irrep1) + cnt1,_index->getNstart(irrep2) + cnt2, _index->getNstart(irrep3) + cnt3, _index->getNstart(irrep4) + cnt4 );

                        char trans = 'T';
                        char notra = 'N';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET !!!

                        int rightdim = linsize2 * linsize3 * linsize4; //(ijkl) -> (ajkl)
                        double * Umx = _unitary->getBlock(irrep1);
                        dgemm_(&notra, &notra, &linsize1, &rightdim, &linsize1, &alpha, Umx, &linsize1, mem1.get(), &linsize1, &beta, mem2.get(), &linsize1);

                        int leftdim = linsize1 * linsize2 * linsize3; //(ajkl) -> (ajkd)
                        Umx = _unitary->getBlock(irrep4);
                        dgemm_(&notra, &trans, &leftdim, &linsize4, &linsize4, &alpha, mem2.get(), &leftdim, Umx, &linsize4, &beta, mem1.get(), &leftdim);

                        int jump1 = linsize1 * linsize2 * linsize3; //(ajkd) -> (ajcd)
                        int jump2 = linsize1 * linsize2 * linsize3;
                        leftdim   = linsize1 * linsize2;
                        Umx = _unitary->getBlock(irrep3);
                        for (int bla=0; bla<linsize4; bla++)
                            dgemm_(&notra, &trans, &leftdim, &linsize3, &linsize3, &alpha, mem1.get()+jump1*bla, &leftdim, Umx, &linsize3, &beta, mem2.get()+jump2*bla, &leftdim);

                        jump2    = linsize1 * linsize2;
                        jump1    = linsize1 * linsize2;
                        rightdim = linsize3 * linsize4;
                        Umx = _unitary->getBlock(irrep2);
                        for (int bla=0; bla<rightdim; bla++)
                            dgemm_(&notra, &trans, &linsize1, &linsize2, &linsize2, &alpha, mem2.get()+jump2*bla, &linsize1, Umx, &linsize2, &beta, mem1.get()+jump1*bla, &linsize1);

                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        HamCI.setVmat(_index->getNstart(irrep1) + cnt1,_index->getNstart(irrep2) + cnt2, _index->getNstart(irrep3) + cnt3,_index->getNstart(irrep4) + cnt4, mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ] );

                    } //end if the problem has orbitals from all 4 selected irreps
                } // end if irrep 4 >= irrep2
            }// end run irrep3
        } // end run irrep2
}

/**
 * This method fills the OneBodyMatrixElements array with
 * the 1-body matrix elements so we can rotate them
 */
void OrbitalTransform::buildOneBodyMatrixElements()
{
    if(!QmatrixWork.size() || !OneBodyMatrixElements.size())
    {
        QmatrixWork.resize(numberOfIrreps);
        OneBodyMatrixElements.resize(numberOfIrreps);

        for (int irrep=0; irrep<numberOfIrreps; irrep++)
        {
            const int size = _index->getNORB(irrep) * _index->getNORB(irrep);
            QmatrixWork[irrep].reset(new double[size]);
            OneBodyMatrixElements[irrep].reset(new double[size]);
        }	
    }

    //#pragma omp parallel for schedule(dynamic)
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        for (int row=0; row<linsize; row++)
        {
            const int HamIndex1 = _index->getNstart(irrep) + row;
            //set the diagonal elements.
            OneBodyMatrixElements[irrep][row*(1+linsize)] = _hamorig->getTmat(HamIndex1,HamIndex1);
            for (int col=row+1; col<linsize; col++)
            {
                const int HamIndex2 = _index->getNstart(irrep) + col;
                // remember blas and lapack are in fortran and there it is convenient to run first over columnsand then over rows.
                OneBodyMatrixElements[irrep][row + linsize * col] = _hamorig->getTmat(HamIndex1,HamIndex2); 
                OneBodyMatrixElements[irrep][col + linsize * row] = OneBodyMatrixElements[irrep][row + linsize * col];
            }
        }
    }
    rotate_old_to_new(OneBodyMatrixElements.data());
}

/**
 * Calculates the rotation of the 1-body matrix
 * T' = Q * T * Q^T
 * It uses the OneBodyMatrixElements array as start and endpoint
 * point and uses the Qmatrixwork as temp storage
 */
void OrbitalTransform::rotate_old_to_new(std::unique_ptr<double []> * matrix)
{
    //#pragma omp parallel for schedule(dynamic)
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        int linsize  = get_norb(irrep);
        if(linsize > 1)
        {
            double * Umx = _unitary->getBlock(irrep);
            double alpha = 1.0;
            double beta  = 0.0;
            char trans   = 'T';
            char notrans = 'N';
            // Qmatrixwork = Umx * OneBodyMatrixElements
            dgemm_(&notrans,&notrans,&linsize,&linsize,&linsize,&alpha,Umx,&linsize,matrix[irrep].get(),&linsize,&beta,QmatrixWork[irrep].get(),&linsize);
            // OneBodyMatrixElements = Qmatrixwork * Umx^T
            dgemm_(&notrans,&trans,  &linsize,&linsize,&linsize,&alpha,QmatrixWork[irrep].get(),&linsize,Umx,&linsize,&beta,matrix[irrep].get(),&linsize);
        }
    }
}

double OrbitalTransform::TmatRotated(const int index1, const int index2) const
{
    if(!OneBodyMatrixElements.size())
    {
        std::cerr << "First build the OneBodyMatrixElements!" << std::endl;
        exit(EXIT_FAILURE);
    }

    const int irrep1 = _hamorig->getOrbitalIrrep(index1);
    const int irrep2 = _hamorig->getOrbitalIrrep(index2);

    if (irrep1 != irrep2)
        //From now on: both irreps are the same.
        return 0.0;

    int shift = _index->getNstart(irrep1);

    return OneBodyMatrixElements[irrep1][index1-shift + _index->getNORB(irrep1) * (index2-shift)];
}

void OrbitalTransform::fillConstAndTmat(Hamiltonian& Ham) const
{
    //Constant part of the energy
    double value = _hamorig->getEconst();
    Ham.setEconst(value);

    //One-body terms: diagonal in the irreps
    for (int irrep=0; irrep<numberOfIrreps; irrep++)
    {
        const int linsize = _index->getNORB(irrep);
        for (int cnt1=0; cnt1<linsize; cnt1++)
            for (int cnt2=cnt1; cnt2<linsize; cnt2++)
            {
                int shift = _index->getNstart(irrep);
                Ham.setTmat( cnt1 + shift , cnt2 + shift , TmatRotated(cnt1 + shift , cnt2 + shift) );
            }
    }
}

double OrbitalTransform::get_difference_orig(Hamiltonian& hamin) const
{
    return _hamorig->get_difference(&hamin);
}

void OrbitalTransform::CheckDeviationFromUnitary() const
{
    _unitary->CheckDeviationFromUnitary();
}

void OrbitalTransform::set_unitary(UnitaryMatrix unit)
{
    _unitary.reset(new UnitaryMatrix(unit));
}

void OrbitalTransform::rotate_active_space_vectors(double *eigenvecss)
{
    _unitary->rotate_active_space_vectors(eigenvecss , mem1.get());
}

double OrbitalTransform::get_norb(int irrep) const
{
    return _index->getNORB(irrep);
}

void OrbitalTransform::update_unitary(double * change)
{
    _unitary->updateUnitary(mem1.get(),mem2.get(),change,1);//First 1, sets multiply.
}

/**
 * Update _ham_rot in place with a jacobi rotation between orbital k and l over
 * an angle of theta. k and l should be in the same irrep.
 * @param k the first orbital
 * @param l the second orbital
 * @param theta the angle to rotation over
 * @return an rotated Hamiltonian
 */
void OrbitalTransform::do_jacobi_rotation(int k, int l, double theta, Hamiltonian * _ham_rot)
{
    if(!_ham_rot)
        _ham_rot = new Hamiltonian(*_hamorig);

    assert(_ham_rot->getOrbitalIrrep(k) == _ham_rot->getOrbitalIrrep(l)) , "Rotate only in the same irrep. If you remain stubborn, redo the calculation for the matrixelements in c1 symmetry, so you can rotate the chosen orbitals.";

    const OptIndex index = _hamorig->get_index_object();

    const int irrep = _index->get_irrep_each_orbital()[k];
    const int linsize = _index->getNORB(irrep);
    const int shift = _index->getNstart(irrep);
    // the relative index in the irrep
    const int k2 = k - shift;
    const int l2 = l - shift;

    const double cos = std::cos(theta);
    const double sin = std::sin(theta);

    // the one particle integrals

    // first copy element that we are gonna overwrite
    const double tmpkk = _ham_rot->getTmat(k,k);
    const double tmpll = _ham_rot->getTmat(l,l);
    const double tmpkl = _ham_rot->getTmat(k,l);

    double tmp;
    tmp = cos*cos*tmpkk+sin*sin*tmpll-2*cos*sin*tmpkl;
    _ham_rot->setTmat(k,k,tmp);

    tmp = cos*cos*tmpll+sin*sin*tmpkk+2*cos*sin*tmpkl;
    _ham_rot->setTmat(l,l,tmp);

    tmp = tmpkl*(cos*cos-sin*sin)+cos*sin*(tmpkk-tmpll);
    _ham_rot->setTmat(k,l,tmp);

    for(int a=shift;a<linsize+shift;a++)
    {
        if(a == k || a == l)
            continue;

        const double tmpk = _ham_rot->getTmat(k,a); 
        const double tmpl = _ham_rot->getTmat(l,a); 

        tmp = cos * tmpk - sin * tmpl;
        _ham_rot->setTmat(a,k,tmp);

        tmp = cos * tmpl + sin * tmpk;
        _ham_rot->setTmat(a,l,tmp);
    }


    // the two particle elements
    // this is code from Sebastian (CheMPS2) adapted to 
    // Jacobi rotations in one irrep
    // Two-body terms --> use eightfold permutation symmetry in the irreps :-)
    for (int irrep1 = 0; irrep1<numberOfIrreps; irrep1++)
        for (int irrep2 = irrep1; irrep2<numberOfIrreps; irrep2++)
        {
            const int productSymm = SymmInfo.directProd(irrep1, irrep2);
            for (int irrep3 = irrep1; irrep3<numberOfIrreps; irrep3++)
            {
                const int irrep4 = SymmInfo.directProd(productSymm, irrep3);
                // Generated all possible combinations of allowed irreps
                if(irrep4>=irrep2 && (irrep4==irrep || irrep3==irrep || irrep2==irrep || irrep1==irrep))
                {
                    const int linsize1 = index.getNORB(irrep1);
                    const int linsize2 = index.getNORB(irrep2);
                    const int linsize3 = index.getNORB(irrep3);
                    const int linsize4 = index.getNORB(irrep4);

                    if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0))
                    {
                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ]
                                            = _ham_rot->getVmat(_index->getNstart(irrep1) + cnt1,_index->getNstart(irrep2) + cnt2, _index->getNstart(irrep3) + cnt3, _index->getNstart(irrep4) + cnt4 );

                        int rightdim = linsize2 * linsize3 * linsize4;
                        // (ijkl) -> (ajkl)
                        if(irrep1 == irrep)
                            for(int d=0;d<rightdim;d++)
                            {
                                const double tmpk = mem1[d*linsize1+k2];
                                const double tmpl = mem1[d*linsize1+l2];

                                mem1[d*linsize1+k2] = cos*tmpk-sin*tmpl;
                                mem1[d*linsize1+l2] = cos*tmpl+sin*tmpk;
                            }


                        int leftdim = linsize1 * linsize2 * linsize3;
                        // (ajkl) -> (ajkd)
                        if(irrep4 == irrep)
                            for(int d=0;d<leftdim;d++)
                            {
                                const double tmpk = mem1[k2*leftdim+d];
                                const double tmpl = mem1[l2*leftdim+d];

                                mem1[k2*leftdim+d] = cos*tmpk-sin*tmpl;
                                mem1[l2*leftdim+d] = cos*tmpl+sin*tmpk;
                            }

                        int jump = linsize1 * linsize2 * linsize3;
                        leftdim = linsize1 * linsize2;
                        // (ajkd) -> (ajcd)
                        if(irrep3 == irrep)
                            for (int strideidx=0;strideidx<linsize4;strideidx++)
                                for(int d=0;d<leftdim;d++)
                                {
                                    const double tmpk = mem1[k2*leftdim+d+jump*strideidx];
                                    const double tmpl = mem1[l2*leftdim+d+jump*strideidx];

                                    mem1[k2*leftdim+d+jump*strideidx] = cos*tmpk-sin*tmpl;
                                    mem1[l2*leftdim+d+jump*strideidx] = cos*tmpl+sin*tmpk;
                                }


                        jump = linsize1 * linsize2;
                        rightdim = linsize3 * linsize4;
                        // (ajcd) -> (abcd)
                        if(irrep2 == irrep)
                            for (int strideidx=0;strideidx<rightdim;strideidx++)
                                for(int d=0;d<linsize1;d++)
                                {
                                    const double tmpk = mem1[k2*linsize1+d+jump*strideidx];
                                    const double tmpl = mem1[l2*linsize1+d+jump*strideidx];

                                    mem1[k2*linsize1+d+jump*strideidx] = cos*tmpk-sin*tmpl;
                                    mem1[l2*linsize1+d+jump*strideidx] = cos*tmpl+sin*tmpk;
                                }

                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        _ham_rot->setVmat(_index->getNstart(irrep1) + cnt1,_index->getNstart(irrep2) + cnt2, _index->getNstart(irrep3) + cnt3,_index->getNstart(irrep4) + cnt4, mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ] );

                    } //end if the problem has orbitals from all 4 selected irreps
                } // end if irrep 4 >= irrep2
            }// end run irrep3
        } // end run irrep2

    _unitary->jacobi_rotation(_index->get_irrep_each_orbital()[k], k,l , theta); //to keep track of the total unitary in orbitaltransform.
    //return _ham_rot;
}

/* vim: set ts=4 sw=4 expandtab :*/
