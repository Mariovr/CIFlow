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
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>


#include "Permutator.h"
#include "Hamiltonian.h"
#include "CIMethod.h"
#include "CIDens.h"
#include "UnitaryMatrix.h"
#include "Options.h"
#include "scpp_assert.h"

#ifdef __OMP__
#include <omp.h>
#endif 

using namespace std;

inline std::ostream& operator << (std::ostream& os, const Vector2d & m) 
{
    for( unsigned r =0; r<m.num_rows(); ++r ) {
	for( unsigned c=0; c<m.num_cols(); ++c ) {
            os << m(r,c);
	    if( c + 1 < m.num_cols() )
	        os << "\t";
        }
        os << "\n";
    }
    return os;
} 
CIDens::CIDens(CIMethod * cim)
{
    _cim = cim;
    _norb = _cim->get_l();
    //Only allocate memory when explicitely asked to construct the density matrices.
    _one_dens.resize(0);
    _two_dens.resize(0);
}


double CIDens::get_two_rdm(int spin , int orb1 , int orb2, int orb3 , int orb4, const vector<Vector2d > & twordm) const
{
    if(spin == 1)
        return twordm[1](orb1 *dim() +orb2 , orb3*dim() + orb4);
    else
    {
        if(orb1 == orb2 || orb3 == orb4)
        {
            return 0.;
        }
        int sign = 1.;
        if(orb1 > orb2)
        {
            int temp = orb2;
            orb2 = orb1;
            orb1 = temp;
            sign *= -1.;
        }
        if(orb3 > orb4)
        {
            int temp = orb3;
            orb3 = orb4;
            orb4 = temp;
            sign *= -1.;
        }

        return sign* twordm[spin](orb1 *dim() + orb2 - (orb1 +1) *(orb1+2)/2  , orb3*dim() + orb4 - (orb3 +1) *(orb3+2)/2  );
    }
}

void CIDens::add_two_rdm(int spin , int orb1 , int orb2 , int orb3 , int orb4, double value, vector<Vector2d > & twordm) 
{
    if(spin == 1)
        twordm[1](orb1 *dim() +orb2 , orb3*dim() + orb4) += value ;
    else
    {
        if(orb1 == orb2 || orb3 == orb4)
        {
            return ;
        }
        if(orb1 > orb2)
        {
            int temp = orb2;
            orb2 = orb1;
            orb1 = temp;
            value *= -1.;
        }
        if(orb3 > orb4)
        {
            int temp = orb3;
            orb3 = orb4;
            orb4 = temp;
            value *= -1.;
        }
        //cout << "first " << orb1 *dim() + orb2 - (orb1 +1) *(orb1+2)/2 << std::endl;
        //cout << "second : "<<  orb3*dim() + orb4 - (orb3 +1) *(orb3+2)/2<< std::endl;
        twordm[spin](orb1 *dim() + orb2 - (orb1 +1) *(orb1+2)/2  , orb3*dim() + orb4 - (orb3 +1) *(orb3+2)/2  ) += value;
    }

}

//We only save the unique elements onerdm is symmetric so we keep the upper triangular matrix:
double CIDens::get_one_rdm(int spin , int orb1 , int orb2, const valarray<valarray<double>> & onerdm)const
{
    if (orb1 > orb2)
    {
        int temp = orb2;
        orb2 = orb1;
        orb1 = temp;
    }
    return onerdm[spin][(orb1 *  dim()) + orb2 - (orb1 * (orb1+1) /2.)];
}

void CIDens::add_one_rdm(int spin , int orb1 , int orb2, double value, valarray<valarray<double>> & onerdm)
{
    //Watch out in this implementation it is not necessary because the only contribution to the offdiagonal elements comes from a function that doesn't use the arrange_order function so orb1 and orb2 is always in correct order. If you use this function to fill in a context where this is not the case uncomment the following (also do it when the SCPP test assert fails).
    //if (orb1 > orb2)
    //{
        //int temp = orb2;
        //orb2 = orb1;
        //orb1 = temp;
    //}
    SCPP_TEST_ASSERT(orb1 <= orb2, "Error in the creation of the one rdm: Make sure orb1 < orb2, orb1 = " << orb1 << " orb2 = " << orb2);
    onerdm[spin][(orb1 *  dim()) + orb2 - (orb1 * (orb1+1) /2.)] += value;
}

void CIDens::allocate_one_memory(valarray<valarray<double>> & onerdm)
{
    //Two sort of 1rdm, one for the up and one for the down electrons, because for all the relevant determinants in a given part of the Hilbert Space -> the number of up and down electrons stays constant.
    onerdm = valarray<valarray<double>>(2);
    for (int i = 0; i < 2; i++) 
    {
        onerdm[i] = valarray<double> (0., dim()* ( dim()+1)/2. );
    }
}

void CIDens::allocate_two_memory(vector<Vector2d > & twordm )
{
    // Four sorts of relevant 2rdm: 0 -> up, up, up ,up ,1 -> up, down, up, down, 2 -> down, up, down, up,
    // 3-> down, down, down, down (where the first two are for creation and the last two for annihilation.) ALL other combinations are irrelevant because when contracted with the matrixelements they give zero because they integrate down with up spinorbitals or they dont keep the number of ups and downs constant, REMARK a_i+ a_j+ a_k a_l -> 2rdm(ij , lk) because 2rdm corresponds to physicist notation of inproduct.
    // In memory we keep only the 0,1,3 cases , because the 2 is completely equal to the 1. 
    // Because other notation of 2rdm = (pq,rs) with p < q , r < s (this for spin orbitals and we assume up < down for corresponding spatial orbitals (first L ups than L downs) )
    // So for case 1 we keep all elements because p < q and r < s always (u,d,u,d) , and for 0, and 3 we keep only the elements with p<q , r<s and the other ones are generated by changing fermion sign.
    twordm= vector<Vector2d >(3);
    SCPP_TEST_ASSERT(twordm.size() == 3, "Error: _two_dens size is not 3 before allocation of new memory. size " <<twordm.size());

    for (int a = 0; a < 3; a++) 
    {
        if (a != 1)
            twordm[a] = Vector2d((dim() -1/2) * dim() , (dim() -1/2) * dim()) ;
        else
            twordm[a] = Vector2d(dim()* dim() , dim() * dim()) ;
    }
}

void CIDens::construct_density(bool twordm)
{
    build_parallel(twordm);
    #ifdef _DEBUG
        test_invariants(twordm);
    #endif
}

void CIDens::reset_1rdm()
{
    _one_dens.resize(0);
}

void CIDens::reset_2rdm()
{
    _two_dens.resize(0);
}

matrix CIDens::spin_summed_1rdm()const {
    matrix spinsummed = matrix( dim(), dim());
    for(int i = 0 ; i < 2 ; i++){
        for(int j = 0 ; j <  dim(); j++){
            for(int k = 0 ; k <  dim(); k++){
                spinsummed(j , k ) = get_one_rdm(i,j,k, _one_dens);
            }
        }    
    }
    return spinsummed;
}

void CIDens::test_invariants(bool twordm)
{
    double som = 0;
    int N = _cim->gNup() + _cim->gNdown();
    for(int a = 0 ; a < 2 ; a++)
    {
        for( int i = 0; i <  dim(); i++)
        {
            som += get_one_rdm(a,i,i, _one_dens);
        }
    }
    SCPP_TEST_ASSERT(fabs(som -N ) < 1e-11 , "Error in density matrices: Trace of the 1rdm : " << som << " is not equal to the sum of Nup: " << _cim->gNup() << " and Ndown:  " << _cim->gNdown()<< std::endl );

    print_one_dens(cout);
    print_two_dens(cout);

    if (twordm)
    {
        som = 0;
        //partial trace over the 2rdm(PQ,PQ) with the SPIN ORBITALS (P<Q) (we assume ups are smaller than down spin orbitals). so iup jup jup iup and id,jd,jd,id -> only offdiag with i < j
        //and iu,jd,iu,jd all elements. And no d,u,d,u elements because we assumed we ordered first the ups and then the downs
        for( int i = 0; i <  dim(); i++)
        {
            for( int j = 0; j <  dim(); j++)
            {
                som += get_two_rdm(1,i,j,i,j);
            }
            for( int j = i+1; j <  dim(); j++)
            {
                som += get_two_rdm(0,i,j,i,j);
                som += get_two_rdm(2,i,j,i,j);
            }
        }
        SCPP_TEST_ASSERT( fabs(som - 1/2. * N *(N-1)) < 1e-11, "Error in density matrices: Trace of the 2rdm: " << som << " is not equal to the number of possible pair formings: 1/2*N *(N-1),  " << 1/2. * N * (N-1));

        SCPP_TEST_ASSERT(fabs(get_dens_energy() - _cim->get_ci_energy() ) < 1e-11, setprecision(14) << "Error in density matrices: Energy of state: " <<  _cim->get_ci_energy() << " is not equal to energy calculated from the density matrices.: " << get_dens_energy() << std::endl );

        if(_cim->get_name() == "FCI")
        {
            double spinsq = get_spin_squared();
            SCPP_TEST_ASSERT(spinsq - floor(fabs(spinsq))  < 1e-9 || fabs(4*(spinsq - floor(fabs(spinsq)) ) - 3) < 1e-9, setprecision(14) << "Error in density matrices: Spin squared is not an integer number, or has .75 as decimal numbers:  spin squared: " <<  spinsq << std::endl );
        }
    }

}

unsigned CIDens::dim() const
{
    return _norb;
}

double CIDens::get_seniority() const
{
    double som = 0;
    for(int i = 0 ; i < 2 ; i++)
    {
        for(int j = 0 ; j < dim() ; j++)
	{
            som += get_one_rdm(i,j,j, _one_dens);
	}
    }
    for (int i = 0; i < dim(); i++)
	{
        som -=  2.*get_two_rdm(1,i,i,i,i,_two_dens); //factor 2 because of contr. of u , d, u,d and d,u,d,u
    }
    return som;
}

double CIDens::get_spin_squared() const
{
    double spsm = 0;
    for (int k = 0; k < dim(); k++) 
    {
        for (int i = 0; i < dim(); i++)
        {
            spsm -= get_two_rdm(1,k,i,i,k,_two_dens);
        }
    }
    //cout << "sz: " << _cim->get_sz()<< std::endl;
    return spsm - _cim->get_sz() + _cim->gNup() +pow(_cim->get_sz(),2) ;
}

double CIDens::get_spin() const
{
    return -1+std::sqrt(1+4*get_spin_squared() )/2.;
}

void CIDens::get_no(matrix &occupations, matrix &NO )const
{
    matrix spinsummed = spin_summed_1rdm();
    for(int irrep = 0 ; irrep < _cim->get_ham()->getNirreps() ; irrep++){
        int nstart =_cim->get_ham()->getNstart(irrep);
        int norb = _cim->get_ham()->getNORB(irrep);
        //handles the case if there are no orbitals associated to a given irrep.
        if( norb > 0){
            matrix NO_part = matrix(norb ,norb );
            matrix occ_part = matrix(norb ,1 );
            for(int i = 0 ; i < norb; i++){
                for(int j = i ; j < norb; j++){
                    NO_part(i,j) = spinsummed(i+nstart,j+nstart);
                    if(j!= i)
                        NO_part(j,i) = spinsummed(j+nstart,i+nstart);
                }
            }
            NO_part.diagonalize(occ_part,NO_part);
            if(CIMethod_debug)
	    {
                NO_part.Print();
		occ_part.Print();
            }
            for(int i = nstart ; i < norb+nstart; i++){
                occupations(i,0) = occ_part(nstart+norb-1 -i , 0);
                for(int j = nstart ; j < norb+nstart ; j++){
                    NO(j,i) = NO_part(j-nstart ,nstart+norb-1-i);
                }
            }
            //spinsummed.diagonalize(occupations, NO);
            //NO contains the natural orbitals and occupations the corresponding occupations.
            if(CIMethod_debug){
                std::cout << "#The occupations of the Natural Orbitals:"  <<std::endl; 
                occupations.Print();
                std::cout << "#The occupations of the Natural Orbitals:"  <<std::endl; 
		NO.Print();
                std::cout << "#Check if the Natural Orbitals are unitary"   <<std::endl;
                matrix unit = matrix(dim() , dim());
                unit.prod(NO , NO , 'N', 'T');
                unit.Print();
            }
        }
    }
}

void CIDens::build_parallel(bool twordm)
{
    #ifdef __OMP__
        int num_t = omp_get_max_threads();
    #else
        int num_t = 1;
    #endif
    unsigned long long num_elems = (_cim->get_dim()*1ull*(_cim->get_dim()+1ull))/2.;
    unsigned long long size_part = num_elems/num_t + 1;

    // every thread should process the lines between i and i+1
    // with i the thread number
    std::vector<unsigned int> workload(num_t+1);
    workload.front() = 0;
    workload.back() = _cim->get_dim();

    for(int i=1;i<num_t;i++)
    {
        auto num_lines = workload[i-1];
        unsigned long long num_elems = 0;

        while(num_elems < size_part)
           num_elems += _cim->get_dim() - num_lines++;

        if(num_lines > _cim->get_dim() )
           num_lines = _cim->get_dim();

        workload[i] = num_lines;
    }
    vector< valarray<valarray<double>> > parts_one(num_t);
    vector<vector<Vector2d >  > parts_two(num_t);
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
        allocate_one_memory(parts_one[me]);
        if(twordm)
            allocate_two_memory(parts_two[me]);

		construct_CI_one_dens(workload[me] , workload[me+1], parts_one[me] );
        if(twordm)
            construct_CI_two_dens(workload[me] , workload[me+1], parts_two[me] );
    }// End parallel

    add_from_vector_one(parts_one); //REMARK because we only filled the upper diagonal => _mat should be a SparseMatrix_CRS_Sym type.
    if (twordm)
        add_from_vector_two(parts_two); //REMARK because we only filled the upper diagonal => _mat should be a SparseMatrix_CRS_Sym type.
}

void CIDens::add_from_vector_one(vector< valarray<valarray<double>> > & parts)
{
    allocate_one_memory(_one_dens);
    //We add all the parts of the one rdm to the result in _one_dens
    for(int i = 0 ; i < parts.size() ; i++ )
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < dim(); k++)
            {
                for (int l = k; l < dim(); l++)
                {
                    double value = get_one_rdm(j,k,l, parts[i]);
                    add_one_rdm(j,k,l, value , _one_dens );
                }
            }
        }
    }
    //print_one_dens(cout);
}

void CIDens::add_from_vector_two(vector<vector<Vector2d >  > & parts)
{
    //We add all the parts of the one rdm to the result in parts[0], after that we use the moveconstructor of vector to move everything in _two_dens.
    allocate_two_memory(_two_dens);
    for(int p = 0 ; p < parts.size() ; p++ )
    {
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < dim(); i++){
                int k;
                if (j == 1)
                    k = 0;
                else
                    k = i+1;
                    
                for (; k < dim(); k++){
                    for (int l = 0; l < dim(); l++) 
                    {
                        int m;
                        if (j == 1)
                            m = 0;
                        else
                            m = l+1;
                        for (; m < dim(); m++) 
                        {
                            if( j ==1 )
                            {
                                if( i*dim() + k <=  m*dim() +  l)
                                {
                                     double value = get_two_rdm(j,i,k,m,l, parts[p]);
                                     add_two_rdm(j,i,k,m,l, value , _two_dens );
                                }
                            }
                            else{
                                if( i *dim() + k - (i+1) *(i+2)/2.  <= l*dim() + m- ( l+1) *( l+2)/2. )
                                {
                                    double value = get_two_rdm(j,i,k,m,l, parts[p]);
                                    add_two_rdm(j,i,k,m,l, value , _two_dens );
                                }

                            }
                        }//end for m
                    }//End for l
                }//End for k
            }   
        }
    }
}


void CIDens::fill_one_dense(TYPE bitstring1 , TYPE bitstring2, TYPE same ,  double contr, int up_or_down,valarray<valarray<double>> & onerdm , vector<Vector2d > & twordm ){
    int * orbs = _cim->get_orbs(bitstring1 ^ bitstring2 , 2);
    int * otherorbs;
    int * sameorbs;
    if(up_or_down == 2){
        otherorbs = _cim->get_orbs(bitstring1 & bitstring2, _cim->gNup() -1 );
        sameorbs = _cim->get_orbs(same, _cim->gNdown() );
    }    
    else if(up_or_down == 3){
        otherorbs = _cim->get_orbs(bitstring1 & bitstring2, _cim->gNdown() -1 );
        sameorbs = _cim->get_orbs(same, _cim->gNup() );
    }
    int diff = _cim->get_nstacked(orbs, bitstring1 & bitstring2);
    if (diff % 2 != 0 )
        contr *= -1;
    if(up_or_down<2)
    {
        add_one_rdm(up_or_down,orbs[0],orbs[1], contr, onerdm); //make sure one_dense is initialized to zero. Contribution of transposed is gone.
    }
    else
    {
        //up_or_down =2 -> up orbs, up_or_down = 3 -> down orbs
        if(up_or_down == 2){
            for(int i = 0 ; i < _cim->gNup()-1 ; i ++)
            {
                add_two_rdm(0,orbs[0],otherorbs[i],orbs[1],otherorbs[i],  contr,twordm); // all ups
                //add_two_rdm(0,otherorbs[i],orbs[0],otherorbs[i],orbs[1] , contr); // all ups
                //add_two_rdm(0,orbs[0],otherorbs[i],otherorbs[i],orbs[1], -1.*contr); // all ups
                //add_two_rdm(0,otherorbs[i],orbs[0],orbs[1],otherorbs[i],  -1.*contr ); // all ups
            }    
            for(int i = 0 ; i < _cim->gNdown() ; i ++)
            {
                add_two_rdm(1,orbs[0] ,sameorbs[i] ,orbs[1] ,sameorbs[i] , contr,twordm); //sameorbs[i] are downs.
            }
            //CONTRIUBTION of the transposed.
            /*
            for(int i = 0 ; i < _cim->gNup()-1 ; i ++)
            {
                add_two_rdm(0,orbs[1],otherorbs[i],orbs[0],otherorbs[i],  contr    ); // all ups
                //add_two_rdm(0,otherorbs[i],orbs[1],otherorbs[i],orbs[0], contr     ); // all ups
                //add_two_rdm(0,otherorbs[i],orbs[1],orbs[0],otherorbs[i], -1.*contr ); // all ups
                //add_two_rdm(0,orbs[1],otherorbs[i],otherorbs[i],orbs[0], -1.*contr ); // all ups
            }    
            for(int i = 0 ; i < _cim->gNdown() ; i ++)
            {
                add_two_rdm(1,orbs[1],sameorbs[i],orbs[0],sameorbs[i], contr) ; //sameorbs[i] are downs.
            }    
            */
        }
        else{
            for(int i = 0 ; i < _cim->gNdown()-1 ; i ++){
                add_two_rdm(2,orbs[0],otherorbs[i],orbs[1],otherorbs[i], contr ,twordm); // all downs
                //add_two_rdm(2,otherorbs[i],orbs[0],otherorbs[i],orbs[1], contr ); // all downs
                //add_two_rdm(2,otherorbs[i],orbs[0],orbs[1],otherorbs[i], -1.*contr ); // all downs
                //add_two_rdm(2,orbs[0],otherorbs[i],otherorbs[i],orbs[1], -1.*contr ); // all downs
            }    
            for(int i = 0 ; i < _cim->gNup() ; i ++){
                //sameorbs[i] ups.
                add_two_rdm(1 ,sameorbs[i] ,orbs[0] ,sameorbs[i] ,orbs[1] , contr ,twordm); //make sure two_dense is initialized to zero.
            }    
            //CONTRIUBTION of the transposed.
            /*
            for(int i = 0 ; i < _cim->gNdown()-1 ; i ++){
                add_two_rdm(2, orbs[1] ,otherorbs[i],orbs[0] ,otherorbs[i], contr); // all downs
                //add_two_rdm(2, otherorbs[i] ,orbs[1],otherorbs[i] ,orbs[0], contr); // all downs
                //add_two_rdm(2, orbs[1] ,otherorbs[i],otherorbs[i] ,orbs[0], -1.*contr); // all downs
                //add_two_rdm(2, otherorbs[i] ,orbs[1],orbs[0] ,otherorbs[i], -1.*contr); // all downs
            }    
            for(int i = 0 ; i < _cim->gNup(); i ++)
            {
                add_two_rdm( 1 ,sameorbs[i] ,orbs[1] ,sameorbs[i] ,orbs[0] , contr ); //make sure two_dense is initialized to zero.
            }     
            */

        }
    }
    delete [] orbs;
    if(up_or_down > 1){
        delete [] otherorbs;
        delete [] sameorbs;
    }    
}

void CIDens::fill_two_dense(TYPE up1 , TYPE up2, double contr, int up_or_down,vector<Vector2d > & twordm){
    int * orbs = _cim->get_orbs(up1^ up2 , 4);
    _cim->arrange_order(orbs , 4 , up1 , up2);
    int checkorbs[2] = {orbs[0] , orbs[1]};
    int diff = _cim->get_nstacked(checkorbs, up1 & up2);
    checkorbs[0] = orbs[2]; checkorbs[1] = orbs[3];
    diff += _cim->get_nstacked(checkorbs,up1 & up2);
    
    if (diff % 2 != 0 ){
        contr *= -1;
    }
    add_two_rdm(up_or_down,orbs[0],orbs[1],orbs[2],orbs[3],  contr ,twordm); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[0],orbs[1],orbs[3],orbs[2],  -1.*contr ); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[1],orbs[0],orbs[2],orbs[3],  -1.*contr ); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[1],orbs[0],orbs[3],orbs[2],  contr ); //make sure two_dense is initialized to zero.
    //Contribution of transposed.
    //add_two_rdm(up_or_down,orbs[2],orbs[3],orbs[0],orbs[1],  contr ); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[2],orbs[3],orbs[1],orbs[0],  -1.* contr ); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[3],orbs[2],orbs[0],orbs[1],  -1.* contr ); //make sure two_dense is initialized to zero.
    //add_two_rdm(up_or_down,orbs[3],orbs[2],orbs[1],orbs[0],  contr ); //make sure two_dense is initialized to zero.
    delete [] orbs;
}

void CIDens::fill_two_dense_upanddown(TYPE up1 , TYPE up2, TYPE down1, TYPE down2, double contr ,vector<Vector2d > & twordm )
{
    int * orbs = _cim->get_orbs(up1 ^ up2 , 2);
    int * downorbs = _cim->get_orbs(down1 ^ down2 , 2);
    _cim->arrange_order(orbs , 2 , up1 , up2);
    _cim->arrange_order(downorbs , 2 , down1 , down2);
    int diff = _cim->get_nstacked(orbs, up1 & up2);
    diff += _cim->get_nstacked(downorbs , down1 & down2);

    if (diff % 2 != 0 /*|| downorbs[0] > orbs[0]*/)
    {
        contr *= -1;
    }

    //TYPE d = 1;
    //SCPP_TEST_ASSERT( ( (up1 & (d<<orbs[0] ) )> 0 )&& ( (up2 & (d<<orbs[1] ) )> 0 ) , "Error one of the different orbitals is not associated to the right up, up1:  orbs: " << orbs[0] << " orbs:  " << orbs[1] << " gives values : " << (up2 & (d <<orbs[1] ) ));
    //SCPP_TEST_ASSERT((down1 & (d<<downorbs[0] ) )> 0 && (down2 & (d<<downorbs[1] ) )> 0 , "Error one of the different orbitals is not associated to the right down, down1:  orbs: " << downorbs[0] << " down2: orbs:  " << downorbs[1] );

    add_two_rdm(1,orbs[0],downorbs[0],orbs[1],downorbs[1], contr, twordm); //make sure two_dense is initialized to zero.
    //Contribution of transposed.,
    //add_two_rdm(1,orbs[1],downorbs[1],orbs[0],downorbs[0], contr); //make sure two_dense is initialized to zero.

    delete [] downorbs;
    delete [] orbs;
}

void CIDens::fill_diagonal(TYPE up1 , TYPE down1, double contr, bool two_rdm,valarray<valarray<double>> & onerdm , vector<Vector2d > & twordm )
{
    int *  orbsup = _cim->get_orbs(up1,_cim->gNup());
    int * orbsdown = _cim->get_orbs(down1,_cim->gNdown());
    //first part of up,up,up,up of the rdms
    for(int i = 0 ; i < _cim->gNup() ; i++){
        if(!two_rdm){
            add_one_rdm(0,orbsup[i],orbsup[i], contr, onerdm);
        }
        else{
            for(int j = i+1 ; j < _cim->gNup() ; j++){
                //REMARK: when i =j it is trivially zero
                add_two_rdm(0,orbsup[i],orbsup[j],orbsup[i],orbsup[j], contr, twordm); //make sure two_dense is initialized to zero
                //add_two_rdm(0,orbsup[j],orbsup[i],orbsup[j],orbsup[i],contr); //make sure two_dense is initialized to zero
                //add_two_rdm(0,orbsup[i],orbsup[j],orbsup[j],orbsup[i], -1.*contr); //make sure two_dense is initialized to zero
                //add_two_rdm(0,orbsup[j],orbsup[i],orbsup[i],orbsup[j],-1.*contr); //make sure two_dense is initialized to zero
            }	    
        }
    }

    //Contribution to the down,down,down,down part of the rdms 
    for(int i = 0 ; i < _cim->gNdown() ; i++){
        if(!two_rdm){
            add_one_rdm(1,orbsdown[i],orbsdown[i], contr, onerdm);
        }
        else{
            for(int j = i+1 ; j < _cim->gNdown() ; j++){
                add_two_rdm(2,orbsdown[i],orbsdown[j],orbsdown[i],orbsdown[j], contr , twordm); //make sure two_dense is initialized to zero
                //add_two_rdm(2,orbsdown[j],orbsdown[i],orbsdown[j],orbsdown[i],contr ); //make sure two_dense is initialized to zero
                //add_two_rdm(2,orbsdown[i],orbsdown[j],orbsdown[j],orbsdown[i], -1.*contr ); //make sure two_dense is initialized to zero
                //add_two_rdm(2,orbsdown[j],orbsdown[i],orbsdown[i],orbsdown[j],-1.*contr ); //make sure two_dense is initialized to zero
            }	    
        }
    }

    //Contribution to the mixing parts.
    if(two_rdm){
        for(int i = 0 ; i < _cim->gNup() ; i++){
            for(int j = 0 ; j < _cim->gNdown() ; j++){
                add_two_rdm(1,orbsup[i],orbsdown[j],orbsup[i],orbsdown[j],  contr,twordm); //make sure two_dense is initialized to zero
            }	    
        }
    }	
    delete [] orbsup;
    delete [] orbsdown;
}

void CIDens::fill_density_matrix(TYPE up1 , TYPE down1 , TYPE up2, TYPE down2, double contr, vector<Vector2d > & twordm){
    if(_cim->_perm->popcount( down1 ^ down2) > 4  || _cim->_perm->popcount(up1 ^ up2) > 4){ 
        return;
    }
    else if (_cim->_perm->popcount( up1 ^ up2) == 4 && _cim->_perm->popcount(down1 ^ down2) == 0){ // two different occupied orbitals
        fill_two_dense(up1, up2, contr, 0, twordm);
    }
    else if(_cim->_perm->popcount(up1 ^ up2) == 2){ // one different occupied orbital-> two posibilities: 1) downs are the same , 2) downs have also one differen occupied orbital.
        if(_cim->_perm->popcount(down1 ^ down2) == 0){
            //1)downs are the same:
            fill_one_dense(up1,up2 , down1 , contr, 2, _one_dens , twordm);
        }//end downs are the same
        else if(_cim->_perm->popcount(down1 ^ down2) == 2){
            //2) downs differ also by own orbital.
            fill_two_dense_upanddown(up1 , up2, down1, down2, contr ,twordm);
        } // ends down and ups differ by one
    } // end ups differ by one.
    else if (_cim->_perm->popcount( down1 ^ down2) == 4 && _cim->_perm->popcount(up1 ^ up2) == 0){ // two different occupied orbitals of the downs
        fill_two_dense(down1, down2, contr, 2,twordm);
    } //end downs differ by two.
    else if(_cim->_perm->popcount(down1 ^ down2) == 2 &&_cim->_perm->popcount(up1 ^ up2) == 0 ){ // one different occupied orbital of downs->  ups are the same (downs and up differ by 1 is already taken into acount above)
        fill_one_dense(down1,down2 , up1 , contr, 3, _one_dens , twordm);
    } // end down differ by one.
    else if(_cim->_perm->popcount( down1 ^ down2) == 0 && _cim->_perm->popcount(up1 ^ up2) == 0){ // diagonal contribution
        fill_diagonal(up1 , down1,contr,1, _one_dens , twordm); //set 1 for 2rdm
    }
}

void CIDens::fill_density_matrix_one(TYPE up1 , TYPE down1 , TYPE up2 , TYPE down2, double contr,valarray<valarray<double>> & onerdm){
    if(_cim->_perm->popcount( up1 ^ up2) > 2 || _cim->_perm->popcount(down1 ^ down2) > 2){
        return; 
    }
    else if (_cim->_perm->popcount( up1 ^ up2) == 2 && _cim->_perm->popcount(down1 ^ down2) == 0){
        fill_one_dense(up1,up2 , down1 , contr, 0, onerdm , _two_dens);
    }
    else if(_cim->_perm->popcount( up1 ^ up2) == 0 && _cim->_perm->popcount(down1 ^ down2) == 2){
        fill_one_dense(down1,down2 , up1 , contr, 1, onerdm , _two_dens);
    }
    else if(_cim->_perm->popcount( up1 ^ up2) == 0 && _cim->_perm->popcount(down1 ^ down2) == 0){
        fill_diagonal(up1, down1, contr,0, onerdm , _two_dens); //set 0  for 1rdm
    }		    
}

double CIDens::get_one_electron_energy(){
    double energy = 0.0;
    for (int i = 0; i < dim(); i++) {
        for (int j = 0; j < dim(); j++) {
            energy += _cim->get_ham()->getTmat(i,j) * get_one_rdm(0,i,j, _one_dens);
            energy += _cim->get_ham()->getTmat(i,j) * get_one_rdm(1,i,j, _one_dens);
        }
    }
    return energy;
}

double CIDens::get_two_electron_energy() {
    double energy = 0.0;
    for (int i = 0; i < dim(); i++) {
        for (int j = 0; j < dim(); j++) {
            for (int k = 0; k < dim(); k++) {
                for (int l = 0; l < dim(); l++) {
                    double vmat = _cim->get_ham()->getVmat(i,j,k,l);
                    energy +=  vmat * get_two_rdm(0,i,j,k,l,_two_dens);
                    energy +=  2*vmat * get_two_rdm(1,i,j,k,l,_two_dens); //for u,d,u,d and d,u,d,u which are equal.
                    energy +=  vmat * get_two_rdm(2,i,j,k,l,_two_dens);
                }
            }
        }
    }
    return 0.5*energy;
}

double CIDens::get_dens_energy() {
    return get_one_electron_energy() + get_two_electron_energy() + _cim->get_econst(); 
}

void CIDens::print_one_dens(std::ostream & os, const valarray<valarray<double>> & onerdm)const{
   SCPP_TEST_ASSERT(_one_dens[0].size() > 0, "Error the size of the onerdm is equal to zero");
   for (int j = 0; j < 2; j++) {
     os <<"#One rdm of the " << j << " spin electrons." <<std::endl;
     int count = 0;
     for (int i = 0; i < dim(); i++) {
	 for (int k = i; k < dim(); k++) {
		 if(fabs(get_one_rdm(j,i,k,onerdm)) >= 1e-12 ){
			 os << i << " " << k << " " << get_one_rdm(j,i,k, onerdm) << std::endl;
			 count ++;
		  }
	  }
     }
     os <<"#We obtained " << count << " non zero elements." << std::endl;
   }
}

void CIDens::print_two_dens(std::ostream & os,const  vector<Vector2d > & twordm)const{
   for (int j = 0; j < 3; j++) {
     os << "#Two rdm case : " << j << "---------------------------------------------" << std::endl;
     int count = 0;
     for (int i = 0; i < dim(); i++){
         int k;
         if (j == 1)
             k = 0;
         else
             k = i+1;
             
	     for (; k < dim(); k++){
		     for (int l = 0; l < dim(); l++) {
                 int m;
                 if (j == 1)
                     m = 0;
                 else
                     m = l+1;
			     for (; m < dim(); m++) {
                     if(fabs( get_two_rdm(j,i,k,l,m,twordm)) > 1e-12  ){
                         os << i << " " << k << " " << l << " " << m << " " << get_two_rdm(j,i,k,l,m,twordm) << std::endl;
                         count ++;
                     }
			     }
	              }
              }
     }
     os <<"#We obtained " << count << " non zero elements." << std::endl;
    }
}

void CIDens::compare_one_dens(CIDens *cid){
   double totverschil = 0;
   for (int j = 0; j < 2; j++) {
       std::cout <<"One rdm of the " << j << " spin electrons." <<std::endl;
       for (int i = 0; i < dim(); i++) {
           for (int k = 0; k < dim(); k++) {
               if(get_one_rdm(j,i,k, _one_dens) != 0 || get_one_rdm(j,i,k, cid->get_one_dm() ) != 0){
		   double verschil = get_one_rdm(j,i,k, _one_dens)- get_one_rdm(j,i,k, cid->get_one_dm() );
		   if(abs(verschil) > 1e-10){
		       std::cout << i << " " << k << " this: " << get_one_rdm(j,i,k, _one_dens) << " other: " << get_one_rdm(j,i,k, cid->get_one_dm() ) << " verschil " << verschil <<std::endl;
                   }
               	   totverschil += pow(verschil,2);
               }
           }
       }
     }
     std::cout <<"total difference is: " << sqrt(totverschil) <<endl;
}

void CIDens::compare_two_dens(CIDens *cid){
   double totverschil = 0;
   for (int j = 0; j < 3; j++) {
     std::cout << "Two rdm case : " << j << "---------------------------------------------" << std::endl;
     for (int i = 0; i < dim(); i++){
	     for (int k = 0; k < dim(); k++){
		     for (int l = 0; l < dim(); l++) {
			     for (int m = 0; m < dim(); m++) {
				 if(get_two_rdm(j,i,k,l,m,_two_dens) != 0 || get_two_rdm(j,i,k,l,m, cid->get_two_dm() ) != 0){
				         double verschil = get_two_rdm(j,i,k,l,m,_two_dens) - get_two_rdm(j,i,k,l,m,cid->get_two_dm() );
		                         if(abs(verschil) > 1e-10){
					     std::cout << i << " " << k << " " << l << " " << m << " this " << get_two_rdm(j,i,k,l,m,_two_dens) << " other " << cid->get_two_rdm(j,i,k,l,m, cid->get_two_dm() ) << " verschil " << verschil <<std::endl;
			                 }
					 totverschil += pow(verschil,2);
				 }
		              }
	              }
              }
        }
    }
    std::cout <<"total difference is: " << sqrt(totverschil) <<endl;
}    

//----------------------------SUBCLASSES----------------------------------------------------------

//---------------------------DensDOCI---------------------------------------------------------//
void DensDOCI::construct_CI_one_dens(unsigned int start , unsigned int end, valarray<valarray<double>> & onerdm ) 
{
    int L = dim();
    Permutator_Bit my_perm(dim());
    TYPE up1 = my_perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start;++idx_begin)
        up1 = my_perm.permutate_bit(up1);
    for(long i = start; i < end; i++) 
    {
        for (int j = 0; j < L; j++) 
        {
            double contr = _cim->get_eigvec(i,0) * _cim->get_eigvec(i,0);
            TYPE shiftbit = 1;
            if(up1 & ( shiftbit << j))
            {
                add_one_rdm(0,j,j,contr, onerdm);
                add_one_rdm(1,j,j,contr, onerdm);
            }
        }
        up1 = my_perm.permutate_bit(up1);
    } 
}

void DensDOCI::construct_CI_two_dens(unsigned int start , unsigned int end, vector<Vector2d > & twordm ) 
{
    vector< vector<int> > vw  ( dim() +1,  vector<int> ( _cim->gNup() + 1 , 0));
    setup_vertex_weights(vw);

    Permutator_Bit my_perm(dim());
    TYPE up1 = my_perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start;++idx_begin)
        up1 = my_perm.permutate_bit(up1);

    for (long i = start; i < end ; i++) 
    {
        double contr = _cim->get_eigvec(i,0) * _cim->get_eigvec(i,0);
        for (int j = 0; j < dim() ; j++) //Loop over the L spatial orbitals.
        {
            // first test the string with itself.
            TYPE shiftbit = 1;
            if(up1 & ( shiftbit << j))
            {
                for (int k = j; k < dim(); k++) 
                {
                    if(up1 & ( shiftbit << k))
                    {
                        if (j != k) 
                        {
                            add_two_rdm(0,j,k,j,k , contr,twordm) ;
                            //add_two_rdm(0,j,k,k,j, -1.*contr);

                            add_two_rdm(2,j,k,j,k, contr,twordm) ;
                            //add_two_rdm(2,j,k,k,j, -1.* contr) ;

                            add_two_rdm(1,k,j,k,j, contr,twordm);
                        }
                        add_two_rdm(1,j,k,j,k, contr,twordm);
                    }
                }

                // Then test all the single excitations you can think off.
                TYPE up_interm = up1 ^ (shiftbit << j);
                for (int k = j + 1; k < dim() ; k++) 
                {
                    if(up1 & ( shiftbit << k))
                    {
                        continue;
                    } 
                    else 
                    {
                        TYPE up2 = up_interm ^ (shiftbit << k);
                        unsigned int index = determine_weight(up2, vw);
                        double contr2 = _cim->get_eigvec(i,0) * _cim->get_eigvec(index,0);
                        add_two_rdm(1,j,j,k,k, contr2,twordm);
                        //add_two_rdm(1,k,k,j,j, contr2, twordm); //Contribution of transposed.
                    } 
                }
            }//End jth orbital is occupied in current determinant.
        }
        up1 = my_perm.permutate_bit(up1);
    }//End for determinants
     
}

void DensDOCI::setup_vertex_weights(vector<vector<int>> & vw) 
{
    /* Set up the vertex_weight vector */
    vw[0][0] = 1; /* Set the first element to one. */
    for (int k = 1; k <= ( _norb  - _cim->gNup()) ; k++) 
    { /* Go down the first column; at most num_orbs-num_elecs can be set to one (must have N elecs in orbitals) */
        vw[k][0] = 1;
    }

    for (int m = 1; m <= _cim->gNup(); m++) 
    { /* Fill the remaining columns */
        for (int k = m; k <= ( _norb  - (_cim->gNup() - m)); k++) 
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

unsigned int DensDOCI::determine_weight(TYPE string , const vector<vector<int> > & vw) 
{
    unsigned int weight = 1;
    int num_elecs_interm = 0;

    for (int k = 0; k < _cim->get_l(); k++) 
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

//---------------------------DensFILE-------------------------------------------------------------//
//REMARK DensFILE can be used for FCI_File and for CI_Big because they have the same permutator object.
void DensFILE::construct_CI_one_dens( unsigned int start , unsigned int end, valarray<valarray<double>> & onerdm){
    Permutator_File my_perm(_cim->get_perm()->get_filename(),dim());
    TYPE up1 = my_perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start;++idx_begin)
        up1 = my_perm.permutate_bit(up1);

    if (CIMethod_debug)
        cout << "Started construct_CI_one_dens in DensFILE" << endl;
    TYPE dets_row[2] = {0};
    TYPE dets_col[2] = {0};
    for(unsigned long i = start ; i < end ; i ++){
        if (CIMethod_debug){
            if( i % 5000 == 0)
        	cout << "line : " << i << " of lines: " << my_perm.get_dim() << endl;
        }
       	_cim->_perm->permutate_bit(dets_row,i);
       	TYPE up1 = dets_row[0];
       	TYPE down1 = dets_row[1];
       	for(unsigned long j = i ; j < _cim->get_dim() ; j ++){
       	    _cim->_perm->permutate_bit(dets_col,j);
       	    TYPE up2 = dets_col[0];
            TYPE down2 = dets_col[1];
	        double contr = _cim->get_eigvec(i,0) * _cim->get_eigvec(j,0);
       	    fill_density_matrix_one(up1,down1, up2, down2,contr, onerdm);
        } // end cols 
    } // end rows 
    if (CIMethod_debug)
        cout << "We constructed the density matrix" << endl;
}

void DensFILE::construct_CI_two_dens(unsigned int start , unsigned int end, vector<Vector2d > & twordm){
    Permutator_File my_perm(_cim->get_perm()->get_filename(),dim());
    TYPE up1 = my_perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start;++idx_begin)
        up1 = my_perm.permutate_bit(up1);

    if (CIMethod_debug)
        cout << "Started construct_CI_two_dens in DensFILE" << endl;
    TYPE dets_row[2] = {0};
    TYPE dets_col[2] = {0};
    for(unsigned long i = start ; i < end ; i ++){
        if (CIMethod_debug){
            if( i % 5000 == 0)
        	cout << "line : " << i << " of lines: " << my_perm.get_dim() << endl;
        }
       	_cim->_perm->permutate_bit(dets_row,i);
       	TYPE up1 = dets_row[0];
       	TYPE down1 = dets_row[1];
       	for(unsigned long j = i ; j < _cim->get_dim() ; j ++){
       	    _cim->_perm->permutate_bit(dets_col,j);
       	    TYPE up2 = dets_col[0];
            TYPE down2 = dets_col[1];
	        double contr = _cim->get_eigvec(i,0) * _cim->get_eigvec(j,0);
       	    fill_density_matrix(up1,down1, up2, down2,contr,twordm);
        } // end cols 
    } // end rows 
    if (CIMethod_debug)
        cout << "We constructed the density matrix" << endl;
}


//------------------------------------------DensFCI------------------------------------------
void DensFCI::construct_CI_one_dens(unsigned int start , unsigned int end, valarray<valarray<double>> & onerdm){
    #ifdef __OMP__
    auto num_t = omp_get_max_threads();
    auto me = omp_get_thread_num();
    cout << "paralel " << me << endl;
    #else 
    auto num_t = 1;
    auto me = 0;
    #endif
    unsigned long long size_part = _cim->gUpDim()/num_t + 1;
    unsigned long long size_part_d = _cim->gDownDim()/num_t + 1;
    unsigned long start_l =  me * size_part;
    unsigned long end_l = start_l + size_part;
    if( end_l > _cim->gUpDim())
        end_l = _cim->gUpDim();

    if (CIMethod_debug)
        cout << "Started construct_CI_one_dens in DensFCI" << endl;

    Permutator_Bit perm(dim() );
    TYPE up1 = perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        up1 = perm.permutate_bit(up1);
    //Contribution of the ups, with the downs equal.
    for(unsigned long i = start_l ; i < end_l ; i ++){
        TYPE up2 = up1;
       	for(unsigned long j = i ; j < _cim->gUpDim() ; j ++){
            TYPE down1 = perm.get_start_int(_cim->gNdown());
       	    for(unsigned long k = 0 ; k < _cim->gDownDim() ; k++){
	        double contr = _cim->get_eigvec(i*_cim->gDownDim()+k,0) * _cim->get_eigvec(_cim->gDownDim() * j + k,0);
       	        fill_density_matrix_one(up1,down1, up2, down1,contr, onerdm);
       	        down1 = perm.permutate_bit(down1);
	    }
       	    up2 = perm.permutate_bit(up2);
        } // end cols 
       	up1 = perm.permutate_bit(up1);
    } // end rows 

    start_l =  me * size_part_d;
    end_l = start_l + size_part_d;
    if( end_l > _cim->gDownDim()-1)
        end_l = _cim->gDownDim()-1;

    TYPE down1 = perm.get_start_int(_cim->gNdown());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        down1 = perm.permutate_bit(down1);
    //Now contribution of the downs, with the ups equal, but with no diagonal elements this time because they are already calculated in the loop above.
    for(unsigned long i = start_l ; i < end_l ; i ++){
	    TYPE down2 = perm.permutate_bit(down1);
       	for(unsigned long j = i+1 ; j < _cim->gDownDim() ; j ++){
            TYPE up1 = perm.get_start_int(_cim->gNup());
       	    for(unsigned long k = 0 ; k < _cim->gUpDim() ; k++){
	        double contr = _cim->get_eigvec(k*_cim->gDownDim()+i,0) * _cim->get_eigvec(_cim->gDownDim() * k + j,0);
       	        fill_density_matrix_one(up1,down1, up1, down2,contr, onerdm);
       	        up1 = perm.permutate_bit(up1);
	    }
       	    down2 = perm.permutate_bit(down2);
        } // end cols 
       	down1 = perm.permutate_bit(down1);
    } // end rows 
    if (CIMethod_debug)
        cout << "We constructed the density matrix" << endl;
}

void DensFCI::construct_CI_two_dens(unsigned int start , unsigned int end, vector<Vector2d > & twordm ){
    if (CIMethod_debug)
        cout << "Started construct_CI_one_dens in DensFCI" << endl;

    #ifdef __OMP__
    auto num_t = omp_get_max_threads();
    auto me = omp_get_thread_num();
    #else 
    auto num_t = 1;
    auto me = 0;
    #endif
    unsigned long long size_part = _cim->gUpDim()/num_t + 1;
    unsigned long start_l =  me * size_part;
    unsigned long end_l = start_l + size_part;
    if( end_l > _cim->gUpDim()-1)
        end_l = _cim->gUpDim()-1;

    Permutator_Bit perm(dim() );
    TYPE up1 = perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        up1 = perm.permutate_bit(up1);

    //Handles everything except when the ups are equal.
    for(unsigned long i = start_l ; i < end_l ; i ++){
	    TYPE up2 = perm.permutate_bit(up1);
       	for(unsigned long j = i+1 ; j < _cim->gUpDim() ; j ++){
            TYPE down1 = perm.get_start_int(_cim->gNdown());
       	    for(unsigned long k = 0 ; k < _cim->gDownDim() ; k++){
                TYPE down2 = perm.get_start_int(_cim->gNdown());
                for(unsigned long l = 0 ; l < _cim->gDownDim() ; l++){
                    double contr = _cim->get_eigvec(i*_cim->gDownDim()+k,0) * _cim->get_eigvec(_cim->gDownDim() * j + l,0);
                    fill_density_matrix(up1,down1, up2, down2,contr,twordm);
                    down2 = perm.permutate_bit(down2);
                }    
                down1 = perm.permutate_bit(down1);
	        }
       	    up2 = perm.permutate_bit(up2);
        } // end cols 
       	up1 = perm.permutate_bit(up1);
    } // end rows 

    up1 = perm.get_start_int(_cim->gNup());
    for(auto idx_begin=0;idx_begin< start_l;++idx_begin)
        up1 = perm.permutate_bit(up1);

//This handles the case when the ups are equal.
    for(unsigned long i = start_l ; i < end_l ; i ++){
        TYPE down1 = perm.get_start_int(_cim->gNup());
        for(unsigned long k = 0 ; k < _cim->gDownDim() ; k++){
            TYPE down2 = down1;
            for(unsigned long l = k; l < _cim->gDownDim() ; l++){
                double contr = _cim->get_eigvec(i*_cim->gDownDim()+k,0) * _cim->get_eigvec(_cim->gDownDim() * i + l,0);
                fill_density_matrix(up1,down1, up1, down2,contr,twordm);
                down2 = perm.permutate_bit(down2);
            }    
            down1 = perm.permutate_bit(down1);
        }
        up1 = perm.permutate_bit(up1);
    } 

    if (CIMethod_debug)
        cout << "We constructed the density matrix" << endl;
}

