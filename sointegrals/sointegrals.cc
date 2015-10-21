#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libmints/sointegral_twobody.h>
#include <libpsio/psio.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <memory>
#include <iomanip>
#include <fstream>

#include <hdf5.h>

#include "../include/Irreps.h"
#include "../include/Lapack.h"
#include "../include/Hamiltonian.h"

INIT_PLUGIN

// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) { \
    if(status < 0) \
    outfile->Printf("%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
} 

using namespace boost;

namespace psi{ namespace sointegrals{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "SOINTEGRALS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_bool("PRINT_INTEGRALS", true);
        /*- Whether to compute two-electron integrals -*/
        options.add_bool("DO_TEI", true);
        // save to a HDF5 file
        options.add_bool("SAVEHDF5", false);
        options.add_str_i("FILENAME", "atomicintegrals");
        options.add_bool("SAVEOVERLAP", true);
        options.add_str_i("OVERLAPNAME", "overlap.txt");
    }

    return true;
}

class ERIPrinter
{
public:

    ERIPrinter() { count = 0; }

    shared_ptr<Hamiltonian> Ham;

    unsigned int count;

    // Our functor...the nice thing about using a C++ functor is that the
    // code here is inlined by the compiler.
    void operator() (int pabs, int qabs, int rabs, int sabs,
                     int pirrep, int pso,
                     int qirrep, int qso,
                     int rirrep, int rso,
                     int sirrep, int sso,
                     double value)
    {
        outfile->Printf("%1d %1d %1d %1d %16.48f \n", 
			pabs, qabs, rabs, sabs, value);
        Ham->setVmat(pabs, rabs, qabs, sabs, value);

        count++;
    }
};

extern "C" PsiReturnType
sointegrals(Options &options)
{
    bool print = options.get_bool("PRINT_INTEGRALS");
    bool doTei = options.get_bool("DO_TEI"); //Transforms to symmetrical orthogonal atomic orbitals.
    bool savehdf5 = options.get_bool("SAVEHDF5");//Save in hdf5 ?
    bool saveoverlap = options.get_bool("SAVEOVERLAP");
    std::string filename = options.get_str("FILENAME");
    std::string overlapname = options.get_str("OVERLAPNAME");
    boost::algorithm::to_lower(filename);

    if(options.get_str("S_ORTHOGONALIZATION") != "SYMMETRIC")
    {
        outfile->Printf("This will only work with symmetric orthogonalisation: AO == OSO\n");
        return Failure;
    }

    shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Form basis object:
    // Create a basis set parser object.
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    // Construct a new basis set.
    shared_ptr<BasisSet> aoBasis = BasisSet::construct(parser, molecule, "BASIS");

    // The integral factory oversees the creation of integral objects
    shared_ptr<IntegralFactory> integral(new IntegralFactory
            (aoBasis, aoBasis, aoBasis, aoBasis));
    
    // N.B. This should be called after the basis has been built, because
    // the geometry has not been
    // fully initialized until this time.
    molecule->print();
    // The basis set is also created from the information stored in the
    // checkpoint file
    // ...it needs to raid the checkpoint file to find those dimensions.
    // Create an SOBasis object using the basis set and integral factory.
    shared_ptr<SOBasisSet> soBasis(new SOBasisSet(aoBasis, integral));

    // Obtain block dimensions from the SO basis
    const Dimension dimension = soBasis->dimension();

    // The matrix factory can create matrices of the correct dimensions...
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(dimension, dimension);

    //shared_ptr<OneBodyAOInt> aOBI(integral->ao_overlap());
    //SharedMatrix sMatao(factory->create_matrix("Overlapao"));
    //aOBI->compute(sMatao);
    //sMatao->print();
    // Form the one-electron integral objects from the integral factory
    shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());
    shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));
 
    const int nmo = dimension.sum();
    const int nirrep   = dimension.n();
    const double NuclRepulsion =  molecule->nuclear_repulsion_energy();

    int nelectrons = 0;
    for(int i=0;i<molecule->natom();i++)
        nelectrons += molecule->true_atomic_number(i);

    nelectrons -= molecule->molecular_charge();
    int nalpha = nelectrons/2;
    int nbeta = nelectrons/2;
    if(nelectrons % 2 != 0)
    {
	    nalpha +=1 ;
    }
    //Deduces the group, to create the Hamiltonian object.
    int SyGroup = 0;
    bool stopFindGN = false;
    std::string SymmLabel =  molecule->sym_label();
    do {
        if (SymmLabel.compare(Irreps::getGroupName(SyGroup))==0)
            stopFindGN = true;
        else
            SyGroup += 1;
    } 
    while ((!stopFindGN) && (SyGroup<42));

    outfile->Printf("If anything went wrong: Is %s equal to %s?\n", SymmLabel.c_str(), (Irreps::getGroupName(SyGroup)).c_str());
    
    std::vector<int> OrbIrreps;
    OrbIrreps.reserve(nmo);
    for (int h=0; h<nirrep; ++h)
        for (int i=0; i<dimension[h]; ++i)
            OrbIrreps.push_back(h);

    // Compute the one electron integrals, telling each object where to
    // store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    shared_ptr<Hamiltonian> Ham(new Hamiltonian(nmo, SyGroup, OrbIrreps.data(), nalpha , nbeta, NuclRepulsion));
    Ham->_overlap = std::vector<double>(Ham->getL() *Ham->getL() ,0.);
    if(saveoverlap)
    {
        std::ofstream ofs (overlapname, std::ofstream::out);
        //saves transformation from symmetric orthogonal orbitals to MO in HDF5 format.
        for (int irrep=0; irrep< nirrep; irrep++)
        {
            int norb = dimension[irrep];
        
            if(norb > 0)
            {
        	std::stringstream irrepname;
        	irrepname << "irrep_" << irrep;
        
        	ofs << irrepname.str() << std::endl;
        	ofs << std::setprecision(15);
        	for(int ll = 0 ; ll < norb ; ll++)
        	{
        	    for(int aa = 0 ; aa < norb ; aa++)
        	    {
        		ofs << sMat->get(irrep,  aa, ll ) << "    "; //We transpose because psi4 saves in columns and unitarymatrix saves in rows.
			sMat->print() ;
			Ham->set_overlap(irrep, aa, ll, sMat->get(irrep, ll , aa) );
        	    }
        	    ofs << std::endl;
       	        }
            }
	}
	ofs.close();
    }

    //Here we print the pure SO integrals which are not orthonormal, (but pointgroup symmetry adapted).
    outfile->Printf("****  Molecular Integrals For DOCI Start Here \n");
    outfile->Printf("Nalpha = %1d \n", nalpha);
    outfile->Printf("Nbeta = %1d \n", nbeta);
    outfile->Printf("Symmetry Label = %s \n", SymmLabel.c_str());
    outfile->Printf("Nirreps = %1d \n", nirrep);
    outfile->Printf("Nuclear Repulsion Energy = %16.48f \n", NuclRepulsion);
    /*  
    outfile->Printf("Dimension of Irreps = ");
    for (int h=0; h<nirrep; ++h)
        outfile->Printf("%2d ", dimension[h]);
    outfile->Printf("\n");*/
    outfile->Printf("Number Of SO Orbitals = %2d \n", nmo);
    outfile->Printf("Irreps Of SO Orbitals = \n");
    for (int h=0; h<nirrep; ++h)
        for (int i=0; i<dimension[h]; ++i)
            outfile->Printf("%2d ", h);
    outfile->Printf("\n");
    outfile->Printf("\n");//there is no DOCC line in SOintegrals.
    outfile->Printf("\n"); //There is no SOCC line in SO integrals.

    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);

     outfile->Printf("CIFlowOverlap: \n");
     int nTot = 0;
     for (int irrep=0; irrep< nirrep; irrep++)
     {
         int norb = dimension[irrep];
         if(norb > 0)
         {
             for(int ll = 0 ; ll < norb ; ll++)
             {
                 for(int aa = 0 ; aa < norb ; aa++)
                 {
     		         outfile->Printf("%1d %1d %16.48f \n",nTot+ll , nTot+aa  , sMat->get(irrep,  ll, aa )) ; 
                 }
             }
             nTot += norb;
         }
     }
    if(print)
    {
        //tMat->print();
        //vMat->print();
        //hMat->print();
    }
    outfile->Printf("*** OEI SO\n");

    int count = 0;
    for (int h=0; h<nirrep; ++h)
    {
        for (int i=0; i<dimension[h]; ++i)
            for (int j=i; j<dimension[h]; ++j)
            {
                outfile->Printf("%1d %1d %16.48f \n", count+i, count+j, hMat->get(h,i,j));
                Ham->setTmat(count+i, count+j, hMat->get(h,i,j));
            }

        count += dimension[h];
    }
    // 1. Obtain an object that knows how to compute two-electron AO
    // integrals.
    shared_ptr<TwoBodyAOInt> tb(integral->eri());
    
    // 2. Create an object that knows how to convert any two-body AO
    // integral to SO.
    shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral));
    
    // 3. Find out how many SO shells we have.
    int nsoshell = soBasis->nshell();

    // 4. We to create an instance of our ERIPrinter
    ERIPrinter printer;
    printer.Ham = Ham;

    outfile->Printf("*** TEI\n");
    
    // 5. Create an SOShellCombintationsIterator to step through the
    // necessary combinations
    SOShellCombinationsIterator shellIter(soBasis, soBasis, soBasis, soBasis);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        // 6. Call the TwoBodySOInt object to compute integrals giving
        // it the
        // instance to our functor.
        eri->compute_shell(shellIter, printer);
    }
    outfile->Printf("****  HF Energy = -10000000000000 \n"); 
    outfile->Printf("****  SO Integrals For DOCI End Here \n");
    outfile->Printf("number of TEI: %d\n", printer.count);
    if(savehdf5)
        Ham->save(filename);

    //This saves the symmetrically orthogonalized atomic orbital matrixelements.
    if(doTei){

    // Construct Shalf
    SharedMatrix eigvec= factory->create_shared_matrix("L");
    SharedMatrix temp= factory->create_shared_matrix("Temp");
    SharedMatrix temp2= factory->create_shared_matrix("Temp2");
    SharedVector eigval(factory->create_vector());

    sMat->diagonalize(eigvec, eigval);//vectors are in column format.

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    int *dimpi = eigval->dimpi();//dimension per irrep.
    double min_S = fabs(eigval->get(0,0));
    for (int h=0; h<nirrep; ++h)
        for (int i=0; i<dimpi[h]; ++i)
        {
            if (min_S > eigval->get(h,i))//eigvals are always positive but watch out for to small eigenvalues with point to linear dependencies.
                min_S = eigval->get(h,i);
            double scale = 1.0 / sqrt(eigval->get(h, i));
            eigval->set(h, i, scale);
        }

    outfile->Printf("Lowest eigenvalue of overlap S = %14.10E\n", min_S);

    if(min_S < options.get_double("S_TOLERANCE") )
        outfile->Printf("WARNING: Min value of overlap below treshold!!!!\n");


    // Create a vector matrix from the converted eigenvalues
    temp2->set_diagonal(eigval);

    temp->gemm(false, true, 1.0, temp2, eigvec, 0.0);
    sMat->gemm(false , false, 1.0, eigvec, temp, 0.0);

    temp->gemm(false, false, 1.0, hMat, sMat, 0.0);
    hMat->gemm(true , false, 1.0, sMat, temp, 0.0);

    //hMat->print();
    //sets the transformed h.
    int count = 0;
    for (int h=0; h<nirrep; ++h)
    {
        for (int i=0; i<dimension[h]; ++i)
            for (int j=i; j<dimension[h]; ++j)
            {
                outfile->Printf("%1d %1d %16.48f \n", count+i, count+j, hMat->get(h,i,j));
                Ham->setTmat(count+i, count+j, hMat->get(h,i,j));
            }

        count += dimension[h];
    }

    shared_ptr<Hamiltonian> Ham2(new Hamiltonian(*Ham));

    OptIndex index = Ham2->get_index_object();
    Irreps SymmInfo;
    SymmInfo.setGroup(Ham2->getNGroup());

    //Create the memory for the orbital transformations.
    unsigned long long maxlinsize = 0;
    for (int irrep=0; irrep< index.getNirreps(); irrep++)
    {
        unsigned int linsize_irrep = index.getNORB(irrep);
        if (linsize_irrep > maxlinsize)  
            maxlinsize  = linsize_irrep;
    }

    //Determine the blocksize for the 2-body transformation
    auto& maxBlockSize = maxlinsize;
    //Allocate 2-body rotation memory: One array is approx (maxBlockSize/273.0)^4 * 42 GiB --> [maxBlockSize=100 --> 750 MB]
    auto maxBSpower4 = maxBlockSize * maxBlockSize * maxBlockSize * maxBlockSize; //Note that 273**4 overfloats the 32 bit integer!!!
    auto sizeWorkmem1 = std::max( std::max( maxBSpower4 , 3*maxlinsize*maxlinsize ) , 1uLL * nmo * nmo ); //For (2-body tfo , updateUnitary, calcNOON)
    auto sizeWorkmem2 = std::max( std::max( maxBSpower4 , 2*maxlinsize*maxlinsize ) , nmo*(nmo + 1uLL) ); //For (2-body tfo, updateUnitary and rotate_to_active_space, rotate2DMand1DM)
    std::unique_ptr<double []> mem1(new double[sizeWorkmem1]);
    std::unique_ptr<double []> mem2(new double[sizeWorkmem2]);


    //Two-body terms --> use eightfold permutation symmetry in the irreps :-)
    for (int irrep1 = 0; irrep1<nirrep; irrep1++)
        for (int irrep2 = irrep1; irrep2<nirrep; irrep2++)
        {
            const int productSymm = SymmInfo.directProd(irrep1,irrep2);
            for (int irrep3 = irrep1; irrep3<nirrep; irrep3++)
            {
                const int irrep4 = SymmInfo.directProd(productSymm,irrep3);
                // Generated all possible combinations of allowed irreps
                if (irrep4>=irrep2)
                {
                    int linsize1 = index.getNORB(irrep1);
                    int linsize2 =index.getNORB(irrep2);
                    int linsize3 =index.getNORB(irrep3);
                    int linsize4 =index.getNORB(irrep4);

                    if ((linsize1>0) && (linsize2>0) && (linsize3>0) && (linsize4>0))
                    {
                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ]
                                            = Ham->getVmat(index.getNstart(irrep1) + cnt1,index.getNstart(irrep2) + cnt2, index.getNstart(irrep3) + cnt3, index.getNstart(irrep4) + cnt4 );

                        char trans = 'T';
                        char notra = 'N';
                        double alpha = 1.0;
                        double beta  = 0.0; //SET !!!

                        int rightdim = linsize2 * linsize3 * linsize4; //(ijkl) -> (ajkl)
                        double * Umx = sMat->pointer(irrep1)[0];
                            //_unitary->getBlock(irrep1);
                        dgemm_(&trans, &notra, &linsize1, &rightdim, &linsize1, &alpha, Umx, &linsize1, mem1.get(), &linsize1, &beta, mem2.get(), &linsize1);

                        int leftdim = linsize1 * linsize2 * linsize3; //(ajkl) -> (ajkd)
                        Umx = sMat->pointer(irrep4)[0];
                        dgemm_(&notra, &notra, &leftdim, &linsize4, &linsize4, &alpha, mem2.get(), &leftdim, Umx, &linsize4, &beta, mem1.get(), &leftdim);

                        int jump1 = linsize1 * linsize2 * linsize3; //(ajkd) -> (ajcd)
                        int jump2 = linsize1 * linsize2 * linsize3;
                        leftdim   = linsize1 * linsize2;
                        Umx = sMat->pointer(irrep3)[0];
                        for (int bla=0; bla<linsize4; bla++)
                            dgemm_(&notra, &notra, &leftdim, &linsize3, &linsize3, &alpha, mem1.get()+jump1*bla, &leftdim, Umx, &linsize3, &beta, mem2.get()+jump2*bla, &leftdim);

                        jump2    = linsize1 * linsize2;
                        jump1    = linsize1 * linsize2;
                        rightdim = linsize3 * linsize4;
                        Umx = sMat->pointer(irrep2)[0];
                        for (int bla=0; bla<rightdim; bla++)
                            dgemm_(&notra, &notra, &linsize1, &linsize2, &linsize2, &alpha, mem2.get()+jump2*bla, &linsize1, Umx, &linsize2, &beta, mem1.get()+jump1*bla, &linsize1);

                        for (int cnt1=0; cnt1<linsize1; cnt1++)
                            for (int cnt2=0; cnt2<linsize2; cnt2++)
                                for (int cnt3=0; cnt3<linsize3; cnt3++)
                                    for (int cnt4=0; cnt4<linsize4; cnt4++)
                                        Ham2->setVmat(index.getNstart(irrep1) + cnt1,index.getNstart(irrep2) + cnt2, index.getNstart(irrep3) + cnt3,index.getNstart(irrep4) + cnt4, mem1[cnt1 + linsize1 * ( cnt2 + linsize2 * (cnt3 + linsize3 * cnt4) ) ] );

                    } //end if the problem has orbitals from all 4 selected irreps
                } // end if irrep 4 >= irrep2
            }// end run irrep3
        } // end run irrep2

    //Set overlap
     int nTot = 0;
     for (int irrep=0; irrep< nirrep; irrep++)
     {
         int norb = dimension[irrep];
         if(norb > 0)
         {
             for(int ll = 0 ; ll < norb ; ll++)
             {
                 for(int aa = 0 ; aa < norb ; aa++)
                 {
     		         Ham2->set_overlap(irrep,ll , aa  , sMat->get(irrep,  ll, aa )) ; 
                 }
             }
             nTot += norb;
         }
     }


        int lastindex = filename.find_last_of("."); 
        string rawname = filename.substr(0, lastindex); 
        Ham2->save(rawname+"orthon");
        //Ham2->save_file(rawname+"orthon");

    }//end doTei transform to symmetric orthogonal.

    return Success;
    }

}} // End namespaces
