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
#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <utility>

#include <../include/MyHDF5.h>
// This allows us to be lazy in getting the spaces in DPD calls
#define ID(x) ints.DPD_ID(x)

INIT_PLUGIN
// macro to help check return status of HDF5 functions
#define HDF5_STATUS_CHECK(status) { \
    if(status < 0) \
    outfile->Printf("%s:%d: Error with HDF5. status code=%d\n", __FILE__, __LINE__, status); \
} 

namespace psi{ namespace mointegrals{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "MOINTEGRALS"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        // If true: saves the unitary transformation from SO to MO integrals in a textfile 
	// with name U_MO_FILENAMEunit_ao_to_mo.txt, and from the symmetric orthogonal atomic 
	// orbitals to MO in hdf5 with U_MO_FILENAMEorthon.h5
        options.add_bool("SAVE_UNITARIES", false); 
        options.add_str_i("U_MO_FILENAME", "unitary-mo"); //base filename for the unitaries.
        options.add_str_i("FILENAME", "psioutput.dat"); //base filename for the unitaries.
        // save to a HDF5 file; implement for big files see(sointegrals.)
        //options.add_bool("SAVEHDF5", false);
        //options.add_str_i("HDF5_FILENAME", "integrals.h5");
    }

    return true;
}


extern "C" PsiReturnType
mointegrals(Options &options)
{
    /*
     * This plugin shows a simple way of obtaining MO basis integrals, directly from a DPD buffer.  It is also
     * possible to generate integrals with labels (IWL) formatted files, but that's not shown here.
     */
    int print = options.get_int("PRINT");
    bool save_unitaries = options.get_bool("SAVE_UNITARIES"); 
    //bool savehdf5 = options.get_bool("SAVEHDF5");
    std::string filename = options.get_str("FILENAME");
    std::string mofilename = options.get_str("U_MO_FILENAME");
    to_lower(filename);
    to_lower(mofilename);
    FILE * integralfile;
    integralfile=fopen(filename.c_str(), "w");
   
    // Grab the global (default) PSIO object, for file I/O
    boost::shared_ptr<PSIO> psio(_default_psio_lib_);

    // Now we want the reference (SCF) wavefunction
    boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();
    if(!wfn) throw PSIEXCEPTION("SCF has not been run yet!");
    
    /*MoldenWriter mollie(wfn);
    mollie.write("moldenwriteroutput.txt");*/ //Please call the MoldenWriter from the psi4 input file from now on.
    
    // Quickly check that there are no open shell orbitals here...
    int nirrep  = wfn->nirrep();
    
    // For now, we'll just transform for closed shells and generate all integrals.  For more elaborate use of the
    // LibTrans object, check out the plugin_mp2 example.
    std::vector<boost::shared_ptr<MOSpace> > spaces;
    spaces.push_back(MOSpace::all);
    IntegralTransform ints(wfn, spaces, IntegralTransform::Restricted);
    ints.transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
    // Use the IntegralTransform object's DPD instance, for convenience
    dpd_set_default(ints.get_dpd_id());
    
    //Readin the MO OEI in moOei & print everything
    int nalpha    = Process::environment.wavefunction()->nalpha();
    int nbeta     = Process::environment.wavefunction()->nbeta();
    int nmo       = Process::environment.wavefunction()->nmo();
    int nIrreps   = Process::environment.wavefunction()->nirrep();
    int *orbspi   = Process::environment.wavefunction()->nmopi();
    int *docc     = Process::environment.wavefunction()->doccpi();
    int *socc     = Process::environment.wavefunction()->soccpi();
    
    int nTriMo = nmo * (nmo + 1) / 2;
    double *temp = new double[nTriMo];
    Matrix moOei("MO OEI", nIrreps, orbspi, orbspi);
    IWL::read_one(psio.get(), PSIF_OEI, PSIF_MO_OEI, temp, nTriMo, 0, 0, filename);
    moOei.set(temp);
    
    moOei.print();
    

    fprintf(integralfile, "****  Molecular Integrals For DOCI Start Here \n");
    //get permutation array to put integral indices in standard order so HF determinant corresponds to 0..01..1|0..01..1
    /*
    int ii = 0;
    std::vector< std::pair<double,int> > V;
    for(std::vector<double>::iterator it = wfn->epsilon_a()->begin() ;it != wfn->epsilon_a()->end()  ;it++){
	    std::pair<double,int> P = std::make_pair(*it,ii);
	    V.push_back(P);
	    ii += 1;
    }
    std::sort(V.begin(),V.end());
    int permute_index[V.size()];
    for(int l = 0 ; l < V.size() ; l++){
	    permute_index[V[l].second] = l;
    }
    */ //add permute_index[] to the OEI, and TEI indices before writing to file and to irreporder
    /*  
    fprintf(integralfile,"MOE alpha = " );
    int i = 0;
    for(std::vector<double>::iterator it = wfn->epsilon_a()->begin() ;it != wfn->epsilon_a()->end()  ;it++){
	    fprintf(integralfile,"%lf ", *it);
	    fprintf(integralfile," %d %d %d",i,permute_index[i] , V[i].second );
	    i++;
    }
    fprintf(integralfile,"\n");*/
    std::string SymmLabel =  Process::environment.molecule()->sym_label();
    fprintf(integralfile, "Nalpha = %1d \n", nalpha);
    fprintf(integralfile, "Nbeta = %1d \n", nbeta);
    fprintf(integralfile, "Symmetry Label = ");
    fprintf(integralfile, SymmLabel.c_str());
    fprintf(integralfile, " \n");
    fprintf(integralfile, "Nirreps = %1d \n", nirrep);
    double NuclRepulsion =  Process::environment.molecule()->nuclear_repulsion_energy();
    fprintf(integralfile, "Nuclear Repulsion Energy = %16.48f \n", NuclRepulsion);
    fprintf(integralfile, "Number Of Molecular Orbitals = %2d \n", nmo);
    fprintf(integralfile, "Irreps Of Molecular Orbitals = \n");
     /*
    std::vector<int> irrep_order(V.size(),0);
    int aa = 0; 
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
	  irrep_order[permutate_index[aa]] = h;
	  aa +=1 ;
       }
    }
    for(int l = 0 ; l < V.size() ; l++){
          fprintf(integralfile, "%1d ",irrep_order[l]);
    }
    */
    for (int h=0; h<nirrep; ++h){
       for (int cnt=0; cnt<moOei.rowspi(h); ++cnt){
          fprintf(integralfile, "%1d ",h);
       }
    }
    fprintf(integralfile, "\n");
    fprintf(integralfile, "DOCC = ");
    for (int h=0; h<nirrep; h++){
       fprintf(integralfile, "%2d ",docc[h]);
    }
    fprintf(integralfile, "\n");
    fprintf(integralfile, "SOCC = ");
    for (int h=0; h<nirrep; h++){
       fprintf(integralfile, "%2d ",socc[h]);
    }
    fprintf(integralfile, "\n");
    
    double EnergyHF = 0.0;
    EnergyHF += NuclRepulsion;

    /*REMARK we adapted sointegrals.cc to save the Hamiltonian in symmetrically orthogonalized so orbitals so if we want the unitary from that basis to the MO:
     * We have to create a unitary that first inverts the transformation from the symmetrically orthogonalized orbitals to the so orbitals and then multiply in the Ca coefficient matrix to get the MO.
     * We resort to this because the Hamiltonian objects can only save orthogonal orbitals. 
     */
     // contains the transformation from Symmetrized atomic Orbitals -> MO (THIS IS NOT FROM Symmetric orthogonalized orbitals to MO)
     boost::shared_ptr<MatrixFactory> factory = Process::environment.wavefunction()->matrix_factory();
     const Dimension &dimension = moOei.rowspi();
     SharedMatrix ca = factory->create_shared_matrix("Transformationao->mo: ");
     ca = Process::environment.wavefunction()->Ca();
     fprintf(integralfile, "CIFlowTransformation: \n");
     //saves transformation from symmetric orthogonal orbitals to MO in HDF5 format.
     for (int irrep=0; irrep< nirrep; irrep++)
     {
         int norb = dimension[irrep];
         std::stringstream irrepname;
         irrepname << "irrep_" << irrep;

         fprintf(integralfile ,"%s\n" , irrepname.str().c_str() );
         if(norb > 0)
         {
             for(int ll = 0 ; ll < norb ; ll++)
             {
                 for(int aa = 0 ; aa < norb ; aa++)
                 {
     		         fprintf(integralfile , "%.15f    " , ca->get(irrep,  aa, ll )) ; //We transpose because psi4 saves in columns and unitarymatrix saves in rows.
                 }
     		         fprintf(integralfile , "\n" ) ; //We transpose because psi4 saves in columns and unitarymatrix saves in rows.
             }
         }
     }
     ca->print();
     
     SharedMatrix S= factory->create_shared_matrix("OverlapMatrixCIFlow: ");
     S = Process::environment.wavefunction()->S();
     fprintf(integralfile, "CIFlowOverlap: \n");
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
     		         fprintf(integralfile , "%1d %1d %16.48f \n",nTot+ll , nTot+aa  , S->get(irrep,  ll, aa )) ; 
                 }
             }
             nTot += norb;
         }
     }
     S->print();
    if(save_unitaries)
    {
        // Construct Shalf
        SharedMatrix eigvec= factory->create_shared_matrix("L");
        SharedMatrix temp= factory->create_shared_matrix("Temp");
        SharedMatrix temp2= factory->create_shared_matrix("Temp2");
        SharedVector eigval(factory->create_vector());

        S->diagonalize(eigvec, eigval);

        // Convert the eigenvales to sqrt(eigenvalues)
        int *dimpi = eigval->dimpi();
        double min_S = fabs(eigval->get(0,0));
        for (int h=0; h<nIrreps; ++h)
            for (int i=0; i<dimpi[h]; ++i)
            {
                if (min_S > eigval->get(h,i))
                    min_S = eigval->get(h,i);
                //double scale = 1.0 / sqrt(eigval->get(h, i));
                double scale = sqrt(eigval->get(h, i));
                eigval->set(h, i, scale);
            }

        fprintf(integralfile, "Lowest eigenvalue of overlap S = %14.10E\n", min_S);

        if(min_S < options.get_double("S_TOLERANCE") )
            fprintf(integralfile, "WARNING: Min value of overlap below treshold!!!!\n");


        // Create a vector matrix from the converted eigenvalues
        temp2->set_diagonal(eigval);

        temp->gemm(false, true, 1.0, temp2, eigvec, 0.0);
        S->gemm( false , false, 1.0, eigvec, temp, 0.0);

        // order of multi. is changed to have the output compatible to UnitaryMatrix class.
        // (columns/rows interchanged)
	// REMARK: in unitarymatrix world we need Ca^+ S^+ phi (because unitary matrix saves in rows and Ca and S saves in columns.)
	//if we do S Ca (and psi4 saves row major and unitary matrix column major)
	//we can just copy the memory and read in with pointers. because the transposition is automatically, because of the different memory layout.
        temp->gemm(false, false, 1.0, S, ca, 0.0);

        fprintf(integralfile, "temp:\n");

        std::ofstream ofs (mofilename+ "unit_ao_to_mo.txt", std::ofstream::out);
	mofilename += "orthon.h5";
        hid_t file_id = H5Fcreate(mofilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        hid_t group_id = H5Gcreate(file_id, "/Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        const Dimension &ca_dims = ca->rowspi();


	    //saves transformation from symmetric orthogonal orbitals to MO in HDF5 format.
        for (int irrep=0; irrep< nirrep; irrep++)
        {
            int norb = ca_dims[irrep];

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
		        ofs << ca->get(irrep,  aa, ll ) << "    "; //We transpose because psi4 saves in columns and unitarymatrix saves in rows.
		    }
		    ofs << std::endl;
		}

                double *p_data = temp->pointer(irrep)[0];

                hsize_t dimarray      = norb * norb;
                hid_t dataspace_id    = H5Screate_simple(1, &dimarray, NULL);
                hid_t dataset_id      = H5Dcreate(group_id, irrepname.str().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, p_data);

                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);
            }
        }
        H5Gclose(group_id);
        H5Fclose(file_id);
        ofs.close();
    } //End save_unitaries.
    
    fprintf(integralfile, "****  MO OEI \n");
    
    nTot = 0;
    for(int h = 0; h < nirrep; ++h){
       for(int cnt = 0; cnt < moOei.rowspi(h); ++cnt){
          for (int cnt2 = 0; cnt2 < moOei.colspi(h); ++cnt2){
             fprintf(integralfile, "%1d %1d %16.48f \n",
                               nTot+cnt,nTot+cnt2, moOei[h][cnt][cnt2]); //index is added to have the right order
          }
          if (cnt <docc[h]) EnergyHF += 2*moOei[h][cnt][cnt];
          else{
             if (cnt <socc[h]+docc[h]) EnergyHF += moOei[h][cnt][cnt];
          }
       }
       nTot += moOei.rowspi(h);
    }

    fprintf(integralfile, "****  MO TEI \n");

    /*
     * Now, loop over the DPD buffer, printing the integrals
     */
    dpdbuf4 K;
    psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    // To only process the permutationally unique integrals, change the ID("[A,A]") to ID("[A>=A]+")
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[A,A]"), ID("[A,A]"),
                  ID("[A>=A]+"), ID("[A>=A]+"), 0, "MO Ints (AA|AA)");
    for(int h = 0; h < nirrep; ++h){
        global_dpd_->buf4_mat_irrep_init(&K, h);
        global_dpd_->buf4_mat_irrep_rd(&K, h);
        for(int pq = 0; pq < K.params->rowtot[h]; ++pq){
            int p = K.params->roworb[h][pq][0];
            int q = K.params->roworb[h][pq][1];
            int psym = K.params->psym[p];
            int qsym = K.params->qsym[q];
            int prel = p - K.params->poff[psym];
            int qrel = q - K.params->qoff[qsym];
            for(int rs = 0; rs < K.params->coltot[h]; ++rs){
                int r = K.params->colorb[h][rs][0];
                int s = K.params->colorb[h][rs][1];
                int rsym = K.params->rsym[r];
                int ssym = K.params->ssym[s];
                int rrel = r - K.params->roff[rsym];
                int srel = s - K.params->soff[ssym];
                // Print out the absolute orbital numbers, the relative (within irrep)
                // numbers, the symmetries, and the integral itself
                fprintf(integralfile, "%1d %1d %1d %1d %16.48f \n",
                                 p,q,r ,  s,  K.matrix[h][pq][rs]);
                if ((p==q) && (r==s)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF += 2*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF += 0.5*K.matrix[h][pq][rs];
                }
                if ((p==s) && (r==q)){
                   if ((prel <docc[psym]) && (rrel < docc[rsym])) EnergyHF -= K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel < docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel <docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                   if ((prel>=docc[psym]) && (prel < socc[psym] + docc[psym]) && (rrel>= docc[rsym]) && (rrel < socc[rsym] + docc[rsym])) EnergyHF -= 0.5*K.matrix[h][pq][rs];
                }
            }
        }
        global_dpd_->buf4_mat_irrep_close(&K, h);
    }
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
    fprintf(integralfile, "****  HF Energy = %16.48f \n", EnergyHF);
    fprintf(integralfile, "****  Molecular Integrals For DOCI End Here \n");


    fclose(integralfile);
    return Success;

}

}} // End Namespaces
