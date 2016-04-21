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

#include <stdlib.h> /* srand, rand*/
#include <iostream>
#include <fstream>
#include <time.h> /* time*/
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <algorithm>
#include <chrono>

#include "Hamiltonian.h"
#include "CIMethod.h"
#include "CIDens.h"
#include "Options.h"
#include "Output.h"
#include "SimulatedAnnealing.h"
#include "UnitaryMatrix.h"
#include "OrbitalTransform.h"
#include "LocalMinimizer.h"
#include "HamConstruct.h"
#include "scpp_assert.h"
#include "Properties.h"

int num = 0;
using namespace std;

double e_in_ham(Hamiltonian * ham, CIMethod *cim){
    cim->set_ham(ham ,true); //ham of cim is optimal basis of min_sen.
    cim->solve();
    std::cout  << "energy in new basis: " << cim->get_ci_energy() << std::endl;
    return cim->get_ci_energy();
}

double get_no_energy(CIMethod * cim , const CIDens & cid, bool print_no)
{
    //This method calculates the energy of the CIMethod in the natural orbital basis defined by the density matrices in cid. REMARK: after this function returns cim contains the transformed Hamiltonian.
    Hamiltonian * ham = new Hamiltonian(*cim->get_ham()); //Save old ham to set back.
    OrbitalTransform orbtrans = OrbitalTransform(cim->get_ham());
    matrix NO = matrix(cim->get_l() , cim->get_l());
    matrix occupations = matrix(cim->get_l() ,1);
    cid.get_no( occupations, NO);
    if(CIMethod_debug){
        NO.Print();
        occupations.Print();
    }	    
    OptIndex opt = ham->get_index_object();
    UnitaryMatrix unit(NO, opt);
    ofstream file("unitary_naturalorbs" + cid.get_cim()->get_name() + cim->get_ham()->get_short_filename()+ ".dat" , std::ofstream::out);
    unit.print_unitary(file);
    file.close();
    //unit.print_unitary();
    //unit.CheckDeviationFromUnitary();
    //orbtrans->set_unitary(unit);
    if(print_no)
    {
        std::ofstream file;
	    file.open("naturaloccupations.dat", fstream::app);
	    file << ham->get_short_filename() << "    " ;
	    for(int i = 0 ; i < ham->getL() ; i++)
	    {
	        file << occupations[i] << "    ";
	    }
	    file << std::endl;
	    file.close();
    }
    orbtrans.rotate_active_space_vectors(NO.getpointer());
    orbtrans.fillHamCI(*ham);
    cim->set_ham(ham);
    cim->solve();
    std::cout  << "energy in natural orbitals: " << cim->get_ci_energy() << std::endl;
    return cim->get_ci_energy();
}

Hamiltonian * get_mmin_ham( CIMethod* cimmin, bool diis)
{
    Iterative_Subotnik * min_sen;
    diis = 0;
    if(diis)
        min_sen = new Iterative_Subotnik_DIIS(cimmin  , 1e-7, 0);
    else
        min_sen = new Iterative_Subotnik( cimmin, 1e-7, 0);
    min_sen->optimize();
    min_sen->save_unitary("mmin");
    Hamiltonian * ham = min_sen->get_ham();
    delete min_sen; 
    return ham;
}


double get_fci_mmin_energy(CIMethod * cim , bool diis)
{
    FCI fci{cim->get_ham() };
    //fci.set_sparsity(1e-14);
    Hamiltonian * mminham = get_mmin_ham(&fci, diis);
    double fcimmine = e_in_ham(mminham , cim );
    std::cout  << "E in fci mmin basis is : " << fcimmine << std::endl;
    return fcimmine;
}

double get_cisd_mmin_energy(CIMethod * cim , bool diis)
{
    std::cout  << "Make sure the cisd determinants for this problem are writen in the file cisddeterminants.dat : " << std::endl;
    FCI_File cisd{cim->get_ham(), "cisddeterminants.dat"};
    Hamiltonian * mminham = get_mmin_ham(&cisd, diis);
    double cisdmmine = e_in_ham(mminham , cim );
    std::cout  << "E in cisd mmin basis is : " << cisdmmine << std::endl;
    return cisdmmine;
}

double get_mmin_doci_energy(CIMethod * cim , bool diis)
{
    DOCI doci { cim->get_ham() } ;
    Hamiltonian * mminham = get_mmin_ham(cim, diis);
    double docimmine = e_in_ham(mminham , &doci );
    std::cout  << "doci mmin E is : " << docimmine << std::endl;
    doci.reset_output("dociinmmin" + cim->get_perm_info() );
    doci.print_output();
    doci.print_rdm();
    return docimmine;
}

double get_fci_no_energy(CIMethod * cim)
{
    //This method calculates the energy of the CIMethod cim in the natural orbital basis of FCI.
    std::cout << "We are going to calculate the CIMethod: " << " in the natural orbital basis of the FCI solution, we start by generating the FCI solution." << std::endl;
    FCI fci {cim->get_ham() };
    fci.solve();
    //std::cout << "FCI total energy: " << fci->get_ci_energy() << std::endl;
    std::unique_ptr<CIDens> densmatfci { new DensFCI(&fci)} ;
    densmatfci->construct_density(false);
    double cim_fcino_energy = get_no_energy(cim , *densmatfci, 1);
    return cim_fcino_energy;
}

void orb_opt(string sort , CIMethod * cim, bool printoutput, bool hamsave)
{
    Hamiltonian * ham = cim->get_ham();
    if(sort == "sim")
    {
        SimulatedAnnealing orbopt {cim  , 0.4 , 0.99 , 1.4 , 0.999};
        orbopt.run_multiple(5, true); //put true when you want to save the optimal unitary transformation.
        std::cout  << "Optimized energy with simulated annealing: " << orbopt.get_opt_ci_energy() << std::endl;
    }
    else if(sort == "local")
    {
        SCPP_TEST_ASSERT(cim->get_name() == "DOCI" , "Error: Local optimizer only works for DOCI at the moment, the cim type is: " << cim->get_name() );
        //handles mcscf orbital optimization.
        LocalMinimizer locmin {cim};
        double e_min = locmin.optimize();
        std::cout  << "Optimized energy with Local Minimizer for DOCI: " << e_min << std::endl;
    }
    else if(sort == "line")
    {
        LineSearch linsearch(cim);
        linsearch.optimize();
    }    
    else if(sort == "fmmin")
    {
       	get_fci_mmin_energy(cim, 1);//Standard we put diis for iterative subotnik on.
    }
    else if (sort == "sdmmin")
    {
       	get_cisd_mmin_energy(cim,1);
    }
    else if(sort == "mmind")
    {
        get_mmin_doci_energy(cim ,1);
    }
    else if(sort == "hmmin")
    {
        Hamiltonian * mminham = get_mmin_ham( cim, 1);
        double mmine = e_in_ham(mminham , cim );
        //std::cout  << "energy in own mmin basis is : " << mmine << std::endl;
    }
    else if(sort == "loadham")
    {
    	//string hamname = "matrixelements_otherbasis/ham" + cim->get_ham()->get_short_filename()+"FCI_Filecisddeterminantsmmin.h5" ; 
    	//string hamname = "hamiltonians/ham" + cim->get_ham()->get_short_filename()+"FCI_Filecisddeterminantsmmin.h5" ; 
    	string hamname = "ham" + cim->get_ham()->get_short_filename() + "FCI_Filecisddeterminantsmmin.h5";
    	//string hamname = "matrixelements_otherbasis/ham" + cim->get_ham()->get_short_filename() + "FCImmin.h5";
	    cout << hamname << endl;
	    cim->load_ham(hamname);
        cim->solve();
        std::cout  << "energy with loaded Hamiltonian: " << cim->get_ci_energy() << std::endl;
    }
    else if(sort == "unit")
    {
        string unitname = "unitarymatrix.dat";
        OptIndex opt = cim->get_ham()->get_index_object();
        //UnitaryMatrix unit(opt , unitname);
        UnitaryMatrix unit(opt );
        //unit.set_random_unitary({{0,1} , {1,2} , {0,2}} );
        unit.set_random_unitary();
        OrbitalTransform orbtrans = OrbitalTransform(cim->get_ham());
        Hamiltonian * ham = new Hamiltonian(*cim->get_ham()); //Save old ham to set back.
        orbtrans.set_unitary(unit);
        orbtrans.fillHamCI(*ham);
        cim->set_ham(ham);
        cim->solve();
        std::cout  << "energy with random transformation of hamiltonian: " << cim->get_ci_energy() << std::endl;
    }
    else if(sort == "fno")
    {
        get_fci_no_energy(cim);
    }
    else if(sort == "no")
    {
        std::unique_ptr<CIDens> densmat = cim->get_density();
        densmat->construct_density(false); //False for the construction and allocation of two rdm, because we only need 1rdm.
        get_no_energy(cim , *densmat ,0);
    }
    if( hamsave)
    {
        string hamname =  cim->get_name() +sort + to_string(num);
        cim->get_ham()->save_file(hamname);
        num += 1;
    }
    if(printoutput)
    {
        cim->reset_output(sort);
        cim->print_output({ "shannon" , "spin_squared" , "mulliken", "seniority"}, 0 ,false);
        cim->print_rdm(0, true);
        //cim->print_ham();
        //DensDOCI densmat {cim};
        //densmat.construct_density(false);
        //std::vector<int> orbs {0,1,2,3,4};
        //std::cout << "Mulliken charges first atom: " << densmat.get_mulliken(orbs) << endl;//Make sure the Hamiltonian contains the overlap of the ao, and the transformation from ao to current matrixelements.
    }
    cim->delete_ham();
    cim->set_ham(ham); //Set CIMethod back to the original hamiltonian.
}


int main ( int argc, char ** argv){
    // General settings and declarations.
    cout.precision(15);
    string method;
    string matrixelements;
    vector<pair< vector<string> , vector<string> > > methods(0);

    //Provide matrixelements.
    struct stat stFileInfo;
    int intStat = stat(matrixelements.c_str(),&stFileInfo);
    while(intStat != 0){
        cout << "Give the path to the matrix elements:\t";
        cin >> matrixelements;
        intStat = stat(matrixelements.c_str(),&stFileInfo);
    }	

    //Provide methods.
    while(1){
        pair<vector<string> , vector<string> > calculation {{} , {} };
        cout << "Give the Hilbertspaces for which you want to calculate the CI energy : " << endl;
        cout << "Possible options are: doci, fci, file, big" << endl;
        cout << "To finish providing methods type : end" <<endl;
        cin >> method;
        transform(method.begin(), method.end(), method.begin(), ::tolower);
        //start method providing.
        if (method == "doci" || method == "fci" || method == "file" || method == "big"){
            calculation.first.push_back(method);
        }
        else if(method == "end"){
            break;
        }
        else{
            continue; //Ask again for a valid method.
        }
        if(method == "file" || method == "big"){
            intStat = stat(method.c_str(),&stFileInfo); //intStat is already defined for the matrixelements. (To check if the file exists)
            while(intStat != 0){
                cout << "Which file contains the determinants?"<< endl;
                cin >> method;
                intStat = stat(method.c_str(),&stFileInfo);
            }
            calculation.first.push_back(method);
        }
        cout << "Orbital Optimization? " << endl;
        cout << "Wanna do orbital optimization, or select different Hamiltonian?" <<endl;
        while(1){
            cout << "Possible options are: sim (simmulated annealing), local, fno, loadham, ..., stop with endm" << endl;
            cin >> method;
            transform(method.begin(), method.end(), method.begin(), ::tolower);
            cout << method;
            if (method == "sim" || method == "mcscf" || method == "line" || method == "fmmin" || method == "fno" || method == "mmind" || method == "hmmin" || method == "loadham" || method == "no" || method == "sdmmin"  || method == "local" || method == "unit" ){
                calculation.second.push_back(method);
            }//end if
            else if(method == "endm"){
                break;
            }
        }//end while orbital optimization.
        methods.push_back(calculation) ; 
    }//end method providing.

    //Create objects to do the work, and execute the work.
    HamConstruct hamconst;
    std::shared_ptr<Hamiltonian> ham = hamconst.generate_ham(matrixelements);
    cout << "The group was found to be " << Irreps::getGroupName(ham->getNGroup()) << endl;
    std::cout << "HF energy (lowest energy of a single determinant in the selection): " << ham->get_hf_energy() << std::endl;
    for (vector<pair<vector<string> , vector<string> > >::iterator cimethod  = methods.begin(); cimethod != methods.end() ; cimethod ++){
        if ((*cimethod).first[0] == "doci" || (*cimethod).first[0] == "fci"){
            if((*cimethod).first[0] == "doci"){
                DOCI * cim = new DOCI(ham.get());
                cim->solve();
                std::cout  << "DOCI energy: " << cim->get_ci_energy() << std::endl;
                //std::vector<int> orbs {0,1,2,3,4};
                //std::cout << "Mulliken charges first atom: " << cim->get_mulliken(orbs) << endl;//Make sure the Hamiltonian contains the overlap of the ao, and the transformation from ao to current matrixelements.
                //std::cout << "Spin Squared: " << cim->get_spin_squared() << std::endl;
                cim->print_output({ "shannon" , "spin_squared" , "mulliken"}, 0 , false);
                //cim->print_output();
                //cim->print_ham();
                cim->print_rdm(0,true);
                //auto begin = std::chrono::high_resolution_clock::now();
                //std::cout <<"start timing." <<std::endl;


                //DensDOCI densmat {cim};
                //densmat.construct_density(true);

                //auto end = std::chrono::high_resolution_clock::now();
                //std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() / std::pow(10,9 ) << "seconds" << std::endl;
                //densmat.print_one_dens(std::cout);
                //densmat.print_two_dens(std::cout);
                //std::cout << "DOCI-DENS one-electron energy: " << densmat.get_one_electron_energy() << std::endl;
                //std::cout << "DOCI-DENS two-electron energy: " << densmat.get_two_electron_energy() << std::endl;
                //std::cout << "DOCI-DENS total energy: " << densmat.get_dens_energy() << std::endl;

                //Take care of the orbital optimizations.
                if( (*cimethod).second.size() != 0)
                    for(auto sort : (*cimethod).second )
                        orb_opt(sort , cim ,  1,  1);
                        //orb_opt(sort, cim ,  0,  0);
                delete cim;
            }
            if((*cimethod).first[0] == "fci"){
                FCI *  cim5 = new FCI(ham.get());
                int neigval = 1;
                cim5->solve(neigval);
                //cim5->check_hermiticity();
                std::cout  << "FCI energy:" << cim5->get_ci_energy() << " Per site: " << cim5->get_ci_energy() / (double) ham->getL() <<std::endl;
                for(int i = 0 ; i < neigval ; i ++)
                {
                  //std::cout << "FCI energy: " << cim5->get_ci_energy(i) << endl;//Make sure the Hamiltonian contains the overlap of the ao, and the transformation from ao to current matrixelements.
                  std::vector<int> orbs {0,1,2,3,4};
                  std::cout << "Mulliken charges first atom: " << cim5->get_mulliken(orbs, i ) << endl;//Make sure the Hamiltonian contains the overlap of the ao, and the transformation from ao to current matrixelements.
                  //std::cout << "Spin Squared: " << cim5->get_spin_squared(i) << std::endl;
                }
                cim5->print_output({  "spin_squared", "mulliken", "seniority", "shannon" , "parenergy"}, 0 , false);
                //cim5->print_output({  "spin_squared"}, 0 ,false);
                //Properties prop { cim5->get_eigvec(0)  , cim5->get_l() , cim5->get_perm() };
                //cout << "shannon entropy "  << prop.shannon_ic() << std::endl;
                //cim5->print_output({ "shannon" , "mulliken", "spin_squared"}, 0 , false);
                //cim5->print_output();

                cim5->print_rdm(0 ,true); //0 -> state , 1-> 2rdm
                //cim5->get_ham()->print_overlap(std::cout);
                //cim5->get_ham()->get_unitary()->print_unitary(std::cout);
                //cim5->print_ham();

                //auto begin = std::chrono::high_resolution_clock::now();
                //std::cout <<"start timing." <<std::endl;
                //DensFCI densmatfci {cim5};
                //densmatfci.construct_density();

                //auto end = std::chrono::high_resolution_clock::now();
                //std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() / std::pow(10,9 ) << "seconds" << std::endl;

                //densmatfci.transform_to_ao(true, false); //options: 2rdm , revert
                //ofstream file(ham->get_short_filename() + "aodens");
                //densmatfci.print_one_dens(std::cout);
                //densmatfci.print_two_dens(file);
                //file.close();
                //std::cout << "FCI-DENS one-electron energy: " << densmatfci.get_one_electron_energy() << std::endl;
                //std::cout << "FCI-DENS two-electron energy: " << densmatfci.get_two_electron_energy() << std::endl;
                //std::cout << "FCI-DENS total energy: " << densmatfci.get_dens_energy() << std::endl;
                if( (*cimethod).second.size() != 0)
                    for(auto sort : (*cimethod).second )
                        orb_opt(sort , cim5 ,  1,  1);
                        //orb_opt(sort, cim5 ,  0,  0);
                delete cim5;
            }

        }// end doci or fci
        else if((*cimethod).first[0] == "file" || (*cimethod).first[0] == "big"){
            SCPP_ASSERT((*cimethod).first.size() == 2 , "The size of the first element of the pair that defines the calculation is not equal to 2 for a filemethod (method , detfile), it is: " << (*cimethod).first.size() << endl)
            string sort = (*cimethod).first[0];
            string detfile = (*cimethod).first[1];
            if(sort == "file"){
                FCI_File *  cim2 = new FCI_File(ham.get(), detfile);
                cim2->solve();
                std::cout  << "CIenergy " << sort << " of all the dets in " << detfile << " :  " << cim2->get_ci_energy() << std::endl;
                //std::cout << "Spin squared: " << cim2->get_spin_squared() << std::endl;
                //Properties prop { cim2->get_eigvec(0)  , cim2->get_l() , cim2->get_perm() };
                //cout << cim2->get_name() + cim2->get_ham()->get_short_filename() << ": shannon entropy: "  << prop.shannon_ic() << std::endl;
            	//cim2->print_output({"shannon", "spin_squared" ,"mulliken" });
            	cim2->print_output({"shannon", "spin_squared" , "seniority" });
                //DensFILE densmatfile {cim2};
                //densmatfile.construct_density();
                //std::vector<int> orbs {0,1,2,3,4};
                //std::cout << "Mulliken charges first atom: " << densmatfile.get_mulliken(orbs) << endl;//Make sure the Hamiltonian contains the overlap of the ao, and the transformation from ao to current matrixelements.
	            //cim2->print_ham();
	            //cim2->print_rdm();
		        //cim2->print_dets();

                //densmatfile.print_one_dens(std::cout);
                //std::cout << "FILE-DENS one-electron energy: " << densmatfile.get_one_electron_energy() << std::endl;
                //densmatfile.get_NO();
                //densmatfile.print_two_dens(std::cout);
                //densmatfile.compare_two_dens(densmat);
                //densmatfile.compare_one_dens(densmat);
                //std::cout << "FILE-DENS two-electron energy: " << densmatfile.get_two_electron_energy() << std::endl;
                //std::cout << "FILE-DENS total energy: " << densmatfile.get_dens_energy() << std::endl;
                if( (*cimethod).second.size() != 0)
                    for(auto sort : (*cimethod).second )
                        orb_opt(sort , cim2 ,  1,  1);
                        //orb_opt(sort, cim2 ,  0,  0);

                delete cim2;
            }    
            else{
                CI_Big *  cim3 = new CI_Big(ham.get(), (*cimethod).first[1]);
                cim3->solve();
                std::cout  << "CIenergy " << sort << " of all the dets in " << (*cimethod).first[1] << " :  " << cim3->get_ci_energy() << std::endl;
                //cim3->print_output();
                //cim3->print_rdm();
                //cim3->print_ham();
                if( (*cimethod).second.size() != 0)
                    for(auto sort : (*cimethod).second )
                        orb_opt(sort , cim3 ,  1,  1);
                        //orb_opt(sort, cim3 ,  0,  0);
                delete cim3;
            }
        } // End if cimethods that use permutator_file
    } //End for loop over methods.
    return 0;
}		/* ----------  end of function main  ---------- */
