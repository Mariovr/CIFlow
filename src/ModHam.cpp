#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string> 
#include <cmath>

#include "Hamiltonian.h"
#include "ModHam.h"

using namespace std;

ModHam::ModHam(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options): Hamiltonian(L , 0 , nup , ndown, Econst)
{
    _params = params;
    _options = options;
}

ModHam::ModHam(const Hamiltonian & Ham, std::vector<double> params , std::vector<std::string> options  ): Hamiltonian(Ham.getL() , Ham.getNGroup() , Ham.get_orb2irrep() , Ham.getnup() , Ham.getndown(), 0.)
{
    _params = params;
    _options = options;
}

ModHam::ModHam(std::string filename):Hamiltonian(filename)
{
    std::string line;
    _params = {};
    _options = {};
    std::string start = "###Model Hamiltonian options:";
    std::ifstream inputfile(filename.c_str());   
    do
    {
        getline(inputfile , line);
    }while(line.find(start) == std::string::npos);
    //int pos = line.find(":");
    //std::stringstream options(line.substr(pos+1 , std::string::npos));
    //while(options)
    //{
//
    //}

    getline(inputfile,line);
    int pos;
    while( line.find("#") == std::string::npos )
    {
       _options.push_back(line);
       getline(inputfile,line);
    }
    while( getline(inputfile,line) ) 
    {
       _params.push_back(atof( line.substr(0, std::string::npos).c_str()) );
    }
    inputfile.close();
    _modham = true;
}

std::string ModHam::get_info() const
{
    std::string info;
    std::ostringstream infof;
    infof << "###Model Hamiltonian options:    " ;
    for (int i = 0 ; i < _options.size() ; i ++)
    {
        infof<< std::endl << _options[i] << "    " ; 
    }
    infof << std::endl ;
    infof << "###Model Hamiltonian parameters:     " ;
    for (int i = 0 ; i < _params.size() ; i ++)
    {
        infof<< std::endl << _params[i] << "    " ; 
    }
    infof << std::endl ;
    info += infof.str();
    return info;
}

void ModHam::save_file(const std::string & filename)
{
    Hamiltonian::save_file(filename);
    std::string filenames = "ham" +  get_short_filename() + filename + ".dat";
    std::ofstream outfile(filenames.c_str(), std::ios_base::app );
    std::string output = get_info();
    outfile << output <<std::endl;
    outfile.close();
}

//-----------------------------Subclasses-----------------------------------------------

Hub1d::Hub1d(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options): ModHam( L ,  nup,  ndown, Econst,  params ,  options)
{
    construct_ham(params , options);
}

Hub1d::Hub1d(std::string filename) : ModHam(filename)
{
    construct_ham(_params , _options);
}

Hub1d::Hub1d(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  ): ModHam(Ham, params , options  )
{
    construct_ham(params , options);
}



//Params: first t, then U (both positive, t is automatically set negative, info determines the boundary conditions)
void Hub1d::construct_ham(std::vector<double> params, std::vector<string> info)
{
    if (info[0] == "reset")
        this->set_zero();
    if(info[1] == "pos")
    {
        // one-particle integrals
        for(int i=0;i<(getL()-1);i++)
        {
           setTmat(i, i+1, -1.*params[0]);
        }

        if(info[2] == "per")
        {
            // periodic boundary condition
            setTmat(0, getL()-1, -1.*params[0]);
        }
        else if(info[2] == "aper")
        {
            //anti-periodic boundary condition
            setTmat(0, getL()-1, 1.*params[0]);
        }
        else if(info[2] == "none")
        {
            //No boundary conditions. (sites on a straigth line)
        }
        else
        {
            std::cerr << "This sort of conditions do not exist please provide: periodic -> 0 , anti-periodic -> 1 , none -> 2"  << std::endl;
            exit(1);
        }

        // two-particle integrals
        for(int i=0;i< getL();i++)
           setVmat(i, i, i, i, params[1]); 
    }
    else if(info[1] == "mom" )
    {
      // Don't forget: you cannot rotate states of different momentum!
      // one-particle integrals
      for(int i=0;i<getL();i++)
         setTmat(i, i, -2*params[0]*std::cos(2*M_PI/(1.0*getL())*i));

      // two-particle integrals
      for(int k1=0;k1<getL();k1++)
         for(int k2=0;k2<getL();k2++)
            for(int k3=0;k3<getL();k3++)
               for(int k4=0;k4<getL();k4++)
                  if((k1+k2) == (k3+k4))
                     setVmat(k1,k2,k3,k4, params[1]*1.0/getL());
    }
}

std::string Hub1d::get_info() const
{
    std::string info = Hamiltonian::get_info();
    info += "###Model Hamiltonian Type: Hub1d\n";
    info += ModHam::get_info();
    return info;
}

RedBCS::RedBCS(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options): ModHam( L ,  nup,  ndown, Econst,  params ,  options)
{
    construct_ham(params , options);
    _modham = true;
}

RedBCS::RedBCS(std::string filename) : ModHam(filename)
{
    construct_ham(_params , _options);
    _modham = true;
}

RedBCS::RedBCS(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  ): ModHam(Ham, params , options  )
{
    construct_ham(params , options);
    _modham = true;
}

//Params: First the L epsilon , then g
void RedBCS::construct_ham(std::vector<double> params, std::vector<std::string> info)
{

    if (info[0] == "reset")
	    this->set_zero();
    for(int i = 0 ; i < this->getL() ; i++)
    {
	    this->setTmat(i,i, params[i]);
    }
    for(int i = 0 ; i < this->getL() ; i++)
	    for(int j = 0 ; j < this->getL() ; j++)
	    {
		    this->setVmat(i,i,j,j,params[getL()]);	
	    }
}

std::string RedBCS::get_info() const
{
    std::string info = Hamiltonian::get_info();
    info += "###Model Hamiltonian Type: RedBCS\n";
    info += ModHam::get_info();
    return info;
}

//void RedBCS::save_file(const std::string & filename)
//{
//    Hamiltonian::save_file(filename);
//    const string filenames = "ham" +  get_short_filename() + filename + ".dat";
//    std::ofstream outfile(filenames.c_str() ,std::ios_base::app );
//    outfile << "###Model Hamiltonian Type: RedBCS" <<std::endl;
//    cout << "###Model Hamiltonian Type: RedBCS" <<std::endl;
//    outfile.close();
//    ModHam::save_file(filenames);
//}

FacInt::FacInt(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options): ModHam( L ,  nup,  ndown, Econst,  params ,  options)
{
    construct_ham(params , options);
    _modham = true;
}

FacInt::FacInt(std::string filename) : ModHam(filename)
{
    construct_ham(_params , _options);
    _modham = true;
}

FacInt::FacInt(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  ): ModHam(Ham, params , options  )
{
    construct_ham(params , options);
    _modham = true;
}

//First the L epsilon, then Chi, then eta
void FacInt::construct_ham(std::vector<double> params, std::vector<std::string> info)
{
    
    double som = 0;
    for(int i = 0 ; i< getL() ; i++)
    {
        som -= params[getL()+1] * params[i] *params[i]*1/2. ; //factor 2. only if degeneracies are 2. (which is expected.)
    }
    setEconst(getEconst() + som);  //getEconst() already contains the non-seniority zero part (which needs to be set by user in advance in the constructor of ModHam), this adds the seniority zero part of Econst
    
    if (info[0] == "reset")
        this->set_zero();
    for(int i = 0 ; i < this->getL() ; i++)
    {
        this->setTmat(i,i, params[i] * params[i] * params[getL()+1]/2.);
    }
    for(int i = 0 ; i < this->getL() ; i++)
	    for(int j = i ; j < this->getL() ; j++)
	    {
		    this->setVmat(i,i,j,j,params[getL()] * params[i] * params[j]);	
	    }
}

std::string FacInt::get_info() const
{
    std::string info = Hamiltonian::get_info();
    info += "###Model Hamiltonian Type: FacInt\n";
    info += ModHam::get_info();
    return info;
}

