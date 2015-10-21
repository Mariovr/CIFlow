#include <string>
#include <fstream>
#include <iostream>

#include "Hamiltonian.h"
#include "ModHam.h"
#include "HamConstruct.h"

using namespace std;

std::string HamConstruct::get_type(std::string filename)
{
    std::string line;
    std::string start = "###Model Hamiltonian Type: ";
    std::ifstream inputfile(filename.c_str()) ;   
    do
    {
        getline(inputfile , line);
    }while(line.find(start) == std::string::npos && !inputfile.eof() );
    int pos = line.find(": ");
    std::string type;
    if (pos != std::string::npos)
        type = line.substr(pos+2, std::string::npos);
    else
        type = "chem";
    inputfile.close();
    return type;
}

std::shared_ptr<Hamiltonian> HamConstruct::generate_ham(std::string filename)
{
    _filename = filename;
    int last_index = filename.find_last_of('.');
    std::shared_ptr<Hamiltonian> ham;
    //Check if inputfile is not in modham format.
    std::string type = get_type(filename);
    cout << type << " is found" << endl;
    if(type.find("chem") != std::string::npos)
    {
        ham.reset( new Hamiltonian(filename) );
        cout << "Generated Chemical Hamiltonian" <<endl;
    }
    else if( type.find("Hub1d") != std::string::npos)
    {
        ham.reset(new Hub1d(filename) );
        cout << "Generated Hub1d Hamiltonian" <<endl;
    }
    else if( type.find("RedBCS") != std::string::npos)
    {
        ham.reset(new RedBCS(filename) );
        cout << "Generated RedBCS Hamiltonian" <<endl;
    }
    else if( type.find("FacInt") != std::string::npos)
    {
        ham.reset(new FacInt(filename) );
        cout << "Generated FacInt Hamiltonian" <<endl;
    }
    else if (type.find("Constrained_DM") != std::string::npos)
    {
        ham.reset(new Constrained_DM(filename) );
        cout << "Generated Constrained_DM Hamiltonian" <<endl;
    }
    else
    {
        cerr << "Type was not found." << endl;
        exit(1);
    }
    return ham;
}



