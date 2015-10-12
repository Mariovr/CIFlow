#ifndef __ModHam__
#define _ModHam_

#include "Hamiltonian.h"

class ModHam : public Hamiltonian
{

    public:
    	ModHam(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options);
    	ModHam(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  );
        ModHam(std::string filename);
	    virtual ~ModHam(){}
	    virtual void construct_ham(std::vector<double> params, std::vector<std::string> info) = 0;
	    virtual std::string get_info() const ;
        void save_file(const std::string & filename = "");
        double get_params(int index){return _params[index];}
        std::string get_options(int index){return _options[index];}
        void set_params(int index, double val){_params[index] = val;}
        void set_options(int index, std::string opt){_options[index] = opt;}

    protected:
        std::vector<double> _params;
        std::vector<std::string> _options;

    
};


//Options: reset, {pos,mom} , {per, aper, none}
class Hub1d : public ModHam
{
    public:

        Hub1d(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options);
        Hub1d(std::string filename );
        Hub1d(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  );
        void construct_ham(std::vector<double> params, std::vector<std::string> info);
        std::string get_info() const;
    private:

};


//Options: reset
class RedBCS: public ModHam
{
    public:
        RedBCS(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options);
        RedBCS(std::string filename ) ;
        RedBCS(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  );
        void construct_ham(std::vector<double> params, std::vector<std::string> info);
        std::string get_info() const;
    private:

};


//Options: reset
class FacInt: public ModHam
{
    public:
        FacInt(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options);
        FacInt(std::string filename);
        FacInt(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  );
        void construct_ham(std::vector<double> params, std::vector<std::string> info);
        std::string get_info() const;
    private:
};

class Constrained_DM : public ModHam
{
    public:
        Constrained_DM(const int L , const int nup, const int ndown, double Econst, std::vector<double> params , std::vector<string> options);
        Constrained_DM(std::string filename);
        Constrained_DM(const Hamiltonian & Ham, std::vector<double> params , std::vector<string> options  );
        void construct_ham(std::vector<double> params, std::vector<std::string> info);
        std::string get_info() const;
    private:

};




#endif //End define __ModHam__
