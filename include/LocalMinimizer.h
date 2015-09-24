#ifndef __LocalMinimizer__
#define __LocalMinimizer__

#include <memory>
#include <random>
#include <vector>
#include <string>

#include "SimulatedAnnealing.h"

class Hamiltonian;
class UnitaryMatrix;
class OrbitalTransform;
class CIMethod;
class CIDens;


/*
 *Kindly provided by ward poelmans
 */

class LocalMinimizer: OrbitalOptimization
{
    public:
        LocalMinimizer(CIMethod * cim);

        LocalMinimizer(CIMethod && cim);

        virtual ~LocalMinimizer() = default;

        double optimize(bool dist_choice);
        double optimize(){return optimize(false); }

        void calc_energy();

        std::vector<std::tuple<int,int,double,double>> scan_orbitals();

        int choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &);

        double get_conv_crit() const;

        void set_conv_crit(double);

        void set_conv_steps(int);

        //std::unique_ptr<CIDens> get_density();

    private:
        //! criteria for convergence of the minimizer
        double _conv_crit;

        //! number of steps in convergence area
        int _conv_steps;

        //! our pseudo-random generator
        std::mt19937 _mt;

        //! only rotation within these irreps (if not empty)
        std::vector<int> _allow_irreps;

};

#endif
