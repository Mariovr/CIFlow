#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include <hdf5.h>
#include <sstream>
#include <chrono>

#include "LocalMinimizer.h"
#include "Hamiltonian.h"
#include "OrbitalTransform.h"
#include "UnitaryMatrix.h"
#include "Irreps.h"
#include "CIDens.h"
#include "CIMethod.h"



/*
 *File based on the Local optimization for exact DOCI from Ward Poelmans. 
 */

LocalMinimizer::LocalMinimizer(CIMethod * cim): OrbitalOptimization(cim)
{
   _conv_crit = 1e-6;
   _conv_steps = 25;

   std::random_device rd;
   _mt = std::mt19937(rd());

   // expect to be comma seperated list of allowed irreps
   char *irreps_env = getenv("DOCI_ALLOWED_IRREPS");
   if(irreps_env && strlen(irreps_env) > 0)
   {
      std::string irreps_string = irreps_env;
      const std::string delim = ",";
      Irreps syminfo(_cim->get_ham()->getNGroup());

      auto start = 0U;
      auto end = irreps_string.find(delim);
      try
      { 
         while (true)
         {
            auto elem = irreps_string.substr(start, end - start);
            if(elem.empty())
               break;

            int cur_irrep = std::stoi(elem);

            if(cur_irrep >= 0 && cur_irrep < syminfo.getNumberOfIrreps())
               _allow_irreps.push_back(cur_irrep);

            start = end + delim.length();

            if(end >= std::string::npos)
               break;

            end = irreps_string.find(delim, start);
         }
      } catch (std::exception& e) {
         std::cout << "Invalid value in DOCI_ALLOWED_IRREPS" << e.what() << std::endl;
      }

      std::sort(_allow_irreps.begin(), _allow_irreps.end());
      std::cout << "Allowed irreps: ";
      for(auto &elem: _allow_irreps)
         std::cout << elem << " ";
      std::cout << std::endl;
   }
}

std::vector< std::tuple<int,int,double,double> > LocalMinimizer::scan_orbitals()
{
   auto start = std::chrono::high_resolution_clock::now();
   Hamiltonian & ham = *_optham;

   std::function<double(int,int)> getT = [&ham] (int a, int b) -> double { return ham.getTmat(a,b); };
   std::function<double(int,int,int,int)> getV = [&ham]  (int a, int b, int c, int d) -> double { return ham.getVmat(a,b,c,d); };

   std::vector< std::tuple<int,int,double,double> > pos_rotations;
   // worst case: c1 symmetry
   pos_rotations.reserve(_cim->get_l()*(_cim->get_l()-1)/2);

   std::unique_ptr<CIDens> rdm { _cim->get_density() } ;// { dynamic_cast<DensDOCI>(dd) };
   rdm->construct_density(true);

   for(int k_in=0;k_in<_cim->get_l();k_in++)
      for(int l_in=k_in+1;l_in<_cim->get_l();l_in++)
         if(_optham->getOrbitalIrrep(k_in) == _optham->getOrbitalIrrep(l_in))
         {
            if(!_allow_irreps.empty() && std::find(_allow_irreps.begin(), _allow_irreps.end(), _optham->getOrbitalIrrep(k_in)) == _allow_irreps.end() )
               continue;

            auto found = rdm->find_min_angle(k_in,l_in,0.3,getT,getV);

            if(!found.second)
               // we hit a maximum
               found = rdm->find_min_angle(k_in,l_in,0.01,getT,getV);

            if(!found.second)
               // we're still stuck in a maximum, skip this!
               continue;

            // skip angles larger than Pi/2
            if(fabs(found.first)>M_PI/2.0)
               continue;

            double new_en = rdm->calc_rotate(k_in,l_in,found.first,getT,getV);

            assert(found.second && "Shit, maximum!");

            pos_rotations.push_back(std::make_tuple(k_in,l_in,found.first,new_en));
         }

   auto end = std::chrono::high_resolution_clock::now();

   //std::cout << "Orbital scanning took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   assert(pos_rotations.size()>0);

   return pos_rotations;
}

/**
 * Do the local minimization
 * @param dist_choice if set to true, we use choose_orbitals to choose
 * which pair of orbitals to use (instead of the lowest one)
 */
double LocalMinimizer::optimize(bool dist_choice)
{
   int converged = 0;
   _opt_energy = 0.;

   auto start = std::chrono::high_resolution_clock::now();

   std::pair<int,int> prev_pair(0,0);

   int iters = 1;

   auto& ham2 = *_optham;

   while(converged<_conv_steps)
   {
      auto list_rots = scan_orbitals();

      std::sort(list_rots.begin(), list_rots.end(),
            [](const std::tuple<int,int,double,double> & a, const std::tuple<int,int,double,double> & b) -> bool
            {
            return std::get<3>(a) < std::get<3>(b);
            });

      //for(auto& elem: list_rots)
         //std::cout << std::get<0>(elem) << "\t" << std::get<1>(elem) << "\t" << std::get<3>(elem)+ham2.getEconst() << "\t" << std::get<2>(elem) << std::endl;

      int idx = 0;
      std::pair<int,int> tmp;

      if(dist_choice)
      {
         idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = choose_orbitalpair(list_rots);

         tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

         if(tmp==prev_pair)
            idx = 0;
      }

      tmp = std::make_pair(std::get<0>(list_rots[idx]), std::get<1>(list_rots[idx]));

      // don't do the same pair twice in a row
      if(tmp==prev_pair)
         idx++;

      const auto& new_rot = list_rots[idx];
      prev_pair = std::make_pair(std::get<0>(new_rot), std::get<1>(new_rot));

      if(dist_choice)
         std::cout << iters << " (" << converged << ") Chosen: " << idx << std::endl;

      assert(ham2.getOrbitalIrrep(std::get<0>(new_rot)) == ham2.getOrbitalIrrep(std::get<1>(new_rot)));
      // do Jacobi rotation twice: once for the Hamiltonian data and once for the Unitary Matrix
      _orbtrans->do_jacobi_rotation(std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot), & ham2);
      //_orbtrans->get_unitary().jacobi_rotation(ham2.getOrbitalIrrep(std::get<0>(new_rot)), std::get<0>(new_rot), std::get<1>(new_rot), std::get<2>(new_rot)); //Not necessary is automatically done.

      update_cim();

      
      std::stringstream h5_name;

      if(iters%10==0)
      {
         h5_name << "unitary-" << iters ;
         _opt_unitary.reset(new UnitaryMatrix(_orbtrans->get_unitary()));
         save_unitary(h5_name.str() , true);
      }

      if(iters%50==0)
      {
         h5_name.str("");
         h5_name << "ham-" << iters << ".h5";
         ham2.save(h5_name.str());
      }

      if(fabs(_opt_energy- _cim->get_ci_energy() )<_conv_crit)
         converged++;

      //std::cout << iters << " (" << converged << ")\tRotation between " << std::get<0>(new_rot) << "  " << std::get<1>(new_rot) << " over " << std::get<2>(new_rot) << " E_rot = " << std::get<3>(new_rot)+ham2.getEconst() << "  E = " <<  _cim->get_ci_energy()  << "\t" << fabs(_opt_energy- _cim->get_ci_energy() ) << std::endl;

      _opt_energy= _cim->get_ci_energy();

      iters++;

      if(iters>1000)
      {
         std::cout << "Done 1000 steps, quiting..." << std::endl;
         break;
      }
   }

   auto end = std::chrono::high_resolution_clock::now();

   std::cout << "Minimization took: " << std::fixed << std::chrono::duration_cast<std::chrono::duration<double,std::ratio<1>>>(end-start).count() << " s" << std::endl;

   _opt_unitary.reset(new UnitaryMatrix(_orbtrans->get_unitary()));
   save_unitary("local", false);
   //ham2.save_file("local");

   return _cim->get_ci_energy();
}

double LocalMinimizer::get_conv_crit() const
{
   return _conv_crit;
}

void LocalMinimizer::set_conv_crit(double crit)
{
   _conv_crit = crit;
}

void LocalMinimizer::set_conv_steps(int steps)
{
   _conv_steps = steps;
}

/**
 * Choose a pair of orbitals to rotate over, according to the distribution of their relative
 * energy change.
 * @param orbs the list returned by scan_orbitals()
 * @return the index of the pair of orbitals in orbs
 */
int LocalMinimizer::choose_orbitalpair(std::vector<std::tuple<int,int,double,double>> &orbs)
{
   std::uniform_real_distribution<double> dist(0, 1);

   const double choice = dist(_mt);

   double norm = 0;

   for(auto &orb_pair: orbs)
      norm += (_opt_energy - std::get<3>(orb_pair));

   double cum = 0;
   for(int i=0;i<orbs.size();i++)
   {
      cum += (_opt_energy - std::get<3>(orbs[i]))/norm;
      if(choice < cum)
         return i;
   }

   assert(0 && "Should never ever be reached!");
   return -1;
}

/*  
std::unique_ptr< CIDens> LocalMinimizer::get_density()
{
   std::unique_ptr<CIDens> cid  = _cim->get_density() ;
   //cid->construct_density(true);
   return cid;
}*/

/* vim: set ts=3 sw=3 expandtab :*/
