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
#ifndef __OrbitalTransform__
#define __OrbitalTransform__

#include <assert.h>
#include <memory>
#include <vector>
#include "Irreps.h"

class Hamiltonian;
class OptIndex;
class UnitaryMatrix;

class OrbitalTransform
{
    public :
        OrbitalTransform(Hamiltonian * ham);
        OrbitalTransform(const OrbitalTransform &);
        virtual ~OrbitalTransform() = default;

        void fillHamCI(Hamiltonian& HamCI);	
        void fillConstAndTmat(Hamiltonian& Ham) const;
        void buildOneBodyMatrixElements();
        void set_unitary(UnitaryMatrix unit);
        Hamiltonian* get_ham() { return _hamorig; }
        void set_ham(Hamiltonian * ham);

        double TmatRotated(const int index1, const int index2) const;
        double get_norb(int irrep) const;
        double get_difference_orig(Hamiltonian &hamin) const;
        void CheckDeviationFromUnitary() const;
        void update_unitary(double *step);

        UnitaryMatrix & get_unitary();

        void do_jacobi_rotation(int k, int l, double theta, Hamiltonian * ham);

        //!multiplies new eigenbasis eigenvecs from the left (basis in rows) in the current unitary, make sure the new basis is in the columns because rotate_active space transposes it to rows automatically.
        void rotate_active_space_vectors(double *eigenvecs);

    private:
        void rotate_old_to_new(std::unique_ptr<double []> * matrix);

        //! the orginal hamiltonian
        Hamiltonian *  _hamorig;
        //! The rotation to perfrom on _hamorig to get the current hamiltonian
        std::unique_ptr<UnitaryMatrix>  _unitary;

        std::unique_ptr<OptIndex> _index;

        int numberOfIrreps;
        Irreps SymmInfo;

        //! Some memory to do the one body work, allocated when needed
        std::vector< std::unique_ptr<double []> > QmatrixWork;
        // The one body matrix elements, used for rotations, only allocated when needed
        std::vector< std::unique_ptr<double []> > OneBodyMatrixElements;

        std::unique_ptr<double []> mem1;
        std::unique_ptr<double []> mem2;

};

#endif

/* vim: set ts=4 sw=4 expandtab :*/
