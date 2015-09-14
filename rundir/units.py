# -*- coding: utf-8 -*-
# CIFlow is a very flexible configuration interaction program
# Copyright (C) Ghent University 2014-2015
#
# This file is part of CIFlow.
#
# CIFlow is developed by Mario Van Raemdonck <mario.vanraemdonck@ugent.be>
# a member of the Ghent Quantum Chemistry Group (Ghent University).
# See also : http://www.quantum.ugent.be
#
# At this moment CIFlow is not yet distributed.
# However this might change in the future in the hope that
# it will be useful to someone.
#
# For now you have to ask the main author for permission.
#
#--
# -*- coding: utf-8 -*-

'''Physicochemical constants in atomic units

   These are the physical constants defined in this module (in atomic units):

'''
boltzmann = 3.1668154051341965e-06
avogadro = 6.0221415e23
lightspeed = 137.03599975303575
planck = 6.2831853071795864769

'''Conversion from and to atomic units

   Internally Horton always uses atomic units. Atomic units are consistent,
   similar to the SI unit system: one does not need conversion factors in the
   middle of a computation. This choice facilitates the programming and reduces
   accidental bugs.

   References for the conversion values:

   * B. J. Mohr and B. N. Taylor,
     CODATA recommended values of the fundamental physical
     constants: 1998, Rev. Mod. Phys. 72(2), 351 (2000)
   * The NIST Reference on Constants, Units, and Uncertainty
     (http://physics.nist.gov/cuu/Constants/index.html)
   * 1 calorie = 4.184 Joules

   **Conventions followed by this module:**

   Let foo be is the value of an external unit in internal (atomic) units. The
   way to use this unit is as follows: ``5*foo`` litterally means `five times
   foo`. The result of this operation is a floating point number for this value
   in atomic units.

   **Examples:**

   If you want to have a distance of five angstrom in internal units:
   ``5*angstrom``.

   If you want to convert a length of 5 internal units to angstrom:
   ``5/angstrom``.

   **Remarks:**

   It is highly recommended to perform unit conversions only when data is read
   from the input or data is written to the output. Do not perform any unit conversion
   in other parts of the program.  

   An often recurring question is how to convert a frequency in internal units
   to a spectroscopic wavenumber in inverse centimeters. This is how it can be
   done::

     >>> from units import centimeter, lightspeed
     >>> invcm = lightspeed/centimeter
     >>> freq = 0.00320232
     >>> print freq/invcm

   These are the conversion constants defined in this module:
'''


# *** Generic ***
au = 1.0


# *** Charge ***

coulomb = 1.0/1.602176462e-19

# Mol

mol = avogadro

# *** Mass ***

kilogram = 1.0/9.10938188e-31

gram = 1.0e-3*kilogram
miligram = 1.0e-6*kilogram
unified = 1.0e-3*kilogram/mol
amu = unified

# *** Length ***

meter = 1.0/0.5291772083e-10

decimeter = 1.0e-1*meter
centimeter = 1.0e-2*meter
milimeter = 1.0e-3*meter
micrometer = 1.0e-6*meter
nanometer = 1.0e-9*meter
angstrom = 1.0e-10*meter
picometer = 1.0e-12*meter

# *** Volume ***

liter = decimeter**3

# *** Energy ***

joule = 1/4.35974381e-18

calorie = 4.184*joule
kjmol = 1.0e3*joule/mol
kcalmol = 1.0e3*calorie/mol
electronvolt = (1.0/coulomb)*joule
rydberg = 0.5

# *** Force ***

newton = joule/meter

# *** Angles ***

deg = 0.017453292519943295
rad = 1.0

# *** Time ***

second = 1/2.418884326500e-17

nanosecond = 1e-9*second
femtosecond = 1e-15*second
picosecond = 1e-12*second

# *** Frequency ***

hertz = 1/second

# *** Pressure ***

pascal = newton/meter**2
bar = 100000*pascal
atm = 1.01325*bar

# *** Temperature ***

kelvin = 1.0

# *** Dipole ***

debye = 0.39343031369146675 # = 1e-21*coulomb*meter**2/second/lightspeed

# *** Current ***

ampere = coulomb/second



