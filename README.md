CIFlow: Very flexible configuration interaction program
===============================================

Prerequisites
-------------

### 1. Boost.Test


This project uses the Boost.Test testing framework. We advise you to
check the results of the unit tests before using a particular binary 
in production runs.

### 2. Mac OS

Enable Macports static libraries by adding the 
following lines to your ~/.profile

```
# Enable static libraries
export C_INCLUDE_PATH=$C_INCLUDE_PATH:$PORTS_HOME/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:$PORTS_HOME/include
export LIBRARY_PATH=$LIBRARY_PATH:$PORTS_HOME/lib
export DYLD_FALLBACK_LIBRARY_PATH=$DYLD_FALLBACK_LIBRARY_PATH:$PORTS_HOME/lib
```

### 3. Psi4 (optional)

CIFlow has a Hamiltonian object which is able to read in matrix elements
from a plugin to
[Psi4, Ab initio quantum chemistry](http://www.psicode.org),
which works on version psi4.0b5 and higher.

Install Psi4 (with plugin option) and set the relevant
directories in the ./mointegrals/Makefile.

Build procedure
---------------

### 1. Laptop

The CIFlow project is packaged in a library. Several driver routines can be
used to accommodate computational projects.

You can build the library, the driver routines, the documentation, the tests
or all these components at once by giving the respective commands in the
base directory:

   > make library

   > make program

   > make doc

   > make tests

   > make all

To view the Doxygen manual, point your favorite webbrowser to ./doc/html/index.html.

### 2. Cluster

To begin work on the Ghent HPC infrastructure, please clone the repository
into $VSC_HOME/devel/CIFlow.

To compile on the cluster, please use the scripts which are located in the 'scr'
directory. These scripts submit the compilation to the respective debug
queues of the HPC infrastructure.

You can either move to the cluster you want to use,

   > module swap cluster/$CLUSTER

and run the cluster_compilation script, or use the full_cluster_compilation script,
which submits compilation instructions to all clusters which are supported.

All resulting binaries should be located in your $VSC_HOME/devel/CIFlow/bin directory,
and should contain the name of the cluster for which they were compiled.

Always check your CIFlow.e files for any errors reported by Boost.Test.

Usage
-----

CIFlow can be used in interactive mode. Run the binary and it asks for the necessary input.
To obtain some inspiration to use the functions defined in the CIFlow library check src/Main.cpp
and the examples in the tests/src directory.

### 1. Laptop specifics

### 2. Cluster specifics

Please load the following modules before using CIFlow:

   > module load HDF5/1.8.14-ictce-7.1.2-serial

   > arpack-ng/3.1.5-ictce-7.1.2-mt 

   > Boost/1.57.0-intel-2015a-Python-2.7.9

Contributors
------------

CIFlow is developed and maintained by Mario Van Raemdonck at the GQCG.

This version contains contributions of the following people from the GQCG and the Center for Molecular Modeling:

Guillaume Acke,
Ward Poelmans, 
https://github.com/SebWouters/CheMPS2  

Developing
----------

If you want to contribute to CIFlow, take into account that: 

It is perfectly ok to work in your own branch and do there whatever you like. It is also ok to upload these branches, such that others can see and review your code. As soon as one of your branches is mature enough, use git rebase to apply your patches to the master branch.
But always ask me to review your patches before they go into the official master branch.

The purpose of the tests is to keep the program working as we add more features. When a patch breaks previous tests, it is simply unacceptable.

When you find a bug, first write a test that fails because of the bug. Then fix the bug. Commit both test and fix in one patch.

This is a great way of fixing bugs and making tests. It also guarantees that the bug you found will not be reintroduced accidentally in later versions of CIFlow.

Documentation and clean coding standards are encouraged but not necessary.
