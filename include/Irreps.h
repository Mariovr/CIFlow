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
#ifndef IRREPS_H
#define IRREPS_H

#include <string>
using std::string;

/** Irreps class.
    \author Sebastian Wouters <sebastianwouters@gmail.com>
    \date February 7, 2013
    
    Class containing the irrep conventions and multiplication tables. The program requires Abelian point groups with real character tables, with hence \f$I_{\alpha} \otimes I_{\alpha} = I_{trivial}\f$.
  
    \section irreps_conv Irrep conventions
    
    The same conventions as in Psi4 (beta3) are used. For convenience, they are listed below:\n
    \f[
    \begin{tabular}{|l|cccccccc|}
    \hline
             & 0  & 1   & 2   & 3   & 4  & 5   & 6   & 7   \\
    \hline
      0: c1  & A  &     &     &     &    &     &     &     \\
      1: ci  & Ag & Au  &     &     &    &     &     &     \\
      2: c2  & A  & B   &     &     &    &     &     &     \\
      3: cs  & A' & A'' &     &     &    &     &     &     \\
      4: d2  & A  & B1  & B2  & B3  &    &     &     &     \\
      5: c2v & A1 & A2  & B1  & B2  &    &     &     &     \\
      6: c2h & Ag & Bg  & Au  & Bu  &    &     &     &     \\
      7: d2h & Ag & B1g & B2g & B3g & Au & B1u & B2u & B3u \\
    \hline
    \end{tabular}
    \f]*/
class Irreps{

   public:
   
      //! Constructor
      Irreps();
      
      //! Constructor 2
      /** \param nGroup The group number (0 <= nGroup <= 7; else isActivated remains false) */
      Irreps(const int nGroup);
      
      //! Set the group
      /** \param nGroup Number from 0 to 7 (7 included)
          \return Validity of the group number (error returns false) */
      bool setGroup(const int nGroup);
      
      //! Whether the group number is already activated
      /** \return Whether the group number is already activated */
      bool getIsActivated() const;
      
      //! Get the group number
      /** \return The group number (-1 means not activated) */
      int getGroupNumber() const;
      
      //! Get the name of the group
      /** \return The group name ("error" means not activated) */
      string getGroupName() const;
      
      //! Get the name of the group corresponding to nGroup
      /** \param nGroup Group number
          \return The group name corresponding to nGroup */
      static string getGroupName(const int nGroup);
      
      //! Get the number of irreps for the currently activated group
      /** \return The number of irreps for the currently activated group (-1 means not activated) */
      int getNumberOfIrreps() const;
      
      //! Get the name of the irrep with number irrepNumber of the activated group. The irrep with number 0 is always the trivial irrep.
      /** \param irrepNumber The irrep number
          \return The irrep name (not activated returns "error1"; wrong number returns "error2") */
      string getIrrepName(const int irrepNumber) const;
      
      //! Get the irrep number of the trivial irrep (always 0)
      /** \return The trivial irrep number */
      static int getTrivialIrrep();
      
      //! Get the direct product of the irreps with numbers n1 and n2 for the currently activated group
      /** \param n1 The number of the first irrep
          \param n2 The number of the second irrep
          \return The number of the direct product (-1 means not activated; -2 means n1 or n2 out of bound) */
      int directProd(const int n1, const int n2) const;
      
      //! Print all info contained in this class
      static void printAll();
   
   private:
   
      //whether a relevant group number is currently set
      bool isActivated;
      
      //the currently set group number (check isActivated)
      int groupNumber;
      
      //static member functions containing the group names, the number of irreps, the names of the irreps, and the multiplication tables
      static string getGroupNamePrivate(const int nGroup);
      static int getNumberOfIrrepsPrivate(const int nGroup);
      static string getIrrepNamePrivate(const int nGroup, const int nIrrep);
      static int directProdPrivate(const int nGroup, const int n1, const int n2);

};
#endif

