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
#include <stdlib.h>
#include <iostream>
#include <string>

#include "Irreps.h"

using std::string;
using std::cout;
using std::endl;

Irreps::Irreps(){

   isActivated = false;

}

Irreps::Irreps(const int nGroup){

   if ((nGroup >= 0) && (nGroup <= 7)){
      isActivated = true;
      groupNumber = nGroup;
   } else {
      isActivated = false;
   }

}

bool Irreps::setGroup(const int nGroup){

   if ((nGroup >= 0) && (nGroup <= 7)){
      isActivated = true;
      groupNumber = nGroup;
   } else {
      isActivated = false;
   }
   
   return isActivated;

}

bool Irreps::getIsActivated() const{

   return isActivated;

}

int Irreps::getGroupNumber() const{

   return isActivated ? groupNumber : -1 ;

}

string Irreps::getGroupName() const{

   return isActivated ? getGroupNamePrivate(groupNumber) : "error" ;

}

string Irreps::getGroupName(const int nGroup){

   return ((nGroup>=0)&&(nGroup<=7)) ? getGroupNamePrivate(nGroup) : "error" ;

}

string Irreps::getGroupNamePrivate(const int nGroup){

   if (nGroup==0) return "c1";
   if (nGroup==1) return "ci";
   if (nGroup==2) return "c2";
   if (nGroup==3) return "cs";
   if (nGroup==4) return "d2";
   if (nGroup==5) return "c2v";
   if (nGroup==6) return "c2h";
   if (nGroup==7) return "d2h";
   return "error";

}

int Irreps::getNumberOfIrreps() const{

   return isActivated ? getNumberOfIrrepsPrivate(groupNumber) : -1 ;

}

int Irreps::getNumberOfIrrepsPrivate(const int nGroup){
   
   if (nGroup == 0) return 1;
   if (nGroup <= 3) return 2;
   if (nGroup <= 6) return 4;
   return 8;

}

string Irreps::getIrrepName(const int irrepNumber) const{

   if (!isActivated) return "error1";
   
   if ((irrepNumber<0) || (irrepNumber >= getNumberOfIrreps())) return "error2";
   
   return getIrrepNamePrivate(groupNumber,irrepNumber);
   
}
   
string Irreps::getIrrepNamePrivate(const int nGroup, const int nIrrep){
   
   if (nGroup == 0){
      if (nIrrep == 0) return "A";
   }
   
   if (nGroup == 1){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "Au";
   }
   
   if (nGroup == 2){
      if (nIrrep == 0) return "A";
      if (nIrrep == 1) return "B";
   }
   
   if (nGroup == 3){
      if (nIrrep == 0) return "A'";
      if (nIrrep == 1) return "A''";
   }
   
   if (nGroup == 4){
      if (nIrrep == 0) return "A";
      if (nIrrep == 1) return "B1";
      if (nIrrep == 2) return "B2";
      if (nIrrep == 3) return "B3";
   }
   
   if (nGroup == 5){
      if (nIrrep == 0) return "A1";
      if (nIrrep == 1) return "A2";
      if (nIrrep == 2) return "B1";
      if (nIrrep == 3) return "B2";
   }
   
   if (nGroup == 6){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "Bg";
      if (nIrrep == 2) return "Au";
      if (nIrrep == 3) return "Bu";
   }
   
   if (nGroup == 7){
      if (nIrrep == 0) return "Ag";
      if (nIrrep == 1) return "B1g";
      if (nIrrep == 2) return "B2g";
      if (nIrrep == 3) return "B3g";
      if (nIrrep == 4) return "Au";
      if (nIrrep == 5) return "B1u";
      if (nIrrep == 6) return "B2u";
      if (nIrrep == 7) return "B3u";
   }
   
   return "error2";

}

int Irreps::getTrivialIrrep(){

   return 0 ;

}

int Irreps::directProd(const int n1, const int n2) const{

   if (!isActivated) return -1;
   if ((n1<0) || (n1>=getNumberOfIrreps()) || (n2<0) || (n2>=getNumberOfIrreps())) return -2;
   
   return directProdPrivate(groupNumber,n1,n2);

}

int Irreps::directProdPrivate(const int nGroup, const int n1, const int n2){

   if (n1 == n2) return 0;
   if (n1 > n2) return directProdPrivate(nGroup, n2, n1);
   
   //n1<n2 remains
   if (nGroup <= 3){ //ci, cs or c2 (c1 only n1=n2 possible)
      return 1; //Off-diag terms are always the non-trivial group
   }
   
   if (nGroup <= 6){ //d2, c2v, c2h (with the conventional order, both diagonals are symmetry axes of the multiplication table)
      if (n1 == 0) return n2;
      if (n1 == 1){
         if (n2 == 2) return 3;
         if (n2 == 3) return 2;
      }
      return 1;
   }
   
   if (nGroup == 7){ //d2h (with the conventional order, both diagonals are symmetry axes of the multiplication table)
      if (n1 == 0) return n2;
      if (n1 == 1){
         if (n2 == 2) return 3;
         if (n2 == 3) return 2;
         if (n2 == 4) return 5;
         if (n2 == 5) return 4;
         if (n2 == 6) return 7;
         if (n2 == 7) return 6;
      }
      if (n1 == 2){
         if (n2 == 3) return 1;
         if (n2 == 4) return 6;
         if (n2 == 5) return 7;
         if (n2 == 6) return 4;
         if (n2 == 7) return 5;
      }
      if (n1 == 3){
         if (n2 == 4) return 7;
         if (n2 == 5) return 6;
         if (n2 == 6) return 5;
         if (n2 == 7) return 4;
      }
      if (n1 == 4){
         if (n2 == 5) return 1;
         if (n2 == 6) return 2;
         if (n2 == 7) return 3;
      }
      if (n1 == 5){
         if (n2 == 6) return 3;
         if (n2 == 7) return 2;
      }
      return 1;  
   }

   return -1;

}

void Irreps::printAll(){

   for (int cnt=0; cnt<8; cnt++){
      cout << "######################################################" << endl;
      cout << "Name = " << getGroupNamePrivate(cnt) << endl;
      cout << "nIrreps = " << getNumberOfIrrepsPrivate(cnt) << endl;
      cout << "Multiplication table :" << endl;
      for (int cnt2=-1; cnt2<getNumberOfIrrepsPrivate(cnt); cnt2++){
         for (int cnt3=-1; cnt3<getNumberOfIrrepsPrivate(cnt); cnt3++){
            if ((cnt2 == -1) && (cnt3 == -1)) cout << "\t";
            if ((cnt2 == -1) && (cnt3 >= 0 )) cout << getIrrepNamePrivate(cnt,cnt3) << "\t";
            if ((cnt3 == -1) && (cnt2 >= 0 )) cout << getIrrepNamePrivate(cnt,cnt2) << "\t";
            if ((cnt3 >= 0) && (cnt2 >= 0 )) cout << getIrrepNamePrivate(cnt,directProdPrivate(cnt,cnt2,cnt3)) << "\t";
         }
         cout << endl;
      }
   }
   cout << "######################################################" << endl;
   
}


