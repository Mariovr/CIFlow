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
#! /usr/bin/env python 
from string import maketrans
import re
import os
import numpy as np
import sys
from itertools import permutations
from math import sqrt
import math

import detwrite as dw
import read_psi as rp
import filecollector as fc
import copy
if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    import plotfunctions as pf

# Using the generator pattern (an iterable)
class Bitstring(object):
    def __init__(self, norb, nalpha):
        self.norb = norb
        self.nalpha = nalpha
        from scipy.special import binom
        self.dim = binom(self.norb , self.nalpha)
        self.reset()
        print 'A Bitstring generator of ', nalpha , ' electrons in ', norb , ' orbitals is constructed with dimension : ' , self.dim , ' and first element: ' , self.current , ' in string notation: ' , "{0:b}".format(self.current)

    def __iter__(self):
        return self

    # Python 3 compatibility
    def __next__(self):
        return self.next()

    def reset(self):
        self.current = self.start_bit(self.nalpha)
        self.pos = 0

    def start_bit( self, nalpha):
        return 2**nalpha -1

    #    00111000     # 56
    #    01001000     # twos complement, -56
    # &= 00001000
    # Determines lowest set bit
    def lowestSet(self, int_type):
        low = (int_type & -int_type)
        lowBit = -1
        while (low):
            low >>= 1
            lowBit += 1
        return(lowBit)

    def next(self):
        if self.pos  == self.dim:
            raise StopIteration()
        else:
            v= self.current
            t = v | ( v - 1); # t gets v's least significant 0 bits set to 1
            #print "{0:b}".format(t)
            self.current = (t+1) # Next set to 1 the most significant bit to change, and set to 0 the least significant ones
            #print "{0:b}".format(self.current)
            self.current = self.current | ((((~t) & (-~t) ) - 1) >> ( self.lowestSet(v)+ 1 )); #finally add the ones before the changed bit to the upmost right and in between zeros
            self.pos += 1 
            return  v

class Reader(object):
    """
    Class to read in outputfiles from CIFlow.
    """
    def __init__(self, filename , regexp = None , read_ham = False):
        self.filename = filename
        self.regexp = regexp
        self.header = self.read_header()
        if read_ham:
            self.readham()
            self.make_matrix()
            self.diagonalize()

    def readham(self):
        with open(self.filename,"r") as file:
            text = file.read()
        self.data = re.findall(self.regexp,text)

    def make_matrix(self):
        dim = self.get_dim()
        self.matrix = np.zeros((dim,dim))
        self.fill_ham(dim = dim)

    def diagonalize(self):
        from scipy.linalg import eigh #see: http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.linalg.eigh.html#scipy.linalg.eigh for more options
        self.eigval, self.eigvec = eigh(self.matrix)

    def print_matrix_data(self):
        print self.matrix
        print np.shape(self.matrix)

    def print_data(self):
        for dat in self.data:
            print dat[0] , '   ' , dat[1]

    def maketranslation(self,first = "dD",after = "ee"): 
        #use on a string: d = "How would you do that", d.translate(transtable) , print d -> "How woule you eo that"
        return maketrans(first, after)

class CIFlow_Reader(Reader):
    def __init__(self, filename , regexp = None , read_ham = False):
        super(CIFlow_Reader,self).__init__(filename , regexp = regexp, read_ham = read_ham)
        self.get_translate_table()

        if read_ham:
            self.get_evals()
    
    def get_evals(self):
        print '#energies of CI matrix:'
        for i in self.eigval[:4]:
            print i+self.header['constE']

    def read_header(self):
        dict = {}
        with open(self.filename,"r") as file:
            for line in file:
                if '#nup' in line:
                    dict['nup'] = int(line.split()[-1])
                if '#ndown' in line:
                    dict['ndown']= int(line.split()[-1])
                if '#norbs' in line:
                    dict['norbs'] = int(line.split()[-1])
                if 'groundstate energy' in line:
                    dict['energy'] = float(line.split()[-1])
                if 'dimension' in line:
                    dict['dim'] = int(line.split()[-1])
                if 'constant energy' in line:
                    dict['constE'] = float(line.split()[-1])
                if 'CIMETHOD' in line:
                    dict['CIMETHOD'] = line.split()[-1]
        return dict

    def get_psi_input(self, savename = None, print_string = False):
        start =False
        psiinput = ""
        with open(self.filename,"r") as file:
            for line in file:
                if '-------------------------' in line:
                    if start:
                        break
                    else:
                        start = True
                    continue
                if start:
                    psiinput += line
        
        if savename:
            with open(savename , 'w') as file:
                file.write(psiinput)

        if print_string:
            print psiinput

    def get_dim(self):
        return int(re.match(r'(\d+),',self.data[-1][0]).group(1))+1

    def fill_ham(self, dim = -1):
        for dat in self.data:
            self.matrix[int(re.match(r'(\d+),',dat[0]).group(1)), int(re.match(r'\d+,\s+(\d+)',dat[0]).group(1))] = float(dat[1])
    
    def read_rdm(self, twordm = True):
        self.ordm = [np.zeros((self.header['norbs'] ,self.header['norbs']) , dtype=np.float ) for i in range(2) ]
        if twordm:
            self.trdm = [np.zeros((self.header['norbs'] ,self.header['norbs'], self.header['norbs'] ,self.header['norbs']) , dtype=np.float ) for i in range(3) ] 
        reg1rdm = r'^(\d+)\s(\d+)\s([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)'  #onerdm
        reg2rdm = r'^(\d+)\s(\d+)\s(\d+)\s(\d+)\s([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)'  #twordm
        startrdm = -1 
        with open(self.filename,"r") as file:
            for line in file:
                if 'rdm of the' in line:
                    startrdm += 1
                elif 'rdm case' in line:
                    startrdm += 1
                if startrdm in [0,1]:
                    match = re.search(reg1rdm , line)
                    if match:
                        self.ordm[startrdm][int(match.group(1)) , int(match.group(2)) ] = float(match.group(3))
                        self.ordm[startrdm][int(match.group(2)), int(match.group(1)) ] = float(match.group(3))
                if startrdm >= 2:
                    if twordm:
                        match = re.search(reg2rdm , line)
                        if match:#watch out 8 fold permutation symmetry
                            #print match.group(0)
                            self.trdm[startrdm-2][int(match.group(1)) , int(match.group(2)) , int(match.group(3)) , int(match.group(4)) ] = float(match.group(5))
                            if startrdm-2 != 1: #only when alpha, beta, alpha,beta block we dont have this symmetry (if nalpha != nbeta)
                                self.trdm[startrdm-2][int(match.group(2)) ,int(match.group(1)) ,  int(match.group(4)) ,int(match.group(3))  ] =  float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(3)) , int(match.group(4)) , int(match.group(1)) , int(match.group(2)) ] = float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(4)) ,int(match.group(3)) , int(match.group(2)),  int(match.group(1))  ] = float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(3)) , int(match.group(4)) , int(match.group(2)) ,int(match.group(1))  ] = -1.* float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(1)) , int(match.group(2)) , int(match.group(4)), int(match.group(3))  ] = -1.*float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(2)) ,int(match.group(1)) ,  int(match.group(3)) , int(match.group(4)) ] = -1. * float(match.group(5))
                                self.trdm[startrdm-2][int(match.group(4)) ,int(match.group(3)) ,  int(match.group(1)) , int(match.group(2)) ] = -1.*float(match.group(5))
                                #print self.trdm[startrdm-2][int(match.group(1)) , int(match.group(2)) , int(match.group(3)) , int(match.group(4)) ] 
              
    
    def print_rdm(self, twordm = True, filename = 'densitymatrix.dat'):
        dim = self.header['norbs']
        with open(filename , 'w') as dfile:
            dfile.write('#One rdm of the alpha electrons:\n')
            for i in range(dim):
                for j in range(dim):
                    dfile.write('%d\t%d\t%.14f\n' %(i,j , self.ordm[0][i,j]))
            dfile.write('#One rdm of the beta electrons:\n')
            for i in range(dim):
                for j in range(dim):
                    dfile.write('%d\t%d\t%.14f\n' %(i,j , self.ordm[1][i,j]))
            dfile.write('#Two rdm of the alpha electrons (up up up up):\n')
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for l in range(dim):
                            dfile.write('%d\t%d\t%d\t%d\t%.14f\n' %(i,j, k ,l ,  self.trdm[0][i,j,k,l]))
            dfile.write('#Two rdm of the mixed spin electrons (up down up down):\n')
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for l in range(dim):
                            dfile.write('%d\t%d\t%d\t%d\t%.14f\n' %(i,j, k ,l ,  self.trdm[1][i,j,k,l]))
            dfile.write('#Two rdm of the beta electrons (up down up down):\n')
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for l in range(dim):
                            dfile.write('%d\t%d\t%d\t%d\t%.14f\n' %(i,j, k ,l ,  self.trdm[2][i,j,k,l]))

    def read_matrix(self, filename, comment = '#'):
        #return  np.loadtxt(self.filename, dtype='float', comments= '#')
        reg = r'(([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)\s*)+'
        matrix = []
        with open(filename , 'r') as file:
            for line in file:
                if line[0] != comment and re.search(reg,line) and not re.search(r'[a-zA-Z]', line[0]):
                    row = map(float , line.split())
                    matrix.append(row)

        return np.array(matrix)

    def get_translate_table(self):
        """
        This function associates bitstring to the matrix indices
        """
        self.groundstate = []
        if self.header['CIMETHOD'] == 'DOCI':
            reg = r'^(\d+)\s+(\d+)\s+([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)' 
        else:
            reg = r'^(\d+)\s+(\d+)\s+(\d+)\s+([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)' 
        self.mapdict = {}
        with open(self.filename,"r") as file:
            index = 0 
            for line in file:
                match = re.search(reg,line)
                if match:
                    if self.header['CIMETHOD'] == 'DOCI':
                        self.mapdict[match.group(2)+ '|' + match.group(2)] = int(match.group(1))
                        self.groundstate.append(float(match.group(3)))
                    else:
                        self.mapdict[match.group(2)+'|'+match.group(3) ] = int(match.group(1))
                        self.groundstate.append(float(match.group(4)))
                    assert(index == int(match.group(1)) ), str(index) + 'not found ' + 'WARNING' + match.group(1)
                    index += 1
                if "#Significant determinants" in line:
                    break
        assert(len(self.mapdict) == self.header['dim']), 'length dict: ' + str(len(self.mapdict)) +'header: '+ str(self.header['dim'])+ self.filename
        assert(len(self.groundstate) == self.header['dim']), 'length groundstate: ' + str(len(self.mapdict)) +'header: '+ str(self.header['dim']) + self.filename
        self.groundstate = np.array(self.groundstate)  
        #if self.groundstate != None:
            #print 'We have read in the groundstate vec: ', self.groundstate

    def get_mat_element(self,fbitstring , sbitstring):
        #the bitstrings should have the form 000111|000111 the up spin part is separated by the down spin part (with a |) the length of each part is norb for FCI and FCI_FILE, for DOCI just the up spin configuration is oke.
        return self.matrix[self.mapdict[fbitstring],self.mapdict[sbitstring]]

    def get_groundstate(self, bitstring):
        #Give the bitstring which corresponds to one basisvector of the basis on which the groundstate is projected.
        return self.groundstate[self.mapdict[bitstring]]

    def get_energies(self):
        with open(self.filename , 'r') as file:
            data = file.read()
        regexp = r'^([\-+]?\d+\.\d+[eEdD]?[\-+]?\d*)\s*$'
        elist = re.findall( regexp, data, re.MULTILINE)
        self.energies = map(float , elist)
        return self.energies

    def get_properties(self):
        properties = {}
        with open(self.filename , 'r') as file:
            bc = 0
            for line in file:
                #print line
                if '#' == line[0] or '-' == line[0] or '\n' == line[0] or ' ' == line[0]:
                    regexp = r'^#Property\s*eigvec\s*:\s*\d+\s*(\w+)\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d*)\s*$'
                    match = re.search(regexp , line)
                    if match:
                        properties[match.group(1)] = match.group(2)
                        print 'we extracted property:' , match.group(1) , ' with value : '  , match.group(2)
                else:
                    bc += 1
                    if bc > 15: #take into account model hamiltonians their parameters don't start with comment
                        break
        self.properties = properties
        return self.properties

    def get_max_det(self):
        max_val = ('' , 0.)
        for i, val in enumerate( self.groundstate):
            if abs(val) > max_val[1]:
                max_val = (i ,abs(val) )
        for key,value in self.mapdict.iteritems(): 
            if value == max_val[0]:
                max_val = (key , max_val[1])
                break
        return max_val

    def get_selected_coef(self, detlist):
        #print self.mapdict
        coeflist = [ self.groundstate[self.mapdict[det]]**2 for det in detlist]
        som = 0
        if coeflist:
            som = reduce(lambda x , y : x+y,  coeflist)
        return som

    #all the dets in output_cif should exist in self, otherwise comment out the lines under except and pass the exception.
    def calc_overlap(self, output_cif):
        som = 0
        for key, value in output_cif.mapdict.items():
            try:
                som += self.groundstate[self.mapdict[key]] * output_cif.groundstate[value]
            except Exception as e:
                #print self.mapdict
                #print key
                #print e 
                #sys.exit(1)
                pass
        return abs(som)

    def calc_1rdm(self):
        dim = self.header['norbs']
        t2index = [np.zeros((dim,dim)) for i in range(2)]
        upbitgenerator = Bitstring(dim,self.header['nup']) #Creates bitstring generator for 5 orbitals and 3 electrons
        downbitgenerator =Bitstring(dim,self.header['ndown']) #Creates bitstring generator for 5 orbitals and 3 electrons 
        for upbit in bitgenerator:
            print 'We are at position: ' , bitgenerator.pos-1 , ' with bitstring: ', "{0:b}".format(bit)

    def calc_2rdm(self):
        dim = self.header['norbs']
        t4index = [np.zeros((dim,dim, dim, dim)) for i in range(3)]

    def construct_cum(self, ordm , trdm):
        dim = self.header['norbs']
        cumul = [np.zeros((dim,dim, dim, dim)) for i in range(3)]
        hfcor = [np.zeros((dim,dim, dim, dim)) for i in range(3)]
        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        cumul[0][i,j,k,l] = 1./2.*trdm[0][i , j , k , l] -1./2.* (ordm[0][i,k]* ordm[0][j , l]  - ordm[0][i, l ] * ordm[0][ j , k])
                        cumul[1][i,j,k,l] = 1./2.*trdm[1][i , j , k , l]-1./2. *(ordm[0][i,k]* ordm[1][j , l]  )
                        cumul[2][i,j,k,l] = 1./2.*trdm[2][i , j , k , l]-1./2. *(ordm[1][i,k]* ordm[1][j , l]  - ordm[1][i, l ] * ordm[1][ j , k])
                        hfcor[0][i,j,k,l] = 1./2.* (ordm[0][i,k]* ordm[0][j , l]  - ordm[0][i, l ] * ordm[0][ j , k])
                        hfcor[1][i,j,k,l] = 1./2.* (ordm[0][i,k]* ordm[1][j , l]  )
                        hfcor[2][i,j,k,l] = 1./2.* (ordm[1][i,k]* ordm[1][j , l]  - ordm[1][i, l ] * ordm[1][ j , k])


        return hfcor , cumul

    def project_herm_sym(self, orblist):
        dim = self.header['norbs']
        for i in range(dim):
            for j in range(dim):
                if i not in orblist and j not in orblist:
                    self.ordm[0][i,j] = 0
                    self.ordm[1][i,j] = 0
                if i in orblist and j not in orblist:
                    self.ordm[0][i,j] = self.ordm[0][i,j] * 1./2. 
                    self.ordm[1][i,j] = self.ordm[1][i,j] * 1./2. 
                #if j in orblist and i not in orblist:
                #    self.ordm[0][i,j] = 0
                #    self.ordm[1][i,j] = 0
                if j in orblist and i not in orblist:
                    self.ordm[0][i,j] /= 2.
                    self.ordm[1][i,j] /= 2.
                for k in range(dim):
                    for l in range(dim):
                        if i not in orblist and j not in orblist and k not in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i in orblist and j not in orblist and k in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] =  self.trdm[0][i , k , j , l] * 1./4.
                            self.trdm[1][i , k , j , l] =  self.trdm[1][i , k , j , l] * 1./4.
                            self.trdm[2][i , k , j , l] =  self.trdm[2][i , k , j , l] * 1./4.
                        if i in orblist and j not in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] /=4.
                            self.trdm[1][i , k , j , l] /=4.
                            self.trdm[2][i , k , j , l] /=4.
                        if i not in orblist and j in orblist and k in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] /=4.
                            self.trdm[1][i , k , j , l] /=4.
                            self.trdm[2][i , k , j , l] /=4.
                        if i not in orblist and j in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] /=4.
                            self.trdm[1][i , k , j , l] /=4.
                            self.trdm[2][i , k , j , l] /=4.
                        if i in orblist and j not in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] =0.
                            self.trdm[1][i , k , j , l] =0.
                            self.trdm[2][i , k , j , l] =0.
                        if i not in orblist and j in orblist and k in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] =0.
                            self.trdm[1][i , k , j , l] =0.
                            self.trdm[2][i , k , j , l] =0.
                        if i not in orblist and j in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] =0.
                            self.trdm[1][i , k , j , l] =0.
                            self.trdm[2][i , k , j , l] =0.
                        if i in orblist and j in orblist and k not in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] =0
                            self.trdm[1][i , k , j , l] =0
                            self.trdm[2][i , k , j , l] =0
                        if i not in orblist and j not in orblist and k in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] =0
                            self.trdm[1][i , k , j , l] =0
                            self.trdm[2][i , k , j , l] =0
                        if i not in orblist and j not in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i not in orblist and j not in orblist and k in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0

                        if i not in orblist and j in orblist and k not in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i in orblist and j not in orblist and k not in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i in orblist and j in orblist and k in orblist and l not in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i not in orblist and j in orblist and k in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i in orblist and j not in orblist and k in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0
                        if i in orblist and j in orblist and k not in orblist and l in orblist:
                            self.trdm[0][i , k , j , l] = 0
                            self.trdm[1][i , k , j , l] = 0
                            self.trdm[2][i , k , j , l] = 0

    def project_herm(self, orblist):
        dim = self.header['norbs']
        ordm = [np.zeros((dim,dim)) for i in range(3)]
        trdm = [np.zeros((dim,dim, dim, dim)) for i in range(3)]
        for i in range(dim):
            for j in range(dim):
                if i in orblist and j in orblist:
                    ordm[0][i,j] = self.ordm[0][i,j]
                    ordm[1][i,j] = self.ordm[1][i,j]
                #if i in orblist and j not in orblist:
                #    ordm[0][i,j] = self.ordm[0][i,j]  
                #    ordm[1][i,j] = self.ordm[1][i,j]  
                for k in range(dim):
                    for l in range(dim):
                        if i in orblist and j in orblist and k in orblist and l in orblist:
                            trdm[0][i , k , j , l] = self.trdm[0][i,k,j,l]
                            trdm[1][i , k , j , l] = self.trdm[1][i,k,j,l]
                            trdm[2][i , k , j , l] = self.trdm[2][i,k,j,l]
                        if i in orblist and k in orblist:
                            trdm[0][i , k , j , l] =  self.trdm[0][i , k , j , l] 
                            trdm[1][i , k , j , l] =  self.trdm[1][i , k , j , l] 
                            trdm[2][i , k , j , l] =  self.trdm[2][i , k , j , l]
        self.ordm = ordm
        self.trdm =  trdm

    def project(self, orblist):
        pass

    def set_rdm(self , ordm , trdm):
        self.ordm = copy.deepcopy(ordm)
        self.trdm = copy.deepcopy(trdm)

    def transform_rdms(self,unitary, twordm = True):
        """
        The basis to which you want to transform the rdms should be in the columns, the current basis in the rows.
        """
        dim = self.header['norbs']
        t4index = None
        if twordm:
            t4index = [np.zeros((dim,dim, dim, dim)) for i in range(3)]

            for i in range(3):
                temp = np.zeros((dim,dim, dim, dim))
                temp2 = np.zeros((dim,dim, dim, dim))
                temp3 = np.zeros((dim,dim, dim, dim))
                for p in range(0,dim):
                    for mu in range(0,dim):
                        temp[p,:,:,:] += unitary[mu,p]*self.trdm[i][mu,:,:,:]
                    for q in range(0,dim):
                        for nu in range(0,dim):
                            temp2[p,q,:,:] += unitary[nu,q]*temp[p,nu,:,:]
                        for r in range(0,dim):
                            for lam in range(0,dim):
                                temp3[p,q,r,:] += unitary[lam,r]*temp2[p,q,lam,:]
                            for s in range(0,dim):
                                for sig in range(0,dim):
                                    t4index[i][p,q,r,s] += unitary[sig,s]*temp3[p,q,r,sig]
        
        t2index = [np.dot(self.ordm[0], unitary),   np.dot(self.ordm[1], unitary)   ]
        t2index = [ np.dot(unitary.T , t2index[0]) ,  np.dot(unitary.T , t2index[1]) ]

        return t2index , t4index 

    def read_ward_rdm(self, wardrdmname):
        """
        REMARK: Wards density matrices only keep unique elements, so for rdm(pq,rs) this is when p < q and r < s (and ups are 0-(L-1) , downs are L-(2L-1)).
        """
        L = self.header['norbs'] #shorthand for notation.
        t4index = [np.zeros((L,L,L,L)) for i in range(4)] 
        regexp = r'\d+\s+\d+\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)'
        with open(wardrdmname, 'r') as file:
            data = file.read()
            listmatch = re.findall(regexp , data)
            for match in listmatch:
                match = [ int(match[i]) for i in range(4)] + [float(match[4])]
                if match[0]< L and match[1] < L and match[2]< L and match[3] < L:
                    t4index[0][match[0] % L, match[1] % L, match[2] % L, match[3]% L ] = match[4]
                if match[0]< L and match[1] >= L and match[2]< L and match[3] >= L:
                    t4index[1][match[0] % L, match[1] % L, match[2] % L, match[3]% L ] = match[4]  
                if match[0]>= L and match[1] >= L and match[2]>= L and match[3] >= L:
                    t4index[2][match[0] % L, match[1] % L, match[2] % L, match[3]% L ] = match[4]
        return t4index


def test_density():
    cifres = CIFlow_Reader('psioutputoutputfci.dat')
    cifres2 = CIFlow_Reader('testfci.dat')
    cifres.read_rdm()
    cifres2.read_rdm()
    rdmw = cifres.read_ward_rdm('mario2.txt')
    L = cifres.header['norbs']
    w = 2
    for w in range(3):
        for i in range(L):
            for j in range(i,L):
                for k in range(L):
                    for l in range(k,L):
                        #if abs(cifres.trdm[w][i,j,k,l] -rdmw[w][i,j,k,l] ) > 1e-8:
                            #print abs(cifres.trdm[w][i,j,k,l] -rdmw[w][i,j,k,l] ) 
                        if abs(cifres.trdm[w][i,j,k,l] -cifres2.trdm[w][i,j,k,l]) > 1e-9 :
                            print abs(cifres.trdm[w][i,j,k,l] - cifres2.trdm[w][i,j,k,l] ) 
                            print w, " " , i , " " , j, " " , k , " " , l , " " , cifres.trdm[w][i,j,k,l] , ' and nopar' , cifres2.trdm[w][i,j,k,l]
                    #if cifres.trdm[0][i,j,k,l] > 1e-10:
                        #print i , " " , j, " " , k , " " , l , " " , cifres.trdm[0][i,j,k,l] , ' and ward ' , rdmw[0][i,j,k,l]


def overlap_analysis():
    outputname = 'overlapanalysis.dat'
    regexham = r'\s+\((\d+,\s*\d+)\)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' #to extract the Hamiltonian. not necessary here.
    numpair = 2
    #root = './results/h2o_6-31g_docipart/'
    root = './results/beh2_sto_3g_symcomp/'
    fname = 'output_files/'

    os.chdir(root)
    #runlist = map( lambda x : x/100. , range(60 , 172 , 1) ) #h2o
    #runlist = map( lambda x : x/100. , range(80 , 278 , 1) ) #beh2
    runlist = [0.86374047590, 1.02249365060, 1.18124682530, 1.34000000000, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170 ] #BeH2 (linear)
    #runlist = [x*math.pi/180. for x in range(55,60) + [59.9]] #, 60]  ] #angle between two triangles in benzene.
    fileinfo = lambda x: (float(re.search(r'[\w\d\-]*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)\.mout[\w\d_]*output[\w_\s]+\.dat' , x).group(1)), 0 if 'loadham' in x else 1 if  'sim' in x else 2 if 'doci' in x else 3 )#if 'CISDP' in x else 3) #len(str(re.search(r'[\w\d\-]*[\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*\.out([\w])*output[\w_]+\.dat' , x).group(1)) )) #first sort by first function then by second function.
    collector = fc.File_Collector(fname , r'output[\w_\s]+.dat', sortfunction = fileinfo)
    with open(outputname , 'w') as to_file:
        #to_file.write('#R\t%s\n' %'\t'.join(['P'+str(i) for i in range(1, numpair +1) ]) )
        to_file.write('#R\t%s\t%s\n' %( "<doci_sim|fci>" , "<doci|fci>" ) )
        for index in range(0,len(collector.plotfiles) , numpair+1):
            overlaplist = []
            cifdoci = CIFlow_Reader(collector.plotfiles[index], regexp = regexham )
            for i in range(numpair):
                cifp = CIFlow_Reader(collector.plotfiles[index+i+1], regexp = regexham )
                print 'overlap with' ,collector.plotfiles[index+i+1] , 'and ' , collector.plotfiles[index]
                overlaplist.append(cifdoci.calc_overlap(cifp))

            print overlaplist
            to_file.write("%f\t%s\n" %(runlist[index/(numpair+1)], '\t'.join(map(str,overlaplist) ))  )

    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        pf.plot_data(outputname, ylabel = 'overlaps', xlabel = 'R', name = 'overlapanalyse')

def overlap_analysis2():
    outputname = 'overlapanalysisdocicisdmmin.dat'
    regexham = r'\s+\((\d+,\s*\d+)\)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' #to extract the Hamiltonian. not necessary here.
    numoverlap = 2
    #root = './results/h2o_6-31g_docipart/'
    root = './results/cisdc1mminbeh2/' 
    fname = 'output_files/'
    os.chdir(root)
    #runlist = map( lambda x : x/100. , range(60 , 172 , 1) ) #h2o
    #runlist = map( lambda x : x/100. , range(80 , 278 , 1) ) #beh2
    runlist = [0.86374047590, 0.9,1.02249365060, 1.1 , 1.18124682530, 1.2, 1.34000000000, 1.38, 1.4, 1.45, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350,2.2 ,  2.29251904820, 2.35 , 2.45127222290, 2.5 , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0 ,  3.08628492170, 3.1 , 3.2, 3.3 ]#, 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)
    runlist = [0.7, 0.86374047590, 0.9,1.02249365060, 1.1 , 1.18124682530, 1.2, 1.34000000000, 1.38, 1.4, 1.45, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350,2.2 ,  2.29251904820, 2.35 , 2.45127222290, 2.5 , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0 ,  3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)
    #runlist = np.append( np.append( np.arange(0.6 , 1.1, 0.03) , np.arange(1.2 , 2. , 0.05)) , np.arange(2.,6., 0.4) )
     
    #runlist = [ 2.5 , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0, 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)
    #runlist = [x*math.pi/180. for x in range(55,60) + [59.9]] #, 60]  ] #angle between two triangles in benzene.
    fileinfo = lambda x: (float(re.search(r'psi[-\d_\w]*([\-+]?\d+[\.,]+\d+[eEDd]?[\-+]?\d*).*output.*\.dat' , x).group(1)), 0 if ('loadham' not in x and 'mmin' not in x and 'fci' in x) else 1 if  ('sim'  not in x and 'loadham' not in x and 'doci' in x and 'no' not in x and 'mmin' not in x) else 2 if ( 'loadham'  in x and 'fci' in x) else 3 if ( ('loadham' in x or 'mmin' in x)  and 'doci' in x) else 4 if ('no' in x ) else 5) #if 'CISDP' in x else 3) #len(str(re.search(r'[\w\d\-]*[\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*\.out([\w])*output[\w_]+\.dat' , x).group(1)) )) #first sort by first function then by second function.
    collector = fc.File_Collector(fname , r'output[\w_\s]+.dat', sortfunction = fileinfo ,notsearch = '(mmin)|(cisd)' , filterf = lambda x : fileinfo(x)[0] > -1 and fileinfo(x)[0] < 10)
    with open(outputname , 'w') as to_file:
        #to_file.write('#R\t%s\n' %'\t'.join(['P'+str(i) for i in range(1, numoverlap*2 ) ]) )
        to_file.write('#R\t%s\t%s\n' %( "<linhyb_cisdsenmin|fci>", "<linhyb_cisdsenmin|fci>" ) )
        for index in range(0,len(collector.plotfiles) , numoverlap* 2):
            overlaplist = []
            ciffci = CIFlow_Reader(collector.plotfiles[index], regexp = regexham )
            cifdoci = CIFlow_Reader(collector.plotfiles[index+1], regexp = regexham )
            ciffcisim = CIFlow_Reader(collector.plotfiles[index+2], regexp = regexham )
            cifdocisim = CIFlow_Reader(collector.plotfiles[index+3], regexp = regexham )
            overlap =ciffci.calc_overlap(cifdoci)
            print 'overlap with' ,collector.plotfiles[index+1] , 'and ' , collector.plotfiles[index]  
            overlaplist.append(ciffci.calc_overlap(cifdoci))
            overlaplist.append(ciffcisim.calc_overlap(cifdocisim))
            print overlaplist
            to_file.write("%f\t%s\n" %(runlist[index/(numoverlap*2)], '\t'.join(map(str,overlaplist) ))  )

    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        pf.plot_data(outputname, ylabel = 'overlaps', xlabel = 'R', name = 'overlapanalyse')

def wavefunction_analysis():
    outputname = 'wavefunctionanalysis2.dat'
    regexham = r'\s+\((\d+,\s*\d+)\)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' #to extract the Hamiltonian.
    root = './results/beh2_docipart/'
    #root = './results/h2o_6-31g_docipart/'
    fname = 'output_files/'
    matrixelements = 'matrix_elements/psi0_sto-3g0.85.out'

    os.chdir(root)
    psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1, read_ints = False)
    runlist = map( lambda x : x/100. , range(80 , 278 , 1) ) 
    #runlist = map( lambda x : x/100. , range(60 , 172 , 1) )  #h2o
    fileinfo = lambda x: float(re.search(r'[\w\d\-]*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)\.outoutput[\w_]+\.dat' , x).group(1))
    collector = fc.File_Collector(fname , r'outputdoci.dat', sortfunction = fileinfo)
    with open(outputname , 'w') as to_file:
        to_file.write('#R\t%s\n' %'\t'.join(['P'+str(i) for i in range(0,psir.values['nalpha']+1 ) ]) )
        for index, file in enumerate(collector.plotfiles):
            cif = CIFlow_Reader(file, regexp = regexham )
            overlaplist = []
            for ex in range(0,cif.header['nup']+1 ):
                ciplist = dw.ci_excitation(cif.header['nup'] , cif.header['ndown'] , cif.header['norbs'] , [ex] , get_reference = lambda x , y , z :psir.get_hf_orbs(), pairex = True)
                if ex == 0:
                    ciplist = [det[0]+'|' + det[1] for det in ciplist[0:1] ]
                else:
                    ciplist = [det[0]+'|' + det[1]  for det in ciplist[1:]] #remove HF det which is already included in the CI(P) method.
                overlaplist.append(cif.get_selected_coef(ciplist))
            print overlaplist
            print 'conservation of probability: ' , sum(overlaplist)
            to_file.write("%f\t%s\n" %(runlist[index], '\t'.join(map(str,overlaplist) ))  )

    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        pf.plot_data(outputname, ylabel = 'partitions', xlabel = 'R', name = 'partitionanalyse')
    
def wavefunction_analysis_fci():
    outputname = 'wavefunctionanalysisfciii.dat'
    root = './results/phdsenebeccpvdzfcimminc1/'
    fname = '.'
    os.chdir(root)

    fileinfo = lambda x: float(re.search(r'[\w\d\-]*([\-+]?\d*[\.,]?\d*[eEDd]?[\-+]?\d*)outputfci[-\d\w_]+\.dat' , x).group(1))
    fileinfo = lambda x : len(x) 
    collector = fc.File_Collector(fname , r'outputfci.*.dat', sortfunction = fileinfo)

    with open(outputname , 'w') as to_file:
        to_file.write('#R\t%s\n' %'\t'.join(['sen'+str(i) for i in range(0,5, 2 ) ]) )
        for index, file in enumerate(collector.plotfiles):
            print file
            cif = CIFlow_Reader(file)
            overlaplist = []
            for ex in range(0,6 ,2):
                ciplist = dw.cimain(cif.header['nup'] , cif.header['ndown'] , cif.header['norbs'] , [[],[]],[ex ] ,  write = False)
                ciplist = [ up + '|' + down for    up , down in ciplist] 
                overlaplist.append(cif.get_selected_coef(ciplist))
            print overlaplist
            print 'conservation of probability: ' , sum(overlaplist)
            to_file.write("%s\t%s\n" %(file , '\t'.join(map(str,overlaplist) ))  )

    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        pf.plot_data(outputname, ylabel = 'partitions', xlabel = 'R', name = 'partitionanalyse')

def wavefunction_analysis_fcibeh2():
    outputname = 'wavefunctionanalysisfcihmmin.dat'
    root = './results/phdsenbeh26-31gmminwfc1/'
    fname = 'output_files'
    fname = '.'
    os.chdir(root)

    fileinfo = lambda x: float(re.search(r'[\w\d\-]*([\-+]?\d+[\.,]+\d+[eEDd]?[\-+]?\d*)[\d\w_]*outputfci[-\d\w_]*\.dat' , x).group(1))
    fileinfo2 = lambda x : len(x)
    sortf = lambda x: ( fileinfo(x) , fileinfo2(x) )
    #fileinfo = lambda x : len(x)
    collector = fc.File_Collector(fname , r'outputfci.*.dat', sortfunction = sortf)

    numsorts = 2
    
    with open(outputname , 'w') as to_file:
        to_file.write('#R\t%s\n' %'\t'.join(['sen'+str(i) for i in range(0,5, 2 ) ]) )
        #for index, file in enumerate(collector.plotfiles):
        #for index in range(0, len(collector.plotfiles), 3 ) :
        for index in range(0, len(collector.plotfiles), numsorts ) :
            ciflist = [CIFlow_Reader(collector.plotfiles[index+i])  for i in range(numsorts)  ] 
            overlaplist = []
            for i in range(numsorts):
                print 'work with ' , collector.plotfiles[index+i]
                print fileinfo(collector.plotfiles[index+i]) ,fileinfo2(collector.plotfiles[index+i]) 
            for cif in ciflist:
                for ex in range(0,8 ,2):
                    ciplist = dw.cimain(cif.header['nup'] , cif.header['ndown'] , cif.header['norbs'] , [[],[]],[ex ] ,  write = False)
                    ciplist = [ up + '|' + down for    up , down in ciplist] 
                    overlaplist.append(cif.get_selected_coef(ciplist))
            print overlaplist
            print 'conservation of probability: ' , sum(overlaplist)
            to_file.write("%s\t%s\n" %(fileinfo(collector.plotfiles[index]) , '\t'.join(map(str,overlaplist) ))  )

    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        pf.plot_data(outputname, ylabel = 'partitions', xlabel = 'R', name = 'partitionanalyse')

def main():
    """
    Use this function to try out some new functionality of the above classes.
    """
    regexham = r'\s+\((\d+,\s*\d+)\)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' #to extract the Hamiltonian.
    root = '.'
    #fname = 'output_files/'
    ciffci = CIFlow_Reader('testfci.dat', regexp = regexham , read_ham= True)
    ciffcipar = CIFlow_Reader( 'psi0_output10outputfci.dat', regexp = regexham , read_ham = True)
    #print ciffci.calc_overlap(cifdoci)
    #print e.get_groundstate('00000000000011|00000000000011') 

    psir = rp.PsiReader('psi0_output10.dat', isbig = False, numorbs = -1 , read_ints = False)

    detlist = dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [range(1,psir.values['nalpha']+psir.values['nbeta']), []], [] , fname = 'determinants.dat' ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0, write = False) #CISDDOCI
    count = 0
    for det in detlist:
        for det2 in detlist:
            #+ because the eigenvectors have already a different phasefactor of 1.
            if abs(ciffci.get_mat_element(det[0]+'|'+det[1], det2[0]+'|'+det2[1]) - ciffcipar.get_mat_element(det[0]+'|'+det[1], det2[0]+'|'+det2[1]) ) > 1e-10 :
                print 'difference in hamiltonian row: ' , det[0]+'|'+det[1] , " col: " , det2[0]+'|'+det2[1] , 'fci: ', ciffci.get_mat_element(det[0]+'|'+det[1], det2[0]+'|'+det2[1]) , 'fciaddres: ' , ciffcipar.get_mat_element(det[0]+'|'+det[1],det2[0]+'|'+det2[1]) 
                count += 1
    print 'There were ' , count , ' different elements'


def mulliken():
    rootdir = './results/mullikenbeh2rangemulliken/output_files'
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-_\w\d]*\.[m]?dat' , x).group(1))
    collector = fc.File_Collector(rootdir , r'psi0_[-\w\d.]*outputfci[\w_]*.dat', sortfunction =fileinfo)
    #rootdir = './results/noconstraineddmhermextrema100angstromloop7/output_files'
    #fileinfo = lambda x: float(re.search(r'run[\=]*([-]?\d+\.\d+[e\-]*\d*)outputfci\.dat' , x).group(1))
    #collector = fc.File_Collector(rootdir , r'Constrained_DM.+run=([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)outputfci\.dat', sortfunction = fileinfo)

    rootdirmat = './results/mullikenbeh2rangemulliken/matrixelements'
    #fileinfomat = lambda x: float(re.search(r'run[\=]*([-]?\d+\.\d+[e\-]*\d*)\.mod' , x).group(1))
    fileinfo2 = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-_\w\d]*\.[m]?out' , x).group(1))
    collectormat = fc.File_Collector(rootdirmat , r'psi0_[-\d\w]+([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)\.mout', sortfunction = fileinfo2)

    with open(os.path.join(os.path.dirname(rootdirmat ) ,  "mullikendata.dat") , "w" ) as file:
        file.write("#Lambda\tFCI_E\tMulliken_A\tMulliken_rest\n")
        for num , outfile in enumerate(collector.plotfiles):
        #for outfile in plotfiles:
            #print 'starting with: ' , outfile
            ciffci = CIFlow_Reader( outfile)
            ciffci.read_rdm(twordm=False)
            
            print 'matrixelements : ' , collectormat.plotfiles[num]
            probinfo = rp.PsiReader(collectormat.plotfiles[num], isbig = False, numorbs = -1 , read_ints = False)
            unitary2 = probinfo.unit
            overlap2 = probinfo.overlap

            unitary = np.zeros((len(ciffci.ordm), len(ciffci.ordm)))
            overlap = np.zeros((len(ciffci.ordm), len(ciffci.ordm)))

            start = 0
            for i in unitary2:
                if i != None:
                    end = len(i)
                    #print start, end , i
                    #print np.array(i)
                    unitary[start:start+end ,start:start+end] = i
                    start += end

            start = 0
            for i in overlap2:
                if i != None:
                    end = len(i)
                    overlap[start:start+end ,start:start+end] = i
                    start += end

            #print ciffci.ordm
            t2index , t4index = ciffci.transform_rdms(unitary, twordm = False)
            ciffci.ordm = t2index
            #print 'One rdm in MO basis (after backtransforming to check the seniority non-zero contribution)' , ciffci.ordm 

            #t2index , t4index = cifdoci.transform_rdms(unitary.T, twordm = False) #for psi4 unitaries which save the 
            #print unitary
            #print overlap
            density = np.dot(t2index, overlap)
            mulliken_n = 2*sum([density[i,i] for i in range(9)] ) 
            mulliken_o = 2*sum([density[i,i] for i in range(9,13)] ) 
            print 'outputfile: ', outfile, 'mulliken Be : ' , mulliken_n  , 'mulliken_h2 : '  , mulliken_o , ' energy: ' , ciffci.header['energy']
            file.write("%f\t%f\t%f\t%f\n" % (fileinfo(outfile), ciffci.header['energy'] , mulliken_n , mulliken_o  ) )

def back_transform():
    rootdir = './results/hechangerep/'
    rootdir = './results/hechangerepstar631g/'
    rootdir = './results/hechangerepstar631g/'
    collector = fc.File_Collector(rootdir , r'[-\w\d.]*output[\w_]*.dat', notsearch = 'rdm_ao' , sortfunction = None)
    for outfile in collector.plotfiles:
        cifread = CIFlow_Reader(outfile)
        cifread.read_rdm(twordm=False)
        #print cifread.ordm
        #print len(cifread.ordm)
        unitary = cifread.read_matrix(os.path.join(rootdir , os.path.join('matrixelements', 'unitary-mounit_ao_to_mo.txt') ) )
        #print unitary
        #print len(unitary)
        t2index , t4index = cifread.transform_rdms(unitary , twordm = False)
        np.savetxt(outfile +'rdm_ao.dat' , t2index)
        
def test_transform():
    cifread = CIFlow_Reader('psioutputoutputfci.dat')
    cifread.read_rdm(twordm= True)
    print cifread.trdm[1][6][5][6][5]
    print cifread.trdm[1][1][2][1][1]
    print cifread.trdm[1][4][4][4][4]
    #print cifread.ordm
    #print len(cifread.ordm)
    unitary = cifread.read_matrix('unitaryao_to_mo.txt' )
    #print unitary
    #print len(unitary)
    t2index , t4index = cifread.transform_rdms(unitary , twordm = True)
    np.savetxt('rdm_ao.dat' , t2index)
    print t4index[1][6][5][6][5]
    print t4index[1][1][2][1][1]
    print t4index[1][4][4][4][4]

def test_max_det():
    #cifread = CIFlow_Reader('psioutputoutputfci.dat')
    cifread = CIFlow_Reader('hamatomicintegralsorthonoutputfci.dat')
    maxdet = cifread.get_max_det()
    print maxdet

def extract_properties():
    #rootdir = './results/shannonebeh2smallvirtual6/output_files'
    rootdir = './results/nisystemnoconstrainedwithwfnatom/'
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)\.[m]?dat[\d]*' , x).group(1))
    collector = fc.File_Collector(rootdir , r'[-\w\d.]*outputfci[\d\w_\.]+.dat[\d]*', notsearch = 'rdm_ao' , sortfunction = fileinfo)
    data = []
    for outfile in collector.plotfiles:
        cifread = CIFlow_Reader(outfile)
        cifread.get_properties()
        #print cifread.properties
        data.append(cifread.properties['parenergy'])
    os.chdir(rootdir)
    with open('propertiesfci.dat', 'w') as file:
        file.write('\n'.join(map(str,data) )  )

def energy_rdm():
    #nameham = "hamnoplussto-3gpatrick10.0new.out"
    #namewf = "hamnoplussto-3gpatrick10.0newoutputfci.dat"
    nameham = "5psioutput.dat"
    namewf =  "5psioutputoutputfci.dat"
    readermo = rp.PsiReader(nameham, read_ints = True)
    namehamn = "5natomhammogost.dat"
    namehamo = "5oatomhammogost.dat"
    #nameham = "ni_system.dat"
    #namewf = "ni_systemoutputfci.dat"
    rootdir = '.'

    readern = rp.PsiReader(namehamn, read_ints = True)
    readero = rp.PsiReader(namehamo, read_ints = True)
    cifread = CIFlow_Reader(namewf)

    cifread.read_rdm()
    cifread.print_rdm()

    print 'Before manipulation rdm trace(ordm) : ' , sum(cifread.ordm[0].diagonal() + cifread.ordm[1].diagonal())

    e_n = readermo.calc_energy( cifread.ordm , cifread.trdm )
    print ' e totaal:   ',  e_n

    import numpy
    unitary = readermo.get_unitary().T
    t2index , t4index = cifread.transform_rdms(unitary , twordm = True)

    cifread.set_rdm(t2index  , t4index  )
    print 'after manipulation 1 rdm trace(ordm) in atom without overlap: ' , sum(cifread.ordm[0].diagonal() + cifread.ordm[1].diagonal())
    cifread.project_herm([0,1,2,3,4])
    t2index2 , t4index2 = cifread.transform_rdms( numpy.linalg.inv(readern.get_unitary().T) , twordm = True)
    cifread.set_rdm(t2index2 , t4index2)
    print 'after projection  manipulation and transformation to mo trace(ordm) : ' , sum(cifread.ordm[0].diagonal() + cifread.ordm[1].diagonal())
    
    cifread.print_rdm(filename = 'projecteddensitymatrix.dat')


    e_n = readern.calc_energy( cifread.ordm , cifread.trdm )
    print ' e totaal:   ',  e_n

    #cifread2 = CIFlow_Reader(namewf)
    #cifread2.read_rdm()

    #L = cifread2.header['norbs']
    #for w in range(3):
    #    for i in range(L):
    #        for j in range(i,L):
    #            for k in range(L):
    #                for l in range(k,L):
    #                    if abs(cifread2.trdm[w][i,j,k,l] -cifread.trdm[w][i,j,k,l]) > 1e-9 :
    #                        print abs(cifread.trdm[w][i,j,k,l] - cifread2.trdm[w][i,j,k,l] ) 
    #                        print w, " " , i , " " , j, " " , k , " " , l , " " , cifread.trdm[w][i,j,k,l] , ' and nopar' , cifread2.trdm[w][i,j,k,l]

    #print 'after manipulation rdm trace(ordm) : ' , sum(cifread.ordm[0].diagonal() + cifread.ordm[1].diagonal())

    #t2index , t4index = cifread.transform_rdms(unitary , twordm = True)
     
    unitary = numpy.linalg.inv(readern.get_unitary() )
    t2 , t4 =  readern.transform_integrals(unitary)
    readern.matrix_to_list(t2 , t4)

    unitary = numpy.linalg.inv(readero.get_unitary() )
    t2 , t4 =  readero.transform_integrals(unitary)
    readero.matrix_to_list(t2 , t4)

    e_nn = readern.calc_energy( t2index , t4index , orblist = [0,1,2,3,4]  , nucrepbool = False)
    e_no = readero.calc_energy( t2index , t4index , orblist = [5,6,7,8,9]  , nucrepbool = False, startrdm = 5)
    print 'e n atom : ' , e_nn , ' e o atom : ' , e_no , ' e interactie :   ',  e_n - e_nn - e_no , ' e all : ' , e_n

    #reader.create_output(fname = 'newoutputao.dat')

def energy_atom():
    rootdir = './results/5bohrnoplusconstrainednatomdd/output_files/'
    rootdir = './results/5bohrnoplusconstrainednatomddmostartplusdiisoffpsi/output_files/'
    rootdir = './results/eqbohrnopluseqpsioutputdatfciconstrainedeqbohr/output_files/'
    rootdir = './results/4bohrnoplus4psioutputdatfciconstrained4bohr/output_files/'
    namemoham = './4psioutput.dat' 
    #namemoham = './hamnoplussto-3gpatrick6.0new.out' 
    readermo = rp.PsiReader(namemoham, read_ints = True) #for the  unitary transformation
    nameham = "4natomhammogost.dat"
    namehamo = "4oatomhammogost.dat"

    #nameconstrained = './results/10bohrnoplusconstrainednatomddmostartplusdiisoffpsi/noconstrainedm_.dat'
    #constrained_data = np.loadtxt(nameconstrained)
    #print constrained_data
    #rootdir = './results/5bohrnoplusnatom/'
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)\.[m]?dat[\d]*' , x).group(1))
    collector = fc.File_Collector(rootdir , r'[-\w\d.]*outputfci[\d\w_\.]+.dat[\d]*', notsearch = 'rdm_ao' , sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= 0 and fileinfo(x) < 10000)

    data = []
    for row , outfile in enumerate(collector.plotfiles ):
        #print outfile , 'we are at', fileinfo(outfile)
        #print constrained_data[row , 0]

        reader = rp.PsiReader(nameham, read_ints = True)
        readero = rp.PsiReader(namehamo, read_ints = True)
        #reader.create_constrained_vals(constrained_data[row , 3], constrained_data[row , 4], [0,1,2,3,4])
        cifread = CIFlow_Reader(outfile)
        cifread.read_rdm()

        import numpy
        unitary = readermo.get_unitary().T
        #unitary = numpy.linalg.inv(reader.get_unitary() ).T
        #unitary = rp.generate_random_unitary(reader.values['norb'])
        print unitary
        t2index , t4index = cifread.transform_rdms(unitary , twordm = True)

        #hartree , cum = cifread.construct_cum(t2index , t4index )

        e_tot = readermo.calc_energy( cifread.ordm , cifread.trdm   )
        print 'etot ' , e_tot
         

        #np.savetxt('rdm_ao.dat' , t2index)
        unitary = numpy.linalg.inv(reader.get_unitary() )
        t2 , t4 =  reader.transform_integrals(unitary)
        reader.matrix_to_list(t2 , t4)
        unitaryo = numpy.linalg.inv(readero.get_unitary() )
        t2o , t4o =  readero.transform_integrals(unitaryo)
        readero.matrix_to_list(t2o , t4o)
        #e_all = reader.calc_energy( t2index , t4index , orblist = [0,1,2,3,4,5,6,7,8,9]  )
        #e_nn = reader.calc_atom_e( t2index , hartree , cum, orblist = [0,1,2,3,4])
        #e_no = readero.calc_atom_e( t2index ,hartree , cum , orblist = [5,6,7,8,9])
        e_nn = reader.calc_energy( t2index , t4index, orblist = [0,1,2,3,4])
        e_no = readero.calc_energy( t2index ,t4index, orblist = [5,6,7,8,9])
        #e_nn = reader.calc_energy( t2index , t4index , orblist = [0,1,2,3,4] ,nucrepbool = False  )
        #e_no = reader.calc_energy( t2index , t4index , orblist = [5,6,7,8,9] ,nucrepbool = False  )
        #print 'e n atom : ' , e_all - -73.4435839702339  - (e_all - e_nn - e_no )/2. , ' e o atom : ' , e_no , ' e interactie :   ',  e_all - e_nn - e_no
        print 'e n check: ' , e_nn  , 'e o : check' , e_no , ' e int: ' , e_tot - e_nn - e_no,  'e tot check '  , e_tot 
        data.append( (fileinfo(outfile) , e_nn , e_no , e_tot-e_nn-e_no ,    e_tot )  )
    with open(os.path.join(rootdir, 'n_atom_e2ghost5.dat') , 'w' ) as file:
        file.write( '\n'.join([ str(dat[0]) + '\t' + str(dat[1]) + '\t' + str(dat[2]) + '\t' + str(dat[3]) + '\t' + str(dat[4])  for dat in data   ]  ) )

def energy_decomp():
    rootdir = './results/nisystemnoconstrainedwithwfnatom/'
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)\.[m]?dat[\d]*' , x).group(1))
    collector = fc.File_Collector(rootdir , r'[-\w\d.]*outputfci[\d\w_\.]+.dat[\d]*', notsearch = 'rdm_ao' , sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= 0 and fileinfo(x) < 10000)
    data = []
    reader = rp.PsiReader("ni_system.dat", read_ints = True)
    #print np.array(reader.unit[0]).T

    for outfile in collector.plotfiles:
        cifread = CIFlow_Reader(outfile)
        cifread.read_rdm()
        #print cifread.ordm
        print cifread.trdm[0][0 , 0, 0, 0]
        e_n = reader.calc_energy( cifread.ordm , cifread.trdm   )
        print outfile ,'  ' , e_n
        data.append( (fileinfo(outfile),  e_n ) )

    os.chdir(rootdir)
    with open('energiedecompfci.dat', 'w') as file:
        file.write('\n'.join( [ str(tup[0]) + '\t' +  str(tup[1])  for tup in data  ]  )  )


if __name__ == "__main__":
    #main()
    #wavefunction_analysis()
    #wavefunction_analysis_fci()
    #wavefunction_analysis_fcibeh2()
    #overlap_analysis()
    #overlap_analysis2()
    #mulliken()
    #test_transform()
    #test_density()
    #back_transform()
    #test_max_det()
    #extract_properties()
    #energy_decomp()
    #energy_rdm()
    energy_atom()
