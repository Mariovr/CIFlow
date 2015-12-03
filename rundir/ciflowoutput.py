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
if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    import plotfunctions as pf

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
        self.ordm = np.zeros((self.header['norbs'] ,self.header['norbs']) , dtype=np.float )
        if twordm:
            self.trdm = [np.zeros((self.header['norbs'] ,self.header['norbs'], self.header['norbs'] ,self.header['norbs']) , dtype=np.float ) for i in range(3) ] 
        reg1rdm = r'^(\d+)\s(\d+)\s([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)'  #onerdm
        reg2rdm = r'^(\d+)\s(\d+)\s(\d+)\s(\d+)\s([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)'  #twordm
        startrdm = 0 
        with open(self.filename,"r") as file:
            for line in file:
                if 'THE 1RDM' in line:
                    startrdm = 1
                elif 'rdm case' in line:
                    startrdm += 1
                if startrdm == 1:
                    match = re.search(reg1rdm , line)
                    if match:
                        self.ordm[int(match.group(1)) , int(match.group(2)) ] = float(match.group(3))
                        self.ordm[int(match.group(2)), int(match.group(1)) ] = float(match.group(3))
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

    def get_max_det(self):
        max_val = ('' , 0.)
        for i, val in enumerate( self.groundstate):
            if val > max_val[1]:
                max_val = (i ,val)
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
        
        t2index = np.dot(self.ordm, unitary)
        t2index = np.dot(unitary.T , t2index)

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
    outputname = 'wavefunctionanalysisfcino.dat'
    regexham = r'\s+\((\d+,\s*\d+)\)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' #to extract the Hamiltonian.
    root = './results/buttestfci/'
    fname = 'output_files/'
    matrixelements = 'matrixelements/psi0_cc-pvdz2.40.mout'

    os.chdir(root)
    psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1, read_ints = False)
    runlist = [2.5]
    import numpy as np
    runlist = list(np.arange(1,2.5,0.2)) + [2.5] + list(np.arange(2.6,4,0.2))
    #runlist = map( lambda x : x/100. , range(60 , 172 , 1) )  #h2o
    fileinfo = lambda x: float(re.search(r'[\w\d\-]*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)output[\w_]+\.dat' , x).group(1))
    collector = fc.File_Collector(fname , r'outputfcino.dat', sortfunction = fileinfo)
    with open(outputname , 'w') as to_file:
        to_file.write('#R\t%s\n' %'\t'.join(['sen'+str(i) for i in range(0,psir.values['nalpha']*2+1,2 ) ]) )
        for index, file in enumerate(collector.plotfiles):
            cif = CIFlow_Reader(file, regexp = regexham )
            overlaplist = []
            for ex in range(0,cif.header['nup']*2+1 ,2):
                ciplist = dw.cimain(cif.header['nup'] , cif.header['ndown'] , cif.header['norbs'] , [[],[]], [ex] ,ref = lambda x , y , z :psir.get_hf_orbs(),  write = False)
                ciplist = [det[0]+'|' + det[1]  for det in ciplist] #remove HF det which is already included in the CI(P) method.
                overlaplist.append(cif.get_selected_coef(ciplist))
            print overlaplist
            print 'conservation of probability: ' , sum(overlaplist)
            to_file.write("%f\t%s\n" %(runlist[index], '\t'.join(map(str,overlaplist) ))  )

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
    #rootdir = './results/nopluspatrickham/'
    #collector = fc.File_Collector(rootdir , r'hamnoplus[-\w\d.]*output[\w_]*.dat', sortfunction = None)
    #collector = fc.File_Collector(rootdir , r'psi[-\w\d.]*output[\w_]*.dat', sortfunction = None)
    plotfiles = ['hamatomicintegralsorthonoutputfci.dat']
    rootdir = '.'
    #for outfile in collector.plotfiles:
    for outfile in plotfiles:
        #print 'starting with: ' , outfile
        cifdoci = CIFlow_Reader( outfile)
        cifdoci.read_rdm(twordm=False)
        #print cifdoci.ordm
        if 'fmmin' and 'sim' in outfile: #Transform the RDMs back from the OO basis to the MO basis.
            if 'sim' in outfile:
                unitary2= cifdoci.read_matrix(os.path.join(rootdir, 'unitarysim.dat') ) 
                unitary2 = unitary2.T #Because unitarysim.dat is generated with the read_psi unitary reader (of .h5 file to .dat file ) which transposes implicitly and saves new basis in columns. Instead of the rows such as UnitaryMatrix.
            else:
                unitary2= cifdoci.read_matrix(os.path.join(rootdir, 'unitary_hamnoplussto-3gpatrick100.0permbitFCI.dat') )
            #print unitary2
            t2index , t4index = cifdoci.transform_rdms(unitary2, twordm = False)
            cifdoci.ordm = t2index
            #print 'One rdm in MO basis (after backtransforming to check the seniority non-zero contribution)' , cifdoci.ordm 

        #unitary = cifdoci.read_matrix(os.path.join(rootdir, 'no_100.0_PB.mo') )
        unitary = cifdoci.read_matrix(os.path.join(rootdir, 'unitary.txt') )

        t2index , t4index = cifdoci.transform_rdms(unitary, twordm = False)
        #overlap = cifdoci.read_matrix(os.path.join(rootdir , 'no_100.0_PB.ove') )
        #unitary = cifdoci.read_matrix(os.path.join(rootdir, 'unitary.txt') ) #MO is row index, AO is column index
        overlap = cifdoci.read_matrix(os.path.join(rootdir , 'overlap.txt') )
        #t2index , t4index = cifdoci.transform_rdms(unitary.T, twordm = False) #for psi4 unitaries which save the 
        #print unitary
        #print overlap
        density = np.dot(t2index, overlap)
        mulliken_n = 2*sum([density[i,i] for i in range(5)] ) 
        mulliken_o = 2*sum([density[i,i] for i in range(5,10)] ) 
        print 'outputfile: ', outfile, 'mulliken n : ' , mulliken_n  , 'mulliken_o : '  , mulliken_o 

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

if __name__ == "__main__":
    #main()
    #wavefunction_analysis()
    #wavefunction_analysis_fci()
    #overlap_analysis()
    #overlap_analysis2()
    #mulliken()
    #test_transform()
    #test_density()
    #back_transform()
    test_max_det()

