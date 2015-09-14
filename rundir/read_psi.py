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
#! /usr/bin/env python 
from string import maketrans
import re
import numpy as np
import sys , os
import mmap
import copy
import shutil
from math import sqrt

#sys.path.append('/home/mario/Dropbox/Doctoraat/CIFlow/rundir')
import detwrite as dw

if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    import h5py

def deleteFromMmap(filename , start,end):
    """
    deletes the bytes between start and end from the file with name filename using mmap (in memory copy of file).
    """
    size = os.stat(filename).st_size
    f = open(filename, 'r+') #r+ opens the file for both reading and writing
    mmapob = mmap.mmap(f.fileno(), size )#, access=mmap.ACCESS_READ)
    length = end - start
    size = len(mmapob)
    newsize = size - length

    mmapob.move(start,end,size-end)
    mmapob.flush()
    mmapob.close()
    f.truncate(newsize)
    f.close()

if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    #WARNING implicitely transposes, SAVES new basis in colums
    class Unitary_Matrix(object):
        def __init__(self,fname):
            self.f = h5py.File(fname,"r")
            self.unitary = [None]*8 #biggest point group in psi4 is D2h, which has 8 irreps.
            self.numorb = 0
            for name in self.f['/Data']:
                irrep_num= int(re.search(r'(\d+)',name).group())
                linsize = sqrt(len(self.f['/Data/' + name][:]))
                self.numorb += linsize
                self.unitary[irrep_num] = np.array(self.f['/Data/' + name][:] ) #transforms implicitly because unitary matrix saves first columns (fast changing index), but when we use numpy.shapeit puts data in rows and then cuts. This is oke, because unitary matrix saves the eigenbasis where we need to transform to in the rows, and this python class will keep them than in the columns. This is consistent with new eigenbasis obtained by diagonalizing a 2 index matrix.
                self.unitary[irrep_num].shape = (linsize , linsize)
                
        def get_unitary(self):
            #At this moment we haven't supported point group symmetry
            unit = np.zeros((self.numorb , self.numorb))
            pos = 0
            for irrep in range(len(self.unitary)):
                if self.unitary[irrep] != None:
                    dim = self.unitary[irrep].shape[0] 
                    unit[pos:pos+dim , pos:pos+dim] = self.unitary[irrep][:,:]
                    pos += dim
            return unit        

        def print_unitary(self):
            for irrep in range(len(self.unitary)):
                if self.unitary[irrep] != None:
                    print 'irrep: ' , irrep
                    print self.unitary[irrep]

        def print_structure(self):
            def p(x):
                print x
            self.f.visit(p)



class CIFlow_Input(object):
    def __init__(self, nalpha , nbeta , nucnuc , norb , sym , nirrep , irreplist , oeilist, teilist , hfenergy, docclist , socclist):
        self.nalpha = nalpha 
        self.nbeta = nbeta 
        self.nucnuc = nucnuc 
        self.norb = norb 
        self.sym = sym
        self.nirrep = nirrep
        self.irreplist = irreplist
        self.oeilist = oeilist
        self.teilist = teilist
        self.hfenergy = hfenergy
        self.docclist = docclist
        self.socclist = socclist
        self.create_strings()

    def create_strings(self):
        self.irrepstring = ' '.join(map(str,self.irreplist))
        self.doccstring = ' '.join(map(str,self.docclist))
        self.soccstring = ' '.join(map(str,self.socclist))
        self.oeistring = '\n'.join(['%d %d %.25f '%(i,j,val)      for i, j ,val in self.oeilist])
        self.teistring = '\n'.join(['%d %d %d %d %.25f '%(i,j,k,l,val) for i, j ,k,l,val in self.teilist])

    def __str__(self):
        text = """****  Molecular Integrals For DOCI Start Here 
Nalpha = %(nalpha)d 
Nbeta = %(nbeta)d 
Symmetry Label = %(sym)s 
Nirreps = %(nirrep)d 
Nuclear Repulsion Energy = %(nucnuc).20f 
Number Of Molecular Orbitals =  %(norb)d 
Irreps Of Molecular Orbitals = 
%(irrepstring)s 
DOCC =  %(doccstring)s #this line is ignored
SOCC =  %(soccstring)s #this line is ignored 
****  MO OEI 
%(oeistring)s
****  MO TEI 
%(teistring)s
****  HF Energy = %(hfenergy).20f 
****  Molecular Integrals For DOCI End Here 
"""%(self.__dict__)
        return text

    def write_file(self,fname = 'output.dat'):
        with open(fname , 'w') as mofile:
            mofile.write(str(self))

class Reader(object):
    """
    class to read in data files 
    """
    def __init__(self, filename ):
        self.filename = filename
        self.read()

    def read(self):
        size = os.stat(self.filename).st_size
        f = open(self.filename)
        self.data = mmap.mmap(f.fileno(), size, access=mmap.ACCESS_READ)
        f.close()

    def process_element(self, regexp , text , func):
        return func(re.search(regexp,text).group(1))

    def process_list(self, regexp , text , func):
        return func(re.findall(regexp,text))

    def maketranslation(self,first = "dD",after = "ee"): 
        #use on a string: d = "How would you do that", d.translate(transtable) , print d -> "How woule you eo that"
        return maketrans(first, after)

class PsiReader(Reader):
    def __init__(self, filename , isbig = False , numorbs = -1 , read_ints = True, valdict = {}, intdict = {}):
        self.values = valdict
        self.ints= intdict
        self.transtables = {'f2c': self.maketranslation()}
        if filename:
            super(PsiReader,self).__init__(filename)
            self.regexps = {'nalpha': (r'Nalpha\s*=\s*(\d+)', int) , 'nbeta' : (r'Nbeta\s*=\s*(\d+)', int), 'norb' : (r'Number Of \w+ Orbitals\s*=\s*(\d+)', int) , 'nirrep' : (r'Nirreps\s*=\s*(\d+)' ,int ) , 'sym' : (r'(?m)^Symmetry\s*Label\s*=\s*([\d\w]*)' , lambda x: x), 'nucrep' : (r'(?m)^Nuclear Repulsion Energy =\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', float), 'irreplist' : (r'Irreps Of \w+ Orbitals\s*=\s*\n([\d\s]+)\n', lambda x : map(int,x.split())), 'hfenergy' : ( '(?m)^\*+\s*HF\s*Energy\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', float), 'DOCC' : ( r'(?m)^DOCC\s+=\s+([\d\s]+)', lambda x: map(int , x.split() ) ) , 'SOCC' : ( r'(?m)^SOCC\s+=\s+([\d\s]+)', lambda x : map(int, x.split() ) ) }

            self.regint = { 'mo4index' : (r'(?im)^(\d+)\s(\d+)\s(\d+)\s(\d+)\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', self.tei_func) , 'mo2index' : (r'(?im)^(\d+)\s*(\d+)\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', self.oei_func )}

            #self.create_star_data()
            self.extract_values()
            print self.values

            if isbig:
                self.delete_lines(numorbs) #numorbs is the number of orbitals you wanna keep
                self.values['irreplist'] = self.values['irreplist'][:numorbs]
                self.values['norb'] = numorbs

            if read_ints:
                self.extract_ints()
            self.data.close()
            
    def extract_values(self):
        for key, value in self.regexps.iteritems():
            try:
                self.values[key] = self.process_element(value[0] , self.data , value[1])
            except:
                print 'We couldnt extract : ' , key
                print 'we continue'

    def extract_ints(self):
        for key, value in self.regint.iteritems():
           self.ints[key] = self.process_list(value[0] , self.data , value[1])

    def extract_energies(self):       
        self.energies = self.process_list(r'(\w+)\s*total\s*energy\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)' , self.data,lambda x: [(element[0] + '_%d' %index, element[1]) for index , element in enumerate(x) ])
        return self.energies

    def create_star_data(self):
        """
        This cuts away all the extra psi output, so only the relevant part for ciflow remains.
        """
        deleteFromMmap(self.filename , 0 , re.search(r'(?im)[*]+\s*Molecular[\s\w]+Start',self.data).start(0))
        #deleteFromMmap(self.filename , re.search(r'(?im)[*]+\s*Molecular[\s\w]+End' , self.data).end(0),  os.stat(self.filename).st_size) #to cut the end to but mostly this has little use because mointegrals is mostly at the end of the psi output
        self.data.close()
        self.read()
        #self.data = self.data[re.search(r'(?im)[*]+\s*Molecular[\s\w]+Start',self.data).start(0) : re.search(r'(?im)[*]+\s*Molecular[\s\w]+End' , self.data).end(0)]
    
    def add_index(self,num ):
        """
        Adds the value num to all the indices of the electron repulsion integrals and the one electron integrals.
        """
        for regname in self.ints.keys():
            num_indices = int(re.search(r'\d+',regname).group())
            for index , line in enumerate(self.ints[regname]):
                self.ints[regname][index] = map(lambda x : x+num , line[0:num_indices]) + [line[num_indices]]

    def swap_orbitals(self,orb1 , orb2):
        """
        Swaps the index of orb1 to orb2 and vice versa.
        We need to run three times over all the integrals to avoid overwriting.
        REMARK: WATCH OUT FOR ORBITALTRANSFORMATIONS AFTER SWAPPING THE ORBITALS IF SYMMETRY IS KEPT.
        because current code in ciflow depends on the fact that the orbitals are first ordered by symmetrysector.
        If you swap orbitals from different symmetry sectors this will not be the case.
        """
        #we change orb1 to a new index, not used for a orbital
        for regname in self.ints.keys():
            num_indices = int(re.search(r'\d+',regname).group())
            for index in range(len(self.ints[regname])):
                for i in range(num_indices):
                    if self.ints[regname][index][i] == orb1:
                        self.ints[regname][index][i] = orb1 + self.values['norb']

        for regname in self.ints.keys():
            num_indices = int(re.search(r'\d+',regname).group())
            for index in  range(len(self.ints[regname])):
                for i in range(num_indices):
                    if self.ints[regname][index][i] == orb2:
                        self.ints[regname][index][i] = orb1 

        for regname in self.ints.keys():
            num_indices = int(re.search(r'\d+',regname).group())
            for index in range(len(self.ints[regname])):
                for i in range(num_indices):
                    if self.ints[regname][index][i] == orb1+self.values['norb']:
                        self.ints[regname][index][i] = orb2

        save = self.values['irreplist'][orb2]  
        self.values['irreplist'][orb2] = self.values['irreplist'][orb1] 
        self.values['irreplist'][orb1] = save

    def list_to_matrix(self):
        dim = self.values['norb']
        temp4in = np.zeros((dim,dim, dim, dim))
        temp2in = np.zeros((dim,dim))
        for index in range(len(self.ints['mo2index'])):
            temp2in[self.ints['mo2index'][index][0] , self.ints['mo2index'][index][1]] = self.ints['mo2index'][index][2]

        for index in range(len(self.ints['mo4index'])):
            temp4in[self.ints['mo4index'][index][0] , self.ints['mo4index'][index][1],self.ints['mo4index'][index][2] , self.ints['mo4index'][index][3]] = self.ints['mo4index'][index][4]

        return temp2in , temp4in

    def matrix_to_list(self,t2 , t4):
        dim = t2.shape[0] #all dimensions are equal to the number of orbitals.
        print dim
        self.ints['mo2index'] = []
        self.ints['mo4index'] = []
        for i in range(dim):
            for j in range(dim):
                if t2[i,j] != 0.:
                    self.ints['mo2index'].append( [i, j , t2[i,j] ] )

        for i in range(dim):
            for j in range(dim):
                for k in range(dim):
                    for l in range(dim):
                        if t4[i,j,k,l] != 0.:
                            self.ints['mo4index'].append ( [i, j, k, l , t4[i,j,k,l] ] )

    def transform_integrals(self,unitary):
        dim = self.values['norb']
        temp = np.zeros((dim,dim, dim, dim))
        temp2 = np.zeros((dim,dim, dim, dim))
        temp3 = np.zeros((dim,dim, dim, dim))
        t4index = np.zeros((dim,dim, dim, dim))

        index2, index4 = self.list_to_matrix()
        
        for p in range(0,dim):
            for mu in range(0,dim):
                temp[p,:,:,:] += unitary[mu,p]*index4[mu,:,:,:]
            for q in range(0,dim):
                for nu in range(0,dim):
                    temp2[p,q,:,:] += unitary[nu,q]*temp[p,nu,:,:]
                for r in range(0,dim):
                    for lam in range(0,dim):
                        temp3[p,q,r,:] += unitary[lam,r]*temp2[p,q,lam,:]
                    for s in range(0,dim):
                        for sig in range(0,dim):
                            t4index[p,q,r,s] += unitary[sig,s]*temp3[p,q,r,sig]
        
        t2index = np.dot(index2, unitary)
        t2index = np.dot(unitary.T , t2index)

        return t2index , t4index 


    def insert(self,extraoei , extratei):
        """
        adds molecular integrals at the end (to simulate the energy of atoms which don't interact) (as in infinte distance)
        (see main)
        """
        self.ints['mo4index'] += extratei
        self.ints['mo2index'] += extraoei

    def add_patrick_mo(self,num = -1, one = 'OEI' , two = 'TEI'):
        with open(one, 'r') as tofile:
            oeilines = tofile.read()
        with open(two, 'r') as tofile:
            teilines = tofile.read()

        oei = r'\s+(\d+)\s+(\d+)\s+([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'
        tei = r'\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'

        oeilist = re.findall(oei , oeilines)
        teilist = re.findall(tei , teilines)
        
        self.ints['mo2index'] = self.oei_func(oeilist)
        self.ints['mo4index'] = self.tei_func(teilist)
        self.add_index(num)

    def add_diego_mo(self , fname = 'h2o-psi-cisd-1_8-mario.out'): #num adds a number to the integral indices but is not necessary because diego uses psi3 to generate his integrals and it uses the same numbering as psi4.
        with open(fname, 'r') as tofile:
            oeilines = tofile.read()

        oei = r'h\((\d+)\)\((\d+)\)\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'
        tei = r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*\([\s\d]+\)\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'

        oeilist = re.findall(oei , oeilines)
        teilist = re.findall(tei , oeilines)
        self.ints['mo2index'] = self.oei_func(oeilist)
        self.ints['mo4index'] = self.tei_func(teilist)

    def oei_func(self,oeilist):
        for index,match in enumerate(oeilist):
            oeilist[index] = []
            for i in range(2):
                oeilist[index].append(int(match[i])) 
            oeilist[index].append(float(match[2].translate(self.transtables['f2c'])))
        return oeilist

    def tei_func(self,teilist):
        for index,match in enumerate(teilist):
            teilist[index] = []
            for i in range(4):
                teilist[index].append(int(match[i])) 
            teilist[index].append(float(match[4].translate(self.transtables['f2c'])))
        return teilist

    def keep_orbs(self, num):
        """
        Sets the number of orbitals equal to num.
        Deletes all integrals with index number above or equal to num (indices start with zero so we have then num orbitals).
        """
        for regname in self.ints.keys():
            num_indices = int(re.search(r'\d+',regname).group())
            for i in range(num_indices):
                self.ints[regname] = filter(lambda x : x[i] < num , self.ints[regname])
        self.values['irreplist'] = self.values['irreplist'][:num]
        self.values['norb'] = num

    def delete_lines(self,num , newfile = 'outputnew.dat'):
        """
        For big files (cant keep the integrals in memory) we need to delete the lines referring to orbitals we throw away this makes the processing easier.
        """
        with open(newfile , 'w') as newf:
            for line in iter(self.data.readline , ''):
                good = True
                for key, value in self.regint.iteritems():
                    num_indices = int(re.search(r'\d+',key).group())
                    match = re.search(value[0] , line) 
                    if match:
                        for i in range(1,num_indices+1):
                            if int(match.group(i)) >= num:
                                good = False
                                break
                if good:
                    newf.write(line)

        self.data.close()
        #self.filename = newfile
        shutil.copy(newfile,self.filename)
        self.read()

    def create_ni_system(self, outname = 'outputsc.dat'):
        oeint =list(self.ints['mo2index'])
        teint =list(self.ints['mo4index'])
        self.add_index(self.values['norb'])
        self.insert(oeint, teint)
        self.values['norb'] *= 2
        self.values['nalpha'] *= 2
        self.values['nbeta'] *= 2
        self.values['nucrep'] *= 2
        self.values['sym'] = 'c1'
        self.values['irreplist'] = [0] *self.values['norb']
        self.create_output(outname)

    def change_repulsion(self):
        self.ints['mo4index'] = [[integrals[0] ,integrals[1] ,integrals[2] , integrals[3], -1.*integrals[4] ] for integrals in self.ints['mo4index']]

    def create_output(self , fname = 'output.dat'):
        ciflow_input = CIFlow_Input(self.values['nalpha'], self.values['nbeta'], self.values['nucrep'], self.values['norb'] , self.values['sym'] , self.values['nirrep'] , self.values['irreplist'], self.ints['mo2index'], self.ints['mo4index'], self.values['hfenergy'], self.values['DOCC'] , self.values['SOCC'])
        ciflow_input.write_file(fname = fname)

    def get_hf_orbs(self, frozen = None, virtual = None):    
        #returns determinants for detwrite cimain, if frozen and virtual is set those are removed,only the orbitals of which the occupation can change in detwrite are retained.
        #those are later explicitly added by cimain in detwrite if the appropriate frozen, virtual and active lists are provided.
        docc = copy.copy(self.values['DOCC']) #shallow copy because the lists only contain integers.
        irreplist = copy.copy(self.values['irreplist'])
        norb = self.values['norb']
        if frozen != None:
            for i, numf in enumerate(frozen):
                docc[i] -= numf
                norb -= numf
                for j in range(numf):
                    irreplist.remove(i)

        if virtual != None:
            for i, numf in enumerate(virtual):
                for j in range(numf):
                    irreplist.remove(i)
                norb -= numf

        socc = copy.copy( self.values['SOCC'])
        hfs = [[0]* norb, [0] * norb]
        for index , irrep in enumerate(irreplist):
            if docc[irrep] != 0:
                hfs[0][index]  = 1
                hfs[1][index] = 1
                docc[irrep] -= 1
            if socc[irrep] != 0:
                hfs[0][index]  = 1
                socc[irrep] -= 1

        
        hfs[0].reverse() 
        hfs[1].reverse() 
        return ( ''.join(map(str,hfs[0]))  , ''.join(map(str, hfs[1]) )  )

    def set_irrep_info(self):
        self.occurences = np.bincount(np.array(self.values['irreplist']))
        self.cum = [0]
        for i in range(1,len(self.occurences)+1):
            self.cum.append(self.cum[i-1] + self.occurences[i-1])
        print self.occurences
        print self.cum

    def set_active_space(self, irreplist , exlist , senlist, detfile ):
        """
        Generates a determinantfile in a given active space determined by irreps. Remark at this moment the active space should be concatenated without gaps.
        """
        self.set_irrep_info() #to set the cumulative array of orbital indices where new irreps start.
        nup,  ndown, norb  = 0, 0 ,0
        activehf = [[],[]]
        for irrep in irreplist:
            nup += self.values['DOCC'][irrep]
            ndown += self.values['DOCC'][irrep]
            nup += self.values['SOCC'][irrep]
            norb += self.occurences[irrep]
            activehf[0] += [1]*self.values['DOCC'][irrep] + [0] * (self.occurences[irrep] - self.values['DOCC'][irrep])
            activehf[1] += [1]*self.values['DOCC'][irrep]  + [0] * (self.occurences[irrep] - self.values['DOCC'][irrep])

        activehf[0].reverse()
        activehf[1].reverse()
        activehf = ( ''.join(map(str,activehf[0]))  , ''.join(map(str, activehf[1]) )  )
            
        upstring , downstring = self.get_hf_orbs() #to determine the frozen orbitals
        frozen = upstring[-self.cum[irreplist[0]] :] #we presume that if there is an odd number of electrons, it is in the activespace.
        virtual = upstring[:-1*(self.cum[irreplist[-1] + 1 ] ) ]
        print 'frozen' , frozen , 'virtual' , virtual

        dw.cimain(nup,ndown ,norb ,exlist , senlist ,fname = detfile ,ref = [lambda x , y, z : activehf ], add_frozen = 0, frozenstring = frozen , virtualstring = virtual) 

def test_unitary():
    """
    Calculate orbital optimized unitary matrix with ciflow for a CI method, read it in with this script, read in the matrixelements, transform them with the unitary and output a new integral file.    Perform the CI calculation on the new integral file, you should obtain the same orbital optimized energy.
    """
    rootdir = './results/nopluspatrickham/'
    unitname = os.path.join(rootdir , 'unitarydocisimnopluspatrick100.h5')  #'unitary_output.datDOCI.h5'
    d = print_unitary(unitname)
    #reader = PsiReader('output.dat', isbig = False, numorbs = None, read_ints = True)
    #reader = PsiReader('./results/h2o_testsunit/output_run/psi0.75sto-3g.out', isbig = False, numorbs = None, read_ints = True)
    #t2,t4 = reader.transform_integrals(d.get_unitary())
    #reader.matrix_to_list(t2 , t4)
    #reader.create_output(fname = 'newoutput.dat')

def print_unitary(fname):
    d = Unitary_Matrix(fname)
    d.print_unitary()
    return d

def benzene():
    reader = PsiReader('results/benzene_deformation/output_run_50/psi0_6-31g1.05.out', isbig = False, numorbs = None, read_ints = False)
    reader.set_active_space([1,2] , [] , [0] , 'detfile.dat')

class HDF5Reader(object):
    def __init__(self,fname):
        self.f = h5py.File(fname,"r")
        self.print_structure()

        datanames = {'nucrep' : 'Econst','nalpha' :'nup', 'nbeta': 'ndown', 'norb': 'L' , 'sym' :'nGroup' , 'irreplist' : 'orb2irrep' }
        self.extract_values(datanames)
        self.extract_ints()
        print self.values
        print self.ints
            
    def print_structure(self):
        def p(x):
            print x
        self.f.visit(p)
        #for name in self.f['/Data']:
            #print name
            #print self.f['/Data/'+name][:]

    def extract_values(self,datanames):
        self.values = {}
        for key, value in datanames.iteritems():
            try:
                self.values[key] = self.f['Data/'+value][:]
            except:
                print 'We couldnt extract : ' , key
                print 'we continue'

        self.set_irrep_info()

    def extract_ints(self):
        self.ints = { } 
        self.ints['mo2index'] = []
        for it , numorb in enumerate(self.occurences):
            if numorb > 0:
                self.ints['mo2index'].append(self.f['/TwoIndex/'+'TwoIndex'+str(it)+'/Matrix elements/'][:])
        self.ints['mo4index'] = self.f['/FourIndex/FourIndexObject/Matrix elements/'][:]    

    def set_irrep_info(self):
        self.occurences = np.bincount(np.array(self.values['irreplist']))
        self.cum = [0]
        for i in range(1,len(self.occurences)+1):
            self.cum.append(self.cum[i-1] + self.occurences[i-1])
        #print self.occurences
        #print self.cum
         

def hdf5_ham():
    fname = 'results/beh2_sto_3g_symcomp/hamiltonians/hampsi0_sto-3g0.86DOCIsim.h5'
    ham = HDF5Reader(fname)

def list_test():
    valdict = {'nalpha': 7 , 'nbeta' : 7 , 'norb' : 10 , 'nirrep' : 1 , 'sym' : 'c1' , 'nucrep' :0.560000000000000053290705182007513940334320068359 , 'irreplist': [0]*10, 'hfenergy' : -75.4675801847 , 'DOCC' : [7], 'SOCC' : [0]}
    dir = 'results/nopluspatrickham'
    reader = PsiReader( "" , isbig = False, numorbs = None, read_ints = True, valdict = valdict)
    reader.add_patrick_mo(one = os.path.join(dir, 'no_100.0_PB.one') , two =  os.path.join(dir , 'no_100.0_PB.two') )
    reader.create_output(fname = os.path.join(dir, 'hamnoplussto-3gpatrick100.0.out') )


def main():
    filename = 'output.dat'
    reader = PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
    #print reader.integral_matrix()
    #the unitary underneath is equivalent with swapping orbital with index 0 and 1.
    unitary = np.diag([1]*5)
    unitary[0,1] = 1
    unitary[1,0] = 1
    unitary[0,0] = 0
    unitary[1,1] = 0
    t2 , t4 =  reader.transform_integrals(unitary)
    reader.matrix_to_list(t2 , t4)
    #reader.swap_orbitals(3,0)
    #print reader.extract_energies()
    #reader.create_ni_system()
    #reader.change_repulsion()
    #reader.add_index(-1)
    #reader.keep_orbs(70)
    reader.create_output(fname = 'newoutput.dat')
    #reader.add_diego_mo()
    #reader.add_patrick_mo()

if __name__ == "__main__":
    #benzene()
    test_unitary()
    #print_unitary('./data/unitary-mo.h5')
    #hdf5_ham()
    #list_test()
    #main()
