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
import random
import copy
import filecollector as fc

#sys.path.append('/home/mario/ownCloud/Doctoraat/CIFlowImproved/rundir')
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
            if(fname[-2:] == "h5"):
                self.f = h5py.File(fname,"r")
                self.unitary = [None]*8 #biggest point group in psi4 is D2h, which has 8 irreps.
                self.numorb = 0
                for name in self.f['/Data']:
                    irrep_num= int(re.search(r'(\d+)',name).group())
                    linsize = sqrt(len(self.f['/Data' + name][:]))
                    self.numorb += linsize
                    self.unitary[irrep_num] = np.array(self.f['/Data' + name][:] ) #transforms implicitly because unitary matrix saves first columns (fast changing index), but when we use numpy.shapeit puts data in rows and then cuts. This is oke, because unitary matrix saves the eigenbasis where we need to transform to in the rows, and this python class will keep them than in the columns. This is consistent with new eigenbasis obtained by diagonalizing a 2 index matrix.
                    self.unitary[irrep_num].shape = (linsize , linsize)
            else:
                self.unitary = [None]*8 #biggest point group in psi4 is D2h, which has 8 irreps.
                with open(fname , 'r') as file:
                    for line in file:
                        if "CIFlowTransformation" in line:
                            break
                    self.numorb = 0
                    for line in file:
                        irrep_num= int(re.search(r'(\d+)',line).group())
                        irreplist = []
                        for line in file:
                            try:
                               irreplist.append(np.array(map(float,line.split() ) ) )
                               self.numorb += 1
                            except ValueError:
                                print 'finished executing irrep: ' , irrep_num
                                break
                    self.unitary[irrep_num] = np.transpose(np.array(irreplist)) 

                
        def get_unitary(self):
            #At this moment point group symmetry is not yet supported
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
    def __init__(self, nalpha , nbeta , nucnuc , norb , sym , nirrep , irreplist , oeilist, teilist , hfenergy, docclist , socclist, overlap = None , unit = None):
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
        self.overlap = overlap
        self.unit = unit
        self.create_strings()

    def create_strings(self):
        self.irrepstring = ' '.join(map(str,self.irreplist))
        self.doccstring = ' '.join(map(str,self.docclist))
        self.soccstring = ' '.join(map(str,self.socclist))
        self.oeistring = '\n'.join(['%d %d %.25f '%(i,j,val)      for i, j ,val in self.oeilist])
        self.teistring = '\n'.join(['%d %d %d %d %.25f '%(i,j,k,l,val) for i, j ,k,l,val in self.teilist])
        if self.overlap != None:
            self.overlapstring = ""
            for irrep, data in enumerate(self.overlap):
                if data != None:
                    self.overlapstring += "irrep_"+ str(irrep) + "\n"
                    self.overlapstring+= '\n'.join([ '    '.join(map( lambda x: "%.12f" % x , dat) ) for dat in data  ]  )

            self.transformationstring = ""
            for irrep, data in enumerate(self.unit):
                if data != None:
                    self.transformationstring+= "irrep_"+ str(irrep) + "\n"
                    self.transformationstring+= '\n'.join([ '    '.join(map( lambda x: "%.12f" % x , dat) ) for dat in data  ]  )

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
"""%(self.__dict__)
        if self.overlap != None:
            text += """CIFlowOverlap: 
%(overlapstring)s
CIFlowTransformation: 
%(transformationstring)s
"""%(self.__dict__)
        text+= """****  MO OEI 
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


class ModHam(CIFlow_Input):
    def __init__(self  , nalpha , nbeta , norb , modtype, options , params, matrixelements = None):
        if matrixelements is None :
            super(ModHam, self).__init__(nalpha , nbeta , 0., norb , 'c1' , 1 , [0]*norb, [] , []  , 0., [nalpha], [0] )
            self.options = options
            self.params = params
            self.modtype = modtype
            self.make_mod_strings()
        else:
            self.psireader = PsiReader(matrixelements, isbig = False , numorbs = -1 , read_ints = True, valdict = {}, intdict = {})
            super(ModHam , self).__init__(self.psireader.values['nalpha'],self.psireader.values['nbeta'],self.psireader.values['nucrep'],self.psireader.values['norb'] ,self.psireader.values['sym'] ,self.psireader.values['nirrep'] ,self.psireader.values['irreplist'],self.psireader.ints['mo2index'],self.psireader.ints['mo4index'],self.psireader.values['hfenergy'],self.psireader.values['DOCC'] ,self.psireader.values['SOCC'],self.psireader.overlap ,self.psireader.unit)
            self.options = options
            self.params = params
            self.modtype = modtype
            self.make_mod_strings()

    def make_mod_strings(self):
        self.optionsstring = '\n'.join(self.options)
        self.paramsstring = '\n'.join(map(str , self.params ) )

    def transform_integrals(self, unitary):
        return self.psireader.transform_integrals(unitary)

    def matrix_to_list(self , t2, t4):
        self.psireader.matrix_to_list(t2 , t4)

    def __str__(self):
        text = super(ModHam , self).__str__()
        text += """###Model Hamiltonian Type: %(modtype)s
###Model Hamiltonian options:
%(optionsstring)s
####Params:
%(paramsstring)s
""" %(self.__dict__ )
        return text

class Reader(object):
    """
    class to read in data files 
    """
    def __init__(self, filename ):
        self.filename = filename
        self.read()

    def read(self):
        size = os.stat(self.filename).st_size
        self.f = open(self.filename)
        self.data = mmap.mmap(self.f.fileno(), size, access=mmap.ACCESS_READ)
        self.f.close()

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
            self.regexps = {'nalpha': (r'Nalpha\s*=\s*(\d+)', int) , 'nbeta' : (r'Nbeta\s*=\s*(\d+)', int), 'norb' : (r'Number Of \w+ Orbitals\s*=\s*(\d+)', int) , 'nirrep' : (r'Nirreps\s*=\s*(\d+)' ,int ) , 'sym' : (r'(?m)^Symmetry\s*Label\s*=\s*([\d\w]*)' , lambda x: x), 'nucrep' : (r'(?m)^Nuclear Repulsion Energy =\s*([\-+]?\d+\.?\d*[eEdD]?[\-+]?\d*)', float), 'irreplist' : (r'Irreps Of \w+ Orbitals\s*=\s*\n([\d\s]+)\n', lambda x : map(int,x.split())), 'hfenergy' : ( '(?m)^\*+\s*HF\s*Energy\s*=\s*([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)', float), 'DOCC' : ( r'(?m)^DOCC\s+=\s+([\d\s]+)', lambda x: map(int , x.split() ) ) , 'SOCC' : ( r'(?m)^SOCC\s+=\s+([\d\s]+)', lambda x : map(int, x.split() ) ) }

            self.regint = { 'mo4index' : (r'(?im)^(\d+)\s(\d+)\s(\d+)\s(\d+)\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', self.tei_func) , 'mo2index' : (r'(?im)^(\d+)\s*(\d+)\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)', self.oei_func )}

            #self.create_star_data()
            self.extract_values()
            self.values = copy.deepcopy(self.values)
            self.overlap = self.read_ao_info(filename, "CIFlowOverlap")
            self.unit = self.read_ao_info(filename, "CIFlowTransformation")
            for i in range(len(self.unit) ):
                if self.unit[i] != None:
                    self.unit[i] = np.array(self.unit[i]).T #save in columns the new basis

            if isbig:
                self.delete_lines(numorbs) #numorbs is the number of orbitals you wanna keep
                self.values['irreplist'] = self.values['irreplist'][:numorbs]
                self.values['norb'] = numorbs


            if read_ints:
                self.extract_ints()
                self.ints = copy.deepcopy(self.ints)
            self.data.close()
            
    def read_ao_info(self, fname , startstring):
        overlap = [None] * 8 #8 is max num of irreps in abelian point groups we consider
        with open(fname , 'r') as file:
            for line in file:
                if startstring in line:
                    break
            irrep_num = -1
            irreplist = []
            for line in file:
                if re.search(r'irrep', line ):
                    if irrep_num != -1:
                        overlap[irrep_num] = irreplist
                    irrep_num = int(re.search(r'(\d+)',line).group())
                    irreplist = []
                elif re.search(r'([\-+]?\d+\.*\d*[eEdD]?[\-+]?\d*)', line):
                    irreplist.append(np.array(map(float,line.split() ) ) )
                else:
                    overlap[irrep_num] = np.array(irreplist )
                    #print "finished reading in ao information this is the result:"
                    #print overlap
                    break 

        return overlap

    def get_unitary(self):
        #At this moment point group symmetry is not yet supported
        unit = np.zeros((self.values['norb'], self.values['norb']))
        pos = 0
        for irrep in range(len(self.unit)):
            if self.unit[irrep] != None:
                dim = self.unit[irrep].shape[0] 
                unit[pos:pos+dim , pos:pos+dim] = self.unit[irrep][:,:]
                pos += dim
        return unit        

    def extract_values(self):
        for key, value in self.regexps.iteritems():
            try:
                self.values[key] = self.process_element(value[0] , self.data , value[1])
            except:
                self.values[key] = [None]
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

    def add_patrick_mo(self,num = -1, one = 'OEI' , two = 'TEI', overlapfile = None, unitfile = None):
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

        self.overlap = [None] * 8 #8 is max num of irreps in abelian point groups we consider
        self.unit = [None] * 8

        if overlapfile != None:
            self.overlap[0] = np.loadtxt(overlapfile  , comments = '#')
            self.unit[0] = np.loadtxt(unitfile , comments = '#')

    def add_diego_mo(self , fname = 'h2o-psi-cisd-1_8-mario.out'): #num adds a number to the integral indices but is not necessary because diego uses psi3 to generate his integrals and it uses the same numbering as psi4.
        with open(fname, 'r') as tofile:
            oeilines = tofile.read()

        oei = r'h\((\d+)\)\((\d+)\)\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'
        tei = r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*\([\s\d]+\)\s*=\s*([\-+]?\d+\.\d+[eEdD]?[\-+]?\d+)'

        oeilist = re.findall(oei , oeilines)
        teilist = re.findall(tei , oeilines)
        self.ints['mo2index'] = self.oei_func(oeilist)
        self.ints['mo4index'] = self.tei_func(teilist)

    def add_pierre_mo(self, fname = 'SphInts_2'):
        with open(fname, 'r') as tofile:
            intlines = tofile.read()

        oei = r'([\-+]?\d+\.\d*[eEdD]?[\-+]?\d*)\s*(\d+)\s+(\d+)\s+0\s+0'
        tei = r'([\-+]?\d+\.\d*[eEdD]?[\-+]?\d*)\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)'

        oeilist = [ (a,b,val)    for val , a,b  in re.findall(oei , intlines) ]
        teilist = [(a,b,c,d , val ) for val , a,b,c,d in re.findall(tei ,  intlines) ]
        
        self.ints['mo2index'] = self.oei_func(oeilist)
        self.ints['mo4index'] = self.tei_func(teilist)
        self.add_index(-1)

        self.overlap = [None] * 8 #8 is max num of irreps in abelian point groups we consider
        self.unit = [None] * 8

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


    def create_constrained_vals(self, lagmul , elpop, orblist):
        self.values['nucrep'] -=  lagmul * elpop 
        partdiag = np.zeros((self.values['norb'], self.values['norb']), float )
        for i in orblist:
            partdiag[i,i] =1.
        unit = self.get_unitary().T
        overlap = self.overlap[0]
        part = np.dot(unit , overlap )
        part2 = np.dot(part , partdiag)
        part = np.dot(part2 , unit.T)
        total = (part + part.T)/2. * lagmul #for hermicity and mult with lagmul 
        num = 0
        for  i , j , integral in self.ints['mo2index']:
            self.ints['mo2index'][num][2] += total[i,j]
            num +=1

    def calc_energy(self , rdm1 , rdm2, orblist = [], nucrepbool = True, startrdm = 0, mul_check = True):
        #remark integrals are written in chemical notation, and rdms are in physical notation.
        if orblist == []:
            orblist = range(0, self.values['norb'])

        if mul_check:
            overlap = self.overlap[0]
            print overlap , rdm1[0]
            part = np.dot(overlap, rdm1[0] )
            som = 0
            for orb in orblist:
                som += part[orb,orb]
            print 'Mulliken charge or orblist: ' , orblist , ' is equal to : ' , 2*som

        nucrep = 0 
        if nucrepbool == True:
            nucrep += self.values['nucrep']

        energy1 = 0
        for i , j , integral in self.ints['mo2index']:
            #print i , j ,integral
            if i in orblist and j in orblist: # and j in orblist:
                energy1 += (integral ) * rdm1[0][i,j] 
                energy1 += (integral )  * rdm1[1][i,j]
        #for i , j , integral in self.ints['mo2index']:
        #    #print i , j ,integral
        #    if j in orblist: # and j in orblist:
        #        energy1 += (integral ) * rdm1[0][i,j] 
        #        energy1 += (integral )  * rdm1[1][i,j]
        #energy1 /=2.
        print 'Read psi 1 electron energy : ' , energy1

        energy2 = 0
        for i , j , k , l, integral in self.ints['mo4index']:
            if i in orblist and k in orblist:# and l in orblist and j in orblist :
                energy2 += integral * rdm2[0][i , k ,j  ,l ]
                energy2 += 2*integral * rdm2[1][i , k , j ,l ]
                energy2 += integral * rdm2[2][i , k ,j  ,l ]

        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if i in orblist and l in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if k in orblist and j in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if l in orblist and j in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #energy2 /= 4.

        print 'Read psi 2 electron energy : ' , energy2/2.
        print 'total energy = : ' , energy1 + 1/2.* energy2 + nucrep
        return energy1 + 1/2.*  energy2 + nucrep

    def calc_energy_sym(self , rdm1 , rdm2, orblist = [], nucrepbool = True, startrdm = 0, mul_check = True):
        #remark integrals are written in chemical notation, and rdms are in physical notation.
        if orblist == []:
            orblist = range(0, self.values['norb'])

        if mul_check:
            overlap = self.overlap[0]
            print overlap , rdm1[0]
            part = np.dot(overlap, rdm1[0] )
            som = 0
            for orb in orblist:
                som += part[orb,orb]
            print 'Mulliken charge or orblist: ' , orblist , ' is equal to : ' , 2*som

        nucrep = 0 
        if nucrepbool == True:
            nucrep += self.values['nucrep']

        energy1 = 0
        for i , j , integral in self.ints['mo2index']:
            #print i , j ,integral
            if i in orblist: # and j in orblist:
                energy1 += (integral ) * rdm1[0][i,j] 
                energy1 += (integral )  * rdm1[1][i,j]
        for i , j , integral in self.ints['mo2index']:
            #print i , j ,integral
            if j in orblist and i not in orblist: # and j in orblist:
                energy1 +=  (integral ) * rdm1[0][i,j] 
                energy1 +=  (integral )  * rdm1[1][i,j]
        print 'Read psi 1 electron energy : ' , energy1

        energy2 = 0
        for i , j , k , l, integral in self.ints['mo4index']:
            if i in orblist and k in orblist:# and l in orblist and j in orblist :
                energy2 += integral * rdm2[0][i , k ,j  ,l ]
                energy2 += 2*integral * rdm2[1][i , k , j ,l ]
                energy2 += integral * rdm2[2][i , k ,j  ,l ]

        for i , j , k , l, integral in self.ints['mo4index']:
            if i in orblist and l in orblist and k not in orblist and j not in orblist :
                energy2 += integral * rdm2[0][i , k ,j  ,l ]
                energy2 += 2*integral * rdm2[1][i , k , j ,l ]
                energy2 += integral * rdm2[2][i , k ,j  ,l ]
        for i , j , k , l, integral in self.ints['mo4index']:
            if k in orblist and j in orblist:# and i not in orblist and l not in orblist :
                energy2 += integral * rdm2[0][i , k ,j  ,l ]
                energy2 += 2*integral * rdm2[1][i , k , j ,l ]
                energy2 += integral * rdm2[2][i , k ,j  ,l ]
        for i , j , k , l, integral in self.ints['mo4index']:
            if l in orblist and j in orblist:# and i not in orblist and k not in orblist :
                energy2 += integral * rdm2[0][i , k ,j  ,l ]
                energy2 += 2*integral * rdm2[1][i , k , j ,l ]
                energy2 += integral * rdm2[2][i , k ,j  ,l ]

        print 'Read psi 2 electron energy : ' , energy2/2.
        print 'total energy = : ' , energy1 + 1/2.* energy2 + nucrep
        return energy1 + 1/2.*  energy2 + nucrep

    def calc_atom_e(self, rdm1 , hartreecor , cum, orblist ):
        #remark integrals are written in chemical notation, and rdms are in physical notation.
        overlap = self.overlap[0]
        print overlap , rdm1[0]
        part = np.dot(overlap, rdm1[0] )
        som = 0
        for orb in orblist:
            som += part[orb,orb]
        print 'Mulliken charge or orblist: ' , orblist , ' is equal to : ' , 2*som

        energy1 = 0
        for i , j , integral in self.ints['mo2index']:
            #print i , j ,integral
            if i in orblist and j in orblist: # and j in orblist:
                energy1 += (integral ) * rdm1[0][i,j] 
                energy1 += (integral )  * rdm1[1][i,j]
        print 'Read psi 1 electron energy : ' , energy1

        energy2 = 0
        for i , j , k , l, integral in self.ints['mo4index']:
            if i in orblist and j in orblist:# and l in orblist and j in orblist :
                energy2 += integral * cum[0][i , k ,j  ,l ]
                energy2 += 2*integral * cum[1][i , k , j ,l ]
                energy2 += integral * cum[2][i , k ,j  ,l ]

        for i , j , k , l, integral in self.ints['mo4index']:
            if i in orblist and k in orblist:# and l in orblist and j in orblist :
                energy2 += integral * hartreecor[0][i , k ,j  ,l ]
                energy2 += 2*integral * hartreecor[1][i , k , j ,l ]
                energy2 += integral * hartreecor[2][i , k ,j  ,l ]

        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if i in orblist and l in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if k in orblist and j in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #for i , j , k , l, integral in self.ints['mo4index']:
        #    if l in orblist and j in orblist:# and l in orblist and j in orblist :
        #        energy2 += integral * rdm2[0][i , k ,j  ,l ]
        #        energy2 += 2*integral * rdm2[1][i , k , j ,l ]
        #        energy2 += integral * rdm2[2][i , k ,j  ,l ]
        #energy2 /= 4.

        print 'Read psi 2 electron energy : ' , energy2/2.
        #print 'total energy = : ' , energy1 + 1/2.* energy2
        #return energy1 + 1/2.*  energy2
        print 'total energy = : ' , energy1 +  energy2
        return energy1 +   energy2

    def create_ham_matrix(self):
        ham1 = np.zeros((self.values['norb'], self.values['norb']), float )
        ham2 = np.zeros((self.values['norb'], self.values['norb'], self.values['norb'], self.values['norb']), float )
        for i , j , integral in self.ints['mo2index']:
            ham1[i,j] = integral
            ham1[j,i] = integral

        for i , j , k , l, integral in self.ints['mo4index']:
            ham2[i,j, k ,l] = integral
            ham2[j ,i,k ,l ] = integral
            ham2[i ,j, l,k ] = integral
            ham2[k ,l, i,j ] = integral
            ham2[l ,k, i,j ] = integral
            ham2[k ,l, j,i ] = integral
            ham2[j,i, l ,k] = integral
            ham2[l,k, j ,i] = integral
        return ham1 , ham2

    def calc_energy_sparse(self , rdm1 , rdm2, orblist = [], nucrepbool = True):
        #remark integrals are written in chemical notation, and rdms are in physical notation.
        nucrep = 0 
        if nucrepbool == True:
            nucrep += self.values['nucrep']

        if orblist == []:
            orblist = range(0, self.values['norb'])

        ham1 , ham2 = self.create_ham_matrix()

        energy1 = 0
        for i in range(0, self.values['norb'] ):
            for j in range(0, self.values['norb'] ):
                if i in orblist and j in orblist:
                    energy1 += ham1[i,j] * rdm1[0][i,j] 
                    energy1 += ham1[i,j]  * rdm1[1][i,j]
        print 'Read psi 1 electron energy : ' , energy1

        energy2 = 0
        for i in range(0, self.values['norb'] ):
            for j in range(0, self.values['norb'] ):
                for k in range(0, self.values['norb'] ):
                    for l in range(0, self.values['norb'] ):
                        if i in orblist and j in orblist and k in orblist and l in orblist :
                            energy2 += ham2[i,j,k,l]* rdm2[0][i , k ,j  ,l ]
                            energy2 += 2*ham2[i,j,k,l]* rdm2[1][i , k , j ,l ]
                            energy2 += ham2[i,j,k,l]* rdm2[2][i , k ,j  ,l ]

        print 'Read psi 2 electron energy : ' , energy2/2.

        print 'total energy = : ' , energy1 + 1/2.* energy2 + nucrep
        return energy1 + 1/2.*  energy2 + nucrep

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

    def add_ni_system(self , psireader, outname = 'ni_system.dat'):
        psireader.add_index(self.values['norb'])
        norb = self.values['norb']
        oeint =list(psireader.ints['mo2index'])
        teint =list(psireader.ints['mo4index'])
        self.insert(oeint, teint)
        self.values['norb'] += psireader.values['norb']
        print self.values['nalpha']
        self.values['nalpha'] += psireader.values['nalpha']
        print self.values['nalpha'] 
        self.values['nbeta'] += psireader.values['nbeta']
        self.values['nucrep'] += psireader.values['nucrep'] 
        self.values['sym'] = 'c1'
        self.values['irreplist'] = [0] *self.values['norb']
        overlap = [None] * 8 #8 is max num of irreps in abelian point groups we consider
        overlap[0] = np.zeros( (self.values['norb'], self.values['norb']) )
        overlap[0][:norb , :norb] = self.overlap[0]
        overlap[0][norb: , norb : ] = psireader.overlap[0]
        unit = [None] * 8
        unit[0] = np.zeros((self.values['norb'], self.values['norb']))
        unit[0][:norb , :norb] = self.unit[0]
        unit[0][norb: , norb:] = psireader.unit[0]
        self.overlap = overlap
        self.unit = unit
        self.create_output(outname)

    def change_repulsion(self, factor = -1.):
        self.ints['mo4index'] = [[integrals[0] ,integrals[1] ,integrals[2] , integrals[3], factor*integrals[4] ] for integrals in self.ints['mo4index']]

    def create_output(self , fname = 'output.dat'):
        ciflow_input = CIFlow_Input(self.values['nalpha'], self.values['nbeta'], self.values['nucrep'], self.values['norb'] , self.values['sym'] , self.values['nirrep'] , self.values['irreplist'], self.ints['mo2index'], self.ints['mo4index'], self.values['hfenergy'], self.values['DOCC'] , self.values['SOCC'], self.overlap , self.unit)
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
        return (self.occurences, self.cum)

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

        dw.cimain(nup,ndown ,norb , [range(1,nup+ndown+1), [] ] , senlist ,fname = detfile ,ref = [lambda x , y, z : activehf ], add_frozen = 0, frozenstring = frozen , virtualstring = virtual) 

class HDF5Reader(object):
    def __init__(self,fname):
        self.f = h5py.File(fname,"r")
        self.print_structure()

        datanames = {'nucrep' : 'Econst','nalpha' :'nup', 'nbeta': 'ndown', 'norb': 'L' , 'sym' :'nGroup' , 'irreplist' : 'orb2irrep' , 'overlap' : 'overlap' }
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
    #fname = 'results/beh2_sto_3g_symcomp/hamiltonians/hampsi0_sto-3g0.86DOCIsim.h5'
    fname = 'hamatomicintegralsorthon.h5'
    ham = HDF5Reader(fname)

def list_test():
    dir = 'Mario_NO'

    fileinfo = lambda x: float(re.search(r'no_([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)_PB\.one' , x).group(1))
    search = r'no.*\.one' 
    hamfiles = fc.File_Collector(dir , search = search ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000000. )

    distancelist = [fileinfo(i) for i in hamfiles.plotfiles]
    print distancelist

    for num , dist in enumerate(distancelist):
        nucrep = 56./dist
        valdict = {'nalpha': 7 , 'nbeta' : 7 , 'norb' : 10 , 'nirrep' : 1 , 'sym' : 'c1' , 'nucrep' : nucrep , 'irreplist': [0]*10, 'hfenergy' : -75.4675801847 , 'DOCC' : [7], 'SOCC' : [0]}
        reader = PsiReader( "" , isbig = False, numorbs = None, read_ints = True, valdict = valdict)
        reader.add_patrick_mo(one = os.path.join(dir, 'no_%.1f_PB.one' % dist) , two =  os.path.join(dir , 'no_%.1f_PB.two' % dist) , overlapfile = os.path.join(dir, 'no_%.1f_PB.ove' % dist), unitfile =os.path.join(dir , 'no_%.1f_PB.mo' % dist)  )
        reader.create_output(fname = 'hamnoplussto-3gpatrick%.1fnew.out' % dist )

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

def reverse_repulsion():
    filename = 'results/hechangerep/matrixelements/psi0_basisset2.00.mout'
    reader = PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
    reader.change_repulsion()
    reader.create_output(fname = 'results/hechangerep/matrixelements/changedrep.dat')

def test_new_format():
    filename = "psioutput.dat"
    filename = "5natomhamao.dat"
    reader = PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
    reader.create_output(fname = 'newpsioutput.dat')


def test_mod_ham():
    nalpha =4 ; nbeta = 4 ; nucnuc = 0 ;  norb  = 8 ;  sym = 'c1' ; nirrep = 1 ; irreplist  = [0]* norb;  hfenergy = 0. ; docclist = [4] ; socclist = [0] ; overlap = None ; unit = None
    modtype = 'Hub1d' ; options = ['notzero' , 'pos', 'none']  ; params = [1. , 4.] 
    hub1d = ModHam(nalpha,nbeta,norb,modtype , options , params)
    hub1d.write_file(fname = 'hub1d.mod')

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
    if fname[-2:] == "h5":
        d.print_structure()
    return d

def benzene():
    reader = PsiReader('results/benzene_deformation/output_run_50/psi0_6-31g1.05.out', isbig = False, numorbs = None, read_ints = False)
    reader.set_active_space([1,2] , [] , [0] , 'detfile.dat')

def generate_random_unitary(norb):
    #Create random unitary matrix based on the QR decomposition (this is valid for every matrix)
    # A = QR with Q unitary and R invertible -> Q = AR^-1
    unit = np.random.uniform(0,1,norb * norb) 
    unit = np.reshape(unit , (norb,norb ))
    q, r = np.linalg.qr(unit)
    print q
    print np.dot(q , np.linalg.inv(q) ) 
    return q

def generate_random_hamiltonian(reader):
    unit = generate_random_unitary(reader.values['norb'])
    t2,t4 = reader.transform_integrals(unit)
    reader.matrix_to_list(t2,t4)
    return reader

def generate_random_hams():
    rootdir = './results/shannon_entropy/beh2'
    filename = './results/shannon_entropy/beh2/psioutput.dat'
    reader = PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
    os.chdir(os.path.join(rootdir , 'random_hamiltonians') )
    for i in range(1000):
        generate_random_hamiltonian(reader)
        name = 'randomhamiltonian' + str(i) +'.dat'
        reader.create_output(fname = name)

def add_ni_system():
    filename1 = "psioutputn+.dat"; filename2 = "psioutputo.dat"
    reader1 = PsiReader(filename1 , isbig = False, numorbs = None, read_ints = True)
    print reader1.values
    reader2 = PsiReader(filename2 , isbig = False, numorbs = None, read_ints = True)
    print reader1.values
    reader1.add_ni_system(reader2)
    print reader1.values

def test_pierre_format():
    filename = "psioutput.dat"
    reader = PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
    reader.add_pierre_mo( "SphInts_2" )
    reader.create_output(fname = 'pierreintegrals.dat')

if __name__ == "__main__":
    #reverse_repulsion()
    #benzene()
    #test_unitary()
    #print_unitary('./tests/data/unitary-moorthon.h5')
    #print_unitary('unitary.txt')
    #test_new_format()
    #test_mod_ham()
    #generate_random_hams()
    #hdf5_ham()
    #list_test()
    #main()
    test_new_format()
    #test_pierre_format()
    #add_ni_system()
