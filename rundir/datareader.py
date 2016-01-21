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
#! /bin/env python
from string import maketrans
import re
import numpy as np

def matchlist(line):
    match = re.search(r'\[(?P<list>[+\-\d\s.eE]*)\]' , line) 
    if match:
        return map(float,match.group('list').split())
    else:
        Exception("No match for a list was found in the following line: %s", line)
      
def matchdict(line):
    match = re.search(r'\{(?P<dict>[+\-\d\s.eE:,]*)\}' , line) 
    if match:
        return dict(map(lambda x: tuple(map(lambda x : int(float(x)),x.split(':'))) ,match.group('dict').split(','))) #little trick at second lambda function to make strings '3.0' transformable to integer 3
    else:
        Exception("No match for a list was found in the following line: %s", line)
  
class Reader(object):
    """
    Abstract data class that reads in data from files
    """
    def __init__(self, filename , comment = '#' ):
        self.filename = filename
        self.comment = comment
        self.kind = None
        self.depvar = {'xas' : 'R' , 'yas' : 'CI energies' }
        self.units = {'x' : '(a.u.)' , 'y' : '(a.u.)'}
        self.read()
  
    def read(self):
        self.read_header()
        self.read_data()

    def get_xlabel(self):    
        return self.depvar['xas'] + ' ' + self.units['x']

    def get_ylabel(self):
        return self.depvar['yas'] + ' '+self.units['y']

    def set_ylabel(self, yasstring):
        self.depvar['yas'] = yasstring

    def convert_au_to_ang(self , x):
        """
        See units.py for how to convert angstroms to atomic units
        """
        meter = 1.0/0.5291772083e-10
        angstrom = 1.0e-10*meter
        return x*angstrom

    def convert_ang_to_au(self, x):
        meter = 1.0/0.5291772083e-10
        angstrom = 1.0e-10*meter
        return x/angstrom
    
    def convert(self, col = 0 , to = "au"):
        if to == "au":
            self.data[:,col] = [self.convert_au_to_ang(x)  for x in self.data[:,col]]
        else:
            self.data[:,col] = [self.convert_ang_to_au(x)  for x in self.data[:,col]]

    def read_header(self):  
        """
        Reads in the header of the outputfiles.
        """
        with open(self.filename, 'r') as file:
            for line in file:
                if  line.startswith(self.comment):
                    self.columnheader = line.split()
                    self.columnheader[0] = self.columnheader[0].strip(self.comment)
                    self.depvar['xas'] = self.columnheader[0]
                    if 'nup' in line:
                        self.nup = int(re.search(r'nup[:=\s]*(\d+)',line).group(1))
                    elif 'ndown' in line:
                        self.ndown = int(re.search(r'ndown[:=\s]*(\d+)',line).group(1))
                    elif 'norbs' in line:
                        self.norbs = int(re.search(r'norbs[:=\s]*(\d+)',line).group(1))
                else:
                    break #header is finished usefull for long outputfiles

        print 'we have read the information from the header of ', self.filename
        print self.columnheader
    
    def read_data(self):
        with open(self.filename,'r') as ref:
            #For more info check: http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html#numpy.genfromtxt  
            self.data = np.genfromtxt(self.filename, dtype='float', comments= self.comment, delimiter=None, skiprows=0, skip_header=0, skip_footer=0, converters=None, missing='', missing_values= 0., filling_values=None, usecols=None, names=None, excludelist=None, deletechars=None, replace_space='_', autostrip=False, case_sensitive=True, defaultfmt='f%i', unpack=None, usemask=False, loose=True, invalid_raise=True)

    def calc_error(self , difcols, absolute = True):
        errorlist = []
        for tuple in difcols:
            if absolute:
                errorlist.append( abs(self.data[:,tuple[0]] - self.data[:, tuple[1]  ]) ) 
            else:
                errorlist.append( self.data[:,tuple[0]] - self.data[:, tuple[1]  ] ) 
        
        return errorlist

    def save_data(self, fname):
        np.savetxt(fname , self.data)
            
  
  

def main():
    fname =   "./h8output.dat"
    d = Reader(fname)
    print d.data
    d.convert(0)
    print d.data
    d.save_data("converted.dat")


if __name__ == "__main__":
    main()
