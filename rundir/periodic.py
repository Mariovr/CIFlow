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
'''Periodic table of elements

   This module contains an object ``periodic`` that can be used as a Pythonic
   periodic table. It can be used as follows::

       >>> import periodic
       >>> periodic['si'].number
       14
       >>> periodic['He'].number
       2
       >>> periodic['h'].symbol
       'H'
       >>> periodic[3].symbol
       'Li'
       >>> periodic['5'].symbol
       'B'
'''


from units import angstrom
import os


__all__ = ['periodic', 'Element', 'Periodic']



class Element(object):
    '''Represents an element from the periodic table.

       The following attributes are supported:

       number
            The atomic number

       symbol
            A string with the symbol of the element.

       cov_radius
            The covalent radius. B. Cordero, V. Gomez, A. E. Platero-Prats, M.
            Reves, J. Echeverria, E. Cremades, F. Barragan, and S. Alvarez,
            Dalton Trans. pp. 2832--2838 (2008), URL
            http://dx.doi.org/10.1039/b801115j

       bs_radius
            The Bragg-Slater radius. J. C. Slater, J. Chem. Phys. 41, 3199
            (1964), URL http://dx.doi.org/10.1063/1.1725697

       vdw_radius
            van der Waals radius. R. S. Rowland and R. Taylor, J. Phys. Chem.
            100, 7384 (1996), URL http://dx.doi.org/10.1021/jp953141+

       wc_radius
            Waber-Cromer radius of the outermost orbital maximum. J. T. Waber
            and D. T. Cromer, J. Chem. Phys. 42, 4116 (1965), URL
            http://dx.doi.org/10.1063/1.1695904
    '''
    def __init__(self, number=None, symbol=None, cov_radius=None, bs_radius=None, vdw_radius=None, wc_radius=None):
        self.number = number
        self.symbol = symbol
        self.cov_radius = cov_radius
        self.bs_radius = bs_radius
        self.vdw_radius = vdw_radius
        self.wc_radius = wc_radius


class Periodic(object):
    '''A periodic table data structure.'''
    def __init__(self, elements):
        '''**Arguments:**

           elements
                A list of :class:`Element` instances.
        '''
        self.elements = elements
        self._lookup = {}
        for element in elements:
            self._lookup[element.number] = element
            self._lookup[element.symbol.lower()] = element

    def __getitem__(self, index):
        '''Get an element from the table based on a flexible index.

           **Argument:**

           index
                This can be either an integer atomic number, a string with the
                elemental symbol (any case), or a string with the atomic number.

           **Returns:** the corresponding :class:`Element` instance
        '''
        result = self._lookup.get(index)
        if result is None and isinstance(index, basestring):
            index = index.strip()
            result = self._lookup.get(index.lower())
            if result is None and index.isdigit():
                result = self._lookup.get(int(index))
                if result is None:
                    raise KeyError('Could not find element %s.' % index)
        return result

def load_periodic():

    convertors = {
        'int': (lambda s: int(s)),
        'float': (lambda s : float(s)),
        'str': (lambda s: s.strip()),
        'angstrom': (lambda s: float(s)*angstrom),
    }

    fn = os.path.join('data' , 'elements.txt')
    nelement = 0
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        with open(fn,'r') as infile:
            rows = {}
            step = 1
            for line in infile:
                line = line[:line.find('#')].strip()
                if len(line) > 0:
                    if step == 1:
                        name, convert = line.split()
                        step = 2
                    elif step == 2:
                        row = []
                        for word in line.split():
                            if word == 'None':
                                row.append(None)
                            else:
                                row.append(convertors[convert](word))
                        nelement = max(nelement, len(row))
                        rows[name] = row
                        step = 1
        elements=[]
        args=[]
        for i in xrange(nelement):
            for name, values in rows.iteritems():
                if len(values) < nelement :
                    for j in range(nelement-len(values)):
                        values.append(None)
                args.append((name,values[i]))
            kwargs = dict(args)
            elements.append(Element(**kwargs))
    else:
        elements = []
        periodic2 = ['' , 'H' , 'He' , 'Li' , 'Be',              'B',              'C',              'N',              'O',              'F',             'Ne',             'Na',             'Mg',             'Al',             'Si',              'P',              'S',             'Cl',             'Ar',              'K',             'Ca',             'Sc',             'Ti',              'V',             'Cr',             'Mn',             'Fe',             'Co',             'Ni',             'Cu',             'Zn',             'Ga',             'Ge',             'As',             'Se',             'Br' ] 

        for i in xrange(len(periodic2)):
            elements.append(Element(number = i , symbol = periodic2[i] ))
            #print 'periodic2[i]' , periodic2[i]


    return Periodic(elements)

periodic = load_periodic()

def main():
    for element in periodic.elements:
        print element.number , ' = ' ,  element.symbol, ' cov rad =  ' , element.cov_radius , ' bs rad =  ' , element.bs_radius, ' vanderwaals rad =  ' , element.vdw_radius, ' wc rad =  ' , element.wc_radius



if __name__ == "__main__":
    main()

