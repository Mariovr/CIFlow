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
import os
import shutil
import subprocess
import re 
import sys
import math

import detwrite as dw
from periodic import periodic
import read_psi as rp
import ciflowoutput as cp
import filecollector as fc

if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    import plotfunctions as pf

def generate_dir(name,infilename = None,args = None, prefix = ''):
    """ generates a seperate dir for all the information of the
    handled problem
    """  
    dirname = "%s" %name
    print ('#making directory for the results: %s' %dirname)
    dir = dirname #name of the directory
    #check if the directory already exists
    if os.path.isdir(dir): #OEPS the directory exists already
        no = 'l'
        while  no != 'y' and no != 'n':
            no = raw_input('remove dir %s: (y,n)' %name).lower()
        if no == 'y':
            shutil.rmtree(dir)
            os.mkdir(dir)
        else:
            print 'the directory you wanted to create does already exist so we copy the files there.'
        
    else:
        os.mkdir(dir)

    if isinstance(infilename,str): 
        shutil.copy(infilename,dir)
    if args is not None:
        for i in args:
            shutil.copy(i,dir)
    print dir        
    if prefix:
        path, file = os.path.split(prefix)
        x = [os.path.join(path , i)  for i in os.listdir(path) if re.search(file, i)]
        for a in x:
            shutil.copy(a,dir)
    os.chdir(dir)
    
def convert_au_to_ang(x):
    """
    See units.py for how to convert angstroms to atomic units
    This function will be removed soon.
    """
    return x*5.2917720859*10**(-11)/(10**(-10))

def input_psi(fname, basis ,name= ' ' , charge_mult = (0,1) , atomlist = [] , positions = [] , units = 'au' , ref = 'rhf' , userbasis = False , functional = None, path_to_plugin = '../mointegrals/mointegrals.so', DOCC = None, sym = None, energies = None, hdf5 = False, guess = 'sad', su = False, basispath = './'):
    '''
    Generate the input for the psi4 program.
    #set reference can be: rhf , rohf, uhf , rks , uks
    '''
    text = """#PYTHONPATH must include directory above plugin directory.
#Define either externally or here, then import plugin.

"""
    molecule = 'molecule %s {\n' %name
    ch_mulst = ' '.join(map(str, charge_mult)) + '\n'
    el_list = []
    for index , atom in enumerate(atomlist):
        el_list.append(periodic[atom].symbol + ' ' + ' '.join(map(str,positions[index])))
    el_list = '\n'.join(el_list)
  
    if sym != None:
        el_list += '\nsymmetry %s\n' %sym
    unitstring = '\nunits %s\n' %units
    endmolecule = '}\n'
    molecule += ch_mulst + el_list + unitstring + endmolecule
  
    if userbasis:
        #molecule += '\nbasis file ./%s%s.gbs\n' %(basispath, basis)
        pass
  
    integraltype = os.path.basename(path_to_plugin) 
    integralname = integraltype[:-3]
    fout = fname[:-4]+ '.mout'

    integraltext = """plugin_load("%(path_to_plugin)s")
set basis %(basis)s
set %(integralname)s print 1
set %(integralname)s filename %(fout)s
set %(integralname)s save_unitaries %(su)s
set reference %(ref)s
set guess %(guess)s
set {
#exact ERI
SCF_TYPE pk
DF_SCF_GUESS false
# strict convergence
E_CONVERGENCE 1e-9
D_CONVERGENCE 1e-9
ints_tolerance 1e-10
}

"""%vars()

    if hdf5:
        integraltext += 'set %s HDF5_FILENAME = %s\n' %(integralname, fname)

    if DOCC != None:
        integraltext += '\nset docc = %s\n' %DOCC
  
  
    if functional:
        integraltext += '\nset dft_functional %s\n' %functional
    
    if 'mo' in integralname:
        endtext = "scf() \n\n"
    else:
        endtext = "\n"
    endtext += """plugin("%s")\n"""%integraltype

    if energies !=None:
        for x in energies:
            endtext += "energy('%s')\n" %x

    psi_input = text + molecule + integraltext + endtext
    f = open(fname , 'w')
    f.write(psi_input)
    f.close()

def create_matrix_elements(elemdir, basissets, runlist, atoms, chmult = (0,1) , moltype = 'dimer', package = 'psi' , units = 'angstrom', path_mo = '../../../../mointegrals/mointegrals.so', DOCC = None, energies = None , sym = None, extrapar= 104.479848, only_input = False, hdf5 = False, guess = 'sad', ref = 'rhf', su = 'False', basispath = './', userbasis = False):
    olddir = os.getcwd()
    generate_dir(elemdir,None)
    filenum = 0 #To read in the matrixelement files in the correct order.
    for basis in basissets:
        for r in runlist:
            if moltype == 'linear':
                positions = [[0,0,0],[0,0,r],[0,0,-1*r]] #for molecules like beh2
            elif moltype == 'dimer':
                positions = [[0,0,0],[0,0,r]] #for molecules like n2
            elif moltype == 'c2v':    
                positions = [[''],[1, r],[1, r, 2 , extrapar]] #for molecules like h2o
            elif moltype == 'atom':
                positions = [[0,0,0]]
            elif moltype == 'benzene':
                positions = benzene(r, extrapar) #REMARK we run over the angle -> r in runlist, and keep the radius of the circle on which all the C atoms are placed constant to extrapar, (extrapar = 1.398 A is standard)
            elif moltype == 'h8':
                positions =  [[0,0,0],[0,0,r],[0,0,2*r],[0,0,3*r],[0,0,4*r],[0,0,5*r], [0,0,6*r] , [0 , 0 , 7*r] ]
            elif moltype == 'h6':
                positions =  [[0,0,0],[0,0,r],[0,0,2*r],[0,0,3*r],[0,0,4*r],[0,0,5*r]]  
            elif moltype == 'but':
                positions = [[0,0,0],[extrapar,0,0], [0,r,0] , [extrapar,r,0]]
            elif moltype == 'square': 
                positions = [[0,0,0 ] , [0,r,0],[0,-r,0],[r,0,0] ,[-r,0,0]]

            else:
                print 'Error: moltype %s is not know' %moltype

            if package == 'psi':    
                fname = "psi%d_" %filenum +basis + "%.2f" %(r) +  ".dat"
                input_psi(fname , basis ,name= 'mol' , charge_mult = chmult , atomlist = atoms , positions = positions , units = units , ref = ref, userbasis = userbasis, path_to_plugin = path_mo , DOCC = DOCC, sym = sym, energies = energies, hdf5= hdf5, guess = guess, su = su, basispath = basispath)
                if not only_input:
                    subprocess.call(["psi4", fname])

        filenum += 1        

    os.chdir(olddir)    

def create_ciflow_input_file(matrixelements , methods , fname = "flow.dat", prin = False):
    with open(fname, 'w') as file:
        if matrixelements == '.':
            matrixelements = matrixelements.lstrip('.')
        file.write(matrixelements+'\n')
        file.write(methods[0]+'\n')
        for cimethod in methods[1:]:
            if( cimethod == "doci" or cimethod == "fci" or cimethod == "file" or cimethod == "big"):
                file.write('endm\n')
            file.write(cimethod+'\n')
        file.write('endm\n')
        file.write('end\n')
        if prin:
            file.write('true\n')

def gkci():
    root = 'Be_exrapolation_gamma' ; fname = 'flow.dat'
    shutil.copy('./ciflow.x', root)#When the matrixelements are already present.
    os.chdir(root) #When the matrixelements are already present.
    dirs = os.listdir(os.getcwd())
    dirs = [dir for dir in dirs if 'tar' not in dir and '.x' not in dir and '.out' not in dir and '.dat' not in dir]
    dirs.sort()
    print dirs
    methods = [('rhf' , None) ]#, ('rks', 'wb97x')]
    for file in dirs:
        #files =  os.listdir(os.path.join(root,dir))  #+ os.listdir('He')
        #files.sort()
        match = re.search(r'(\w+)_([\w\-\d]+)_\d+.txt', file)
        basis = match.group(2)
        if 'V5' in basis:
            el = "Be"
            print el , basis
            name = '%s_%s_gkcioutput.dat' %(el,basis)
            with open(name, 'w') as tofile:
                for method in methods:
                    #input_psi(basis+'input.dat', basis ,name= ' ' , charge_mult = (0,1) , atomlist = [periodic[el].number] , positions = [(0,0,0)] , units = 'au' , ref = method[0] , sym = 'c1' , functional = method[1], path_to_plugin = '../../mointegrals/mointegrals.so')
                    #subprocess.call(["psi4", basis + "input.dat"])
                    tofile.write('#method: %s\tfunctional: %s\n' %(method[0] , str(method[1])))
                    #for file in files:
                    print file
                    tofile.flush()
                    numdet = int(re.search(r'\w+_[\w\-\d]+_(\d+).txt', file).group(1))
                    print 'dets: ' , numdet
                    #shutil.copyfile(os.path.join(os.path.join(root,dir),file) , 'determinants.dat')
                    try:
                        if numdet >= 60000:
                          create_ciflow_input_file(basis+"input.out" , ['big',file , 'none'], fname = fname)
                        else:
                          create_ciflow_input_file(basis+"input.out" , ['file',file , 'none'], fname = fname)
                        doci_en =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(fname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line)
                            doci_en += line
                        match_ci = re.search(r'CIenergy[\w\s\d._/-]+:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                        tofile.write('%s\t%s\n' %(file.strip('.txt') , match_ci.group(1)))
                    except subprocess.CalledProcessError as e:
                        print e.returncode
                        tofile.write('%s\t%s\n' %(file.strip('.txt') ,e.returncode ) )


def scale_basis(namein , nameout , scale, atomold , atomnew):
    with open(namein,'r') as file:
        basisdat = file.read()

    basisdat = re.sub('%s(\s*)0' %atomold  ,r'%s\1 0'%(atomnew) , basisdat) 

    def replfunc(match):
        new = float(match.group(2))* scale
        return match.group(1) + '%.7f' %new +  match.group(3)

    basisdat = re.sub(r'(?m)^(\s*)([-\d.]+)(\s*[-\d.]+)' , replfunc , basisdat)

    with open(nameout, 'w') as file2:
        file2.write(basisdat)
        


def atom_main():
  doci_energies = [] ; hf_energies = [] ; fci_energies = []; cisd_energies= [] ; cisddoci_energies = []
  flowname = "flow.dat"; detfile1= "cisddeterminants.dat" ; detfile2 = "cisddocideterminants.dat"
  methods = ["doci", "file" , detfile1,"file", detfile2]
  basissets = ["isobasis"]# ,"cc-pvdz"]
  name = "iso_electronic_scaled"
  positions = [(0 , 0 , 0)]
  atomlist = range(10,28) #the iso electronic serie

  for basis in basissets:
      with open(name + "_" + basis + ".dat" , 'w') as f:
        f.write('#R\tHF\tDOCI\tCISD\tCISDDOCI\tFCI\n')
        for index , atom in enumerate(atomlist): #fun fact index equals charge here :)
          if index > 0:
              scale_basis('isobasissave.gbs' , 'isobasis.gbs',atom/10. * atom/10., periodic[10].symbol , periodic[atom].symbol)

          input_psi("input.dat", basis ,name= ' ' , charge_mult = (index,1) , atomlist = [atom] , positions = positions , units = 'au' , ref = 'rhf', userbasis = True)
          subprocess.call(["psi4", "input.dat"])

          with open('output.dat' , 'r') as ifile:
              data = ifile.read()

          nup =int(re.search(r'\nNalpha\s*=\s*(\d+)', data).group(1))
          ndown = int(re.search(r'\nNbeta\s*=\s*(\d+)', data).group(1))
          norb = int(re.search(r'\nNumber Of Molecular Orbitals\s*=\s*(\d+)', data).group(1))
          dw.cimain(nup,ndown,norb, [1,2], [],fname = detfile1,ref = [dw.get_hf_det])
          dw.cimain(nup,ndown,norb,[1,2], [0] , fname = detfile2,ref = [dw.get_hf_det])

          create_ciflow_input_file("output.dat" , methods , fname = flowname)
          doci_en =""
          process = subprocess.Popen(["./bin/ciflowi.x"] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
          for line in iter(process.stdout.readline, ''):
              sys.stdout.write(line)
              doci_en += line
          match_doci = re.search(r'DOCI\s*energy\s*:\s*([\-+]?\d*[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
          match_hf = re.search(r'HF\s*energy[\(\)\.\s\w\d]+:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
          #match_fci = re.search(r'CIenergy\s*big\s*[\w\s\d]+\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
          #match_fci = re.search(r'FCI\s*energy\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
          match_cisd = re.search(r'CIenergy[\w\s\d]+cisddeterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
          match_cisddoci = re.search(r'CIenergy[\w\s\d]+cisddocideterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)', doci_en)
          print "The atom: " , periodic[atom].symbol ,  " HF energy == ", match_hf.group(1), " Doci energy = " , match_doci.group(1) , "CISD energy = " , match_cisd.group(1) , "CISDDOCI energy = " , match_cisddoci.group(1)," FCI energy = " , 0 #match_fci.group(1) 
          doci_energies.append(float(match_doci.group(1)) )
          hf_energies.append(float(match_hf.group(1)) )
          #fci_energies.append(float(match_fci.group(1)) )
          fci_energies.append(0.0)
          cisd_energies.append(float(match_cisd.group(1)) )
          cisddoci_energies.append(float(match_cisddoci.group(1)) )
          f.write("%d\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n" %(atom, float(match_hf.group(1)), float(match_doci.group(1)), float(match_cisd.group(1)),float(match_cisddoci.group(1)), 0)) #float(match_fci.group(1))))
  if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
    pf.plot_data(name + "_" + basis + ".dat", xlabel = "atom number (a.u.)", ylabel ="E (Hartree)" )

def mainsize(*args , **kwargs):
  doci_energies = [] ; hf_energies = [] ; fci_energies = []; cisd_energies= [] ; cisddoci_energies = []
  flowname = "flow.dat"; detfile1= "cisddeterminants.dat" ; detfile2 = "cisddocideterminants.dat"
  methods = ["doci", "file" , detfile1,"file", detfile2]
  basissets = ["sto-3g"]# ,"cc-pvdz"]
  name = 'beh2size'
  #runlist = list(np.arange(0.6 , 1.8 , 0.1)) #+ list(np.arange(4,11,1))
  runlist = [0.86374047590, 1.02249365060, 1.18124682530, 1.34000000000, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170 ] #BeH2
  #runlist = [0.6402695456, 0.7461049954, 0.8519404452, 0.957775895, 1.0636113448, 1.1694467946, 1.2752822444, 1.3811176942, 1.486953144, 1.5927885938, 1.6986240436, 1.8044594934, 1.9102949432, 2.016130393, 2.1219658428, 2.2278012926, 2.3336367424 ] #H2O
  #runlist = [0.848, 0.973, 1.098, 1.223, 1.348, 1.473, 1.598, 1.723, 1.848, 1.973] #N2

  nup = 3; ndown = 3 ; norb = 7 ; exlist = [1,2] ; senlist = [] 

  for basis in basissets:
      with open(name + "_" + basis + ".dat" , 'w') as f:
        f.write('#R\tHF\t2DOCI\t2CISD\t2CISDDOCI\t2FCI\tscHF\tscDOCI\tscCISD\tscCISDDOCI\tFCI\n')
        for r in runlist:
            positions = [[0,0,0],[0,0,r],[0,0,-1*r]]
            atoms = [4,1,1]
            input_psi("input.dat", basis ,name= ' ' , charge_mult = (0,1) , atomlist = atoms , positions = positions , units = 'angstrom' , ref = 'rhf', userbasis = False)
            subprocess.call(["psi4", "input.dat"])
            for i in range(2):
                if i == 0:
                    dw.cimain(nup,ndown,norb,exlist,senlist ,fname = detfile1,ref = [dw.get_hf_det])
                    dw.cimain(nup,ndown,norb,[1,2], [0] , fname = detfile2,ref = [dw.get_hf_det])
                    create_ciflow_input_file("output.dat" , methods , fname = flowname)
                    doci_en =""
                    process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line)
                        doci_en += line
                    match_doci = re.search(r'DOCI\s*energy\s*:\s*([\-+]?\d*[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_hf = re.search(r'HF\s*energy[\(\)\.\s\w\d]+:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    #match_fci = re.search(r'FCI\s*energy\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_cisd = re.search(r'CIenergy[\w\s\d]+cisddeterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_cisddoci = re.search(r'CIenergy[\w\s\d]+cisddocideterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)', doci_en)
                    print "The distance: " , r ,  " HF energy == ", match_hf.group(1), " Doci energy = " , 2*match_doci.group(1) , "CISD energy = " ,2* match_cisd.group(1) , "CISDDOCI energy = " ,2* match_cisddoci.group(1)#," FCI energy = " ,2* match_fci.group(1) 
                    #f.write("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f" %(r, 2*float(match_hf.group(1)), 2*float(match_doci.group(1)),2* float(match_cisd.group(1)),2*float(match_cisddoci.group(1)),2* float(match_fci.group(1))))
                    f.write("%.15f\t%.15f\t%.15f\t%.15f\t%.15f" %(r, 2*float(match_hf.group(1)), 2*float(match_doci.group(1)),2* float(match_cisd.group(1)),2*float(match_cisddoci.group(1))))

                if i == 1:
                    reader = pr.PsiReader('output.dat', isbig = False)
                    reader.create_ni_system(outname = 'output.dat')
                    dw.cimain(nup*2,ndown*2,norb*2,[1,2], [] , fname = detfile1,ref = [lambda nup,ndown , norb : (('0'*((norb-nup)/2)+'1' * (nup/2)) *2, ('0'*((norb-ndown)/2)+'1' * (ndown/2) ) *2 )])
                    dw.cimain(nup*2,ndown*2,norb*2,[1,2], [0] , fname = detfile2,ref = [lambda nup,ndown , norb : (('0'*((norb-nup)/2)+'1' * (nup/2)) *2, ('0'*((norb-ndown)/2)+'1' * (ndown/2) ) *2 )])
                    create_ciflow_input_file("output.dat" , methods , fname = flowname)
                    doci_en =""
                    process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line)
                        doci_en += line
                    match_doci = re.search(r'DOCI\s*energy\s*:\s*([\-+]?\d*[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_hf = re.search(r'HF\s*energy[\(\)\.\s\w\d]+:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    #match_fci = re.search(r'FCI\s*energy\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_cisd = re.search(r'CIenergy[\w\s\d]+cisddeterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)' , doci_en)
                    match_cisddoci = re.search(r'CIenergy[\w\s\d]+cisddocideterminants.dat\s*:\s*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)', doci_en)
                    print "The distance: " , r ,  " HF energy == ", match_hf.group(1), " Doci energy = " , 2*match_doci.group(1) , "CISD energy = " ,2* match_cisd.group(1) , "CISDDOCI energy = " ,2* match_cisddoci.group(1)#," FCI energy = " ,2* match_fci.group(1) 
                    #f.write("%.15f\t%.15f\t%.15f\t%.15f\t%.15f\n" %(float(match_hf.group(1)), float(match_doci.group(1)), float(match_cisd.group(1)),float(match_cisddoci.group(1)), float(match_fci.group(1))))
                    f.write("\t%.15f\t%.15f\t%.15f\t%.15f\n" %(float(match_hf.group(1)), float(match_doci.group(1)), float(match_cisd.group(1)),float(match_cisddoci.group(1))))

def process_output(ciflow_out, regexp =r'[\w\s\d.()]*[eE]nergy[\w\s\d\-.()+]*:\s*([\-+]?\d*[\.,]?\d+[eEDd]?[\-+]?\d*)' , func = lambda x : float(x)):
    energies = []
    for line in ciflow_out.split('\n'):
        match = re.search(regexp ,line)
        if match:
            energies.append('%.15f' % func(match.group(1)))
        else:
            #print 'No match'
            pass
    return energies       

def print_output(afh , energies, methods):
    print 'R: ' + str(afh) + ' ||'
    print methods
    print energies
            
def create_header(afh, methods, psienergies, extra = None):
    string = '#' + afh + '    '+ 'HF' + '    '
    itmethods = iter(methods)
    methodstring = []
    calculation = ""
    method = itmethods.next()
    try:
        while True :
            calculation =  method.upper()   
            if method == "file" or  method == "big":
                method = itmethods.next()
                calculation += '(' + method.upper() + ')'
            methodstring.append(calculation)
            method = itmethods.next()    
            while  method != "doci" and method != "fci" and method != "file" and method != "big":
                    methodstring.append(calculation+ '$_{'+method.upper()+'}$')
                    method = itmethods.next()    
    except StopIteration:
        print 'created header'

    if extra:
        for pair in extra:
            methodstring.insert(pair[0] , pair[1])#insert works with index and then value.
    string += '    '.join(methodstring) 
    string += '    ' + '    '.join(psienergies) + '\n'           
    return string        

def main_opt(*args , **kwargs):
    ciflowoutputfile = "ciflowoutput2.txt" 
    flowname = "flow.dat"; detfile1= "determinants1.dat" ; detfile2 = "determinants2.dat" ; detfile3 = "determinants3.dat" ; detfile4 = "determinants4.dat" ; detfile5 = "determinants5.dat" ; detfile6 = "determinants6.dat"; detfile7 = 'frozencisdpibenzene.dat' ; detfile8 = 'cisdpibenzenesigma.dat' ; detfile9= 'CISDirrepref+cisd.dat' ; detfile10 = "CISDT(P).dat" ; detfile11 = "lin_hyb.dat" ;  detfile12 = "ref2.dat"
    detfile14 = 'cisdtdeterminants.dat'
    detfile15 = 'cisdtqdeterminants.dat'
    detfile12 = 'linhybcisdd'
    detfile20 = 'sen0.dat'
    detfile21 = 'sen2.dat'
    detfile22 = 'sen4.dat'
    detfile23 = 'sen6.dat'
    detfile24 = 'sen8.dat'
    detfile25 = 'sen10.dat'
    detfile30 = 'sen0-4.dat'
    detfile100 = 'cisddbar.dat'
    detfile200 = 'cisddbartbar.dat'
    change_ints = False#if you changed the standard psi4 integrals somewhere during the process
    #flowname = "flow.dat"; detfile1= "CIS(P)" ; detfile2 = "CISD(P)" ; detfile3 = "CISDT(P)" ; detfile4 = "CISDTQ(P)"
    #methods = ["doci","none" ] #,  "file", detfile2 , "none"]#, "file" , detfile3, "none" , "file" , detfile4 , "none" ,"file" , detfile5,"none"]#, "file" , detfile6 , "none"] #everything on
    psienergies = [] #provide here a list of energies to be calculated by psi.
    methods = ["file", detfile2 , "file" , detfile3 , "file" , detfile4] + ["file" ,detfile100 , "file" ,detfile200 ]+["doci"  , "local" , "fno" , "sim" , "fmmin" , "file" , detfile21 , "file" , detfile22 , "file",  detfile23 , "fci" ] 
    extra = None
    #methods = [  "file" , detfile1, "sim", "file", detfile2 , "sim" ,"file" ,detfile3 , "sim" , "file", detfile4, "sim", "fci"]
    basissets = ['sto-3g']
    import numpy as np
    runlist = np.arange(0.4 , 8.1, 0.01) 
    atoms = [6 , 8]
    #atoms = [6,6,6,6,6,6,1,1,1,1,1,1] #benzene
    #runlist = [x*math.pi/180. for x in range(55,61) ]#+ [59.5]] #, 60]  ] #angle between two triangles in benzene.
    #runlist = [2.224207488557] 
    name = 'noplustruncmethmulcompnisystem'
    rootdir = 'results/noplusmulnisystemtruncmeth' #relative to the directory that contains this script
    #rootdir = './results/phdsenebeccpvdzfcimminc1' #relative to the directory that contains this script
    exe = 'ciflow.x'
    #elemdir = 'matrixelements'
    #elemdir = 'output_run'
    elemdir = '.'

    #generate_dir(rootdir , exe )
    generate_dir(rootdir , exe , args= ["10000psioutput.dat", exe])
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    #create_matrix_elements(elemdir , basissets, runlist, atoms, chmult = (0,1) , moltype = 'dimer', package = 'psi' , units = 'bohr', path_mo = '../../../../mointegrals/mointegrals.so' , DOCC = None, energies = psienergies, sym = 'c1', hdf5 = False, guess = 'read', extrapar = 0.741, ref = 'rhf', su = False, basispath = '../../../data/basissets/')#, extrapar = 1.398) #extrapar is size for benzene, and angle for c2v, benzene extrapar = 1.398 C-C distance , extrapar = 104.479848 for angle h2o, 0.741 = extrapar for h2 equilibrium geometry in but
    print 'current dir:' , os.getcwd()

    outputfile = open(ciflowoutputfile , 'w')
    #fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.[m]?out' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]+\d+[eEd]?[\-+]?\d*)[-\w\d_]*\.[m]?dat' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'FCIunit(\d+)\.[m]?dat' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'.*-[\w\d]*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)orthon\.h5' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'hamnoplussto-3gpatrick([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)new\.out' , x).group(1))
    fileinfo = lambda x: 0 

    hamfiles = {}
    #search = r'psi.+%s.+mout' 
    search = r'psioutput'
    #search = r'hampsi0_%s.+mmin0.dat'
    #search = r'hampsi0_%s.+local2.dat'
    #search = 'hamnoplus.*%s.+out'
    #search = r'hampsi.+%s.+.FCIhmmin0.dat'
    #search = r'randomham.+\d+\.dat' 
    #search = r'ham.+%s.+out' 
    #search = r'hampsi0.+%s.+smmind0.dat' 
    #search = r'ham.+%s.+FCIunit\d+\.dat' 
    #search = r'psi.+%s.+orthon.h5'
    for basis in basissets:
        #hamfiles[basis] = fc.File_Collector(elemdir , search = search %basis ,notsearch = r'(\.sw\w)|(unitary)',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1. and fileinfo(x) <= 1000. )
        hamfiles[basis] = fc.File_Collector(elemdir , search = search  ,notsearch = r'(\.sw\w)|(unitary)',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1. and fileinfo(x) <= 1000. )
        #hamfiles[basis] = fc.File_Collector( '.', search = search %basis ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) <= 10. and fileinfo(x) >= 2  and fileinfo(x) in runlist  and not '.sw' in x  )
        #hamfiles[basis] = fc.File_Collector(elemdir , search = search ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000. )

    #fileinfo2 = lambda x: float(re.search(r'psi0_sto-3g1.29outputfciunit(\d+).dat' , x).group(1))
    #search = r'psi.+outputfciunit\d+.dat' 
    #outputfilesfci = fc.File_Collector('output_files', search = search  ,notsearch = r'\.sw\w',sortfunction = fileinfo2, filterf =  lambda x : fileinfo2(x) >= 0 and fileinfo2(x) < 10000. )
    #hamfiles[basissets[0]].plotfiles += ['ni_system.dat']
    #print hamfiles[basissets[0]]

    afh = 'R'
    for basis in basissets:
        fname = name + "_" +basis + ".dat"
        psir = rp.PsiReader(hamfiles[basis].plotfiles[0], isbig = False, numorbs = -1 , read_ints = False) #just to get number of electrons and orbitals
        dw.generate_all_ex(psir.values['nalpha'],psir.values['nbeta'],psir.values['norb'], dw.get_hf_det(psir.values['nalpha'],psir.values['nbeta'],psir.values['norb'])  , aname = 'determinants', addfrozen = 0)
        #( psir.get_hf_orbs()[0][:], psir.get_hf_orbs()[1][:] ) 
        #dw.biggest_det_ex(outputfilesfci.plotfiles[index])
        dw.generate_all_sen(psir.values['nalpha'],psir.values['nbeta'],psir.values['norb'], 'sen', addfrozen = 0)
        #reflist = dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'], [1] , [] , pairex = True)
        dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2] , [2]], [] ,fname = detfile100 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] ) #CISD
        dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2] , [2,3]], [] ,fname = detfile200 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] ) #CISD
        #dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2,3,4] , []], [] ,fname = detfile15 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] ) #CISD
        #dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[] ,[]], [0,2] ,fname = detfile2 ,ref = [lambda x , y , z : dw.get_hf_det(x,y,z)]) #SEN0-2
        with open(fname , 'w') as f:
            header = create_header(afh , methods, psienergies, extra = extra)
            print header
            f.write(header)
            for index , matrixelements in enumerate(hamfiles[basis].plotfiles):
                print matrixelements
                if change_ints: #True if you want to change the standard psi4 integrals.    
                    psir = rp.PsiReader(matrixelements, isbig = True, numorbs = 60 , read_ints = True) #for benzene.
                    if '1.05' not in matrixelements:
                        #psir.set_active_space([1,2] , [1,2] , [] , detfile7) #for pisystem of deformed benzene
                        psir.set_active_space([1] , [[1,2],[]] , [] , detfile7) #for pisystem of deformed benzene
                    else:
                        #psir.set_active_space([2,3,4,5] , [1,2] , [] , detfile7) #for pisystem of benzene
                        psir.set_active_space([1] , [[1,2],[2]] , [] , detfile7) #for pisystem of benzene
                    #psir.set_active_space([0,1,2,3] , [1] , [] , detfile3) #for pisystem of deformed benzene
                    #psir.set_active_space([1,2] , [1,2] , [] , detfile4) #for pisystem of deformed benzene
                    #psir.keep_orbs(56)
                    newname = os.path.basename(matrixelements)+ 'new.dat'
                    psir.create_output(newname)
                    matrixelements = newname

                create_ciflow_input_file(matrixelements , methods , fname = flowname)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass
                energies = process_output(ci_flow) 
                mullikencharge = process_output(ci_flow, regexp = "Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)") 
                if psienergies: #True if the list is not empty 
                    psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1 , read_ints = False)
                    energies +=  [x[1] for x in psir.extract_energies()]

                print_output(matrixelements, energies , methods)
                f.write("%.15f\t%s\t%s\n" %(fileinfo(matrixelements), '\t'.join(energies ), '\t'.join(mullikencharge)) )
                if change_ints:
                    os.remove(matrixelements)
  
    shutil.copy('../../htest.py' , '.')
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)
    outputfile.close()
    try:
        generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
        os.chdir('..')
        generate_dir('unitaries', None, prefix = 'unitary_')
        os.chdir('..')
        generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
        os.chdir('..')
    except:
        pass
    
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
 

def prepare_dmrg(name):
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*).*\.out' , x).group(1))
    elemdir = "./results/h2o_ccpvdz"
    basis = "cc-pvdz"
    hamfiles = fc.File_Collector(elemdir , search = r'psi.+%s.+out' %basis ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) > -1. )

    with open(name, 'w') as f:
        for ham in hamfiles.plotfiles:
            f.write("%s\n" % ham)
        f.write("end")    

def benzene(angle,size, ch_dist = 1.1):
    """
    Returns the coordinates of a deformed benzene molecule on a circle of radius size (for benzene size is also the C-C bond length), when the two reversed with equal sides triangles in the davidstar of benzene are rotated to each other with angle: angle , when angle = 60 degrees = pi/3. we get the coordinates of benzene.
    REMARK: provide angle in radials
                    psir.swap_orbitals(16,18)
                    psir.swap_orbitals(29,23)
                    if index == 0:
                        psir.swap_orbitals(36,24)
                        psir.swap_orbitals(41,25)
                        psir.swap_orbitals(42,26)
                        psir.swap_orbitals(44,27)
                        psir.swap_orbitals(45,28)
                        psir.swap_orbitals(48,29)
                    elif index == 1:
                        psir.swap_orbitals(36,24)
                        psir.swap_orbitals(41,25)
                        psir.swap_orbitals(42,26)
                        psir.swap_orbitals(44,27)
                        psir.swap_orbitals(45,28)
                        psir.swap_orbitals(48,29)
                    elif index ==2:
                        psir.swap_orbitals(36,24)
                        psir.swap_orbitals(39,25)
                        psir.swap_orbitals(40,26)
                        psir.swap_orbitals(44,27)
                        psir.swap_orbitals(45,28)
                        psir.swap_orbitals(48,29)
                    elif index ==3:
                        psir.swap_orbitals(36,24)
                        psir.swap_orbitals(40,25)
                        psir.swap_orbitals(41,26)
                        psir.swap_orbitals(44,27)
                        psir.swap_orbitals(45,28)
                        psir.swap_orbitals(48,29)
                    elif index ==4:
                        psir.swap_orbitals(35,24)
                        psir.swap_orbitals(40,25)
                        psir.swap_orbitals(41,26)
                        psir.swap_orbitals(46,27)
                        psir.swap_orbitals(47,28)
                        psir.swap_orbitals(48,29)
                    elif index ==5:
                        psir.swap_orbitals(35,24)
                        psir.swap_orbitals(40,25)
                        psir.swap_orbitals(41,26)
                        psir.swap_orbitals(46,27)
                        psir.swap_orbitals(47,28)
                        psir.swap_orbitals(48,29)

    """
    import numpy as np
    C = [[size * math.cos(angle/2.) , -1.*size*math.sin(angle/2.) , 0 ],[size * math.cos(angle/2.) , size*math.sin(angle/2.) ,0],[-1.*size*math.sin((math.pi/3.-angle)/2.) ,size*math.cos((math.pi/3.-angle)/2.) ,0],[-1.*size*math.sin((math.pi/3.+angle)/2.) ,size*math.cos((math.pi/3.+angle)/2.) ,0],[-1.*size*math.sin((math.pi/3.+angle)/2.) ,-1.*size*math.cos((math.pi/3.+angle)/2.) ,0],[-1.*size*math.sin((math.pi/3.-angle)/2.) ,-1.*size*math.cos((math.pi/3.-angle)/2.) ,0]]
    return C+ list(np.array(C)*(ch_dist+size)/size)
    
def test_benzene():
    import numpy as np
    from numpy import linalg as LA
    test =  benzene(55*math.pi/180.,1.398)
    print LA.norm(np.array(test[1])-np.array(test[2])) , 'has to be approx ' , 1.502
    print LA.norm(np.array(test[0])-np.array(test[1])) , 'has to be approx ' , 1.291

def farnaz(*args , **kwargs):
    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"; psienergies = [] ;  extra = None
    methods = ["doci" , "local" , "doci" , "fmmin" ]
    basissets = ['6-31g++']
    atoms = [2,2]
    name = 'he2varrg'
    rootdir = './results/he2varrg/' #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    runlist = [ 0.7,  0.86374047590, 0.9, 1. , 1.18124682530, 1.2, 1.34000000000, 1.38, 1.4, 1.45, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350,2.2 ,  2.29251904820, 2.35 , 2.45127222290, 2.5 , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0 ,  3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)

    generate_dir(rootdir , exe)
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    create_matrix_elements(elemdir , basissets, runlist, atoms, chmult = (0,1) , moltype = 'dimer', package = 'psi' , units = 'angstrom', path_mo = '../../../../mointegrals/mointegrals.so' , DOCC = None, energies = psienergies, sym = 'c1', hdf5 = False, guess = 'sad', extrapar = 104.479848, ref = 'rhf', userbasis = False, basispath = '../../../data/basissets/', su = True)#, extrapar = 1.398) #extrapar is size for benzene, and angle for c2v, benzene extrapar = 1.398 C-C distance

    outputfile = open(ciflowoutputfile , 'w')
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.[m]?out' , x).group(1))
    hamfiles = {}
    search = r'psi.+%s.+mout' 
    for basis in basissets:
        hamfiles[basis] = fc.File_Collector(elemdir , search = search %basis ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000. )
  
    afh = 'R'
    for basis in basissets:
        fname = name + "_" +basis + ".dat"
        with open(fname , 'w') as f:
            header = create_header(afh , methods, psienergies, extra = extra)
            print header
            f.write(header)
            for index , matrixelements in enumerate(hamfiles[basis].plotfiles):
                print matrixelements

                create_ciflow_input_file(matrixelements , methods , fname = flowname)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass
                energies = process_output(ci_flow) 

                print_output(matrixelements, energies , methods)
                f.write("%.15f\t%s\n" %(runlist[index], '\t'.join(energies )) )
  
    
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)
    outputfile.close()
    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
    generate_dir('unitaries', None, prefix = 'unitary_')
    os.chdir('..')
    generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
    os.chdir('..')
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
 

def farnaz2(*args , **kwargs):
    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"; psienergies = [] ;  extra = None
    methods = ["fci" , "none"]
    basissets = ['basisset']
    atoms = [2]
    name = 'hechangedrepcontinue'
    rootdir = './results/hechangedrepcontinue/' #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    import numpy as np
    runlist = np.arange(-1., 1. , 0.1)

    generate_dir(rootdir , exe)
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    create_matrix_elements(elemdir , basissets, [2.], atoms, chmult = (0,1) , moltype = 'atom', package = 'psi' , units = 'angstrom', path_mo = '../../../../mointegrals/mointegrals.so' , DOCC = None, energies = psienergies, sym = 'c1', hdf5 = False, guess = 'sad', extrapar = 104.479848, ref = 'rhf', userbasis = True, basispath = '../../../data/basissets/', su = True)#, extrapar = 1.398) #extrapar is size for benzene, and angle for c2v, benzene extrapar = 1.398 C-C distance

    outputfile = open(ciflowoutputfile , 'w')
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.[m]?out' , x).group(1))
    hamfiles = {}
    search = r'psi.+%s.+mout' 
    for basis in basissets:
        hamfiles[basis] = fc.File_Collector(elemdir , search = search %basis ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000. )
  
    print hamfiles[basissets[0]].plotfiles

    afh = 'R'
    for basis in basissets:
        fname = name + "_" +basis + ".dat"
        with open(fname , 'w') as f:
            header = create_header(afh , methods, psienergies, extra = extra)
            print header
            f.write(header)
            for fac in runlist:
                psir = rp.PsiReader(hamfiles[basissets[0]].plotfiles[0], isbig = False, numorbs = 60 , read_ints = True) #for benzene.
                psir.change_repulsion(factor = fac)
                newname = 'new'+ str(fac)    +'.dat'
                psir.create_output(newname)
                matrixelements = newname
                print matrixelements

                create_ciflow_input_file(matrixelements , methods , fname = flowname)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass
                energies = process_output(ci_flow) 

                print_output(matrixelements, energies , methods)
                f.write("%.15f\t%s\n" %(fac, '\t'.join(energies )) )
  
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)
    
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
 
    outputfile.close()
    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
    generate_dir('unitaries', None, prefix = 'unitary_')
    os.chdir('..')
    generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
    os.chdir('..')

def create_modham_input(elemdir, nalpha , norb , modtype , options , params , runlist, filename = None):
    olddir = os.getcwd()
    generate_dir(elemdir,None, args = [filename])
    #generate_dir(elemdir,None)
    filenum = 0
    for r in runlist:
        params[-1] = r
        hub1d = rp.ModHam(nalpha, nalpha,norb,modtype , options , params, matrixelements = filename)
        hub1d.write_file(fname = str(filenum) + modtype + 'nalpha' + str(nalpha) + '_' + 'norb' + str(norb)+ '_'  + '_'.join(map(str,params)) +'run=' + str(r)+ '.mod')

        filenum += 1        

    os.chdir(olddir)    


import tarfile
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename+'.tar.gz', "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))



def hubloop():
    nalpharange = [1 , 2, 3]
    optionrange = {2 : ['per', 'none']}
    options = ['notzero' , 'pos', 'per'] 
    for alpha in nalpharange:
        for key, value in optionrange.iteritems():
            for val in value:
                options[key] = val
                hub1d(alpha, options)
            


def hub1d(nalpha, options):
    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"
    methods = ["fci" , "none"]
    name = 'hub1dmpdper'
    rootdir = './results/hub1dmpdnorb6pair' + str(nalpha) + options[2]  + 'U0-100'+'/' #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    import numpy as np
    runlist = range(0,101)

    olddir = os.getcwd()
    generate_dir(rootdir , exe)
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.

    norb  = 6 ;  
    modtype = 'Hub1d' ; params = [1. , 4.] 
    create_modham_input(elemdir, nalpha , norb , modtype , options , params , runlist)


    outputfile = open(ciflowoutputfile , 'w')
    fileinfo = lambda x: float(re.search(r'run[-\w\d=]*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)\.mod' , x).group(1))
    hamfiles = {}
    search = modtype+ '.+mod' 
    hamfiles= fc.File_Collector(elemdir , search = search  ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000. )
  
    print hamfiles.plotfiles

    afh = 'R'
    fname = name + "_" + ".dat"
    with open(fname , 'w') as f:
        header = create_header(afh , methods, [] , extra = None)
        print header
        f.write(header)
        for index , matrixelements in enumerate(hamfiles.plotfiles):
            print matrixelements
            create_ciflow_input_file(matrixelements , methods , fname = flowname )
            ci_flow =""
            if methods:
                try:
                    ci_flow =""
                    process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line) #writes intermediate output to screen.
                        outputfile.write(line)
                        outputfile.flush()
                        ci_flow += line
                except Exception as e:    
                    print 'Ciflow gave the following error', e
                    print ci_flow
                    pass
            energies = process_output(ci_flow) 

            print_output(matrixelements, energies , methods)
            f.write("%.15f\t%s\n" %(runlist[index], '\t'.join(energies )) )
  
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)

    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
    make_tarfile( 'guillaumehub1_na_'+str(nalpha)+'_norb_'+ str(norb) + '_U_' + str(params[1]) + options[2], 'output_files')
    
 
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        #plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
    os.chdir(olddir)


def con_dm_wrap():
    fileinfo = lambda x: float(re.search(r'hamnoplussto-3gpatrick([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)new\.out' , x).group(1))
    search = r'hamnoplussto-3gpatrick.*new\.out' 
    hamfiles = [i for i in os.listdir('ham_patrick') if(  ( re.search(search , i) and fileinfo(i) <= 10. and fileinfo(i) >= 2  and fileinfo(i) * 2 % 2 == 0 and not '.sw' in i  and not fileinfo(i) == 5.)) ]
    hamfiles = sorted(hamfiles , key = fileinfo  )
    hamfiles 
    print hamfiles
    for mat in hamfiles:
        con_dm(mat)

def con_dm(mat):
    matrixelements = "eqpsioutput.dat"
    #matrixelements = "hamnoplussto-3gpatrick100.0.out"
    #matrixelements = "ni_system.dat"
    matrixelements = mat

    detfile1= "cisddeterminants.dat" 
    ciflowoutputfile = "ciflowoutput2.txt" ; flowname = "flow.dat"
    name = 'noconstrainedm3'
    methods = ["fci"]
    methods = ["file" , detfile1]
    #rootdir = os.path.join('./results/' ,  mat + 'noplusconstrainednatomddpsi4cisdt6-31g/')
    #methods = ["doci", "local"]
    rootdir = os.path.join('./results/' ,  mat + 'noplusconstrainednatomddpsi4cisdtq/')
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    import numpy as np
    runlist = list(np.arange(5.,8.01,0.2)) 


    olddir = os.getcwd()
    generate_dir(rootdir , exe,  args = [matrixelements] )
    #matrixelements = os.path.split(matrixelements)[1]
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #shutil.copy(matrixelements, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir)
    outputfile = open(ciflowoutputfile , 'w')
    os.mkdir('output_files')
  
    afh = 'R'
    fname = name + "_" + ".dat"
    #solsave = 2.
    #width = 10.
    with open(fname , 'w') as f:
        header = create_header(afh , methods, [] , extra = None).replace('\n', '') + "\tlambda" +"\tMulliken_A\n"
        print header
        f.write(header)
        psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1 , read_ints = False) #just to get number of electrons and orbitals
        dw.cimain(psir.values['nalpha']-2,psir.values['nbeta']-2 ,psir.values['norb']-2,[[1,2,3,4] , []], [] ,fname = detfile1 ,ref =  [lambda x , y , z : dw.get_hf_det(psir.values['nalpha']-2, psir.values['nbeta']-2,psir.values['norb'] -2)] , add_frozen = 2) #CISD
        #dw.cimain(psir.values['nalpha']-2,psir.values['nbeta']-2 ,psir.values['norb']-2,[[1,2,3] , []], [] ,fname = detfile1 ,ref =  [lambda x , y , z : dw.get_hf_det(psir.values['nalpha']-2, psir.values['nbeta']-2,psir.values['norb']-2 )] , add_frozen = 2) #CISD
        #for index , matrixelements in enumerate(hamfiles.plotfiles):
        for r in runlist:
            norb  = 10 ;  
            nalpha = 7;

            def func(x , pri = False):
                print 'r ' , r , ' x ' , x
                modtype = 'Constrained_DM' ; params = [0. , 1. , 2. , 3. , 4. , 5,6,7,8, r , x]  ; options = [ "lj"]
                condm = rp.ModHam(nalpha, nalpha,norb,modtype , options , params, matrixelements =matrixelements)
                condm.write_file(fname = 'conelements.mod')
                create_ciflow_input_file('conelements.mod', methods , fname = flowname, prin = pri)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass

                energies = process_output(ci_flow) 
                mullikencharge = process_output(ci_flow, regexp = "Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)") 
                return (energies , [ abs(float(m) - r) for m in mullikencharge] )
                #return ( mullikencharge, energies )

            if r <= 6.:
                solsave = gss2(func, 4. , -1., tol=1e-8)
            else: 
                solsave = gss2(func, 1. , -5. , tol=1e-8)
            energies, mullikencharge = func(solsave , pri = True)
            #shutil.copy( 'conelementsoutputdocilocal0.dat' , os.path.join('output_files' , 'conelementsoutputdocilocal' + str(r)+ '.dat') )
            #shutil.copy( 'conelementsoutputfci.dat' , os.path.join('output_files' , 'conelementsoutputfci' + str(r)+ '.dat') )
            shutil.copy( 'conelementscisddeterminantsoutputci_file.dat' , os.path.join('output_files' , 'conelementsoutputcisd' + str(r)+ '.dat') )
            #shutil.copy( 'conelementsoutputdocilocal0.dat' , os.path.join('output_files' , 'conelementsoutputdocilocal' + str(r)+ '.dat') )

            print_output('searchfor:'+str(r), energies , methods)
            f.write("%.15f\t%s\t%f\t%f\n" %(r, '\t'.join(energies ), solsave , mullikencharge[-1] ) ) 
  
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)

    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
 
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        #plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
    os.chdir(olddir)


'''python program for golden section search'''
gr=(math.sqrt(5)-1)/2
def gss(f,a,b,tol=1e-5):
    '''
    golden section search
    to find the minimum of f on [a,b]
    f: a strictly unimodal function on [a,b]

    example:
    >>> f=lambda x:(x-2)**2
    >>> x=gss(f,1,5)
    >>> x
    2.000009644875678

    '''
    c=b-gr*(b-a)
    d=a+gr*(b-a)
    fc=f(c)[1][-1];fd=f(d)[1][-1]
    while abs(c-d)>tol:       
        if fc<fd:
            b=d
            d=c  
            c=b-gr*(b-a)
            fd=f(d)[1][-1];fc=f(c)[1][-1]
        else:
            a=c
            c=d  
            d=a+gr*(b-a)
            fc=f(c)[1][-1];fd=f(d)[1][-1]
    return (b+a)/2 

def gss2(f, a, b, tol=1e-5):
    '''
    golden section search
    to find the minimum of f on [a,b]
    f: a strictly unimodal function on [a,b]

    example:
    >>> f = lambda x: (x-2)**2
    >>> x = gss(f, 1, 5)
    >>> x
    2.000009644875678

    '''
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc=f(c)[1][-1];fd=f(d)[1][-1]
    while abs(c - d) > tol:
        print fc , ' vals ' , fd
        if fc < fd:
            b = d
            c = b - gr * (b - a)
            d = a + gr * (b - a)
            fd = fc 
            fc=f(c)[1][-1]
        else:
            a = c
            c = b - gr * (b - a)
            d = a + gr * (b - a)
            fc = fd
            fd=f(d)[1][-1]

        # we recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop

    return (b + a) / 2

def con_dmloop():
    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"
    detfile1= "cisddeterminants.dat" 
    methods = ["file" , detfile1]
    name = 'noconstrainedm2'
    rootdir = './results/6bohrcisdnoplusconstrained'  #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    import numpy as np
    runlist = list(np.arange(-100 , 1 ,  1.) )

    #matrixelements = "psioutput.dat"
    #matrixelements = "hamnoplussto-3gpatrick100.0.out"
    matrixelements = "6psioutput.dat"

    olddir = os.getcwd()
    generate_dir(rootdir , exe, args = [matrixelements] )
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    outputfile = open(ciflowoutputfile , 'w')

    norb  = 10 ;  
    nalpha = 7;
    modtype = 'Constrained_DM' ; params = [0. , 1. , 2. , 3. , 4. ,  6.8 , 4.]  ; options = [ "lj"]
    create_modham_input(elemdir,nalpha  , norb , modtype , options , params , runlist ,matrixelements)


    fileinfo = lambda x: float(re.search(r'run[\=]*([-]?\d+\.\d+[e\-]*\d*)\.mod' , x).group(1))
    hamfiles = {}
    search = modtype+ '.+mod' 
    hamfiles= fc.File_Collector(elemdir , search = search  ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1000 and fileinfo(x) < 1000. )
  
    afh = 'R'
    fname = name + "_" + ".dat"
    with open(fname , 'w') as f:
        header = create_header(afh , methods, [] , extra = None)
        psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1 , read_ints = False) #just to get number of electrons and orbitals
        dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2] , []], [] ,fname = detfile1 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] ) #CISD
        print header
        f.write(header)
        for index , matrixelements in enumerate(hamfiles.plotfiles):
            create_ciflow_input_file(matrixelements, methods , fname = flowname)
            ci_flow =""
            if methods:
                try:
                    ci_flow =""
                    process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line) #writes intermediate output to screen.
                        outputfile.write(line)
                        outputfile.flush()
                        ci_flow += line
                except Exception as e:    
                    print 'Ciflow gave the following error', e
                    print ci_flow
                    pass

            energies = process_output(ci_flow) 
            mullikencharge = process_output(ci_flow, regexp = "Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)") 
            #mullikencharge = float(re.search("Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)", ci_flow).group(1) )

            print_output('R:'+str(runlist[index]), energies , methods)
            f.write("%.15f\t%s\t%s\n" %(fileinfo(matrixelements), '\t'.join(energies ) , '\t'.join(mullikencharge ) ) )
  
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)

    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
 
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        #plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
    os.chdir(olddir)

def correct_condm():
    dirname = 'results/noconstraineddmkeepwfhamnoplussto-3gpatrick5.0new.outfci'
    shutil.copy('ciflow.x', dirname)#When the matrixelements are already present.
    matrixelements = "hamnoplussto-3gpatrick5.0new.out"
    #shutil.copy('ni_system.dat', dirname)#When the matrixelements are already present.
    os.chdir(dirname)

    with open('noconstrainedm_.dat' , 'r') as file:
        data = file.read()
        lambdalist = re.findall(r'\s(-?\d*\.\d+)\s*$', data, re.M)

    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"
    methods = ["fci" ]
    name = 'noconstrainedmmulliken'
    elemdir = 'matrixelements'
    import numpy as np
    mullikenlist = list(np.arange(5. , 8.01 ,  0.01) )

    assert(len(mullikenlist) == len(lambdalist) )

    runlist = zip(mullikenlist , lambdalist)
    print runlist


    outputfile = open(ciflowoutputfile , 'w')

    afh = 'lambda'
    fname = name + "_" + ".dat"
    with open(fname , 'w') as f:
        header = create_header(afh , methods, [] , extra = None)
        print header
        f.write(header)
        for  mul ,  lam in runlist:
            norb  = 10 ;  
            nalpha = 7;
            modtype = 'Constrained_DM' ; params = [0. , 1. , 2. , 3. , 4. ,  mul , lam]  ; options = [ "lj"]
            condm = rp.ModHam(nalpha, nalpha,norb,modtype , options , params, matrixelements =matrixelements)
            condm.write_file(fname = 'conelements.mod')
            create_ciflow_input_file('conelements.mod', methods , fname = flowname)
            ci_flow =""
            if methods:
                try:
                    ci_flow =""
                    process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line) #writes intermediate output to screen.
                        outputfile.write(line)
                        outputfile.flush()
                        ci_flow += line
                except Exception as e:    
                    print 'Ciflow gave the following error', e
                    print ci_flow
                    pass

            energies = process_output(ci_flow) 
            mullikencharge = process_output(ci_flow, regexp = "Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)") 
            #mullikencharge = float(re.search("Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)", ci_flow).group(1) )

            print_output('R:'+str(mul)+'\t' + str(lam), energies , methods)
            f.write("%.15f\t%s\t%f\t%s\n" %(float(mul) , '\t'.join(energies ) , float(lam) , '\t'.join(mullikencharge ) ) )
  
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)

    generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
 
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        #plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )
    os.chdir(olddir)


def main_benzene(*args , **kwargs):
    ciflowoutputfile = "ciflowoutput.txt" 
    flowname = "flow.dat"
    change_ints = True #if you changed the standard psi4 integrals somewhere during the process
    extra = None
    detfile1 = "pisystemfci"
    psienergies = []
    methods = [  "file" , detfile1, "hmmin"]

    basissets = ['6-31g']
    runlist = [x*math.pi/180. for x in range(55,61) ] #+ [59.5]] #, 60]  ] #angle between two triangles in benzene.

    atoms = [6,6,6,6,6,6,1,1,1,1,1,1] #benzene

    name = 'benzenefcipisystem'
    rootdir = './results/benzenepisystemdeformfci2sad' #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'

    generate_dir(rootdir , exe , prefix = None)
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    create_matrix_elements(elemdir , basissets, runlist, atoms, chmult = (0,1) , moltype = 'benzene', package = 'psi' , units = 'angstrom', path_mo = '../../../../mointegrals/mointegrals.so' , DOCC = None, energies = psienergies, sym = None, hdf5 = False, guess = 'sad', extrapar = 1.398 , ref = 'rhf', su = False, basispath = '../../../data/basissets/')#, extrapar = 1.398) #extrapar is size for benzene, and angle for c2v (for h2o extrapar = 104.479848) , benzene extrapar = 1.398 C-C distance

    outputfile = open(ciflowoutputfile , 'w')
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.[m]?out' , x).group(1))

    hamfiles = {}
    search = r'psi.+%s.+mout' 
    for basis in basissets:
        hamfiles[basis] = fc.File_Collector(elemdir , search = search %basis ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1. and fileinfo(x) < 10000. )

    afh = 'R'
    for basis in basissets:
        fname = name + "_" +basis + ".dat"
        psir = rp.PsiReader(hamfiles[basis].plotfiles[0], isbig = False, numorbs = -1 , read_ints = False) #just to get number of electrons and orbitals
        #dw.cimain(7,7, 10,[[1,2,3,4,5,6],[]], [] , fname = detfile6  ,ref = [ lambda x , y , z : dw.get_hf_det(x,y ,z)]) #hyb_lin
        with open(fname , 'w') as f:
            header = create_header(afh , methods, psienergies, extra = extra)
            print header
            f.write(header)
            for index , matrixelements in enumerate(hamfiles[basis].plotfiles):
                print matrixelements
                if change_ints: #True if you want to change the standard psi4 integrals.    
                    psir = rp.PsiReader(matrixelements, isbig = True, numorbs = 60 , read_ints = True) #for benzene.
                    if '1.05' not in matrixelements:
                        psir.set_active_space([1,2] , [] , [] , detfile1) #for pisystem of deformed benzene
                    else:
                        psir.set_active_space([2,3,4,5] , [] , [] , detfile1) #for pisystem of benzene
                    #psir.keep_orbs(56)
                    newname = os.path.basename(matrixelements)+ 'new.dat'
                    psir.create_output(newname)
                    matrixelements = newname

                create_ciflow_input_file(matrixelements , methods , fname = flowname)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow.x" ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass
                energies = process_output(ci_flow) 
                mullikencharge = process_output(ci_flow, regexp = "Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)") 
                if psienergies: #True if the list is not empty 
                    psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1 , read_ints = False)
                    energies +=  [x[1] for x in psir.extract_energies()]

                print_output(matrixelements, energies , methods)
                f.write("%.15f\t%s\t%s\n" %(runlist[index], '\t'.join(energies ), '\t'.join(mullikencharge)) )
                if change_ints:
                    os.remove(matrixelements)
  
    shutil.copy('../../htest.py' , '.')
    os.remove(exe) #can be handy to keep it in the dir, to exactly reproduce results later. (but be warned its big.)
    outputfile.close()
    try:
        generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
        os.chdir('..')
        generate_dir('unitaries', None, prefix = 'unitary_')
        os.chdir('..')
        generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
        os.chdir('..')
    except:
        pass
    
    if (os.getenv('VSC_INSTITUTE_LOCAL') != 'gent'):
        plotter = pf.Plot_Files(fname)
        plotter.data[0].depvar['yas'] = 'Energy'  #change the future y-axis label 
        plotter.data[0].depvar['xas'] = '$R$'  #is normally set to the column header 0
        plotter.data[0].units['x'] = r'(\AA)'
        plotter.data[0].units['y'] = r'(E$_h$)'
        plotter.generate_plot(xlimg = None , ylimg =None , exname = '' , prefix = True, titel =  name, name = fname, depcol = 0, ylist = None )


if __name__ == "__main__":
    #hubloop()
    #hub1d()
    #main_opt()
    #main_benzene()
    #correct_condm()
    #print create_header("R" , ["fci" , "doci" , "local" , "file" , "dfddf" , "local" , "mmin" , "file" ,"dfdfdf" ], [] )
    #con_dm("hamhamnoplussto-3gpatrick10.0newDOCIsim0.dat")
    #con_dm("hamnoplussto-3gpatrick5.0new.out")
    #con_dm("hamnoplussto-3gpatrick6.0new.out")
    #con_dm("hamnoplussto-3gpatrick5.0new.out")
    #con_dm("hampsiham002.00cnminorthon.dat")
    #con_dm("hamnoplussto-3gpatrick2.0new.out")
    #con_dm("n2_2bohr.dat")
    con_dm("psioutput.dat")
    #con_dmloop()
    #con_dm_wrap()
    #farnaz()
    #farnaz2()
    #print benzene(math.pi/3 , 1.398)
    #test_benzene()
    #prepare_dmrg('dmrginput.dat')
    #input_psi('input.dat', 'sto-3g' ,name= ' ' , charge_mult = (0,1) , atomlist = [4] , positions = [[0,0,0]] , units = 'au' , ref = 'rhf' , userbasis = False , functional = None, path_to_plugin = '../mointegrals/mointegrals.so', DOCC = None, sym = None, energies = ['cisd', 'ccsd' , 'fci'])
    #mainsize()
    #atom_main()
    #gkci()
    #scale_basis('isobasis.gbs' , 100, 'Be' , 'Ar')
