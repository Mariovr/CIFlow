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

def create_ciflow_input_file(matrixelements , methods , fname = "flow.dat"):
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
    change_ints = False#if you changed the standard psi4 integrals somewhere during the process
    #flowname = "flow.dat"; detfile1= "CIS(P)" ; detfile2 = "CISD(P)" ; detfile3 = "CISDT(P)" ; detfile4 = "CISDTQ(P)"
    #methods = ["doci","none" ] #,  "file", detfile2 , "none"]#, "file" , detfile3, "none" , "file" , detfile4 , "none" ,"file" , detfile5,"none"]#, "file" , detfile6 , "none"] #everything on
    psienergies = [] #provide here a list of energies to be calculated by psi.
    methods = ["fci", "fno" , "fci" , "hmmin" , "file" , detfile1 , "fno" ,"file" , detfile1 , "fmmin" , "file" ,detfile14 , "fmmin" , "file" ,detfile14 ,"fno" , "file" , detfile15 ,"fmmin" ,"file",  detfile15 , "fno"] #, "doci" , "fmmin"]
    methods = [  "file", detfile2 , "hmmin","file" , detfile3 ] + ["doci"  , "local" , "file" , detfile21, "file" , detfile22 ]
    methods = [  "file", detfile2 ] 
    #methods =  ["doci"  , "local" , "file" , detfile21, "file" , detfile22 ]
    #extra = [(1,"DOCI_(MMin)")]
    extra = None
    #methods = [  "file" , detfile1, "sim", "file", detfile2 , "sim" ,"file" ,detfile3 , "sim" , "file", detfile4, "sim", "fci"]
    basissets = ['6-31g']
    import numpy as np
    runlist = np.arange(1. , 8.1, 0.1) 
    atoms = [4,1,1] #beh2
    name = 'phdsenhierbeh26-c1symmetry'
    rootdir = './results/phdsenhierbeh26-31gc1symmetry' #relative to the directory that contains this script
    exe = 'ciflow.x'
    elemdir = 'matrixelements'
    #elemdir = 'matrixelements_otherbasis'
    #elemdir = 'random_hamiltonians'
    #elemdir = 'output_run'
    #elemdir = '.'

    generate_dir(rootdir , exe  )
    #generate_dir(rootdir , exe , args= ["ni_system.dat"], prefix = r'./ham_patrick/hamnoplus')
    #generate_dir(rootdir , exe )
    #generate_dir(rootdir , exe , args= ["ni_system.dat"])
    #shutil.copy(exe, rootdir)#When the matrixelements are already present.
    #os.chdir(rootdir) #When the matrixelements are already present.
    create_matrix_elements(elemdir , basissets, runlist, atoms, chmult = (0,1) , moltype = 'linear', package = 'psi' , units = 'bohr', path_mo = '../../../../mointegrals/mointegrals.so' , DOCC = None, energies = psienergies, sym = 'c1', hdf5 = False, guess = 'read', extrapar = 0.741, ref = 'rhf', su = False, basispath = '../../../data/basissets/')#, extrapar = 1.398) #extrapar is size for benzene, and angle for c2v, benzene extrapar = 1.398 C-C distance , extrapar = 104.479848 for angle h2o, 0.741 = extrapar for h2 equilibrium geometry in but

    outputfile = open(ciflowoutputfile , 'w')
    fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.[m]?out' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]+\d+[eEd]?[\-+]?\d*)[-\w\d_]*\.[m]?dat' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'FCIunit(\d+)\.[m]?dat' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'.*-[\w\d]*([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)orthon\.h5' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)[-\w\d]*\.out' , x).group(1))
    #fileinfo = lambda x: float(re.search(r'hamnoplussto-3gpatrick([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*)new\.out' , x).group(1))

    hamfiles = {}
    search = r'psi.+%s.+mout' 
    #search = r'hampsi0_%s.+sdmmin3.dat'
    #search = r'hampsi0_%s.+local2.dat'
    #search = 'hamnoplus.*%s.+out'
    #search = r'hampsi.+%s.+.FCIhmmin0.dat'
    #search = r'randomham.+\d+\.dat' 
    #search = r'ham.+%s.+out' 
    #search = r'hampsi0.+%s.+smmind0.dat' 
    #search = r'ham.+%s.+FCIunit\d+\.dat' 
    #search = r'psi.+%s.+orthon.h5'
    for basis in basissets:
        hamfiles[basis] = fc.File_Collector(elemdir , search = search %basis ,notsearch = r'(\.sw\w)|(unitary)',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -2. and fileinfo(x) <= 1000. )
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
        dw.generate_all_ex(psir.values['nalpha'],psir.values['nbeta'],psir.values['norb'], ( psir.get_hf_orbs()[0][:], psir.get_hf_orbs()[1][:] )  , aname = 'determinants', addfrozen = 0)
        #dw.biggest_det_ex(outputfilesfci.plotfiles[index])
        dw.generate_all_sen(psir.values['nalpha'],psir.values['nbeta'],psir.values['norb'], 'sen', addfrozen = 0)
        #reflist = dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'], [1] , [] , pairex = True)
        #dw.cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2] , []], [] ,fname = detfile1 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] ) #CISD
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
 

import tarfile
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename+'.tar.gz', "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

if __name__ == "__main__":
    main_opt()
    #print create_header("R" , ["fci" , "doci" , "local" , "file" , "dfddf" , "local" , "mmin" , "file" ,"dfdfdf" ], [] )
    #input_psi('input.dat', 'sto-3g' ,name= ' ' , charge_mult = (0,1) , atomlist = [4] , positions = [[0,0,0]] , units = 'au' , ref = 'rhf' , userbasis = False , functional = None, path_to_plugin = '../mointegrals/mointegrals.so', DOCC = None, sym = None, energies = ['cisd', 'ccsd' , 'fci'])
    #scale_basis('isobasis.gbs' , 100, 'Be' , 'Ar')
