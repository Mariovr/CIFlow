#!/usr/bin/env python
#
#PBS -N golettcon_dmtest4bohrbohrfci
#PBS -o golettcon_dmwf4bohrborhfci.file2
#PBS -e golettcon_dmwf4bohrborhfci.file2
#PBS -l walltime=71:59:59
#PBS -l nodes=1:ppn=16 ##PBS -q long ##possible ques are: debug (15m,1h), bshort(1h,12h), short (1h,12h), long (48h,72h), default-> short
##############################PBS -l vmem=60gb
#PBS -q long
#PBS -m abe  ##no mail put -m n 

#
## change to gent gaussian group
import os , sys,shutil,datetime 
os.system('module load matplotlib')
os.system('module load numpy')
os.system('module load h5py')
sys.path.append(os.path.join(os.getenv('HOME'), 'devel/CIFlow/rundir/'))
sys.path.append(os.path.join(os.getenv('HOME'), 'devel/CIFlow/rundir/data/'))
import subprocess
import re 
import detwrite as dw
import read_psi as rp
import htest as ht
import math

periodic = ['' , 'H' , 'He' , 'Li' , 'Be',              'B',              'C',              'N',              'O',              'F',             'Ne',             'Na',             'Mg',             'Al',             'Si',              'P',              'S',             'Cl',             'Ar',              'K',             'Ca',             'Sc',             'Ti',              'V',             'Cr',             'Mn',             'Fe',             'Co',             'Ni',             'Cu',             'Zn',             'Ga',             'Ge',             'As',             'Se',             'Br' ]

class File_Collector(object):
    def __init__(self, rootdir , search , notsearch = '.png' , notdir = 'xyvwa' , filelist = None , sortfunction = None , rev = False, filterf = None):
        if filelist != None:
            self.plotfiles = filelist
        else:
            self.plotfiles = []
        self.sortfunction = sortfunction
        self.readfiles(rootdir , search , notsearch = notsearch , notdir = notdir)
        self.sortplotfiles(rev)
        if filterf != None:
            self.plotfiles = filter(filterf , self.plotfiles)
        print self.plotfiles
  
    def addfiles(self , *args):
        for i in args:
            self.plotfiles.append(i)
    
    def sortplotfiles(self, rev = False):
        if self.sortfunction != None:
            self.plotfiles = sorted(self.plotfiles ,key = self.sortfunction , reverse = rev )
        else:
            print 'No sort function given so the order of the files doesn\'t matter for the figure'
  
    def readfiles(self, dirname , search , notsearch = 'rgvar' , notdir = 'xyvwa'):
        """
        If you want to plot data from a single file use readdata instead, this is a wrapper for readdata if you want to plot data from multiple files
        """
        print('We are in the following directory: %s looking for files that contain %s and not %s' %(dirname, search , notsearch))
        dirlist =  os.listdir(dirname)
        for filep in dirlist:
            filep = os.path.join(dirname,filep) 
            if os.path.islink(filep):
                pass
            elif os.path.isdir(filep):
                m = re.search(notdir , filep)
                if m is None:
                    self.readfiles(filep , search, notsearch = notsearch, notdir = notdir )
            elif os.path.isfile(filep) and '.swp' not in filep: 
                nm = re.search(notsearch, filep)
                m = re.search(search , filep)
                if m is not None and nm is None:
                    self.plotfiles.append(filep)
            else:
                pass

import tarfile
def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename+'.tar.gz', "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

def main():
    name = 'beh2cc-pvdzcisddocidebug'
    hamdirname = "beh2ccpvdzham"
    pokemon = os.environ['VSC_INSTITUTE_CLUSTER']
    print 'we have run this script on cluster: ' , pokemon
    active_space = False
    flowname = "flow.dat"; detfile1 = "cisddeterminants.dat" ; detfile2 = "cisddocideterminants.dat" ; detfile3 = 'CIS(P)' ; detfile4 = 'CISD(P)' ; detfile6 = 'fcifile.dat' ; detfile7 = 'be_ccpv5z96100.txt'
    detfile10= "linhyb23p.dat"
    detfile11 = "linhyb2p.dat"
    #methods = [ "big", detfile6, "loadham"] #everything on#everything on
    methods = [ "file" , detfile1 , "hmmin", "file" ,detfile10 , "loadham"] #everything on#everything on
    methods = ["doci" , "loadham" , "file", detfile3 , "loadham" , "file" , detfile4 , "loadham" ] #everything on#everything on
    #methods = []
    #psienergies = ['fci']
    psienergies = []
    basissets = ['cc-pvdz']
    #runlist = [ 3.4 , 3.5 ,  3.8 , 3.9 , 4. ,  4.5 ,4.8 ,  5., 5.6, 6.] #BeH2 (linear)
    #runlist = [1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 4. , 4.5 , 5.] #BeH2
    runlist = [0.7, 0.9 , 1, 1.1 , 1.2 , 1.3 ,1.4 , 1.5 , 1.6 , 1.7, 1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 ,3.6,3.7,3.8,3.9 , 4. , 4.5 , 5., 6.] #BeH2
    #runlist = [0.6402695456, 0.7 , 0.7461049954, 0.8 , 0.8519404452, 0.9 , 0.957775895, 1. , 1.0636113448, 1.1, 1.1694467946, 1.2 , 1.2752822444, 1.3, 1.3811176942, 1.4, 1.486953144, 1.5, 1.5927885938, 1.6, 1.6986240436 , 1.72 , 1.8 , 1.85 , 1.9 , 1.95 , 2. , 2.05, 2.2278012926, 2.4 , 2.6 , 3. , 3.5 , 4. , 4.5 , 5. ,6.]#, 2.1 , 2.15, 2.2278012926 ]#, 1.8044594934, 1.9102949432, 2.016130393, 2.1219658428, 2.2278012926, 2.3336367424 ] #H2O (c2v)
    #runlist = [2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.2,3.3, 3.5,3.8 , 4.]
    #runlist = [2.32,2.42,2.52,2.61,2.71,2.99,3.09,3.19,3.38,3.57, 3.76,4.05]
    #runlist = [1.9] +  runlist
    #runlist = [1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 4. , 4.5 , 5.] #BeH2
    #runlist = [2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 4. , 4.5 , 5.] #BeH2
    #runlist = []

    arrayid = 2
    #hereunder for array jobs for heavy calculations. But make sure it is always defined
    #arrayid =int(os.getenv('PBS_ARRAYID'))
    #arraystep = 1
    #name += str(arrayid)
    #try:
    #    runlist = runlist[(arrayid-1) * arraystep : arrayid *arraystep]
    #except:
    #    print 'encoutered exception in cut of runlist.'
    #    runlist = runlist[-4:-1]

    #runlist = [1.38, 1.4, 1.45, 1.5, 1.65750634940, 1.81625952410, 1.97501269880] #, 2.13376587350,2.2 ,  2.29251904820, 2.35 , 2.45127222290, 2.5 , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0 ,  3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)
    #runlist = [2.13376587350,2.2 ,  2.29251904820, 2.35 , 2.45127222290, 2.5 ]# , 2.61002539760,2.7 ,  2.76877857230, 2.92753174700, 3.0 ,  3.08628492170, 3.1 , 3.2, 3.3 , 3.4 , 3.5 , 3.6 ,3.7 , 3.8 , 3.9 , 4. , 4.5 , 5. , 6.] #BeH2 (linear)
    #runlist = [ 3.08628492170, 3.1 , 3.2, 3.3 , 3.4 ,  4. ,   5. , 6.] #BeH2 (linear)

    #runlist = [0.86374047590, 1.02249365060, 1.18124682530, 1.34000000000, 1.49875317470, 1.65750634940, 1.81625952410, 1.97501269880, 2.13376587350, 2.29251904820, 2.45127222290, 2.61002539760, 2.76877857230, 2.92753174700, 3.08628492170 ] #BeH2
    #runlist = [0.6402695456, 0.7461049954, 0.8519404452, 0.957775895, 1.0636113448, 1.1694467946, 1.2752822444, 1.3811176942, 1.486953144, 1.5927885938, 1.6986240436 , 2.2278012926 ]#, 1.8044594934, 1.9102949432, 2.016130393, 2.1219658428, 2.2278012926, 2.3336367424 ] #H2O
    #runlist = [ 1.2 , 1.6] #[ 0.7,0.8 ,0.8519404452, 0.957775895, 1.1694467946, 1.2752822444, 1.3811176942]#, 1.8044594934, 1.9102949432, 2.016130393, 2.1219658428, 2.2278012926, 2.3336367424 ] #H2O
    #runlist = [0.848, 0.9 , 0.973, 1.0 , 1.098, 1.223, 1.3, 1.348, 1.4, 1.473,1.5,  1.598, 1.723,1.8 ,  1.848, 1.973] #N2
    #atoms = [7,7] # set DOCC = [3, 0 ,0, 0,0 , 2, 1 ,1] 
    #atoms = [1, 1,1,1,1,1,1,1 ] #h8
    #atoms = [8, 1,1 ]
    atoms = [4,1,1]
    #atoms = [6,6,6,6,6,6,1,1,1,1,1,1] #benzene
    #import math
    #runlist = [x*math.pi/180. for x in range(55,60) + [59.9 , 60]  ] #angle between two triangles in benzene.

    print 'Start to load necessary modules.'
    os.system('module swap cluster/%s' %pokemon)
    os.system('module load PSI/4.0b6-20160201-intel-2016a-mt-Python-2.7.11')
    os.system('module load arpack-ng/3.3.0-intel-2016a')
    os.system('module load matplotlib')
    os.system('module load h5py')
    #os.system('module load Boost/1.55.0-ictce-7.1.2-Python-2.7.8')
    os.system('export MKL_NUM_THREADS=16')
    #os.system('export KMP_STACKSIZE=32m')

    ORIGDIR= os.getenv('HOME')
    mofile = os.path.join(ORIGDIR, 'devel/CIFlow/lib/mointegrals_%s.so' %pokemon)
    exefile = os.path.join(ORIGDIR, 'devel/CIFlow/bin/ciflow_%s.x' %pokemon )
    hamdir = os.path.join(ORIGDIR, hamdirname)
    #det_file = os.path.join(ORIGDIR, detfile10 )
    print 'We submitted the script in: ', ORIGDIR
    print 'Psi 4 plugin is at: ' , mofile

    WORKDIR= os.path.join(os.getenv('VSC_SCRATCH_NODE'),os.getenv('PBS_JOBID'))

    if os.path.isdir(WORKDIR):
        shutil.rmtree(WORKDIR)
    os.mkdir(WORKDIR)
    print 'Workdir is: ', WORKDIR

    shutil.copy(mofile,WORKDIR)
    shutil.copy(exefile,WORKDIR)
    hamdest = os.path.join(WORKDIR, hamdirname)
    shutil.copytree(hamdir,hamdest)
    #shutil.move(det_file,WORKDIR)

    os.chdir(WORKDIR)
    print 'We work on node: ' , os.getenv('HOSTNAME')

    #PSIINPUT = os.path.join(os.getcwd() , 'output_run')
    #print 'We Start to create the outputdir ' ,PSIINPUT
    #if os.path.isdir(PSIINPUT):
    #    shutil.rmtree(PSIINPUT)
    #os.mkdir(PSIINPUT)

    #angle = 104.479848 
    #for basis in basissets:
    #    for r in runlist:
    #        #positions = [[0,0,0],[0,0,r],[0,0,-1*r]] #for molecules like beh2
    #        #positions = [[0,0,0],[0,0,r]] #for molecules like n2
    #        #positions = [[0,0,0],[0,0,r],[0,0,2*r], [0,0,3*r], [0,0,4*r] , [0,0,5*r] , [0,0,6*r], [0,0,7*r]] #For H_8
    #        positions = [[''],[1, r],[1, r, 2 , angle]] #for molecules like h2o
    #        #positions = [[0,0,0]] #for an atom
    #        #positions = ht.benzene(r,  1.398) #REMARK we run over the angle -> r in runlist, and keep the radius of the circle on which all the C atoms are placed constant to extrapar, (extrapar = 1.398 A is standard)
    #        ht.input_psi(os.path.join('output_run' ,"psi%.2f" %(r) +"basis"+basis+  ".dat"), basis ,name= 'mol' , charge_mult = (0,1) , atomlist = atoms , positions = positions , units = 'angstrom' , ref = 'rhf', userbasis = False, path_to_plugin = './mointegrals_%s.so' %pokemon, DOCC = None, sym = None, energies = psienergies)

    print 'Contents of rundir before run are: ', os.listdir(os.getcwd())
    #for file in os.listdir(PSIINPUT):
    #    print 'processing psi input file: ' , file
    #    exfile = os.path.join('output_run',file)
    #    os.system('psi4 ' + exfile)

    #print 'Contents of rundir after run are: ', os.listdir(os.getcwd())
    #print 'Contents of inputdir are: ', os.listdir('./output_run')
    

    hamfiles = File_Collector('.', search = r'psi.+out' ,notsearch = r'\.sw\w',sortfunction = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?).*\.m?out' , x).group(1)) , filterf = lambda x: (float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*).*\.m?out' , x).group(1)) == round(runlist[arrayid],2 ) )  )
    #hamfiles =File_Collector('output_run', search = r'psi.+out' ,notsearch = r'\.sw\w', filterf = lambda x: float(re.search(r'([\-+]?\d+[\.,]?\d+[eEDd]?[\-+]?\d*).*\.m?out' , x).group(1)) > 1.1)
                
    psir = rp.PsiReader(hamfiles.plotfiles[0], isbig = False, numorbs = -1 , read_ints = False)

    #frozen = [1,0,0,0,0,1,0,0]  ; virtual = [1,1,1,1,1,2,0,0]; active = [5,0,2,2,0,4,3,3]
    #assert( sum(frozen) + sum(virtual) + sum(active) == psir.values['norb'])
    #dw.cimain(psir.values['nalpha'] - sum(frozen),psir.values['nbeta'] -sum(frozen)  ,psir.values['norb'] - sum(frozen) - sum(virtual) , [[1,2], []], [],fname = detfile1 ,ref = [lambda x , y , z : psir.get_hf_orbs(  frozen = frozen , virtual = virtual )] , active = active, frozen =  frozen, virtual = virtual ,add_frozen = 0, add_virtual = 0)
    #dw.cimain(psir.values['nalpha'] - sum(frozen),psir.values['nbeta'] -sum(frozen)  ,psir.values['norb'] - sum(frozen) - sum(virtual) , [[1,2], [2,3]], [],fname = detfile10 ,ref = [lambda x , y , z : psir.get_hf_orbs(  frozen = frozen , virtual = virtual )] , active = active, frozen =  frozen, virtual = virtual ,add_frozen = 0, add_virtual = 0)
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'] , [[1,2],[2]], [] , fname = detfile11 ,add_frozen = 0, add_virtual = 0,ref =  [ lambda x , y , z : psir.get_hf_orbs()]   ) #CISDD
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'] , [[1,2],[]], [0] , fname = detfile2 ,add_frozen = 0, add_virtual = 0,ref =  [ lambda x , y , z : psir.get_hf_orbs()]   ) #CISDDOCI
    dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[[], [1]], [] , fname = detfile3 ,add_frozen = 0, ref =  [lambda x , y , z : psir.get_hf_orbs()]) #CIS(P)
    dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[[], [1,2]], [] , fname = detfile4 ,add_frozen = 0, ref =  [lambda x , y , z : psir.get_hf_orbs()]) #CISD(P)
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[[1,2], [] ],[]  , fname = detfile1 ,add_frozen = 0, ref =  [lambda x , y , z : psir.get_hf_orbs()]) #CISD 
    #dw.cimain(psir.values['nalpha']-1,psir.values['nbeta']-1, psir.values['norb']-1,[[1,2,3,4], [] ],[]  , fname = detfile6 ,add_frozen = 1, ref =  [lambda x , y , z : dw.get_hf_det(x,y,z)] ) #FCI
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[range(1,psir.values['nbeta']*2 +1), [] ],[]  , fname = detfile6 ,add_frozen = 0, ref =  [lambda x , y , z : dw.get_hf_det(x,y,z)] ) #FCI
    afh = 'R'
    for basis in basissets:
        fname = "ciflowresults_" + basis + ".dat"
        with open(fname , 'w') as f:
            header = ht.create_header(afh , methods,psienergies)
            print header
            f.write(header)
            for index , matrixelements in enumerate(hamfiles.plotfiles):
                psir = rp.PsiReader(matrixelements, isbig = False, numorbs = -1 , read_ints = False)
                if active_space:
                    psir.set_active_space([0] , [[1,2],[]] , [] , detfile1) #for pisystem of deformed benzene
                    psir.set_active_space([0] , [[1,2], [2]]  , [] , detfile11) #for pisystem of deformed benzene
                ht.create_ciflow_input_file(matrixelements , methods , fname = flowname)
                ci_flow =""
                try:
                    process = subprocess.Popen(["./ciflow_%s.x" %pokemon ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                    for line in iter(process.stdout.readline, ''):
                        sys.stdout.write(line)
                        ci_flow += line
                except Exception as e:    
                    print 'Ciflow gave the following error', e
                    print ci_flow
                    pass

                """
                try:
                    ci_flow = subprocess.check_output(['./ciflow_%s.x' %pokemon ] , stdin =open(flowname , 'r') )
                except Exception as e:    
                    print 'Ciflow gave the following error', e
                    print ci_flow
                    pass
                """
                energies = ht.process_output(ci_flow) + [x[1] for x in psir.extract_energies()]
                ht.print_output(matrixelements, energies , methods)
                f.write("%.15f\t%s\n" %(runlist[index], '\t'.join(energies )) )
  
    os.remove('./ciflow_%s.x' %pokemon)
    os.remove('./mointegrals.so')
    #os.remove(detfile1)
    ht.generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
    ht.generate_dir('unitaries', None, prefix = 'unitary_')
    os.chdir('..')
    ht.generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
    os.chdir('..')

    OUTPUTDIR = os.path.join(ORIGDIR, name+str(arrayid) )
    print 'We Start to create the outputdir ' , OUTPUTDIR
    if os.path.isdir(OUTPUTDIR):
        shutil.rmtree(OUTPUTDIR)

    make_tarfile(OUTPUTDIR , WORKDIR)


    DATADIR = os.path.join(os.getenv('VSC_DATA'), name + str(arrayid) + os.getenv('PBS_JOBID'))  
    make_tarfile(DATADIR, WORKDIR)
    print 'We made a copy of the jobdata to: ', DATADIR    

def create_modham_input(nalpha , norb , modtype , options , params , runlist):
    filenum = 0
    for r in runlist:
        params[1] = r
        hub1d = rp.ModHam(nalpha, nalpha,norb,modtype , options , params)
        hub1d.write_file(fname = str(filenum) + modtype + 'nalpha' + str(nalpha) + '_' + 'norb' + str(norb) + '_'  +options[2]+ '_'.join(map(str,params)) +'run=' + str(r)+ '.mod')

        filenum += 1        


def hub1d():
    pokemon = os.environ['VSC_INSTITUTE_CLUSTER']
    flowname = "flow.dat"
    psienergies = []
    methods = ["fci" , "none"]
    nalpha = 4
    norb = 8
    options = ['notzero' , 'pos', 'per'] 
    name = 'hub1dmpdnorb' + str(norb) + 'pair' + str(nalpha) + options[2]  + 'U0-100' #relative to the directory that contains this script
    import numpy as np
    runlist = range(0,101)
    runlist = [0.1 , 1 , 4, 10 , 100]

    modtype = 'Hub1d' ; params = [1. , 4.] 
  
    arrayid = 2
    #hereunder for array jobs for heavy calculations. But make sure it is always defined
    #arrayid =int(os.getenv('PBS_ARRAYID'))
    #arraystep = 1
    #name += str(arrayid)
    #try:
    #    runlist = runlist[(arrayid-1) * arraystep : arrayid *arraystep]
    #except:
    #    print 'encoutered exception in cut of runlist.'
    #    runlist = runlist[-4:-1]

    print 'Start to load necessary modules.'
    os.system('module load cluster/%s' %pokemon)
    os.system('module load PSI/4.0b6-20160201-intel-2016a-mt-Python-2.7.11')
    os.system('module load arpack-ng/3.3.0-intel-2016a')
    os.system('module load OpenBabel/2.3.2-intel-2014b-Python-2.7.8')
    os.system('module load matplotlib')
    os.system('module load h5py')
    #os.system('module load Boost/1.55.0-ictce-7.1.2-Python-2.7.8')
    os.system('export MKL_NUM_THREADS=16')
    #os.system('export KMP_STACKSIZE=32m')

    ORIGDIR= os.getenv('HOME')
    mofile = os.path.join(ORIGDIR, 'devel/CIFlow/lib/mointegrals_%s.so' %pokemon)
    exefile = os.path.join(ORIGDIR, 'devel/CIFlow/bin/ciflow_%s.x' %pokemon )
    #hamdir = os.path.join(ORIGDIR, hamdirname)
    #det_file = os.path.join(ORIGDIR, detfile10 )
    print 'We submitted the script in: ', ORIGDIR
    print 'Psi 4 plugin is at: ' , mofile

    WORKDIR= os.path.join(os.getenv('VSC_SCRATCH_NODE'),os.getenv('PBS_JOBID'))

    if os.path.isdir(WORKDIR):
        shutil.rmtree(WORKDIR)
    os.mkdir(WORKDIR)
    print 'Workdir is: ', WORKDIR

    shutil.copy(mofile,WORKDIR)
    shutil.copy(exefile,WORKDIR)
    #hamdest = os.path.join(WORKDIR, hamdirname)
    #shutil.copytree(hamdir,hamdest)
    #shutil.move(det_file,WORKDIR)

    os.chdir(WORKDIR)
    print 'We work on node: ' , os.getenv('HOSTNAME')
    create_modham_input(nalpha , norb , modtype , options , params , runlist)

    fileinfo = lambda x: float(re.search(r'run[-\w\d=]*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)\.mod' , x).group(1))
    search = modtype+ '.+mod' 
    hamfiles= File_Collector('.', search = search  ,notsearch = r'\.sw\w',sortfunction = fileinfo, filterf =  lambda x : fileinfo(x) >= -1 and fileinfo(x) < 1000. )

    print 'Contents of rundir before run are: ', os.listdir(os.getcwd())
    
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[[], [1]], [] , fname = detfile3 ,add_frozen = 0, ref =  [lambda x , y , z : psir.get_hf_orbs()]) #CIS(P)
    #dw.cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'],[[], [1,2]], [] , fname = detfile4 ,add_frozen = 0, ref =  [lambda x , y , z : psir.get_hf_orbs()]) #CISD(P)
    afh = 'R'
    fname = "ciflowresults.dat"
    with open(fname , 'w') as f:
        header = ht.create_header(afh , methods,psienergies)
        print header
        f.write(header)
        for index , matrixelements in enumerate(hamfiles.plotfiles):
            ht.create_ciflow_input_file(matrixelements , methods , fname = flowname)
            ci_flow =""
            try:
                process = subprocess.Popen(["./ciflow_%s.x" %pokemon ] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                for line in iter(process.stdout.readline, ''):
                    sys.stdout.write(line)
                    ci_flow += line
            except Exception as e:    
                print 'Ciflow gave the following error', e
                print ci_flow
                pass

            """
            try:
                ci_flow = subprocess.check_output(['./ciflow_%s.x' %pokemon ] , stdin =open(flowname , 'r') )
            except Exception as e:    
                print 'Ciflow gave the following error', e
                print ci_flow
                pass
            """
            energies = ht.process_output(ci_flow) 
            ht.print_output(matrixelements, energies , methods)
            f.write("%.15f\t%s\n" %(runlist[index], '\t'.join(energies )) )
  
    os.remove('./ciflow_%s.x' %pokemon)
    os.remove('./mointegrals.so')
    #os.remove(detfile1)
    ht.generate_dir('output_files', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
    ht.generate_dir('unitaries', None, prefix = 'unitary_')
    os.chdir('..')
    ht.generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
    os.chdir('..')

    OUTPUTDIR = os.path.join(ORIGDIR, name+str(arrayid) )
    print 'We Start to create the outputdir ' , OUTPUTDIR
    if os.path.isdir(OUTPUTDIR):
        shutil.rmtree(OUTPUTDIR)

    make_tarfile(OUTPUTDIR , WORKDIR)


    DATADIR = os.path.join(os.getenv('VSC_DATA'), name + str(arrayid) + os.getenv('PBS_JOBID'))  
    make_tarfile(DATADIR, WORKDIR)
    print 'We made a copy of the jobdata to: ', DATADIR    

def con_dm():
    pokemon = os.environ['VSC_INSTITUTE_CLUSTER']
    ciflowoutputfile = "ciflowoutput.txt" ; flowname = "flow.dat"
    methods = ["doci" , "local"]
    methods = ["fci"]
    name = 'noconstrainedm4borhfci'
    import numpy as np
    runlist = list(np.arange(5.,8.01,0.01)) 

    matrixelements = "psioutput3.dat"

    arrayid = 2
    #hereunder for array jobs for heavy calculations. But make sure it is always defined
    #arrayid =int(os.getenv('PBS_ARRAYID'))
    #arraystep = 1
    #name += str(arrayid)
    #try:
    #    runlist = runlist[(arrayid-1) * arraystep : arrayid *arraystep]
    #except:
    #    print 'encoutered exception in cut of runlist.'
    #    runlist = runlist[-4:-1]

    print 'Start to load necessary modules.'
    os.system('module load cluster/%s' %pokemon)
    os.system('module load PSI/4.0b6-20160201-intel-2016a-mt-Python-2.7.11')
    os.system('module load arpack-ng/3.3.0-intel-2016a')
    os.system('module load OpenBabel/2.3.2-intel-2014b-Python-2.7.8')
    os.system('module load matplotlib')
    os.system('module load h5py')
    #os.system('module load Boost/1.55.0-ictce-7.1.2-Python-2.7.8')
    os.system('export MKL_NUM_THREADS=16')
    #os.system('export KMP_STACKSIZE=32m')

    ORIGDIR= os.getenv('HOME')
    mofile = os.path.join(ORIGDIR, 'devel/CIFlow/lib/mointegrals_%s.so' %pokemon)
    exefile = os.path.join(ORIGDIR, 'devel/CIFlow/bin/ciflow_%s.x' %pokemon )
    #hamdir = os.path.join(ORIGDIR, hamdirname)
    #det_file = os.path.join(ORIGDIR, detfile10 )
    matfile = os.path.join(ORIGDIR,  matrixelements)

    print 'We submitted the script in: ', ORIGDIR
    print 'Psi 4 plugin is at: ' , mofile

    WORKDIR= os.path.join(os.getenv('VSC_SCRATCH_NODE'),os.getenv('PBS_JOBID'))

    if os.path.isdir(WORKDIR):
        shutil.rmtree(WORKDIR)
    os.mkdir(WORKDIR)
    print 'Workdir is: ', WORKDIR

    shutil.copy(mofile,WORKDIR)
    shutil.copy(exefile,WORKDIR)
    shutil.copy(matfile , WORKDIR)
    #hamdest = os.path.join(WORKDIR, hamdirname)
    #shutil.copytree(hamdir,hamdest)
    #shutil.move(det_file,WORKDIR)

    os.chdir(WORKDIR)
    print 'We work on nore: ' , os.getenv('HOSTNAME')

    outputfile = open(ciflowoutputfile , 'w')

    ht.generate_dir('output_files', None)
    os.chdir(WORKDIR)
    matrixelements = os.path.basename(matrixelements)
  
    afh = 'R'
    fname = name + "_" + ".dat"
    with open(fname , 'w') as f:
        header = ht.create_header(afh , methods, [] , extra = None).replace('\n', '') + "\tlambda" +"\tMulliken_A\n"
        print header
        f.write(header)
        for r in runlist:
            norb  = 10 ;  
            nalpha = 7;

            def func(x):
                modtype = 'Constrained_DM' ; params = [0. , 1. , 2. , 3. , 4. ,  r , x]  ; options = [ "lj"]
                condm = rp.ModHam(nalpha, nalpha,norb,modtype , options , params, matrixelements =matrixelements)
                condm.write_file(fname = 'conelements.mod')
                ht.create_ciflow_input_file('conelements.mod', methods , fname = flowname)
                ci_flow =""
                if methods:
                    try:
                        ci_flow =""
                        process = subprocess.Popen(["./ciflow_%s.x" %pokemon] , stdin =open(flowname , 'r'), stdout = subprocess.PIPE)
                        for line in iter(process.stdout.readline, ''):
                            sys.stdout.write(line) #writes intermediate output to screen.
                            outputfile.write(line)
                            outputfile.flush()
                            ci_flow += line
                    except Exception as e:    
                        print 'Ciflow gave the following error', e
                        print ci_flow
                        pass

                print ci_flow
                energies = ht.process_output(ci_flow) 
                mullikencharge = float(re.search("Mulliken[\s\w]+:\s*([\-+]?\d+[\.,]?\d*[eEDd]?[\-+]?\d*)", ci_flow).group(1) )
                return (energies , mullikencharge)

            if r <  7 :
                extremumval = gss(func, -0.3 , 2.5,tol=1e-8)
            elif r < 8.1 :
                extremumval = gss(func,-4 ,0.1,tol=1e-8)
            energies, mullikencharge = func(extremumval)
            shutil.copy( 'conelementsoutputfci.dat' , os.path.join('output_files' , 'conelementsoutputfci' + str(r)+ '.dat') )

            ht.print_output('searchfor:'+str(r), energies , methods)
            f.write("%.15f\t%s\t%f\t%f\n" %(r, '\t'.join(energies ), extremumval , mullikencharge) )
  
    ht.generate_dir('output_files2', None, prefix = r'output[\w_]+.dat')
    os.chdir('..')
 
    os.remove('./ciflow_%s.x' %pokemon)
    os.remove('./mointegrals.so')
    #os.remove(detfile1)
    ht.generate_dir('unitaries', None, prefix = 'unitary_')
    os.chdir('..')
    ht.generate_dir('matrixelements_otherbasis', None, prefix = 'hamp')
    os.chdir('..')
    OUTPUTDIR = os.path.join(ORIGDIR, name+str(arrayid) )
    print 'We Start to create the outputdir ' , OUTPUTDIR
    if os.path.isdir(OUTPUTDIR):
        shutil.rmtree(OUTPUTDIR)

    make_tarfile(OUTPUTDIR , WORKDIR)

    DATADIR = os.path.join(os.getenv('VSC_DATA'), name + str(arrayid) + os.getenv('PBS_JOBID'))  
    make_tarfile(DATADIR, WORKDIR)
    print 'We made a copy of the jobdata to: ', DATADIR    


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
    fc=f(c)[0][-1];fd=f(d)[0][-1]
    while abs(c-d)>tol:       
        if fc<fd:
            b=d
            d=c  
            c=b-gr*(b-a)
            fd=fc;fc=f(c)[0][-1]
        else:
            a=c
            c=d  
            d=a+gr*(b-a)
            fc=fd;fd=f(d)[0][-1]
    return (b+a)/2

if __name__ == "__main__":
    main()
    #hub1d()
    #con_dm()
