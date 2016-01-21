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

from itertools import combinations
import sys
import re
import read_psi as rp
import ciflowoutput as co

def perm_onestring(n,k):
    result = []
    for bits in combinations(range(n), k):
        s = ['0'] * n
        for bit in bits:
            s[bit] = '1'
        result.append(''.join(s))
    return result

def add_excitation( bs, num =1):
    bs = list(bs)
    occ = [pos for pos,i in enumerate(bs) if i == '1']
    emp = [pos for pos,i in enumerate(bs) if i == '0']
    excitations = []
    occcomp = combinations(occ,num)
    for ocsel in occcomp:
        bsc = list(bs)
        for i in ocsel: 
            bsc[i] = '0'
        empcomp = combinations(emp,num)
        for empsel in empcomp:
            for j in empsel:
                bsc[j] = '1'
            excitations.append(''.join(bsc))
            for j in empsel:
                bsc[j] = '0'
    return excitations

def excite_determinant(det, num =1):
    detlist = []
    for numex in range(num+1):
        numa = numex
        numb = num -numex
        exalist = add_excitation(det[0], num =numa)
        exblist = add_excitation(det[1], num =numb)
        detlist += [(exa , exb) for exa in exalist for exb in exblist]
    return detlist

def pair_excitation(detp , num = 1):
    exalist = add_excitation(detp[0], num = num) #We assume even number of electrons, otherwise isolate the odd electron and add it at the end.
    return zip(exalist , exalist)


def get_hf_det(numup,numdown , numorb):
    """
    returns the (0..01..1 , 0..01..1) #nup times 1 and ndown times 1 , to reverse -> add [::-1]
    REMARK: THIS is mostly not the hf determinant, because psi mo integrals are labelled first by irrep and then by energy. Use the get_hf_orbitals form the psireader class instead.
    """
    return ''.join(['0']*(numorb-numup)+['1'] * numup) , ''.join(['0']*(numorb-numdown)+['1'] * numdown  )

def write_to_file(bitstrings , fname = 'determinants.dat', header = '#determinants for the CIMethod\n'):
    #[(00011,00011),....,....]
    text = header + '\n'.join(['|'.join(bit) for bit in bitstrings])
    with open(fname , 'w') as file:
        file.write(text)

def write_games(detlist, fname = 'determinantscisd.dat'):
    import numpy as np
    newlist = [' '.join(map(str,np.array(map(int,list(up[::-1]))) + np.array(map(int,list(down[::-1])) ) ) ) for up,down in detlist]
    newlist = list(set(newlist)) # eliminate dubbles
    numdet = len(newlist)
    text = '%d\n' %numdet + '\n'.join(newlist)
    with open(fname , 'w') as file:
        file.write(text)

def ci_excitation(nup , ndown , norb , exlist , get_reference = get_hf_det, pairex = False):
    """
    Needs a list of the excitations (exlist), the filename where the determinants are going to be written to (fname), the number of up electrons(nup) and the number of down electrons (ndown), the number of orbitals (norb).
    """
    detlist = [ get_reference(nup,ndown,norb)]
    for ex in exlist:
        if pairex:
            detlist += pair_excitation(detlist[0] , num = ex) #only pair excitations cis_p = intersection DOCI and CISD
        else:
            detlist += excite_determinant(detlist[0], num = ex) #normal cis , cisd, cisdt , ... to a reference determinant

    return detlist


def doci(nup , norb):
    """
    Use this only when nup equals ndown, if not use the more general seniority_enhancer
    """
    uploop = perm_onestring(norb,nup)
    downloop = list(uploop)
    return zip(uploop,downloop)

def fci(nup,ndown , norb):
    uploop = perm_onestring(norb,nup)
    downloop = perm_onestring(norb,ndown)
    return [(up,down) for up in uploop for down in downloop]

def seniority_enhancer(nup, ndown , norb , senlist = [0]):
    """
    Generates all the determinants with seniority in seniority list (for example doci is all determinants with seniority equal to zero)
    TODO: test more (seems correct for now)
    """
    sz2 = nup-ndown
    if norb > nup+ndown:
        maxsen = nup + ndown
    else:
        maxsen = 2*norb - nup - ndown
    detlist = []
    for sen in senlist:
        assert(abs(sz2)<= sen ) , 'with nup: %d, ndown: %d, seniority: %d is not possible' %(nup , ndown , sen)
        assert(sen <= maxsen ) , 'with nup: %d, ndown: %d , norb:%d, seniority: %d is not possible' %(nup , ndown ,norb, sen)
        neworb = norb-sen 
        if sz2 >= 0:
            npar = ndown - (sen -sz2)/2
        else:
            npar = nup - (sen - abs(sz2))/2
        newlist = doci(npar, neworb)
        combup = combinations(range(norb),nup-npar)
        for comup in combup:
            downpos = [i for i in range(norb) if i not in comup]
            combdown = combinations(downpos,ndown-npar)
            for comdown in combdown:
                totcomb = [(i, 'u') for i in comup]
                totcomb += [(i, 'd') for i in comdown]
                totcomb = sorted(totcomb, key = lambda x : x[0])
                #print totcomb
                newdet = []
                for det in newlist:
                    up = list(det[0])
                    down = list(det[1])
                    #print up , down
                    for socc in totcomb:
                        #adding the singly occupied electrons
                        if socc[1] == 'u':
                            up.insert(socc[0],'1')
                            down.insert(socc[0] , '0')
                        if socc[1] == 'd':
                            up.insert(socc[0] , '0')
                            down.insert(socc[0] , '1')
                    newdet.append(( ''.join(up), ''.join(down) ))
                #print newdet

                if newdet != []:
                    detlist += newdet
                else:
                    detlist += newlist

    return detlist


"""
Returns a list of reference determinants which spread over the entire norb, in pairs example 5 up, 5down, 10 orbs -> reference list contains the determinant with 5 ups and 5 downs in the lowest 5 spatial single particle levels, and 5 ups and 5 downs in the highest single particle levels.
"""
def get_reference_list(nup , ndown , norb):
    reflist = []
    nbeginu = 0
    nbegind = 0
    nstopu = nup
    nstopd = ndown
    while((nbeginu + nup) < norb and (nbegind + ndown) < norb):
        reflist.append( ( ''.join(['0']*(norb-nstopu)+['1'] * nup+ ['0'] * nbegind) , ''.join(['0']*(norb-nstopd)+['1'] * ndown+ ['0'] * nbeginu ) ) )
        nbeginu += nup
        nbegind += nup
        nstopu += nup
        nstopd += ndown
        if (nbeginu + nup) > norb or (nbegind + ndown) > norb:
            reflist.append( (''.join(['1']*nup +['0']*(norb- nup)) , ''.join(['1'] * ndown+['0']*(norb-ndown)) )  )
    return reflist


def addactvirt(det, active , frozen , virtual):
    det = list(det)
    det = det[::-1]
    pos = 0
    for i in range(len(active) ):
       det[pos:pos] = ['1'] * frozen[i] 
       det[pos +frozen[i]+ active[i] : pos + frozen[i] + active[i] ] = ['0'] * virtual[i]
       pos += frozen[i] + active[i] + virtual[i]
      
    det = det [::-1]
    det = ''.join(det)
    return det 

def insert_frozen_virtual(detlist, active , frozen , virtual):
    detlist = [ (addactvirt(deta,active , frozen , virtual), addactvirt(detb, active , frozen , virtual)    )   for deta , detb in detlist]
    return detlist



def cimain(nup , ndown , norb , exlist, senlist , fname = 'determinants.dat' , exheader = '', ref = [get_hf_det], active = None , virtual = None , frozen = None, add_frozen = 0, add_virtual = 0, frozenstring = None, virtualstring = None, write = True):
    """
    set exlist = [] if you don't want any excitations of the reference determinant, [0] if you only want reference determinant
    set senlist = [] if you don't want any seniority based creation of determinants
    add_frozen is the number of frozen electron pairs in the lowest molecular orbitals.
    """
    detlist = []
    if len(exlist[0]) != 0 or len(exlist[1]) != 0:
        try:
            for reference in ref: 
                if not hasattr(reference ,'__call__'):
                    a = reference[0]
                    b = reference[1]
                    reference = lambda x,y,z : (a, b)

                detlist += ci_excitation(nup , ndown , norb , exlist[0], get_reference = reference , pairex = False)
                detlist += ci_excitation(nup , ndown , norb , exlist[1], get_reference =  reference, pairex = True)

        except TypeError as e:
            print e , 'wrap the reference in a list.'
            sys.exit(1)

    if len(senlist) != 0:
         detlist += seniority_enhancer(nup, ndown , norb , senlist = senlist)

    detlist = list(set(detlist)) # remove double determinants

    if active != None:
        detlist = insert_frozen_virtual(detlist, active , frozen , virtual)
        
    
    if add_frozen != 0:
        detlist = [(detup + '1'* add_frozen, detdown + '1'* add_frozen) for detup , detdown in detlist]
        norb += add_frozen
        nup += add_frozen
        ndown += add_frozen

    if add_virtual != 0:
        detlist = [('0'* add_virtual +detup , '0'* add_virtual +detdown ) for detup , detdown in detlist]
        norb += add_virtual
    
    if frozenstring != None:
        detlist = [(detup + frozenstring, detdown + frozenstring) for detup , detdown in detlist]
    if virtualstring != None:
        detlist = [(virtualstring + detup , virtualstring + detdown ) for detup , detdown in detlist]

    head = exheader+'#determinants for the CIMethod created by exciting a reference determinant created by: %s\n' % ( str(ref) )+'#nup = %d\n#ndown = %d\n#norb = %d\n#excitations = %s' %(nup,ndown,norb, str(exlist)) + '\n#all determinants with following seniorities are included = %s\n#we added %d frozen electron pairs\n' %( str(senlist), add_frozen ) 

    if write:    
        write_to_file(detlist , fname = fname, header = head)
        #write_games(testex)
    return detlist

def main(*args , **kwargs):
    exlist = [[1,2,3,4,5,6] ,[]]; senlist = [] 
    fname = 'psioutput.dat'; detfile1 = 'determinants.dat' ; detfile2 = 'determinants2.dat'

    psir = rp.PsiReader(fname, isbig = False, numorbs = -1 , read_ints = False)
    #cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[[1,2,3,4,5,6], []  ], [],fname = detfile1 ,ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #CISD
    #cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [[1,2], [2] ], [ ] , fname = detfile2 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #CISDDOCI
    #cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [[], [] ], [ 0] , fname = detfile1 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #DOCI
    #cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [[1,2], [] ], [ ] , fname = detfile2 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #CISD
    #cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [[], [] ], [ 0] , fname = detfile1 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #FCI
    cimain(psir.values['nalpha'],psir.values['nbeta'], psir.values['norb'], [range(1, psir.values['nalpha'] + psir.values['nbeta'] ) , [] ], [ ] , fname = detfile1 ,ref =  [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0) #FCI
    #detlist = doci(psir.values['nalpha'],psir.values['norb'])
    #write_to_file(detlist , fname = 'determinantscibig.dat', header = '#determinants for FCI\n')

def num_dets():
    fname = 'psi0_output10.dat' 
    #fname = 'beh20.86ccpvdz.out' 
    psir = rp.PsiReader(fname, isbig = False, numorbs = -1 , read_ints = False)
    detmethods = []
    #detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[range(1, psir.values['nbeta']+psir.values['nalpha']+1), []  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    from scipy import special as sc
    detmethods.append(sc.binom(psir.values['norb'] , psir.values['nalpha']) * sc.binom(psir.values['norb'] , psir.values['nbeta']))
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [] , [1]  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [] , [1,2]  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [] , [1,2,3]  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [] , []  ], [0],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [1,2] , []  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [1,2] , [1,2]  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [1,2] , [1,2,3]  ], [],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods.append( len(cimain(psir.values['nalpha'],psir.values['nbeta'] ,psir.values['norb'],[ [1,2] , []  ], [0],fname = '',ref = [lambda x , y , z : psir.get_hf_orbs()] , add_frozen = 0 , write = False) ) )
    detmethods = [ (value , value / float(detmethods[0]) * 100.)    for value in detmethods]
    with open("numdetsn2.dat" , "w") as file:
        for value in detmethods:
            file.write("%f    %.9f\n" % (value[0] , value[1]))

def generate_all_ex(nup ,ndown , norb , refdet, aname = 'determinants'):
    exlist = [[],[] ] #single ex, pair ex
    maxex = 2*norb-nup-ndown
    maxex = max(maxex , nup+ndown)
    for ex_max in range(1, maxex+1):
        exlist[0] = range(1,ex_max+1) #
        print exlist
        name = aname + str(ex_max) +".dat"
        cimain(nup , ndown , norb ,exlist,  [], fname = name , ref = [lambda x , y , z : refdet ])

def generate_all_sen(nup ,ndown , norb , aname = 'sendeterminants'):
    senlist = [] #single ex, pair ex
    startsen = abs(nup-ndown)
    numpair = min(nup, ndown)
    maxpairbreaking = min(numpair , norb - max(nup , ndown) )
    maxextrasen = 2*maxpairbreaking
    totmaxsen = 2*maxpairbreaking + startsen
    for sen in range(startsen, totmaxsen+1,2):
        senlist = range(startsen,sen +1,2) #
        print senlist
        name = aname + str(sen) +".dat"
        cimain(nup , ndown , norb , [[],[]],  senlist, fname = name , ref = [lambda x , y , z : get_hf_det() ])

def biggest_det_ex(wffile = 'psioutputoutputfci.dat'):
    cifread = co.CIFlow_Reader(wffile)
    maxdet = cifread.get_max_det()
    print 'maxdet = ' , maxdet
    #generate_all_ex(cifread.header['nup'], cifread.header['ndown'], cifread.header['norbs'], get_hf_det(cifread.header['nup'],cifread.header['ndown'],cifread.header['norbs']) )
    generate_all_ex(cifread.header['nup'], cifread.header['ndown'], cifread.header['norbs'], tuple(maxdet[0].split('|')) )

def test_main():
    fname = 'psi0_output10.dat' ; detfile1 = 'determinants.dat' ; detfile2 = 'determinants2.dat'
    psir = rp.PsiReader(fname, isbig = False, numorbs = -1 , read_ints = False)
    #cimain(psir.values['nalpha']-1,psir.values['nbeta']-1 ,psir.values['norb']-1,[[1,2,3,4], [] ], [],fname = detfile1 ,ref = [lambda x , y , z : get_hf_det(x,y,z)] , add_frozen = 1) #CISD
    frozen = [1,0,0,0,0,1,0,0]  ; virtual = [1,1,1,1,1,2,0,0]; active = [5,0,2,2,0,4,3,3]
    assert( sum(frozen) + sum(virtual) + sum(active) == psir.values['norb'])
    cimain(psir.values['nalpha'] - sum(frozen),psir.values['nbeta'] -sum(frozen)  ,psir.values['norb'] - sum(frozen) - sum(virtual) , [[1,2], []], [],fname = detfile2 ,ref = [lambda x , y , z : psir.get_hf_orbs(  frozen = frozen , virtual = virtual )] , active = active, frozen =  frozen, virtual = virtual ,add_frozen = 0, add_virtual = 0) #CISD

if __name__ == "__main__":
    #test_main()
    #reflist = get_reference_list(3,3,10)
    #main()
    #num_dets()
    #print len(fci( 6, 6, 12))
    #cimain(7, 7, 10, [ [1,2,3,4] , [] ], [ ] , fname = "cisddeterminants.dat" ,ref =  [lambda x , y , z :  get_hf_det(7,7,10)] , add_frozen = 0) #FCI
    #biggest_det_ex()
    generate_all_sen(3,5,6)
    generate_all_sen(3,5,10, 'sendetextra')
