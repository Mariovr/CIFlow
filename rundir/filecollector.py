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
# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2014;
# (c) Ghent University, 2014
# -*- coding: utf-8 -*-
#! /usr/bin/env python 

import os
import re

class File_Collector(object):
    def __init__(self, rootdir , search , notsearch = '.png' , notdir = 'xyvwa' , filelist = None , sortfunction = None , rev = False, filterf = None):
        if filelist != None:
            self.plotfiles = filelist
        else:
            self.plotfiles = []
        self.readfiles(rootdir , search , notsearch = notsearch , notdir = notdir)
        self.sortplotfiles(rev, sortfunction)
        if filterf != None:
            self.plotfiles = filter(filterf , self.plotfiles)
        print self.plotfiles
  
    def addfiles(self , *args):
        for i in args:
            self.plotfiles.append(i)
    
    def sortplotfiles(self, rev = False, sortfunction = None):
        if sortfunction != None:
            self.plotfiles = sorted(self.plotfiles ,key = sortfunction , reverse = rev )
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
                print m
                if m is None:
                    self.readfiles(filep , search, notsearch = notsearch, notdir = notdir )
            elif os.path.isfile(filep) and '.swp' not in filep: 
                nm = re.search(notsearch, filep)
                m = re.search(search , filep)
                if m is not None and nm is None:
                    self.plotfiles.append(filep)
            else:
                pass
    
def main(*args , **kwargs):
    pass

if __name__ == "__main__":
    main()

