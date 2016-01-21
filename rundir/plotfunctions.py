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
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2014-2015;
# (c) Ghent University, 2014
#!/usr/bin/env python
import numpy as np
import pylab as pl
import os,sys,shutil
import re
import matplotlib
import math

import datareader as dr
import filecollector as fc

## for Palatino and other serif fonts use:
#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']}) #adjust fonts
#matplotlib.rc('font',**{'family':'serif','serif':['Palatino']})
#matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#matplotlib.rc('text', usetex=True)

def makemovie(name = None):
  # makes a movie from all the .png files in the current directory
  print 'Starting to create a movie, with all the .png files in directory: %s ' %str(os.getcwd())
  if name != None:
    dirname = name
  else:
    dirname = str(os.getcwd())
  command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           dirname+'.avi')

  os.spawnvp(os.P_WAIT, 'mencoder', command)

class Plot_Files(object):
    def __init__(self, name , searchstr = 'plotenergy' , notsearch =r'\.swp|\.png' , regexp = None, sortfunction = None):
        """
        initialisation of the plotting class, we create a figure to which we add one axes
        """
        
        self.data = [] ; self.prefix = []

        if isinstance(name , (list , tuple)):
            self.readdata(name, regexp = regexp ,  substr = None)
        elif os.path.isdir(name):
            self.procesfiles(name , searchstr, notsearch = notsearch , regexp = regexp , sortfunction = sortfunction)
        elif os.path.isfile(name):
            self.readdata([name], regexp = regexp)

        self.fig = pl.figure(1,figsize=(18,14), dpi=100 , frameon = True , facecolor = '0.75' , edgecolor = 'w')
        self.fig.add_subplot('111' , axisbg = 'w' , projection = 'rectilinear') #if you want to add axes on particular place: fig.add_axes([0.15, 0.1, 0.7, 0.3]) where ->  [begin , bottom to start axes , width , height ]
        self.separated = False #if we have a list and need to plot the plots separated
        self.axnum = 0
  
    def add_axes(self,pos = [0.5 , 0.2 , 0.4 , 0.3], axisbg = None , projection = 'rectilinear'):
        self.fig.add_axes(pos , axisbg = axisbg, projection = projection)
  
    def readdata(self, reflist , comment = '#' , regexp = None , substr = None, filename = True, clear_data = True):
        """
        Creates a Read object (see datareader.py) which handles the header and data reading.
        """
        if clear_data:
            self.clear_data()
        if os.path.isfile(str(reflist)):
            reflist = [reflist] #if we work with only one file this wraps it automatically in right format
        for ref in reflist:
            print('start with the collection of data from file %s' %ref)
            if not filename:
              self.prefix.append( os.path.dirname(ref) + '/')
            else:
              self.prefix.append(re.sub('\.dat$' , '' , ref))
            self.data.append(dr.Reader(ref))
  
    def clear_data(self):
        self.data = []

    def adjust_data(self, func):
        for i in range(len(self.data)):
            for j in range(1,len(self.data[i].data[0,:])):
                self.data[i].data[:,j] = func(self.data[i].data[:,j])
            
    def procesfiles(self, dirname , search , notsearch = r'\.sw*|\.png', notdir = 'awfwfr', sortfunction = None , rev = False , regexp = None , substr = None , filelist = None , filename = True):
      filecol = fc.File_Collector(dirname , search , notsearch = notsearch ,filelist = filelist , sortfunction = sortfunction , rev =rev )
      self.readdata(filecol.plotfiles, regexp = regexp ,  substr = substr, filename = filename)
  
    def generate_plot(self, xlimg = None , ylimg =None , exname = '' , prefix = True, titel = None ,name = 'energyplot', depcol = 0, ylist = None , save = True, color = ['','','','','','','','','', '','', '' , '' , '' , ''], exdir = './', samedir = False, datanum = 0, label = None):
      """
      Some nice plots to visualize the data with matplotlib.
      """
      print ('start with the generation of plots')
      self.separated = False
      #if ylist is defined it just plots the columns with those indices. otherwise it plots all the columns
      if not ylist:
        ylist = range(0,self.data[0].data.shape[1]) 
        ylist.remove(depcol)
        

      for i in ylist:
          if label == None:
              label = self.data[datanum].columnheader[i]
          self.plotwrap(depcol, i, 'energy' , name =  exname, xlim = xlimg , ylim = ylimg , prefix = prefix ,save = False, label = label , color = color[i-1], datanum = datanum)
      self.layout(self.data[datanum].get_xlabel() , self.data[datanum].get_ylabel(), tit = titel, xlim = xlimg , ylim = ylimg )
      if save:
          self.savefig(name+exname,  prefix = prefix, exdir = exdir, samedir =  samedir)  
  
    def plotwrap(self, xindex, yindex, yas, name = None, titel = None ,color = 'r' , sort = '' , label = None , xlim = None , ylim = None , prefix = False,save = True, datanum = 0):
      self.fig.axes[self.axnum].plot(self.data[datanum].data[:len(self.data[datanum].data[:,yindex]),xindex],self.data[datanum].data[:,yindex], color+sort , label = label)
      #self.fig.axes[self.axnum].scatter(self.data[datanum].data[:len(self.data[datanum].data[:-3,yindex]),xindex],self.data[datanum].data[:-3,yindex])
      if self.separated == True and save:
        self.layout(self.data[datanum].get_xlabel(), self.data[datanum].get_ylabel(), tit = titel, xlim = xlim , ylim = ylim)
        self.savefig(name, filenum = i , prefix = prefix)
      if self.separated == False and save:
        self.layout(self.data[datanum].get_xlabel() , self.data[datanum].get_ylabel(), tit = titel, xlim = xlim , ylim = ylim)
        self.savefig(name + 'together' , prefix = prefix)
  
    def scatterplot(self ,  xvars , yvars , colorvars , colormap = None ):
        cm = pl.cm.get_cmap(colormap)
        sc = self.fig.axes[self.axnum].scatter(xvars ,yvars, c=colorvars, cmap = cm )
        pl.colorbar(sc)
  
    def normalize_to_groundstate(self):
      print('Warning we normalize all the excited states to the groundstate energy')
      gronddat = self.data[0].data
      for i in range(1,len(self.data)):
        dif = np.shape(gronddat )[0] - np.shape(self.data[i])[0]
        print dif
        if dif < 0 :
          self.data[i].data = self.data[i].data[0:dif ,:] 
        elif dif > 0:
          gronddat = gronddat[: -1.*dif , :]
        print np.shape(gronddat) , np.shape(self.data[i].data)
        self.data[i].data[:,1:3] = self.data[i].data[:,1:3] - gronddat[:,1:3] #we made sure that the data of the groundstateenergy is first in the rgdata list
      del(self.data[0].data, self.prefix[0])
      
    def layout(self ,  xlab , ylab , xlim = None , ylim = None , tit = None , legendhand = None , legendlab = None , legendpos = 0 , finetuning = False , axbg = None , fs = 22, ticksize = 14):
      """
      In this function we finetune some aspects of the axes for all the tuning possibilitys see: http://matplotlib.org/api/axes_api.html
      especially the set functions ;)
      """
      print('We are starting with the layout')
      self.fig.axes[self.axnum].set_xlabel(xlab, fontsize = fs)
      self.fig.axes[self.axnum].set_ylabel(ylab , fontsize = fs)
      if xlim != None:
        self.fig.axes[self.axnum].set_xlim(xlim) #good value for xlim in the case of a xi path is : (2*self.rgeq.energiel[0]-5*(self.rgeq.energiel[1]-self.rgeq.energiel[0]),2*self.rgeq.energiel[-1]+0.5) 
      if ylim != None:
        self.fig.axes[self.axnum].set_ylim(ylim)
      if tit != None:
        self.fig.axes[self.axnum].set_title(tit , fontsize = fs)
      if legendlab != None:
        self.fig.axes[self.axnum].legend(legendhand , legendlab, loc = legendpos)  #if you want to add extra info
      
      #self.fig.axes[self.axnum].ticklabel_format(style='sci', axis='y') #force scientifique notation for y axis
      #self.fig.axes[self.axnum].yaxis.major.formatter.set_powerlimits((0,0))
      self.fig.axes[self.axnum].yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter(useOffset=False))
      #self.fig.axes[self.axnum].yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
      self.fig.axes[self.axnum].xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
      #self.fig.axes[self.axnum].locator_params(tight=True, nbins=3)
      #self.fig.axes[self.axnum].xaxis.major.locator.set_params(nbins=3)


      if self.axnum > 3:
          self.fig.axes[self.axnum].yaxis.major.locator.set_params(nbins=3) 

      for tick in self.fig.axes[self.axnum].xaxis.get_major_ticks():
        tick.label.set_fontsize(ticksize) 
      for tick in self.fig.axes[self.axnum].yaxis.get_major_ticks():
        tick.label.set_fontsize(ticksize) 
    
    
      self.fig.subplots_adjust(hspace=.5)
  
      leg = self.fig.axes[self.axnum].legend(loc = legendpos) #draws the legend on axes[axnum] all the plots that you labeled are now depicted in legend
  
      if axbg != None:
        self.fig.axes[self.axnum].set_axis_bgcolor(axbg)
      """
      if you forgot to add a label to a line with linenumber: lnum you can do: self.fig.axes[self.axnum].lines[lnum].set_label('this is my new label')
      the underneath is the same as : h , l = self.fig.axes[self.axnum].get_legend_handles_labels()
                                      self.fig.axes[self.axnum].legend(h,l)
      """
      if finetuning == True:
        # the matplotlib.patches.Rectangle instance surrounding the legend
        frame  = leg.get_frame()  
        frame.set_facecolor('0.80')    # set the frame face color to light gray
  
        # matplotlib.text.Text instances you can change all properties of labels
        for t in leg.get_texts():
          t.set_fontsize('small')    # the legend text fontsize
  
        # matplotlib.lines.Line2D instances
        for l in leg.get_lines():
          l.set_linewidth(1.5)  # the legend line width
  
    def savefig(self , name , filenum = 0 , samedir = False , prefix = True, exdir = './', typ = 'pdf'):
      """
      After we are satisfied with our figure we save it with this function: dpi = pixels per inch, under a name determined by the savestring function().
      """
      #REMARK watch out with the translation of the dot to nothing when you gave as arguments the current working directory '.' because
      #if you do this it is not possible to save the file in the appropriate place because the folder doesn't exist anymore 
      #because the first . dissapeared you can only remove . from floats or extensions not from current dir (maybe build in check that if the first letter of the filename is a dot then that dot is not removed)
      figname = self.savestring(name , filenum , samedir = samedir , prefix = prefix , exdir = exdir, typ = typ)
      self.fig.savefig(figname , dpi = 300 , facecolor = 'w' , edgecolor = 'w')
      self.fig.clf()
      #self.axnum +=1
      #self.fig.add_subplot('32%i' % (self.axnum+1) , axisbg = 'w' , projection = 'rectilinear')
      self.fig.add_subplot('111' , axisbg = 'w' , projection = 'rectilinear')
  
    def savestring(self , name , filenum , samedir = False , prefix = True, exdir = './', filenumwrite = False, typ = 'png'):
      """
      This function generates the name whereunder the figure is going to be saved
      """
      print exdir
      if filenumwrite:
          name += str(filenum)
      if prefix == True:
        if samedir:
          """
          Making use of some implementation detail of savefig, if we read in files from all different directory's, the prefixes contain the path of those files relative to the rootdirectory. So if you save the file we save it with first the prefix and then the name , so the figures end up in the same directory as the files. If you don't want this behaviour we need to remove the / in the prefixs so fig.savefig will not recognize it as a path so all the figures end up in the current working directory. Remark we only remove the / because if all the figures end up in same dir we need the path information to distinguish them.
          """
          self.prefix = [pre.translate(None , '/.')  for pre  in self.prefix]
        return '%s%s%s.%s' %(exdir , self.prefix[filenum], name , typ)

      else:
        return '%s%s.%s' %(exdir, name, typ)
  
    def writetext(self ,text , pos , axnum = 0, hor = None  ,ver = None , rot = None ,fs =14 , transform = None):
      self.fig.axes[self.axnum].text(pos[0] ,pos[1] ,text , rotation = rot ,horizontalalignment = hor, verticalalignment = ver , fontsize = fs, transform = transform) #, color = 'black', style = 'italic')
  
    def least_sqr_fit(self,x, y):
      """
      Calculates the least square fit of a list of independend variables x and dependend variables y.
      It returns a list of function values of the best fitted straight line, with the given x values as independend variables and also a list with the parameters
      that define the line. It's also possible to fit at the same time multiple datasets with the same xvalues just give y the form [(v1 , v2 , v3) , (v1 , v2 , v3), ... ]
      Where the first tuple consists of the function values of x1 the second of x2 .... , So you get immediately three fitted lines, with the coefficients in a[0][0] , a[0][1]
      , a[0][2] for the first, second and third rico for the three lines same for the bisection point with y axis
      """
      A = np.array([ x, np.ones(len(x))])
      # linearly generated sequence
      a,f,g,h = np.linalg.lstsq(A.T,y) # obtaining the parameters
      print 'de gevonden rechte = %.10f x + %.10f' %(a[0], a[1])
      lined = map(lambda g: a[0]*g +a[1],x) # regression line
      return lined , a

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
    
    def convert(self , col = 0 , to = "au", which = 0):
        if to == "au":
            self.data[which].data[:,col] = [self.convert_au_to_ang(x)  for x in self.data[which].data[:,col]]
        else:
            self.data[which].data[:,col] = [self.convert_ang_to_au(x)  for x in self.data[which].data[:,col]]


def main():
    angstrom = False 
    fname = "results/beh2_sto_3g_symcompc1/beh2fcidocicompc1_sto-3g.dat"
    fname = "results/beh2_sto_3g_symcompc1/overlapanalysisextra6.dat"
    fname = "results/beh2_both/cis_pairsref_sto-3g.dat"
    fname = "results/beh2cisdocic1/ciflowresults_sto-3g.dat"
    fname = './results/c2sym_linhybsdd/c2sym_linhybsdd_sto-3g.dat'
    fname = './results/beh2_testlscisdmminddoci/beh2cisdmmindocitestno_sto-3g.dat'
    fname = './results/beh2_testlscisdmminddoci/beh2cisdmmindocitestno_sto-3g.dat'
    fname = 'beh2linhybresults.dat'
    fname = 'results/beh2linhyb6-31gdata/beh2linhyb6-31gc1act/ciflowresults_6-31g.dat'
    fname = 'results/cisdc1mminbeh2/overlapanalysisdocicisdmmin.dat'
    fname = 'results/h26-31gdoci/h26-31gdoci_6-31g.dat'
    fname = 'results/h2cc-ptzgdoci/h2cc-pqzdocinew_cc-pvtz.dat'
    fname = 'beh2testdatamoihiallopen.dat'
    fname = 'results/beh2localtest/beh2localtest_sto-3g.dat'
    fname = 'results/noconstraineddmhermextrema/noconstrainedm_.dat'
    #fname = 'results/noplusconstraineddmhermextrema100angstrom/noconstrainedm_.dat'
    fname = 'results/noconstraineddmhermextrema100angstromgoodtrans/noconstrainedm_.dat'
    #fname = 'results/noconstraineddmhermextrema0.1angstrom/noconstrainedm_.dat'
    #fname = 'results/noconstraineddmhermextremainfinitedistance/noconstrainedm_.dat'
    fname = 'results/noconstraineddmhermextremainfiniteangstromscan7.98/noconstrainedm_.dat'
    
    plotter = Plot_Files(fname)
    #title = 'NO$^+$ Infinite Angstrom scan different ham (STO-3G)'
    title = 'NO$^+$ Infinite Angstrom scan 7.98 (STO-3G)'
    xlim = None
    ylim = None
    #xlim = (0.5,3.5)
    #ylim = ( -15.9 , -14.5) 

    if angstrom:
        plotter.reader.convert(0) #converts angstrom to au.
        plotter.reader.units['x'] = '(a.u.)'

    plotter.data[0].depvar['xas'] = 'lambda ' #change the future y-axis label 
    plotter.data[0].depvar['yas'] = 'FCI energy' #change the future y-axis label 
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'lambdascanfci')

    #plotter.data[0].depvar['xas'] = 'Mulliken charge on N'#change the future y-axis label 
    #plotter.data[0].depvar['yas'] = 'Extremum lambda'#change the future y-axis label 
    #plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [3], titel = title, name = 'plot' , exname = 'mullikenlambda')
    #plotter.data[0].depvar['yas'] = 'Energy'#change the future y-axis label 
    #plotter.data[0].units['y'] = r'(E$_h$)'
    #plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'mullikenenergy')
    #plotter.data[0].depvar['xas'] = 'Extremum lambda'#change the future y-axis label 
    #plotter.data[0].units['x'] = r'(a.u.)'
    #plotter.generate_plot(depcol = 3 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'lambdaenergy')
    option = None
    #here under is an example to show the flexibility of the plotclass
    if option == 'inset':
        """
        Example of how you need to draw a small inset in a larger plot of a particular
        area of interest with matplotlib
        """
        plotter.reader.depvar['xas'] = r'$\eta$ '  #change the future x-axis label to latex 
        begin =0
        stop = None
        plotter.plotrgvars(begin = begin , stop = stop , name = 'etanul2', save = False)
        begin =  9880
        stop = None
        plotter.rgindex = 5
        plotter.reader.npair = 9
        plotter.add_axes([0.5,0.2,0.3,0.3])
        plotter.fig.axes[1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(5))
        #see also: http://scipy-lectures.github.io/intro/matplotlib/matplotlib.html#ticks
        plotter.plotrgvars(begin = begin , stop =stop, axnum = 1)
        #plotter.separated = True

def togethermulliken():
    fname = 'results/noconstraineddmhermextrema100angstromgoodtrans/noconstrainedm_.dat'
    fname2 = 'results/noconstraineddmhermextremainfinitedistance/noconstrainedm_.dat'
    fname3 = 'results/noconstraineddmhermextrema0.1angstrom/noconstrainedm_.dat'
    fname4 = 'results/noconstraineddmhermextrema/noconstrainedm_.dat'

    plotter = Plot_Files([fname, fname2, fname3 , fname4])
    title = 'NO$^+$ Infinite Angstrom, 100angstrom (STO-3G)'
    xlim = None
    ylim = None

    plotter.data[0].depvar['xas'] = 'Mulliken charge on N'#change the future y-axis label 
    plotter.data[0].depvar['yas'] = 'Lambda at extremum'#change the future y-axis label 
    plotter.data[0].units['y'] = r'(a.u.)'
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [3], titel = title, name = 'plot' , exname = 'mullikenlambdatogether', save = False, label = '100angstrom')
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [3], titel = title, name = 'plot' , exname = 'mullikenlambdatogether', save = False, datanum = 1, label = 'infinite distance')
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [3], titel = title, name = 'plot' , exname = 'mullikenlambdatogether', save = False, datanum = 2, label = '0.1 angstrom')
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [3], titel = title, name = 'plot' , exname = 'mullikenlambdatogether', save = True, datanum = 3, label = '1.225 angstrom')
    plotter.data[0].depvar['yas'] = 'Energy'#change the future y-axis label 
    plotter.data[0].units['y'] = r'(E$_h$)'
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'mullikenenergytogether', save = False, label = '100 angstrom')
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'mullikenenergytogether', save = False, datanum = 1, label = 'infinite distance' )
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'mullikenenergytogether', save = False, datanum =2, label = '0.1 angstrom' )
    plotter.generate_plot(depcol = 0 , xlimg = xlim, ylimg = ylim, ylist = [2], titel = title, name = 'plot' , exname = 'mullikenenergytogether', save = True, datanum = 3, label = '1.225 angstrom')


def shannon_scatter():
    #fname = '/home/mario/DOCI-results/shannon_scattereq/shannon_scattereq_sto-3g.dat'
    fname = './results/shannon_plots/shannontotsenioritybeh2sto-3g.dat'
    #fname = 'results/shannon_scattereq6-31gfcivec/shannon_scattereq6-31gexcitation_6-31g.dat'
    plotter = Plot_Files(fname)
    plotter.data[0].depvar['yas'] = 'Energy'#change the future y-axis label 
    plotter.data[0].depvar['xas'] = '$I_c$'#change the future y-axis label 
    plotter.data[0].units['x'] = r''
    plotter.data[0].units['y'] = r'(E$_h$)'
    xlim = (0. , 10.)
    ylim = (None , None)
    #excitations = ['S','D','T','Q','5']
    seniorities = ['0','2','4','6']
    for i in range(2,5):
    #for i in range(2,7):
        #title = 'CI' +  ''.join(excitations[:i-1])
        title = 'SEN' +  ''.join(seniorities[:i-1])
        #plotter.generate_plot(xlimg = xlim, ylimg = ylim, ylist = [i], titel = title, name = 'plot' , exname = 'excitationscattereqsto-3g', exdir= '', prefix = False)
        plotter.generate_plot(xlimg = xlim, ylimg = ylim, ylist = [i], titel = title, name = 'plot' , exname = 'seniorityscattereqsto-3g', exdir= '', prefix = True)


if __name__ == '__main__':
    #main()  
    #makemovie()
    #shannon_scatter()
    togethermulliken()
