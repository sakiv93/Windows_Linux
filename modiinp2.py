# -*- coding: utf-8 -*-
##
## modiinp.py
## PYTHON-Skript im Paket ABAQUSER 2.0
##
## Skript zum Erweitern einer ABAQUS *.inp um die Befehle zum Erzeugen
## eines *.fil.
##
## Das Skript ist nicht ausgetestet.
##
## Die Verwendung des Skriptes ist freigegeben
## Aenderungen im Skript sind zu kennzeichnen, im Header zu vermerken
## und dem Autor mitzuteilen.  
##
## Stephan Roth
## Stephan.Roth@imfd.tu-freiberg.de
##
## Freiberg, 12.11.2010
##
##########################################################################
##
print '\nPYTHON-Skript modiinp STARTED'
## Modulimports
from string import *
import subprocess
from time import *
import sys
##########################################################################
##
jobname=sys.argv[1]
inputfile=  open(jobname+'.inp', 'r')
outputfile= open(jobname+'.inp_temp', 'w')
outputfile.close()
outputfile= open(jobname+'.inp_temp', 'a')
k_fileformat=0
k_fieldoutput=0
k_nodeout=0
k_oldfil=0
for line in inputfile.readlines():        
     row=split(line)
     k_write=1
     if len(row)!=0 and ( k_oldfil==0 or ( k_oldfil==1 and  \
        row[0][0]=='*' and row[0][1]!='*') ):
          k_oldfil=0
          if len(row)>1 and upper(row[0])=='*OUTPUT,' and \
             upper(row[1][0])=='F':
               k_fieldoutput=1
          elif len(row)>1 and ( (upper(row[0])=='*OUTPUT,' and \
               upper(row[1][0])!='F') or upper(row[0])=='*STEP,' ):
               k_fieldoutput=0
          elif len(row)>1 and upper(row[0])=='*NODE' and \
               upper(row[1][0])=='O' and k_fieldoutput==1:
               k_nodeout=1
               nodeout_list=[]
          elif len(row)!=0 and k_nodeout==1 and row[0][0]!='*':
               nodeout_list.append(line)
          elif len(row)!=0 and k_nodeout==1 and row[0][0]=='*':
               k_nodeout=0
          elif len(row)>1 and upper(row[0])=='*END' and \
               upper(row[1])=='STEP':
               if k_fileformat==0:
                    k_fileformat=1
                    outputfile.writelines( \
                    '*FILE FORMAT, ZERO INCREMENT\n')               
               outputfile.writelines('*EL FILE, FREQUENCY=1\n')               
               outputfile.writelines('\n')               
               outputfile.writelines('SDV\n')               
               outputfile.writelines('*NODE FILE, FREQUENCY=1\n')               
               for nodeout_line in nodeout_list:
                    outputfile.writelines(nodeout_line)               
          elif len(row)>1 and upper(row[0])=='*FILE' and \
               upper(row[1][0])=='F':
               print \
               '   WARNING: FOUND KEYWORD *FILE FORMAT, CHECK INPUTFILE'
               k_write=0
          elif len(row)>1 and ( (upper(row[0])=='*EL' and \
               upper(row[1][0])=='F') or (upper(row[0])=='*NODE' and \
               upper(row[1][0])=='F') ):
               print \
               '   WARNING: FOUND RESULTS FILE REQUEST, CHECK INPUTFILE'
               k_oldfil=1
               k_write=0
     elif k_oldfil==1 and ((len(row)!=0 and ((len(row[0])>1 and \
           row[0][0]!='*' and row[0][1]!='*') or len(row[0])<2)) or \
           len(row)==0):
          k_write=0
     if k_write==1:   
          outputfile.writelines(line)
inputfile.close()
outputfile.close()
subprocess.call('rm '+jobname+'.inp', shell=True) 
subprocess.call('mv '+jobname+'.inp_temp'+' '+jobname+'.inp', shell=True) 
print 'PYTHON-Skript modiinp COMPLETED\n'
##########################################################################
