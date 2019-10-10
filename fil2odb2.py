# -*- coding: utf-8 -*-
##
## fin2odb.py
## PYTHON-Skript im Paket ABAQUSER 2.0
##
## Skript zum Erweitern einer ABAQUS *.odb durch UEL-Angaben:
##   - Elementeigenschaften aus *.info
##   - Elementdefinitionen aus *.inp
##   - Knotenvariable und SDV aus modifiziertem *.fil
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
print 'PYTHON-Skript fil2odb STARTED'
## Modulimports
from abaqusConstants import *
from odbAccess import *
from string import *
from time import *
import sys
import struct
import os
import subprocess
##########################################################################
##
## Subroutine Auslesen des *.inp in eine Liste
def read_inputfile(inpfile):
   A_name=''
   P_name=''
   I_name=''
   P_list=[]
   M_list=[]
   S_list=[]
   meshcounter=-1
   k_newpartorinstance=-1
   nncounter=0
   inputfile=open(inpfile, 'r')
   k_write=0
   for line in inputfile.readlines():        
        row=split(line)
        if len(row)!=0 and \
        (k_write==0  or (k_write!=0 and row[0][0]=='*')):
             k_write=0
             nncounter=nncounter+1
             if (upper(row[0])=='*NODE,') or \
                ((upper(row[0])=='*NODE') and (len(row)==1)):
                  k_write=1
                  if k_newpartorinstance==1:
                       meshcounter=meshcounter+1
                       k_newpartorinstance=0
                       M_list.append([meshcounter,A_name,I_name,P_name, \
                                     [],[]])
                  nset_name='NSETnn'+str(nncounter)
                  for i in range(len(row)-1):
                       param=row[i+1]
                       if upper(param[0])=='N':
                            nset_name=split(split(param,"=")[1],",")[0]
                  M_list[-1][4].append([nset_name,[]])
             elif (upper(row[0])=='*ELEMENT,'):
                  elset_name='ELSETnn'+str(nncounter)
                  for i in range(len(row)-1):
                       param=row[i+1]
                       if upper(param[0])=='E':
                            elset_name=split(split(param,"=")[1],",")[0]
                  for i in range(len(row)-1):
                       param=row[i+1]
                       if (upper(param[0])=='T') and \
                          (split(upper(param),"=")[1][0]=='U'):
                            k_write=2
                            k_morenodes=0
                            eltype=split(split(row[1],"=")[1],",")[0]
                            M_list[-1][5].append([elset_name,eltype,[]])
             elif (upper(row[0])=='*ASSEMBLY,'):
                  A_name=split(split(row[1],"=")[1],",")[0]
                  P_list.append([A_name,[]])
                  P_name=''
                  I_name=''
             elif (upper(row[0])=='*PART,'):
                  P_name=split(split(row[1],"=")[1],",")[0]
                  k_newpartorinstance=1
             elif (upper(row[0])=='*INSTANCE,'):
                  k_write=3
                  for i in range(len(row)-1):
                       parameter=row[i+1]
                       par=split(parameter,"=")[0]
                       parval=split(split(parameter,"=")[1],",")[0]
                       if upper(par)=='NAME':
                            I_name=parval
                       elif upper(par)=='PART':
                            P_name=parval
                  k_newpartorinstance=1
                  P_list[-1][-1].append([I_name,P_name,0, \
                                        [0,0,0,0,0,0,0,0,1,0]])
             elif (upper(row[0])=='*END') and (upper(row[1])=='ASSEMBLY'):
                  A_name=''
             elif (upper(row[0])=='*END') and (upper(row[1])=='INSTANCE'):
                  I_name=''
             elif (upper(row[0])=='*END') and (upper(row[1])=='PART'):
                  P_name=''
             elif (upper(row[0])=='*STEP,'):
                  for i in range(len(row)-1):
                       parameter=row[i+1]
                       par=split(parameter,"=")[0]
                       parval=split(split(parameter,"=")[1],",")[0]
                       if upper(par)=='NAME':
                            Step_name=parval
                  S_list.append(Step_name) 
        elif len(row)!=0 and k_write==1 and row[0][0]!='*':
             M_list[-1][4][-1][1].append([])
             M_list[-1][4][-1][1][-1].append(int(split(row[0],",")[0]))
             for i in range(len(row)-1):
                  M_list[-1][4][-1][1][-1].append( \
                                           float(split(row[i+1],",")[0]))
        elif len(row)!=0 and k_write==2 and row[0][0]!='*':
             if k_morenodes==0:
                  M_list[-1][5][-1][2].append([])
             for i in range(len(row)):
                  M_list[-1][5][-1][2][-1].append( \
                                           int(split(row[i],",")[0]))
             if row[i][-1]==',':
                  k_morenodes=1
             else:
                  k_morenodes=0
        elif len(row)!=0 and k_write==3 and row[0][0]!='*':
             for i in range(len(row)):
                  if len(row)==3:
                       P_list[-1][1][-1][3][i]=float(split(row[i],",")[0])
                  elif len(row)==7:
                       P_list[-1][1][-1][3][i+3]= \
                       float(split(row[i],",")[0])
   inputfile.close()
   for i in range(len(P_list)):
        A_name=P_list[i][0]
        for j in range(len(P_list[i][1])):
             I_name=P_list[i][1][j][0]
             P_name=P_list[i][1][j][1]   
             for k in range(len(M_list)):
                  P_name_M=M_list[k][3]
                  if P_name==P_name_M:
                       P_list[i][1][j][2]=k 
             for k in range(len(M_list)):
                  A_name_M=M_list[k][1]
                  I_name_M=M_list[k][2]
                  if ( (A_name==A_name_M) and (I_name==I_name_M) ):
                       P_list[i][1][j][2]=k 
   return [P_list,M_list,S_list]
##########################################################################
##
## Subroutine zum Auslesen des *.info in eine Liste
def read_infofile(inpfile):
   inputfile= open(inpfile, 'r')
   ueltypecounter=-1
   U_list=[]
   k_write=0
   for line in inputfile.readlines():
        row = split(line)
        if len(row)!=0 and (k_write==0  or (k_write!=0 and len(row[0])>1 \
                       and row[0][0]=='*' and row[0][1]!='*')):
             k_write=0
             if upper(row[0])=='*UEL,':
                  U_list.append([])
                  ueltypecounter=ueltypecounter+1
                  for i in range(5):
                       U_list[ueltypecounter].append([])
                  k_write=1
                  for param in range(len(row)-1):
                       if (param==0) or (param==4):
                            U_list[ueltypecounter][0].append( \
                            split(split(row[param+1],"=")[1],",")[0])
                       else:
                            U_list[ueltypecounter][0].append( \
                            int(split(split(row[param+1],"=")[1],",")[0]))
                  uelid=U_list[ueltypecounter][0][0]                  
             if upper(row[0])=='*SDV,':
                  uelid=split(row[1],"=")[1]                  
                  k_write=3
             if upper(row[0])=='*NV,':
                  uelid=split(row[1],"=")[1]                  
                  k_write=4
             if upper(row[0])=='*DIM':
                  k_write=5
        elif len(row)!=0 and k_write!=0 and row[0][0]!='*':
             for i in range(ueltypecounter+1):
                  if U_list[i][0][0]==uelid:
                       uelnumber=i
             if (k_write==3) or (k_write==4):
                  U_list[uelnumber][k_write].append([])
             for i in range(len(row)):
                  if k_write==2:
                       U_list[uelnumber][k_write].append(int(row[i]))
                  if k_write==1:
                       U_list[uelnumber][k_write].append(int(row[i]))
                       if i==len(row)-1:
                            k_write=2
                  if (k_write==3) or (k_write==4):
                       if i<3:
                            U_list[uelnumber][k_write][-1].append(row[i])
                       else:
                            U_list[uelnumber][k_write][-1].append( \
                                                           int(row[i]))
                  if k_write==5:
                       D=row[i]
   inputfile.close()
   if upper(D)=='THREE_D':
        Dim=THREE_D
   elif upper(D)=='TWO_D_PLANAR':
        Dim=TWO_D_PLANAR
   elif upper(D)=='AXISYMMETRIC':
        Dim=AXISYMMETRIC
   return [U_list,Dim]
##########################################################################
##
## Subroutine zur Erzeugung einer Liste zur Zuordnung v. Komponentenlabels 
## und Invarianten in Abhaengigkeit der Tensorstufe
def GetTensorprop():
   Tensorprop=[]
   Tensorprop.append([])
   Tensorprop[-1].append(SCALAR)
   Tensorprop[-1].append([''])
   Tensorprop[-1].append([])
   Tensorprop.append([])
   Tensorprop[-1].append(VECTOR)
   Tensorprop[-1].append(['1','2','3'])
   Tensorprop[-1].append([MAGNITUDE])
   Tensorprop.append([])
   Tensorprop[-1].append(TENSOR_3D_FULL)
   Tensorprop[-1].append(['11','22','33','12','13','23'])
   Tensorprop[-1].append([MISES, TRESCA, PRESS, INV3, MAX_PRINCIPAL, \
                          MID_PRINCIPAL, MIN_PRINCIPAL])
   Tensorprop.append([])
   Tensorprop[-1].append(TENSOR_3D_SURFACE)
   Tensorprop[-1].append(['11','22','12'])
   Tensorprop[-1].append([MAX_PRINCIPAL, MIN_PRINCIPAL, \
                          MAX_INPLANE_PRINCIPAL, MIN_INPLANE_PRINCIPAL])
   Tensorprop.append([])
   Tensorprop[-1].append(TENSOR_3D_PLANAR)
   Tensorprop[-1].append(['11','22','33','12'])
   Tensorprop[-1].append([MISES, TRESCA, PRESS, INV3, MAX_PRINCIPAL, \
                          MID_PRINCIPAL, MIN_PRINCIPAL, \
                          MAX_INPLANE_PRINCIPAL, MIN_INPLANE_PRINCIPAL, \
                          OUTOFPLANE_PRINCIPAL])
   Tensorprop.append([])
   Tensorprop[-1].append(TENSOR_2D_SURFACE)
   Tensorprop[-1].append(['11','22','12'])
   Tensorprop[-1].append([MAX_PRINCIPAL, MIN_PRINCIPAL, \
                          MAX_INPLANE_PRINCIPAL, MIN_INPLANE_PRINCIPAL])
   Tensorprop.append([])
   Tensorprop[-1].append(TENSOR_2D_PLANAR)
   Tensorprop[-1].append(['11','22','33','12'])
   Tensorprop[-1].append([MISES, TRESCA, PRESS, INV3, MAX_PRINCIPAL, \
                          MID_PRINCIPAL, MIN_PRINCIPAL, \
                          MAX_INPLANE_PRINCIPAL, MIN_INPLANE_PRINCIPAL, \
                          OUTOFPLANE_PRINCIPAL])
   return [Tensorprop]
##########################################################################
##
## Subroutine zur Berechnung der Transformationsmatrix (homogene 
## Koordinaten), siehe Maple-Sheet
def GetTrafomat(T_list):
   Pi=acos(-1)
   Ux=T_list[0]
   Uy=T_list[1]
   Uz=T_list[2]
   Vx=T_list[6]-T_list[3]
   Vy=T_list[7]-T_list[4]
   Vz=T_list[8]-T_list[5]
   P1x=T_list[3]
   P1y=T_list[4]
   P1z=T_list[5]
   P2x=T_list[6]
   P2y=T_list[7]
   P2z=T_list[8]
   phi=T_list[9]*Pi/180
   abs2=sqrt(Vx**2+Vy**2)
   abs3=sqrt(Vx**2+Vy**2+Vz**2)
   if abs2>0:
        T11 = ((Vx*Vz*cos(phi)/(abs2*abs3)+Vy*sin(phi)/abs2)*Vz/abs3+Vx* \
              abs2/abs3**2)*Vx/abs2+(Vx*Vz*sin(phi)/(abs2*abs3)-Vy* \
              cos(phi)/abs2)*(-Vy)/abs2
        T21 = ((Vx*Vz*cos(phi)/(abs2*abs3)+Vy*sin(phi)/abs2)*Vz/abs3+Vx* \
              abs2/abs3**2)*Vy/abs2+(Vx*Vz*sin(phi)/(abs2*abs3)-Vy* \
              cos(phi)/abs2)*Vx/abs2
        T31 = -(Vx*Vz*cos(phi)/(abs2*abs3)+Vy*sin(phi)/abs2)* \
              abs2/abs3+Vx*Vz/abs3**2
        T12 = ((Vy*Vz*cos(phi)/(abs2*abs3)-Vx*sin(phi)/abs2)*Vz/abs3+Vy* \
              abs2/abs3**2)*Vx/abs2+(Vy*Vz*sin(phi)/(abs2*abs3)+Vx* \
              cos(phi)/abs2)*(-Vy)/abs2
        T22 = ((Vy*Vz*cos(phi)/(abs2*abs3)-Vx*sin(phi)/abs2)*Vz/abs3+Vy* \
              abs2/abs3**2)*Vy/abs2+(Vy*Vz*sin(phi)/(abs2*abs3)+Vx* \
              cos(phi)/abs2)*Vx/abs2
        T32 = -(Vy*Vz*cos(phi)/(abs2*abs3)-Vx*sin(phi)/abs2)* \
              abs2/abs3+Vy*Vz/abs3**2
        T13 = (-abs2*cos(phi)*Vz/abs3**2+Vz*abs2/abs3**2)* \
              Vx/abs2-sin(phi)*(-Vy)/abs3
        T23 = (-abs2*cos(phi)*Vz/abs3**2+Vz*abs2/abs3**2)* \
              Vy/abs2-sin(phi)*Vx/abs3
        T33 = abs2**2*cos(phi)/abs3**2+Vz**2/abs3**2
        T41 = -(-(((Vx*Vz*cos(phi)/(abs2*abs3)+Vy*sin(phi)/abs2)*Vz/ \
              abs3+Vx*abs2/abs3**2)*Vx/abs2+(Vx*Vz*sin(phi)/(abs2*abs3)- \
              Vy*cos(phi)/abs2)*(-Vy)/abs2)*P1x-(((Vx*Vz*cos(phi)/(abs2* \
              abs3)+Vy*sin(phi)/abs2)*Vz/abs3+Vx*abs2/abs3**2)*Vy/abs2+ \
              (Vx*Vz*sin(phi)/(abs2*abs3)-Vy*cos(phi)/abs2)*Vx/abs2)* \
              P1y-(-(Vx*Vz*cos(phi)/(abs2*abs3)+Vy*sin(phi)/abs2)*abs2/ \
              abs3+Vx*Vz/abs3**2)*P1z+P1x)
        T42 = -(-(((Vy*Vz*cos(phi)/(abs2*abs3)-Vx*sin(phi)/abs2)*Vz/ \
              abs3+Vy*abs2/abs3**2)*Vx/abs2+(Vy*Vz*sin(phi)/(abs2*abs3)+ \
              Vx*cos(phi)/abs2)*(-Vy)/abs2)*P1x-(((Vy*Vz*cos(phi)/(abs2* \
              abs3)-Vx*sin(phi)/abs2)*Vz/abs3+Vy*abs2/abs3**2)*Vy/abs2+ \
              (Vy*Vz*sin(phi)/(abs2*abs3)+Vx*cos(phi)/abs2)*Vx/abs2)* \
              P1y-(-(Vy*Vz*cos(phi)/(abs2*abs3)-Vx*sin(phi)/abs2)*abs2/ \
              abs3+Vy*Vz/abs3**2)*P1z+P1y)
        T43 = -(-((-abs2*cos(phi)*Vz/abs3**2+Vz*abs2/abs3**2)*Vx/abs2- \
              sin(phi)*(-Vy)/abs3)*P1x-((-abs2*cos(phi)*Vz/abs3**2+Vz* \
              abs2/abs3**2)*Vy/abs2-sin(phi)*Vx/abs3)*P1y-(abs2**2* \
              cos(phi)/abs3**2+Vz**2/abs3**2)*P1z+P1z)
   else:
        T11 = cos(phi)
        T21 = Vz*sin(phi)/abs3
        T31 = 0.
        T12 = -Vz*sin(phi)/abs3
        T22 = cos(phi)
        T32 = 0.
        T13 = 0.
        T23 = 0.
        T33 = 1.
        T41 = -(-Vz*sin(phi)*P1y/abs3-cos(phi)*P1x+P1x)
        T42 = -( Vz*sin(phi)*P1x/abs3-cos(phi)*P1y+P1y)
        T43 = 0.
   T=[]
   T.append([T11,T21,T31])
   T.append([T12,T22,T32])
   T.append([T13,T23,T33])
   T.append([T41+Ux,T42+Uy,T43+Uz])
   return [T]
##########################################################################
##
## Subroutine zum Auslesen des *.fil und Aktualisieren des *.odb
def readfil_updateodb(filfile,odb_name,S_list,P_list,M_list,U_list,Dim):
   A='<8s'
   D='<1d'
   I='<ii'
   #######################################################################
   ## Generieren eines neuen .fil mit Beginn bei 4byte, ab da alle 512byte 
   ## Bloecke ohne Zwischenbyte
   #######################################################################
   ##
   lt = localtime()
   print strftime("      modify fil STARTED  : %H:%M:%S", lt), '...'
   file=open(filfile,'rb')
   filestat=os.stat(filfile)
   filesize=filestat.st_size 
   print '         Read ', filfile, ', ', filesize, ' bytes'
   filfile_red=split(filfile,'.')[0]+'_red.fil'
   file_red=open(filfile_red,'wb')  
   file_red.close()
   file_red=open(filfile_red,'ab')  
   file.seek(4)
   k_EOF=0
   i=0
   while k_EOF==0:
        file.seek(4+i*(8*513))
        file_red.write(file.read(8*512))
        file_red.flush()
        if file.read(1)=='':
             k_EOF=1
             print '         EOF after ', i*(8*512), 'bytes from ', \
                    filfile, '(', filesize, ' bytes)'
        i=i+1
   file_red.close()
   file.close()
   lt = localtime()
   print strftime("      modify fil COMPLETED: %H:%M:%S", lt), '...'
   #######################################################################
   ## Aktualisieren des .odb um Modelldefinitionen (Knoten,Elemente,Parts)
   #######################################################################
   ##
   lt = localtime()
   print strftime("      append model data    STARTED  : %H:%M:%S", lt), \
          '...',
   [Tensorprop] = GetTensorprop()
   Updateodb = openOdb(path=odb_name+'.odb', readOnly=False)
   Part_list=[]
   for i in range(len(M_list)):
        A_name=M_list[i][1]
        I_name=M_list[i][2]
        P_name=M_list[i][3]
        if I_name=='':
             UEL_partname = P_name+'_U'
        else:
             UEL_partname = A_name+'.'+I_name+'.'+P_name+'_U'
        UEL_part = Updateodb.DiscretePart(name=UEL_partname, \
                   embeddedSpace=Dim, type=DEFORMABLE_BODY)
        Part_list.append(UEL_part)
        node_max=0
        for nodeset in M_list[i][4]:
             UEL_part.addNodes(nodeData=nodeset[1],nodeSetName=nodeset[0])
             for node_row in nodeset[1]:
                  if node_row[0]>node_max:
                       node_max=node_row[0] 
        for elemset in M_list[i][5]:
             for j in range(len(U_list)):
                  if U_list[j][0][0]==elemset[1]:
                       elem_id=j
             elemdef=[]
             for k in range(len(elemset[2])):
                  elemdef.append([elemset[2][k][0]])
                  for j in range(len(U_list[elem_id][1])):
                       elemdef[-1].append( \
                                   elemset[2][k][U_list[elem_id][1][j]])
             UEL_part.addElements(elementData=elemdef, \
                                  type=U_list[elem_id][0][4], \
                                  elementSetName=elemset[0])
   for i in range(len(P_list)):
        for j in range(len(P_list[i][1])):
             instname = P_list[i][1][j][0]+'_U'
             [Trafomat]=GetTrafomat(P_list[i][1][j][3])
             instance_new=Updateodb.rootAssembly.Instance(name=instname, \
                          object=Part_list[P_list[i][1][j][2]], \
                          localCoordSystem=Trafomat)
   Updateodb.save()
   lt = localtime()
   print strftime("COMPLETED: %H:%M:%S", lt)
   #######################################################################
   ## inkrementweises Auslesen des reduzierten *.fil  
   ## und Aktualisieren des *.odb
   #######################################################################
   ##
   file=open(filfile_red,'rb')
   nv_typ='dummy'
   k_nostep=1
   byte=0
   file.seek(byte)
   k_EOF=0
   while k_EOF==0:
        file.seek(byte)
        num=struct.unpack(I,file.read(8))[0]
        rc=struct.unpack(I,file.read(8))[0]
	#print 'num',num
	#print 'rc',rc
	if rc==2000:
             file.seek(8,1)
             steptime=struct.unpack(D,file.read(8))[0]
             file.seek(24,1)
             new_step=struct.unpack(I,file.read(8))[0]
             if k_nostep==1:
                  step_number=0
                  k_nostep=0
                  step=new_step
             elif (k_nostep==0) and (step!=new_step): 
                  step_number=step_number+1
                  step=new_step
             IP_list=[]
             N_list=[]
             lt = localtime()
             print '     ', '%8s' % ljust(S_list[step_number],8), \
                   ' %10.8f' % steptime, \
                   strftime("STARTED  : %H:%M:%S", lt), '...',
        elif rc==1911:
             if num==4:
                  k_output=1
             elif num==5:
                  k_output=0
                  file.seek(16,1)
                  elemtyp=struct.unpack(A,file.read(8))[0]
                  for i in range(len(U_list)):
                       if U_list[i][0][0]==elemtyp[:len(U_list[i][0][0])]:
                            elem_id=i
                  num_SDV_var=U_list[elem_id][0][2] 
                  num_ip=U_list[elem_id][0][1]
                  num_sdvpip=U_list[elem_id][0][3]
                  num_SDV=len(U_list[elem_id][3]) 
                  SDV_counter=0
                  k_moreSDV=0
        elif rc==1:
             elem=struct.unpack(I,file.read(8))[0]
        elif rc==5:
	     #print 'k,cou,num',k_moreSDV,SDV_counter,num_SDV_var
             if k_moreSDV==0:
                  SDV_list=[]
             for i in range(num-2):
                  SDV_list.append(struct.unpack(D,file.read(8))[0])
             SDV_counter=SDV_counter+num-2
             if SDV_counter<num_SDV_var:
                  k_moreSDV=1
             else:
                  k_moreSDV=0
		  SDV_counter=0
             if k_moreSDV==0:
                  IP_list.append([elem])
                  for i in range(num_SDV):
                       IP_list[-1].append([])
                       for j in range(num_ip):
                            IP_list[-1][-1].append([])
                            for k in range(len(U_list[elem_id][3][i])-3):
                                 IP_list[-1][-1][-1].append( \
                                 SDV_list[j*num_sdvpip+ \
                                          U_list[elem_id][3][i][k+3]-1])
        elif rc>100 and rc<238:
             if rc==101:
                  nv_typ='U'
             elif rc==102:
                  nv_typ='V'
             elif rc==103:
                  nv_typ='A'
             elif rc==104:
                  nv_typ='RF'
             elif rc==105:
                  nv_typ='EPOT'
             elif rc==107:
                  nv_typ='COORD'
             elif rc==108:
                  nv_typ='POR'
             elif rc==201:
                  nv_typ='NT11'
             nv_pos=len(N_list)
             for i in range(nv_pos): 
                  if N_list[i][0]==nv_typ:
                       nv_pos=i   
             if nv_pos==len(N_list):
                  N_list.append([nv_typ,[],[]])  
             num_NV=num-3
             node=struct.unpack(I,file.read(8))[0]
             N_list[nv_pos][1].append(node)
             N_list[nv_pos][2].append([])
             for i in range(num_NV):
                  N_list[nv_pos][2][-1].append( \
                                        struct.unpack(D,file.read(8))[0])
        elif rc==2001 and k_nostep==0:
             Step=Updateodb.steps[S_list[step_number]]
             Step_frame=Step.getFrame(frameValue=steptime, match=CLOSEST)
             D_list=[]
             for value_row in N_list:
                  node=value_row[1]
                  inputvalues=value_row[2]
                  typ=value_row[0]
                  pos=NODAL
                  D_list.append([typ,'',[],'','',pos,node,inputvalues])
             for elem_row in IP_list:
                  elem=elem_row[0]
                  for part_i in M_list:
                       for eset_i in part_i[5]:
                            for element_i in eset_i[2]:
                                 if element_i[0]==elem:
                                      elem_typ=eset_i[1]
                  for i in range(len(U_list)):
                       if U_list[i][0][0]==elem_typ:
                            elem_id=i
                  num_SDV=len(U_list[elem_id][3]) 
                  for i in range(num_SDV):
                       values=elem_row[i+1]
                       typ=U_list[elem_id][3][i][1]
                       desc=U_list[elem_id][3][i][0]
                       type=U_list[elem_id][3][i][2]
                       for j in range(len(Tensorprop)):
                            if str(Tensorprop[j][0])==type:
                                 Tensorprop_id=j
                                 type=Tensorprop[j][0]
                       clabel=[]
                       for j in range(len(values[0])):
                            clabel.append( \
                                   typ+Tensorprop[Tensorprop_id][1][j])                 
                       label=[elem]
                       inv = Tensorprop[Tensorprop_id][2]
                       inputvalues=[]
                       for j in range(len(U_list[elem_id][2])):
                           inputvalues.append( \
                                       values[U_list[elem_id][2][j]-1])
                       pos=INTEGRATION_POINT
                       k_fd=-1
                       for j in range(len(D_list)):
                            fielddata=D_list[j]
                            if fielddata[0]==typ and fielddata[5]==pos:
                                 k_fd=j
                       if k_fd==-1:
                            k_fd=len(D_list)
                            D_list.append( \
                                   [typ,desc,clabel,type,inv,pos,[],[]])
                       D_list[k_fd][6].append(elem)
                       for j in range(len(inputvalues)):
                            D_list[k_fd][7].append(inputvalues[j])
             for fielddata in D_list:
                  if Step_frame.fieldOutputs.keys().count(fielddata[0])==0:
                       Field_new=Step_frame.FieldOutput( \
                                 name=fielddata[0], \
                                 description=fielddata[1], \
                                 componentLabels=fielddata[2], \
                                 type=fielddata[3], \
                                 validInvariants=fielddata[4] )
                  else:
                       Field_new = Step_frame.fieldOutputs[fielddata[0]]
  		  #print 'f0',fielddata[0]
		  #print 'f1',fielddata[1]
		  #print 'f2',fielddata[2]
		  #print 'f3',fielddata[3]
		  #print 'f4',fielddata[4]
		  #print 'f5',fielddata[5]
		  #print 'f6',fielddata[6]
		  #print 'f7',fielddata[7]



                  Field_new.addData(position=fielddata[5], \
                                    instance=instance_new, \
                                    labels=fielddata[6], \
                                    data=fielddata[7])
             Updateodb.save()
             lt = localtime()
             print strftime("COMPLETED: %H:%M:%S", lt)
        byte=byte+8*num
        file.seek(byte)
        if file.read(1)=='':
             k_EOF=1
             print '      EOF after ', byte, 'bytes'
   file.close()
   subprocess.call('rm '+filfile_red, shell=True)
   Updateodb.close()
   filestat=os.stat(odb_name+'.odb')
   filesize=filestat.st_size 
   print '      Filesize ODB: ', filesize, 'bytes'
##########################################################################
##
## Hauptprogramm
jobname=sys.argv[1]
infofile=sys.argv[2]
inpfile=jobname+'.inp'
filfile=jobname+'.fil'
lt = localtime()
print strftime("   Read INFO     STARTED  : %H:%M:%S", lt), "...",
[U_list,Dim]=read_infofile(infofile)
lt = localtime()
print strftime("COMPLETED: %H:%M:%S", lt)
lt = localtime()
print strftime("   Read INP      STARTED  : %H:%M:%S", lt), "...",
[P_list,M_list,S_list]=read_inputfile(inpfile)
lt = localtime()
print strftime("COMPLETED: %H:%M:%S", lt)
lt = localtime()
print strftime("   Update ODB    STARTED  : %H:%M:%S", lt), "..."
readfil_updateodb(filfile,jobname,S_list,P_list,M_list,U_list,Dim)
lt = localtime()
print strftime("   Update ODB    COMPLETED: %H:%M:%S", lt)
print 'PYTHON-Skript fil2odb COMPLETED'
##########################################################################
