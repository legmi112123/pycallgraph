#!/usr/bin/python

# Therese: 10-Jul-2015
# generate the dmdgp file including all i_i-4 restraints along the
# atom order, the restraints being calculated from a PDB file
# syntaxe: ./generate_dmdgp_i_i-6_ss5_CBO.py 2LXZ_model1.pdb 2lxz_NH2.dmdgp > 2lxz_PDBe_i-i-6_ss5_CBO_NH2.dmdgp

import os,sys,math

filepdb = sys.argv[1]		#takes 2 inputs
filedmdgp = sys.argv[2]

def calcdst(liste_order,i,j,xat,yat,zat): #compute distance between 2 3D points liste_order[i] liste_order[j]

  idx1 = liste_order[i]
  idx2 = liste_order[j]
  #print "idx1",idx1,"idx2",idx2
  dx = xat[idx1] - xat[idx2]
  dx2 = dx*dx
  dy = yat[idx1] - yat[idx2]
  dy2 = dy*dy
  dz = zat[idx1] - zat[idx2]
  dz2 = dz*dz
  dist2 = dx2+dy2+dz2
  dist = math.sqrt(dist2)

  return dist

def calcdstCBO(key1,key2,xat,yat,zat):

  dx = xat[key1] - xat[key2]
  dx2 = dx*dx
  dy = yat[key1] - yat[key2]
  dy2 = dy*dy
  dz = zat[key1] - zat[key2]
  dz2 = dz*dz
  dist2 = dx2+dy2+dz2
  dist = math.sqrt(dist2)

  return dist

def outedge(verti,vertj,dij,liste_order,ii,jj,edgeline): 

  if (verti < vertj) & (dij != 0):
     if not edgeline.has_key(str(verti)+" "+str(vertj)):
        edgeline[str(verti)+" "+str(vertj)] = str(verti)+"\t"+str(vertj)+"\t"+"D\t"+str(dij)+"\t # "+liste_order[ii]+" "+liste_order[jj]
        #print str(verti)+"\t"+str(vertj)+"\t"+"D\t"+str(dij)+"\t # "+liste_order[ii]+" "+liste_order[jj]

  if (verti > vertj) & (dij != 0):
     if not edgeline.has_key(str(vertj)+" "+str(verti)):
        edgeline[str(vertj)+" "+str(verti)] = str(vertj)+"\t"+str(verti)+"\t"+"D\t"+str(dij)+"\t # "+liste_order[ii]+" "+liste_order[jj]
        #print str(vertj)+"\t"+str(verti)+"\t"+"D\t"+str(dij)+"\t # "+liste_order[ii]+" "+liste_order[jj]

  return 

def outedge2(idxi,idxj,dij,keyi,keyj,edgeline): 

  if idxi < idxj:
     if not edgeline.has_key(str(idxi)+" "+str(idxj)):
        edgeline[str(idxi)+" "+str(idxj)] = str(idxi)+"\t"+str(idxj)+"\t"+"D\t"+str(dij)+"\t # "+keyi+" "+keyj
        #print str(idxi)+"\t"+str(idxj)+"\t"+"D\t"+str(dij)+"\t # "+keyi+" "+keyj
  if idxi > idxj:
     if not edgeline.has_key(str(idxj)+" "+str(idxi)):
        edgeline[str(idxj)+" "+str(idxi)] = str(idxj)+"\t"+str(idxi)+"\t"+"D\t"+str(dij)+"\t # "+keyj+" "+keyi
        #print str(idxj)+"\t"+str(idxi)+"\t"+"D\t"+str(dij)+"\t # "+keyj+" "+keyi

  return 

def outedge3(idxi,idxj,dij,keyi,keyj,edgeline): 

  if idxi < idxj:
     edgeline[str(idxi)+" "+str(idxj)] = str(idxi)+"\t"+str(idxj)+"\t"+"D\t"+str(dij)+"\t # "+keyi+" "+keyj
     #print str(idxi)+"\t"+str(idxj)+"\t"+"D\t"+str(dij)+"\t # "+keyi+" "+keyj
  if idxi > idxj:
     edgeline[str(idxj)+" "+str(idxi)] = str(idxj)+"\t"+str(idxi)+"\t"+"D\t"+str(dij)+"\t # "+keyj+" "+keyi
     #print str(idxj)+"\t"+str(idxi)+"\t"+"D\t"+str(dij)+"\t # "+keyj+" "+keyi

  return 

# main 
# read the PDB atom coordinates 
pdb = open(filepdb,"r")
lines = pdb.readlines()
pdb.close()

#ATOM     60  N   LEU A   1     -10.336   6.233   4.888  1.00  0.00           N
#ato = []
#res = []
#nor = []
xat = {}
yat = {}
zat = {}
for line in lines:
    if line[0:4] == "ATOM": 
       nomato = line[12:16]
       #print "nomato avant",nomato,len(nomato)
       ato = line[12:16].strip()
       res = line[17:20].strip()
       nor = line[22:26].strip()
       #print "nomato apres",nomato,len(nomato)
       #ato.append(nomato)
       #res.append(line[17:20])
       #nor.append(int(line[22:26]))
       #print line
       #print "ato",ato,"res",res,"nor",nor
       if (res == "GLY") & (ato == "HA3"):
          ato = "HA1"
       idx = "".join([ato,"-",res,"-",nor])
       #xat.append(float(line[30:38]))
       xat[idx] = float(line[30:38])
       #yat.append(float(line[38:46]))
       yat[idx] = float(line[38:46])
       #zat.append(float(line[46:54]))
       zat[idx] = float(line[46:54])
       #print "idx",idx

# get the maximum residue number (the sequence is supposed to start at 1)
nbres = int(nor)

# read the order list of atoms from dmdgp
dmd = open(filedmdgp,"r")
lines = dmd.readlines()
dmd.close()

atoms = 0
residues = 0
vertices = 0
order = 0
idx = []
nom = {}
resid = {}
numres = {}
index_vert = {}
convert_atom_names = {'C1': 'C', 'CA1': 'CA', 'CB1': 'CB', 'HA1': 'HA', 'HA2': 'HA2', 'N1': 'N', 'O1': 'O', 'O2': 'O2', 'H1': 'H', 'H2': 'H2'}
liste_order = []
liste_vertidx = []
resname = []

for line in lines:
    #print "line16",line[0:16]
    if line[0:14] == "end atom_names":
       atoms = 0
    if line[0:12] == "end residues":
       residues = 0
       # determine the inverse converter array, giving the vertice idx from the
       # atom name, residue name and residue number
       for ii in idx:
         #print "idx",idx
         #print "ii",ii
         #print "nom tot",nom
         #print "nom",nom[str(ii)],ii
         if nom[str(ii)] in convert_atom_names.keys(): 
            namevert = convert_atom_names[nom[str(ii)]]
         else:
            namevert = nom[str(ii)]
         if (namevert == "H") & (int(numres[str(ii)]) == 1): 
            namevert = "H1"
         if (namevert == "O") & (int(numres[str(ii)]) == nbres): 
            #print "coucou",namevert,int(numres[str(ii)])
            namevert = "O1"
         #print "namevert",namevert,"int(numres[str(ii)])",int(numres[str(ii)]),"nbres",nbres
         key = namevert+"-"+resid[str(ii)]+"-"+numres[str(ii)]
         index_vert[key] = ii
         #print "index_vert[key]",index_vert[key],"key",key
        
    if line[0:12] == "end bp_order":
       order = 0
    if line[0:12] == "end vertices":
       vertices = 0
       
    if vertices == 1:
       #print "line vert",line
       ff = line.split()
       idx.append(int(ff[0]))
       #print "idx",idx

    if atoms == 1:
       ff = line.split()
       nom_atome = ff[0]
       for ii in range(1,len(ff)):
           nom[ff[ii]] = ff[0]

    if residues == 1:
       #print line
       ff = line.split()
       resname.append(ff[1])
       for ii in range(2,len(ff)):
           resid[ff[ii]] = ff[1]
           numres[ff[ii]] = ff[0]

    if order == 1:
       ff = line.split()
       idxvert = ff[0]
       #print "idxvert",idxvert,"nom[idxvert]",nom[idxvert]
       liste_vertidx.append(int(idxvert))
       if nom[str(ii)] in convert_atom_names.keys(): 
          #print "nom[str(ii)]",nom[str(ii)],"str(ii)",str(ii)
          #print "nom[idxvert]",nom[idxvert]
          namevert = convert_atom_names[nom[idxvert]]
       else:
          namevert = nom[idxvert]
       if (namevert == "H") & (int(numres[idxvert]) == 1):
          namevert = "H1"
       if (namevert == "O") & (int(numres[idxvert]) == nbres):
          namevert = "O1"
       if (resid[idxvert] == "GLY") & (nom[idxvert] == "HA1"):
          namevert = 'HA1'
       residue = "".join([resid[idxvert],"-",numres[idxvert]]) 
       #print idxvert,namevert,residue
       liste_order.append("".join([namevert,"-",resid[idxvert],"-",numres[idxvert]]))

    if line[0:16] == "begin atom_names":
       atoms = 1
    if line[0:14] == "begin residues":
       residues = 1
    if line[0:14] == "begin bp_order":
       order = 1
    if line[0:14] == "begin vertices":
       vertices = 1

#print "liste_order",liste_order
#print "liste_vertidx",liste_vertidx
#print "len",len(liste_order)

# read again the dmdgp file in order to produce the new one
dmd = open(filedmdgp,"r")
lines = dmd.readlines()
dmd.close()

edges = 0
dihed = 0
edgeline = {}

for line in lines:
    
    #print "edges",edges,"line",line[0:len(line)-1]

    if line[0:9] == "end edges":
       edges = 0
    if line[0:19] == "end dihedral_angles":
       dihed = 0

    if (edges == 1):
      if first == 1:
       # generate all i_i-5 distances 
       for ii in range(0,len(liste_order)-7):
          #i1 = liste_[ii]
          #i2 = liste_vertidx[ii+1]
          #i3 = liste_vertidx[ii+2]
          #i4 = liste_vertidx[ii+3]
          i1 = ii
          i2 = ii+1
          i3 = ii+2
          i4 = ii+3
          i5 = ii+4
          i7 = ii+6
          vert1 = liste_vertidx[i1]
          vert2 = liste_vertidx[i2]
          vert3 = liste_vertidx[i3]
          vert4 = liste_vertidx[i4]
          vert5 = liste_vertidx[i5]
          vert7 = liste_vertidx[i7]
          #print i1,i2,i3,i4
          nom1 = liste_order[i1].split('-')[0]
          nom2 = liste_order[i2].split('-')[0]
          nom3 = liste_order[i3].split('-')[0]
          nom4 = liste_order[i4].split('-')[0]
          nom5 = liste_order[i5].split('-')[0]
          nom7 = liste_order[i7].split('-')[0]
          #print "nom1",nom1,"nom2",nom2,"nom3",nom3,"nom4",nom4,"nom5",nom5,"nom7",nom7
          res1 = int(liste_order[i1].split('-')[2])
          res2 = int(liste_order[i2].split('-')[2])
          #print "res1",res1,"res2",res2
          if (nom1 != "H1") & (nom1 != "H2") & (nom2 != "H1") & (nom2 != "H2") & (nom1 != "O1") & (nom1 != "O2") & (nom2 != "O1") & (nom2 != "O2"): 
            d12 = calcdst(liste_order,i1,i2,xat,yat,zat)
            outedge(vert1,vert2,d12,liste_order,i1,i2,edgeline)
          #print i1,i2,d12,liste_order[i1],liste_order[i2]
          if (nom1 != "H1") & (nom1 != "H2") & (nom3 != "H1") & (nom3 != "H2") & (nom1 != "O1") & (nom1 != "O2") & (nom3 != "O1") & (nom3 != "O2"): 
            d13 = calcdst(liste_order,i1,i3,xat,yat,zat)
            outedge(vert1,vert3,d13,liste_order,i1,i3,edgeline)
          if (nom1 != "H1") & (nom1 != "H2") & (nom4 != "H1") & (nom4 != "H2") & (nom1 != "O1") & (nom1 != "O2") & (nom4 != "O1") & (nom4 != "O2"): 
            d14 = calcdst(liste_order,i1,i4,xat,yat,zat)
            outedge(vert1,vert4,d14,liste_order,i1,i4,edgeline)
          if (nom1 != "H1") & (nom1 != "H2") & (nom5 != "H1") & (nom5 != "H2") & (nom1 != "O1") & (nom1 != "O2") & (nom5 != "O1") & (nom5 != "O2"): 
            d15 = calcdst(liste_order,i1,i5,xat,yat,zat)
            outedge(vert1,vert5,d15,liste_order,i1,i5,edgeline)
          if (nom1 != "H1") & (nom1 != "H2") & (nom7 != "H1") & (nom7 != "H2") & (nom1 != "O1") & (nom1 != "O2") & (nom7 != "O1") & (nom7 != "O2"): 
            d17 = calcdst(liste_order,i1,i7,xat,yat,zat)
            outedge(vert1,vert7,d17,liste_order,i1,i7,edgeline)
          first = 0

       # generate all distances: CBi-Ni, CBi-Ci, CBi-CAi, CBi-HAi
       # generate all distances: Oi-Ci, Oi-CAi, Oi-Ni+1
       for ii in range(0,len(resname)-1):
           #print "index_vert",index_vert
           if resname[ii] != "GLY": 
            keyCB = "CB-"+resname[ii]+"-"+str(ii+1)
            idxCB = index_vert["CB-"+resname[ii]+"-"+str(ii+1)]
           keyCA = "CA-"+resname[ii]+"-"+str(ii+1)
           idxCA = index_vert["CA-"+resname[ii]+"-"+str(ii+1)]
           keyHA = "HA-"+resname[ii]+"-"+str(ii+1)
           idxHA = index_vert["HA-"+resname[ii]+"-"+str(ii+1)]
           keyC = "C-"+resname[ii]+"-"+str(ii+1)
           idxC = index_vert["C-"+resname[ii]+"-"+str(ii+1)]
           keyN = "N-"+resname[ii]+"-"+str(ii+1)
           idxN = index_vert["N-"+resname[ii]+"-"+str(ii+1)]
           keyO = "O-"+resname[ii]+"-"+str(ii+1)
           idxO = index_vert["O-"+resname[ii]+"-"+str(ii+1)]
           if ii < len(resname)-1: 
              residp1 = str(ii+2)
              resnamep1 = resname[ii+1]
              keyNp1 = "N-"+resnamep1+"-"+residp1
              idxNp1 = index_vert["N-"+resnamep1+"-"+residp1]
              keyHp1 = "H-"+resnamep1+"-"+residp1
              idxHp1 = index_vert["H-"+resnamep1+"-"+residp1]
           #
           #print "keyCB",keyCB,"keyCA",keyCA
           #print "idxCB",idxCB,"idxCA",idxCA
           #print "liste_order",liste_order
           if resname[ii] != "GLY": 
            dCBCA = calcdstCBO(keyCB,keyCA,xat,yat,zat)
            outedge2(idxCB,idxCA,dCBCA,keyCB,keyCA,edgeline) 
            #
            dCBHA = calcdstCBO(keyCB,keyHA,xat,yat,zat)
            outedge2(idxCB,idxHA,dCBHA,keyCB,keyHA,edgeline) 
            #
            dCBC = calcdstCBO(keyCB,keyC,xat,yat,zat)
            outedge2(idxCB,idxC,dCBC,keyCB,keyC,edgeline) 
            #
            dCBN = calcdstCBO(keyCB,keyN,xat,yat,zat)
            outedge2(idxCB,idxN,dCBC,keyCB,keyN,edgeline) 
           #
           dOC = calcdstCBO(keyO,keyC,xat,yat,zat)
           outedge2(idxO,idxC,dOC,keyO,keyC,edgeline) 
           #
           dOCA = calcdstCBO(keyO,keyCA,xat,yat,zat)
           outedge2(idxO,idxCA,dOCA,keyO,keyCA,edgeline) 
           if ii < len(resname)-1: 
              dONp1 = calcdstCBO(keyO,keyNp1,xat,yat,zat)
              outedge2(idxO,idxNp1,dONp1,keyO,keyNp1,edgeline) 
              dOHp1 = calcdstCBO(keyO,keyHp1,xat,yat,zat)
              outedge2(idxO,idxHp1,dOHp1,keyO,keyHp1,edgeline) 

       # generates the edge between the NH2 N terminal group and the other atoms of residue 1
       fileNTER = open("NTER.pdb","r")
       lines = fileNTER.readlines()
       fileNTER.close()

       xatNTER = {}
       yatNTER = {}
       zatNTER = {}
       for line in lines:
        if line[0:4] == "ATOM":
          nomato = line[12:16]
       	  #print "nomato avant",nomato,len(nomato)
       	  ato = line[12:16].strip()
       	  res = resname[0]
       	  nor = line[22:26].strip()
       	  #print "nomato apres",nomato,len(nomato)
       	  #ato.append(nomato)
       	  #res.append(line[17:20])
       	  #nor.append(int(line[22:26]))
       	  #print line
       	  #print "ato",ato,"res",res,"nor",nor
       	  if (res == "GLY") & (ato == "HA3"):
          	ato = "HA1"
       	  idx = "".join([ato,"-",res,"-",nor])
       	  xatNTER[idx] = float(line[30:38])
       	  yatNTER[idx] = float(line[38:46])
       	  zatNTER[idx] = float(line[46:54])

       keyH1 = "H1-"+resname[0]+"-1"
       idxH1 = index_vert["H1-"+resname[0]+"-1"]
       keyH2 = "H2-"+resname[0]+"-1"
       idxH2 = index_vert["H2-"+resname[0]+"-1"]
       keyCA = "CA-"+resname[0]+"-1"
       idxCA = index_vert["CA-"+resname[0]+"-1"]
       keyHA = "HA-"+resname[0]+"-1"
       idxHA = index_vert["HA-"+resname[0]+"-1"]
       keyC = "C-"+resname[0]+"-1"
       idxC = index_vert["C-"+resname[0]+"-1"]
       keyN = "N-"+resname[0]+"-1"
       idxN = index_vert["N-"+resname[0]+"-1"]
       keyO = "O-"+resname[0]+"-1"
       idxO = index_vert["O-"+resname[0]+"-1"]

       dNH1 = calcdstCBO(keyN,keyH1,xatNTER,yatNTER,zatNTER)
       outedge3(idxN,idxH1,dNH1,keyN,keyH1,edgeline) 
       dNH2 = calcdstCBO(keyN,keyH2,xatNTER,yatNTER,zatNTER)
       outedge3(idxN,idxH2,dNH2,keyN,keyH2,edgeline) 
       dH1H2 = calcdstCBO(keyH1,keyH2,xatNTER,yatNTER,zatNTER)
       outedge3(idxH1,idxH2,dH1H2,keyH1,keyH2,edgeline) 
       dH1CA = calcdstCBO(keyH1,keyCA,xatNTER,yatNTER,zatNTER)
       outedge3(idxH1,idxCA,dH1CA,keyH1,keyCA,edgeline)
       dH2CA = calcdstCBO(keyH2,keyCA,xatNTER,yatNTER,zatNTER)
       outedge3(idxH2,idxCA,dH2CA,keyH2,keyCA,edgeline)
       dH1HA = calcdstCBO(keyH1,keyHA,xatNTER,yatNTER,zatNTER)
       outedge3(idxH1,idxHA,dH1HA,keyH1,keyHA,edgeline)
       dH2HA = calcdstCBO(keyH2,keyHA,xatNTER,yatNTER,zatNTER)
       outedge3(idxH2,idxHA,dH2HA,keyH2,keyHA,edgeline)
       #
       #edgeline["18 19"] = "18\t19\tD\t1.68487809729\t#  H-1   H2-1"

       # generates the edge between the O2 atom and the other atoms in the last residue
       fileCTER = open("CTER.pdb","r")
       lines = fileCTER.readlines()
       fileCTER.close()

       xatCTER = {}
       yatCTER = {}
       zatCTER = {}
       for line in lines:
        if line[0:4] == "ATOM":
          nomato = line[12:16]
          #print "nomato avant",nomato,len(nomato)
          ato = line[12:16].strip()
          res = resname[len(resname)-1]
          nor = str(len(resname))
          #print "nomato apres",nomato,len(nomato)
          #ato.append(nomato)
          #res.append(line[17:20])
          #nor.append(int(line[22:26]))
          #print line
          #print "ato",ato,"res",res,"nor",nor
          if (res == "GLY") & (ato == "HA3"):
                ato = "HA1"
          idx = "".join([ato,"-",res,"-",nor])
          xatCTER[idx] = float(line[30:38])
          yatCTER[idx] = float(line[38:46])
          zatCTER[idx] = float(line[46:54])

       keyO2 = "O2-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxO2 = index_vert["O2-"+resname[len(resname)-1]+"-"+str(len(resname))]
       keyCA = "CA-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxCA = index_vert["CA-"+resname[len(resname)-1]+"-"+str(len(resname))]
       keyHA = "HA-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxHA = index_vert["HA-"+resname[len(resname)-1]+"-"+str(len(resname))]
       keyC = "C-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxC = index_vert["C-"+resname[len(resname)-1]+"-"+str(len(resname))]
       keyN = "N-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxN = index_vert["N-"+resname[len(resname)-1]+"-"+str(len(resname))]
       keyO1 = "O1-"+resname[len(resname)-1]+"-"+str(len(resname))
       idxO1 = index_vert["O1-"+resname[len(resname)-1]+"-"+str(len(resname))]

       #print "keyN",keyN,"keyO",keyO
       dNO1 = calcdstCBO(keyN,keyO1,xatCTER,yatCTER,zatCTER)
       outedge3(idxN,idxO1,dNO1,keyN,keyO1,edgeline)
       dNO2 = calcdstCBO(keyN,keyO2,xatCTER,yatCTER,zatCTER)
       outedge3(idxN,idxO2,dNO2,keyN,keyO2,edgeline)
       dO1O2 = calcdstCBO(keyO1,keyO2,xatCTER,yatCTER,zatCTER)
       outedge3(idxO1,idxO2,dO1O2,keyO1,keyO2,edgeline)
       dO1CA = calcdstCBO(keyO1,keyCA,xatCTER,yatCTER,zatCTER)
       outedge3(idxO1,idxCA,dO1CA,keyO1,keyCA,edgeline)
       #print "key0",keyO,"keyCA",keyCA,"idxO",idxO,"idxCA",idxCA,"dOCA",dOCA
       # dCCA: required for generate_dmdgp_i_i-6_ss5_CBO.py
       dCCA = calcdstCBO(keyC,keyCA,xatCTER,yatCTER,zatCTER)
       outedge3(idxC,idxCA,dCCA,keyC,keyCA,edgeline)
       dO2CA = calcdstCBO(keyO2,keyCA,xatCTER,yatCTER,zatCTER)
       outedge3(idxO2,idxCA,dO2CA,keyO2,keyCA,edgeline)
       dCO1 = calcdstCBO(keyC,keyO1,xatCTER,yatCTER,zatCTER)
       outedge3(idxC,idxO1,dCO1,keyC,keyO1,edgeline)
       dCO2 = calcdstCBO(keyC,keyO2,xatCTER,yatCTER,zatCTER)
       outedge3(idxC,idxO2,dCO2,keyC,keyO2,edgeline)
       dO1HA = calcdstCBO(keyO1,keyHA,xatCTER,yatCTER,zatCTER)
       outedge3(idxO1,idxHA,dO1HA,keyO1,keyHA,edgeline)
       dCHA = calcdstCBO(keyC,keyHA,xatCTER,yatCTER,zatCTER)
       outedge3(idxC,idxHA,dCHA,keyC,keyHA,edgeline)
       dO2HA = calcdstCBO(keyO2,keyHA,xatCTER,yatCTER,zatCTER)
       outedge3(idxO2,idxHA,dO2HA,keyO2,keyHA,edgeline)

       # generate the restraints for CB if not GLY
       if resname[len(resname)-1] != "GLY":
          keyCB = "CB-"+resname[len(resname)-1]+"-"+str(len(resname))
          idxCB = index_vert["CB-"+resname[len(resname)-1]+"-"+str(len(resname))]
          # generate all distances: CBi-Ni, CBi-Ci, CBi-CAi, CBi-HAi
          dCBCA = calcdstCBO(keyCB,keyCA,xatCTER,yatCTER,zatCTER)
          outedge3(idxCB,idxCA,dCBCA,keyCB,keyCA,edgeline)
          #
          dCBHA = calcdstCBO(keyCB,keyHA,xatCTER,yatCTER,zatCTER)
          outedge3(idxCB,idxHA,dCBHA,keyCB,keyHA,edgeline)
          #
          dCBC = calcdstCBO(keyCB,keyC,xatCTER,yatCTER,zatCTER)
          outedge3(idxCB,idxC,dCBC,keyCB,keyC,edgeline)
          #
          dCBN = calcdstCBO(keyCB,keyN,xatCTER,yatCTER,zatCTER)
          outedge3(idxCB,idxN,dCBC,keyCB,keyN,edgeline)

       file = open("tmp","w")
       for k, v in edgeline.iteritems():
           file.write("%s \n" % v)
       file.close()
       #print "coucou"
       p = os.popen('/bin/sort -n -k 1,1 -k 2,2 tmp',"r")
       #os.system("/bin/sort -n -k 1,1 -k 2,2 tmp")
       while 1:
          linep = p.readline()
          if not linep: break
          print linep[0:len(linep)-1]
       first = 0
       #os.system("/bin/sort -n -k 1,1 -k 2,2 tmp > tmp1")
       #os.system("/bin/ls -lt")
          
    if (edges == 0) & (dihed == 0):
    # not edge part: so, let it as it is     
    # but do not write any dihedral angle restraint
      print line[0:len(line)-1]
      #oo = 0

    if line[0:11] == "begin edges":
       edges = 1
       first = 1

    if line[0:21] == "begin dihedral_angles":
       dihed = 1


