#!/usr/bin/env python3

#AUTEURS: Felix VANDERMEEREN et Yoann Missolo
#Projet de modelisation moleculaire

#Script de calcul des angles et distances au cours d'une trajectoire de DM d'un 
#acide nucleique, par parsing des fichiers pdb associes 

#A-T: 
# * N1 ||| H-N3
# * N6-H6_1 et _2 ||| O4

#G-C:
# * N2-H2_1 et _2 ||| O2
# * N1-H ||| N3
# * O6 ||| H4_1 et _2 -N4




import math as m
import numpy as np
import time

start_time = time.time()



def parse(line): #OK
    """Récupère les coord et le type d'atome, de base et renvoie un dico"""
    
    idAtom = line[12:16]
    base = line[17:20]
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    
    return {'base':base, 'idAtom':idAtom, 'x':x, 'y':y, 'z':z}
    

def getAtom(PDBobj): #OK
    """Repère les atomes d'interet pour les lH (2 pour dist, 3 pour les angles 
    dans un PDB d'ADN et renvoie une liste de dico """
    
    listAtoms_dist, listAtoms_angle = [], []
    
    for line in PDBobj:
        atom = parse(line)
        
        if atom['base'] == "ADE":
            neededAtoms = ["N1", "N6"]
            
            if atom['idAtom'].strip() == "H61" or atom['idAtom'].strip() == "H62":
                listAtoms_angle.append(atom) #Needed pour calul angle
   
            if atom['idAtom'].strip() == "N6":            
                listAtoms_dist.append(atom) #Doublon pour calul dist
                 
                 
        elif atom['base'] == "THY": 
            neededAtoms = ["O4", "N3"]
            
            if atom['idAtom'].strip() == "O4":
                listAtoms_dist.append(atom) #Doublon pour calul dist
                
            if atom['idAtom'].strip() == "H3":
                listAtoms_angle.append(atom) #Needed pour calul angle            
            
            
        elif atom['base'] == "GUA":
            neededAtoms = ["O6", "N1", "N2"]
            
            if atom['idAtom'].strip() == "O6" or atom['idAtom'].strip() == "N2":
                listAtoms_dist.append(atom) #Doublon pour calul dist
                
            if atom['idAtom'].strip() == "H21" \
                                            or atom['idAtom'].strip() == "H22" \
                                            or atom['idAtom'].strip() == "H1":    
                listAtoms_angle.append(atom) #Needed pour calul angle
              

        elif atom['base'] == "CYT":
            neededAtoms = ["O2", "N3", "N4"]
            if atom['idAtom'].strip() == "O2" or atom['idAtom'].strip() == "N4":
                listAtoms_dist.append(atom) #Doublon pour calul dist                
                
            if atom['idAtom'].strip() == "H41" or atom['idAtom'].strip() == "H42":
                listAtoms_angle.append(atom) #Needed pour calul angle


                
        if atom['idAtom'].strip() in neededAtoms:
            listAtoms_dist.append(atom)    
            listAtoms_angle.append(atom)
                    

    return listAtoms_dist, listAtoms_angle
    
    

def calcDist(atom1, atom2):
    """Calcule la dist (en angstrom) entre D et A d'une lH"""
    x1, x2 = atom1['x'], atom2['x']
    y1, y2 = atom1['y'], atom2['y']
    z1, z2 = atom1['z'], atom2['z']
    
    return round(((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5, 2)
   
   
    
def super_calc_dist(slice1, slice2, saut):
    """Calcule toutes les dist (en angstrom) entre les D et A de 2 slices d'une 
    liste de dico d'atomes d'interet"""

    listRes = [0] * saut

    for i in range(saut):
        atom1, atom2 = slice1[i], slice2[i]
        
        listRes[i] = calcDist(atom1, atom2)

    return listRes



def calcAngle(donneur, hydrog, accept):
    """Calcule l'angle (en degres) formé par les 3 atomes d'une lH"""
    
    x_donneur, x_hydrog, x_accept = donneur['x'], hydrog['x'], accept['x']
    y_donneur, y_hydrog, y_accept = donneur['y'], hydrog['y'], accept['y']
    z_donneur, z_hydrog, z_accept = donneur['z'], hydrog['z'], accept['z']
    
    u = [x_donneur-x_hydrog, y_donneur-y_hydrog, z_donneur-z_hydrog]
    v = [x_accept-x_hydrog, y_accept-y_hydrog, z_accept-z_hydrog]
    
    prod_scal = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    norm_u = calcDist(donneur, accept)
    norm_v = calcDist(hydrog, accept)

    if -1 <= prod_scal/(norm_u * norm_v) <= 1: 
        return round(180/m.pi * m.acos(prod_scal/(norm_u * norm_v)))
    else:
        return round(prod_scal/(norm_u * norm_v) *500)
    
       
   
def super_calcAngle(listeAtoms, base):
    """Calcule tous les angles (en degres) entre les D, H et A d'une fusion de 2
    slices d'une liste de dico d'atomes d'interet""" 
    
    if base == "ADE": #Si c'est A--T
        matAngles = np.zeros((3, 3), dtype=object)
        
        matAngles[0, 0:2] = ["A--T_H3",
                        calcAngle(listeAtoms[0], listeAtoms[5], listeAtoms[4])]
        
        matAngles[1, 0:2] = ["A--T_H61",
                        calcAngle(listeAtoms[1], listeAtoms[2], listeAtoms[6])]
        matAngles[2, 0:2] = ["A--T_H62",
                        calcAngle(listeAtoms[1], listeAtoms[3], listeAtoms[6])]
        
        
    elif base == "THY": #Si c'est T--A
        matAngles = np.zeros((3, 3), dtype=object)
        
        matAngles[0, 0:2] = ["T--A_H3", 
                        calcAngle(listeAtoms[0], listeAtoms[1], listeAtoms[3])]
        
        matAngles[1, 0:2] = ["T--A_H61",
                        calcAngle(listeAtoms[2], listeAtoms[5], listeAtoms[4])]
        matAngles[2, 0:2] = ["T--A_H62",
                        calcAngle(listeAtoms[2], listeAtoms[6], listeAtoms[4])]


    elif base == "GUA": #Si c'est G--C
        matAngles = np.zeros((5, 3), dtype=object)
        
        matAngles[0, 0:2] = ["G--C_H21",
                        calcAngle(listeAtoms[0], listeAtoms[1], listeAtoms[6])]
        matAngles[1, 0:2] = ["G--C_H22",
                        calcAngle(listeAtoms[0], listeAtoms[2], listeAtoms[6])]

        matAngles[2, 0:2] = ["G--C_H1",
                        calcAngle(listeAtoms[3], listeAtoms[4], listeAtoms[7])]

        matAngles[3, 0:2] = ["G--C_H41",
                        calcAngle(listeAtoms[5], listeAtoms[9], listeAtoms[8])]
        matAngles[4, 0:2] = ["G--C_H42",
                        calcAngle(listeAtoms[5], listeAtoms[10], listeAtoms[8])]


    elif base == "CYT": #Si c'est C--G
        matAngles = np.zeros((5, 3), dtype=object)
        
        matAngles[0, 0:2] = ["C--G_H21",
                        calcAngle(listeAtoms[0], listeAtoms[6], listeAtoms[5])]
        matAngles[1, 0:2] = ["C--G_H22",
                        calcAngle(listeAtoms[0], listeAtoms[7], listeAtoms[5])]
        
        matAngles[2, 0:2] = ["C--G_H1",
                        calcAngle(listeAtoms[1], listeAtoms[9], listeAtoms[8])]

        matAngles[3, 0:2] = ["C--G_H41",
                        calcAngle(listeAtoms[2], listeAtoms[3], listeAtoms[10])]
        matAngles[4, 0:2] = ["C--G_H42",
                        calcAngle(listeAtoms[2], listeAtoms[4], listeAtoms[10])]
        
        
    return matAngles


    
def arrayToFile(M, fichier, j):
    """Permet d'écrire un matrice dans un fichier texte
    (parce que les fonctions Numpy font pas ce qu'on voulait)"""
        
    for line in M: #On parcourt chaque ligne de la mat
        ligne = "{}{:<8s}    {:>5s}    {:>4s}".format(str(j), line[0],
                                                    str(line[1]), str(line[2]))
                                                    
        fichier.write(ligne + '\n')

  

########
##MAIN##
########

nbFichiers = 2

for i in range(nbFichiers):

    #Fichier d'entree:
    pdb = open("./pdbdna/" + "frame_"+ str(i+1) +"_dna.pdb")
    listAtoms_dist, listAtoms_angle = getAtom(pdb)
    nbAtoms_angle, nbAtoms_dist = len(listAtoms_angle), len(listAtoms_dist)

    #Fichier de sortie:
    out = open("./outS/" + "mat" + str(i+1) +".dat", "w")
    out.write("{:^8s}    {:5s}    {:4s}".format("lH", "angle"+str(i+1),
                                                    "dist"+str(i+1)) + "\n")


    #Parcours des 2 listes en meme temps, pour calcul dist et angles:
    h_angle, b_angle = 0, nbAtoms_angle
    h_dist, b_dist = 0, nbAtoms_dist

    matVal = np.zeros((int(nbAtoms_angle/2), 3), dtype=object)

    saut_T = 3
    saut_A = 4
    saut_C = 5
    saut_G = 6


    i, j = 0, 1
    while h_angle < b_angle or h_dist < b_dist:
        atom = listAtoms_angle[h_angle]
        
        if atom['base'] == "ADE":
            saut_h = saut_A
            saut_b = saut_T
            nbVal = 3
               
               
        elif atom['base'] == "THY":
            saut_h = saut_T
            saut_b = saut_A
            nbVal = 3
            
            
        elif atom['base'] == "GUA":
            saut_h = saut_G
            saut_b = saut_C
            nbVal = 5
        
        
        elif atom['base'] == "CYT":
            saut_h = saut_C
            saut_b = saut_G
            nbVal = 5
     
     
        slice_haut_angle = listAtoms_angle[h_angle : h_angle+saut_h]
        slice_bas_angle = listAtoms_angle[b_angle-saut_b : b_angle]    
        matAngles = super_calcAngle(slice_haut_angle+slice_bas_angle, 
                                                                atom['base'])
        matVal[i : i+nbVal, ] = matAngles
        
        slice_haut_dist = listAtoms_dist[h_dist : h_dist+nbVal]
        slice_bas_dist = listAtoms_dist[b_dist-nbVal : b_dist]
        matVal[i : i+nbVal, 2] = super_calc_dist(slice_haut_dist, slice_bas_dist, nbVal)
        
        arrayToFile(matVal[i : i+nbVal], out, j)
        
        i += nbVal
        j += 1
        h_angle += saut_h
        b_angle -= saut_b
        h_dist += nbVal
        b_dist -= nbVal


    out.close()
    pdb.close()


print("--- %s seconds ---" % (time.time() - start_time))
