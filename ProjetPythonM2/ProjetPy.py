#! /usr/bin/env python3

import re
#Projet Python M2, Parsing du fichier pdb #Module Liaisons Hydrophobes

RESD_Hydro = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ATOM_Dnn = ["N", "NE", "NE2", "NH1", "NH2", "ND1", "ND2", "NZ", "OH", "OG1"]
ATOM_aCC = ["O", "OD", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"] 

#N est donneur de LH, excepté sur la proline qui n'en a pas


pdbname = sys.args[1]
pdbfile = open(pdbname, 'r')

#creer une classe des atomes avec les infos tirées du pdb
#and re.search(infoATOM[3], ATOM_Hydro): pour selectionner les hydrophobes à partir d'une liste

class Atom:
    """
    Ceci est la classe Atome
    """
    
    def __init__(self, Num=0, Nom="CC", Resid="AAA", Chain="A", NumRes=0, 
                    x=0.0, y=0.0, z=0.0):
        """
        definition des objts atome
        """
        self.Num=Num
        self.Nom=Nom
        self.Resid=Resid
        self.Chain=Chain
        self.NumRes=NumRes
        self.x=x
        self.y=y
        self.z=z


def Parse_and_Constructor(pdbfile):
    """Récupère les coord et le type d'atome, de base et renvoie un dico"""
    ListAtom = []
    for line in pdbfile:
        infoATOM = line.split()
        if infoATOM[0] == "ATOM": 
            atom=Atom(infoATOM[1], infoATOM[2], infoATOM[3], infoATOM[4], infoATOM[5], \
                        infoATOM[6], infoATOM[7], infoATOM[8])
            ListAtom.append(atom)

    return ListAtom

def Select_Hydro(ListAtom, RESD_hydro):
    """
    Sélectionne la liste des atomes succeptibles de créer des Liasons Hydrophobes
    """
    List_Hydro = []
    for atome in ListAtom:
        if re.search(self.Resid, RESD_Hydro):
            List_Hydro.append(atome)

    return List_Hydro

def Select_Disulfure(ListAtom):
    """
    Sélectionne la liste des atomes succeptibles de faire des ponts S-S
    """
    List_SS = []
    for atome in ListAtom:
        if self.Nom=="SG":
            List_SS.append(atome)

    return List_SS

def Select_Hbond(ListAtom, ATOM_hydro):
    """
    Sélectionne la liste des atomes donneurs et accepteurs de LH
    """
    List_Hbond = []
    for atome in ListAtom:
        if re.search(self.Nom, ATOM_Dnn) or re.search(self.Nom, ATOM_Acc):
            if self.Resid=="PRO" and self.Nom!="N":
            List_Hbond.append(atome)

    return List_Hydro

def Calcul_Dist:
    return
    
    
