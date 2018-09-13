#! /usr/bin/env python3
import sys
import re

#Projet Python M2, Parsing du fichier pdb #Module Liaisons Hydrophobes

RESD_Hydro = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ATOM_Dnn = ["N", "NE", "NE2", "NH1", "NH2", "ND1", "ND2", "NZ", "OH", "OG1"]
ATOM_Acc = ["O", "OD", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"] 
RESD_IonPos = ["ARG", "LYS", "HIS"]
RESD_IonNeg = ["ASP", "GLU"]
RESD_Aromtq = ["HIS", "PHE", "PRO", "TRP", "TYR"]
ATOM_Cycle = ["ND1", "CE1", "NE2", "CD2", "CG", "N", "CA", "CB", "CD", "CD1" \
                ,"CZ", "CE2", "NE1", "CZ2", "CH2", "CZ3", "CE3"] 

#N est donneur de LH, excepté sur la proline qui n'en a pas
#Pour les LH, faire la distance entre Le N et O moins la distance entre N et H


pdbname = sys.argv[1]
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
        if re.search(atome.Resid, RESD_Hydro):
            List_Hydro.append(atome)

    return List_Hydro

def Select_Disulfure(ListAtom):
    """
    Sélectionne la liste des atomes succeptibles de faire des ponts S-S
    """
    List_SS = []
    for atome in ListAtom:
        if atome.Nom=="SG":
            List_SS.append(atome)

    return List_SS

def Select_Hbond(ListAtom, ATOM_Dnn, ATOM_Acc):
    """
    Sélectionne la liste des atomes donneurs et accepteurs de LH
    """
    List_Hbond = []
    for atome in ListAtom:
        if re.search(atome.Nom, ATOM_Dnn) or re.search(atome.Nom, ATOM_Acc):
            if atome.Resid=="PRO" and atome.Nom=="N":
                continue
            List_Hbond.append(atome)

    return List_Hbond


def Select_ResIRESD_AromtqRESD_Aromtqon(ListAtom, RESD_IonPos, RESD_IonNeg):
    """
    Sélectionne la liste des atomes associés à des residus chargés
    """
    List_ResIon = []
    for atome in ListAtom:
        if re.search(atome.Resid, RESD_IonPos) or re.search(atome.Resid, RESD_IonNeg):
            List_ResIon.append(atome)

    return List_ResIon


def Select_AtmCyc(ListAtom, RESD_Aromtq, ATOM_Cycle):
    """
    Sélectionne la liste des atomes présent dans les cycles des aa aromatiques
    """
    List_AtmCyc = []
    for atome in ListAtom:
        if re.search(atome.Resid, RESD_Aromtq) and re.search(atome.Nom, ATOM_Cycle):
            List_AtmCyc.append(atome)

    return List_AtmCyc






#ListAtom = Parse_and_Constructor(pdbfile)
#for atome in ListAtom:
#    print(atome.Resid)


    
    
