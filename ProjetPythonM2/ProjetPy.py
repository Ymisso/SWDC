#! /usr/bin/env python3
import sys
import Atome as A


RESD_Hydro = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ATOM_Dnn = ["N", "NE", "NE2", "NH1", "NH2", "ND1", "ND2", "NZ", "OH", "OG1"]
ATOM_Acc = ["O", "OD", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH"] 
RESD_IonPos = ["ARG", "LYS", "HIS"]
RESD_IonNeg = ["ASP", "GLU"]
RESD_Aromtq = ["HIS", "PHE", "PRO", "TRP", "TYR"]
ATOM_Cycle = ["ND1", "CE1", "NE2", "CD2", "CG", "N", "CA", "CB", "CD", "CD1" \
                ,"CZ", "CE2", "NE1", "CZ2", "CH2", "CZ3", "CE3"] 

#Pour les LH, faire la distance entre Le N et O moins la distance entre N et H

pdbname = sys.argv[1]
pdbfile = open(pdbname, 'r')

def Select_Hydro(ListAtom, RESD_Hydro):
    """
    Sélectionne la liste des atomes succeptibles de créer des Liasons Hydrophobes
    """
    List_Hydro = []
    for atome in ListAtom:
        if atome.Resid in RESD_Hydro and atome.Nom=="CA":
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
        if atome.Nom in ATOM_Dnn or atome.Nom in ATOM_Acc:
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
        if atome.Resid in RESD_IonPos or atome.Resid in RESD_IonNeg:
            List_ResIon.append(atome)

    return List_ResIon


def Select_AtmCyc(ListAtom, RESD_Aromtq, ATOM_Cycle):
    """
    Sélectionne la liste des atomes présent dans les cycles des aa aromatiques
    """
    List_AtmCyc = []
    for atome in ListAtom:
        if atome.Resid in RESD_Aromtq and atome.Nom in ATOM_Cycle:
            List_AtmCyc.append(atome)

    return List_AtmCyc

###########################################################
# Fonction de calcul


#Hydrophobe
def Calc_Int_Hydro(List_Hydro):
    """
    Calcule les distances entre les carbones Apha des résidus hydrophobes
    """
    Hydrop = []
    for at1 in List_Hydro:
        for at2 in List_Hydro:
            if round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 5.3 \
                and round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 \
                 + (float(at1.z)-float(at2.z))**2)**0.5, 2) > 0:
                Hydrop.append([at1, at2])
    return Hydrop




ListAtom = A.Parse_and_Constructor(pdbfile)
List_Hydro = Select_Hydro(ListAtom, RESD_Hydro)
Hydrophobe = Calc_Int_Hydro(List_Hydro)

#dblon = []
#for Liste in Hydrophobe:
#    if [Liste[0].Num, Liste[1].Num] in dblon:
#        continue 
#    print(Liste[0].Num, Liste[0].Resid, Liste[0].Chain, \
#          Liste[1].Num, Liste[1].Resid, Liste[1].Chain )

#    dblon.append([Liste[0].Num, Liste[1].Num])
#    dblon.append([Liste[1].Num, Liste[0].Num])

##########################################################
#Ponts Disulfures

def Calc_Pts_SS(List_SS):
    """
    Calcule les distances entre les carbones Apha des résidus hydrophobes
    """
    PontSS = []
    for at1 in List_SS:
        for at2 in List_SS:
            if round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2) < 2.2 \
                and round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 \
                 + (float(at1.z)-float(at2.z))**2)**0.5, 2) > 0:

                dist = round(((float(at1.x)-float(at2.x))**2 + (float(at1.y)-float(at2.y))**2 + \
                (float(at1.z)-float(at2.z))**2)**0.5, 2)
                PontSS.append([at1, at2, dist])
    return PontSS

List_SS = Select_Disulfure(ListAtom)
PontSS = Calc_Pts_SS(List_SS)

dblon = []
for Liste in PontSS:
    if [Liste[0].Num, Liste[1].Num] in dblon:
        continue 
    print(Liste[0].Num, Liste[0].Resid, Liste[0].Chain, \
          Liste[1].Num, Liste[1].Resid, Liste[1].Chain, Liste[2] )

    dblon.append([Liste[0].Num, Liste[1].Num])
    dblon.append([Liste[1].Num, Liste[0].Num])









    
    
