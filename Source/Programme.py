#! /usr/bin/env python3
import sys
import Atome as A
import Hydrophobe as Hdp
import PontSS as Pss
import LiaisonsHydg as HB
import Ionique as Iq


RESD_Hydro = ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"]
ATOM_Dnn = ["N", "NE", "NE2", "NH1", "NH2", "ND1", "ND2", "NZ", "OH", "OG1"]
ATOM_Acc = ["O", "OD", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "SG"] 
RESD_IonPos = ["ARG", "LYS", "HIS"]
RESD_IonNeg = ["ASP", "GLU"]
RESD_Aromtq = ["HIS", "PHE", "PRO", "TRP", "TYR"]
ATOM_Cycle = ["ND1", "CE1", "NE2", "CD2", "CG", "N", "CA", "CB", "CD", "CD1" \
                ,"CZ", "CE2", "NE1", "CZ2", "CH2", "CZ3", "CE3"] 


pdbname = sys.argv[1]
pdbfile = open(pdbname, 'r')


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

ListAtom = A.Parse_and_Constructor(pdbfile)
List_Hydro = Hdp.Select_Hydro(ListAtom, RESD_Hydro)
Hydrophobe = Hdp.Calc_Int_Hydro(List_Hydro)

#AFFICHAGE DES INTÉRACTIONS HYDROPHOBES
dblon = []
#Liste qui permet d'éviter qu'une interaction n'apparaisse 2 fois
print("TYPE D'INTÉRACTION : Hydrophobe")
print("\n")
print("{:>6s} {:<6s} {:<6s} | {:>6s} {:<6s} {:<6s}".format("Numero", \
        "Residu", "Chaine", "Numero", "Résidu", "Chaine"))

for Liste in Hydrophobe:
    if [Liste[0].Num, Liste[1].Num] in dblon:
        continue 
    print("{:>6s} {:<6s} {:<6s} | {:>6s} {:<6s} {:<6s}".format(\
          Liste[0].NumRes, Liste[0].Resid, Liste[0].Chain, \
          Liste[1].NumRes, Liste[1].Resid, Liste[1].Chain ))

    dblon.append([Liste[0].Num, Liste[1].Num])
    dblon.append([Liste[1].Num, Liste[0].Num])

print("\n")

##########################################################

List_SS = Pss.Select_Disulfure(ListAtom)
PontSS = Pss.Calc_Pts_SS(List_SS)

#AFFICHAGE DES INTÉRACTIONS PONTS DISULFURES
dblon2 = []
print("TYPE D'INTÉRACTION : Ponts Disulfure S-S")
print("\n")
print("{:>6s} {:<6s} {:<6s} | {:>6s} {:<6s} {:<6s} | {:>8s}".format("Numero",\
          "Residu", "Chaine", \
          "Numero", "Residu", "Chaine", "Distance") )


for Liste in PontSS:
    if [Liste[0].Num, Liste[1].Num] in dblon2:
        continue 
    print("{:>6s} {:<6s} {:<6s} | {:>6s} {:<6s} {:<6s} | {:>8.2f}".format(Liste[0].NumRes,\
          Liste[0].Resid, Liste[0].Chain, \
          Liste[1].NumRes, Liste[1].Resid, Liste[1].Chain, Liste[2]) )

    dblon2.append([Liste[0].Num, Liste[1].Num])
    dblon2.append([Liste[1].Num, Liste[0].Num])

##########################################################

List_HbondD = HB.Select_HbondD(ListAtom, ATOM_Dnn)
List_HbondA = HB.Select_HbondA(ListAtom, ATOM_Acc)
LH = HB.Calcul_LH(List_HbondD, List_HbondA)

print("\n")
#AFFICHAGE DES LIAISONS HYDROGÈNES
dblon3 = []
print("TYPE D'INTÉRACTION : Liaisons hydrogènes")
print("\n")
print("{:>6s} {:<6s} {:<7s} {:<6s} | {:>6s} {:<6s} {:<7s} {:<6s} | {:>12s}".format("Numero",\
          "Residu", "atm dnn", "Chaine", \
          "Numero", "Residu", "atm acc", "Chaine", "Distance D-A"))


for Liste in LH:
    if [Liste[0].Num, Liste[1].Num] in dblon3:
        continue 
    print("{:>6s} {:<6s} {:<7s} {:<6s} | {:>6s} {:<6s} {:<7s} {:<6s} | {:>12.2f}".format(Liste[0].NumRes,\
          Liste[0].Resid, Liste[0].Nom, Liste[0].Chain, \
          Liste[1].NumRes, Liste[1].Resid, Liste[1].Nom, Liste[1].Chain, Liste[2]) )

    dblon3.append([Liste[0].Num, Liste[1].Num])
    dblon3.append([Liste[1].Num, Liste[0].Num])

###########################################################
List_ResIonP = Iq.Select_RESD_AromtqP(ListAtom, RESD_IonPos)
List_ResIonN = Iq.Select_RESD_AromtqN(ListAtom, RESD_IonNeg)
Ionq = Iq.Calcul_IntIon(List_ResIonP, List_ResIonN)

print("\n")
#AFFICHAGE DES INTERACTIONS IONIQUES
dblon4 = []
print("TYPE D'INTÉRACTION : Ionique")
print("\n")
print("{:>6s} {:<6s} {:<6s} | {:>6s} {:<6s} {:<6s} | {:>10s}".format("Numero",\
          "Residu", "Chaine", \
          "Numero", "Residu", "Chaine", "Distance"))


for Liste in Ionq:
    if [Liste[0].Num, Liste[1].Num] in dblon4:
        continue 
    print("{:>6s} {:<6s} {:<7s} {:<6s} | {:>6s} {:<6s} {:<7s} {:<6s} | {:>10.2f}".format(Liste[0].NumRes,\
          Liste[0].Resid, Liste[0].Chain, \
          Liste[1].NumRes, Liste[1].Resid, Liste[1].Chain, Liste[2]) )

    dblon4.append([Liste[0].Num, Liste[1].Num])
    dblon4.append([Liste[1].Num, Liste[0].Num])










    
    
