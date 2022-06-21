import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import defaultdict
from synonymes.GeneOntology import GeneOntology

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID
from synonymes.mirnaID import miRNA, miRNASynonymeTYPE, miRNAPART


if __name__ == '__main__':

    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-x', '--mirna-xls', type=argparse.FileType("r"), required=True, help="input ontology file")
    parser.add_argument('-a', '--mirna-alias', type=argparse.FileType("r"), required=True, help="input ontology file")
    parser.add_argument('-s', '--syn', type=argparse.FileType("w"), required=True, help="output synonym file")
    args = parser.parse_args()

    ent2syns = {}
    for line in args.mirna_alias:
        line = line.strip().split("\t")

        entity = line[0]
        syns = line[1].split(";")
        ent2syns[entity] = syns

    import pandas as pd
    allMirnas = pd.read_excel(args.mirna_xls.name)

    print(allMirnas.columns)
    vAllSyns = {}
    for ri, row in allMirnas.iterrows():
        
        accID = row["Accession"]
        mirnaNamePre = row["ID"]

        if not mirnaNamePre.startswith(("hsa", "mmu")):
            continue

        accMature1ID = row["Mature1_Acc"]
        accMature2ID = row["Mature2_Acc"]

        mirnaNameMat1 = row["Mature1_ID"]
        mirnaNameMat2 = row["Mature2_ID"]

        #print(mirnaNameMat1, pd.isna(mirnaNameMat1))
        #print(mirnaNameMat2, pd.isna(mirnaNameMat2))

        mirPre = miRNA(mirnaNamePre)
        mirMat1 = miRNA(mirnaNameMat1)

        mirMat2 = None
        if not pd.isna(mirnaNameMat2):
            mirMat2 = miRNA(mirnaNameMat2)


        entryName = "-".join(mirnaNamePre.split("-")[1:3]).replace("mir", "miR")
        familySyn = vAllSyns.get(entryName, Synonym(entryName))

        allMirs = [x for x in [mirPre, mirMat1, mirMat2] if not x is None]
        for mirna in allMirs:

            allIDs = set(mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.MIALL] ))
            for syn in allIDs:
                familySyn.addSyn(syn)

        if accID in ent2syns:
            for y in ent2syns[accID]:
                familySyn.addSyn(y)
        if accMature1ID in ent2syns:
            for y in ent2syns[accMature1ID]:
                familySyn.addSyn(y)
        if accMature2ID in ent2syns:
            for y in ent2syns[accMature2ID]:
                familySyn.addSyn(y)

        vAllSyns[entryName] = familySyn


    globalKeywordExcludes = loadExludeWords()
    vPrintSyns = handleCommonExcludeWords([vAllSyns[x] for x in vAllSyns], None, mostCommonCount=None, maxCommonCount=1)
    printToFile(vPrintSyns, args.syn.name, codec='utf8')