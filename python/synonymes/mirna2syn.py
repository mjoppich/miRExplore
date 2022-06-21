import shlex
from collections import Counter, defaultdict

from utils.DataFrame import DataFrame
from database.MIRFamily import MIRFamilyDB
from database.ORGMIRs import ORGMIRDB
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from synonymes.mirnaID import miRNA, miRNASynonymeTYPE, miRNAPART
from utils.idutils import printToFile, dataDir, loadExludeWords



if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-m', '--mirbase', type=argparse.FileType("r"), required=True, help="input mirbase file")
    parser.add_argument('-s', '--syn', type=argparse.FileType("w"), required=True, help="output synonym file")
    args = parser.parse_args()


    mirbase = DataFrame.parseFromFile(dataDir + "/miRExplore/mirnas_mirbase.csv", bConvertTextToNumber=False)
    filename = dataDir + "/miRExplore/miFam.dat"
    familyDB = MIRFamilyDB(filename)

    orgmirs = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

    MIMAT2MIRNA = {}
    MI2MIMAT = defaultdict(set)


    for row in mirbase:

        matureAcc1 = None
        matureID1 = None
        matureAcc2 = None
        matureID2 = None

        MIid = mirnaAccession = row['Accession']


        matureAcc1 = row['Mature1_Acc']
        matureID1 = row['Mature1_ID']

        matureID2 = row['Mature2_ID']
        matureAcc2 = row['Mature2_Acc']

        if not (matureAcc1 == None or matureAcc1 == 'None'):
            MIMAT2MIRNA[matureAcc1] = miRNA(matureID1)
            MI2MIMAT[MIid].add(matureAcc1)

        if not (matureAcc2 == None or matureAcc2 == 'None'):
            MIMAT2MIRNA[matureAcc2] = miRNA(matureID2)
            MI2MIMAT[MIid].add(matureAcc2)


    def makeFamilySynonymes():

        vFamSyns = []

        for family in familyDB:

            for (miID, miName) in family.childMIMATs:

                mirnas = set()

                for mimat in MI2MIMAT[miID]:
                    mirnas.add( MIMAT2MIRNA[mimat] )


            familySyn = Synonym(family.familyID)

            for mirna in mirnas:

                allIDs = set(mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.FAMILY] ))
                for syn in allIDs:
                    familySyn.addSyn(syn)

            vFamSyns.append(familySyn)

        return vFamSyns


    def makeOrgSynonymes():

        dOrgSyns = {}
        addedMirnaIDs = list()

        for mimat in MIMAT2MIRNA:

            mirna = MIMAT2MIRNA[mimat]

            orgmir = orgmirs.mimat2orgmir.get(mimat, None)

            if orgmir == None:
                continue

            synid = orgmir

            if synid in dOrgSyns:
                orgSyn = dOrgSyns[synid]
            else:
                orgSyn = Synonym( synid )

            allsyns = set(mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.MIORG] ))

            for syn in allsyns:
                orgSyn.addSyn(syn)

            dOrgSyns[synid] = orgSyn


        vOrgSyns = [dOrgSyns[x] for x in dOrgSyns]

        return vOrgSyns


    def makeMISynonymes():
        vMISyns = []

        for mi in MI2MIMAT:

            mimats = MI2MIMAT[mi]
            mirnas = [MIMAT2MIRNA[x] for x in mimats]

            orgSyn = Synonym('ORG'+mi)

            for mirna in mirnas:
                allIDs = set(mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.MI] ))

                for syn in allIDs:
                    orgSyn.addSyn(syn)

            vMISyns.append(orgSyn)

        return vMISyns

    def makeMIMATSynonymes():

        mimatSyns = []

        for mimat in MIMAT2MIRNA:

            mirna = MIMAT2MIRNA[mimat]

            mirnaSyn = Synonym(mimat)

            for syn in mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.MIMAT]):
                mirnaSyn.addSyn(syn)

            mimatSyns.append(mirnaSyn)

        return mimatSyns



    synFiles = {}
    synFiles['mirna_families.syn'] = makeFamilySynonymes()
    synFiles['mirna_org.syn'] = makeOrgSynonymes()
    synFiles['mirna_mi.syn'] = makeMISynonymes()
    synFiles['mirna_mimat.syn'] = makeMIMATSynonymes()




    globalKeywordExcludes = loadExludeWords()

    for synFilename in synFiles:
        print(synFilename)

        singleMature = ['mir', 'miR', 'MiR', 'Mir' ,'microRNA', 'MicroRNA', 'micro-RNA', 'Micro-RNA', 'miRNA', 'MiRNA']
        pluralMature = [x+'s' for x in singleMature]
        singleMature.append('MiRNAS')
        singleMature.append('miRNAs')

        allupper = [x.upper() for x in singleMature+pluralMature]

        for x in singleMature+pluralMature+allupper:
            globalKeywordExcludes['generic'].add(x)

        vAllSyns = synFiles[synFilename]

        vPrintSyns = handleCommonExcludeWords(vAllSyns, globalKeywordExcludes, mostCommonCount=66, maxCommonCount=21)

        printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/" + synFilename)
