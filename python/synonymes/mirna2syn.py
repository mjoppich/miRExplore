import shlex
from collections import Counter, defaultdict
from porestat.utils.DataFrame import DataFrame

from database.MIRFamily import MIRFamilyDB
from synonymes.Synonym import Synonym
from synonymes.mirnaID import miRNA, miRNASynonymeTYPE
from utils.idutils import printToFile, dataDir

mirbase = DataFrame.parseFromFile(dataDir + "/miRExplore/mirnas_mirbase.csv", bConvertTextToNumber=False)
filename = dataDir + "/miRExplore/miFam.dat"
familyDB = MIRFamilyDB(filename)

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

    vOrgSyns = []

    for mi in MI2MIMAT:

        mimats = MI2MIMAT[mi]
        mirnas = [MIMAT2MIRNA[x] for x in mimats]

        orgSyn = Synonym('ORG'+mi)

        for mirna in mirnas:
            allIDs = set(mirna.make_strings( miRNA.compositions()[miRNASynonymeTYPE.MIORG] ))

            for syn in allIDs:
                orgSyn.addSyn(syn)

        vOrgSyns.append(orgSyn)

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
#synFiles['mirna_families.syn'] = makeFamilySynonymes()
#synFiles['mirna_org.syn'] = makeOrgSynonymes()
#synFiles['mirna_mi.syn'] = makeMISynonymes()
synFiles['mirna_mimat.syn'] = makeMIMATSynonymes()


for synFilename in synFiles:

    print(synFilename)

    vAllSyns = synFiles[synFilename]

    synCounter = Counter()
    for synonym in vAllSyns:

        for syn in synonym:
            synCounter[syn] += 1

    setCommonWords = set()

    for synWordCount in synCounter.most_common(66):

        if synCounter[synWordCount[0]] <= 2:
            continue

        setCommonWords.add(synWordCount[0])
        print(synWordCount[0] + " " + str(synWordCount[1]))

    vPrintSyns = []

    for synonym in vAllSyns:

        synonym.removeCommonSynonymes(setCommonWords)

        if len(synonym) > 0:
            vPrintSyns.append(synonym)

    printToFile(vPrintSyns, dataDir + "/miRExplore/textmine/synonyms/" + synFilename)
