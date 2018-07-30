from collections import defaultdict

from porestat.utils.DataFrame import DataFrame
from database.MIRFamily import MIRFamilyDB
from neo4j.v1 import GraphDatabase, basic_auth

from database.Neo4JInterface import neo4jInterface
from synonymes.mirnaID import miRNA, miRNAPART
from utils.idutils import dataDir, printToFile

mirbase = DataFrame.parseFromFile(dataDir + "/miRExplore/mirnas_mirbase.csv", bConvertTextToNumber=False)
filename = dataDir + "/miRExplore/miFam.dat"
familyDB = MIRFamilyDB(filename)

print(mirbase.getHeader())

def getOrgMIRNAID(mirnaid, listedIDs, prefix, getID):
    mirna = miRNA(matureID1)

    idPart = getID(mirna)

    if idPart == None:
        print("no ID! " + matureID1)
        return None

    idPart = idPart.upper()

    idpartIdx = -1
    if idPart in listedIDs:
        idpartIdx = listedIDs.index(idPart)
    else:
        idpartIdx = len(listedIDs)
        listedIDs.append(idPart)

    synid = prefix + str(idpartIdx)

    return (synid,listedIDs)


addedOrgMIRNAIDs = list()
addedOrgMIIDs = list()

orgMIRNAID2MIMAT = defaultdict(set)
mimat2MI = {}
mimat2ORGMI = {}

for row in mirbase:

    matureAcc1 = None
    matureID1 = None
    matureAcc2 = None
    matureID2 = None


    mi = row['Accession']
    rowID = row['ID']

    matureID1 = row['Mature1_ID']
    matureAcc1 = row['Mature1_Acc']

    matureID2 = row['Mature2_ID']
    matureAcc2 = row['Mature2_Acc']

    (orgMIR, addedOrgMIRNAIDs) = getOrgMIRNAID(rowID, addedOrgMIRNAIDs, 'ORGMIR',
                                               lambda mirna: mirna.getPart(miRNAPART.ID))
    (orgMI, addedOrgMIIDs) = getOrgMIRNAID(rowID, addedOrgMIIDs, 'ORGMI',
                                           lambda mirna: mirna.getPart(miRNAPART.ID) + mirna.getPart(miRNAPART.PRECURSOR, ""))

    if not (matureID1 == None or matureID1 == 'None'):
        orgMIRNAID2MIMAT[orgMIR].add(matureAcc1)
        mimat2MI[matureAcc1] = mi
        mimat2ORGMI[matureAcc1] = orgMI

    if not (matureID2 == None or matureID2 == 'None'):
        orgMIRNAID2MIMAT[orgMIR].add(matureAcc2)
        mimat2MI[matureAcc2] = mi
        mimat2ORGMI[matureAcc2] = orgMI



allCreatedIDs = ['ORGMIR\tMIMAT\tMI\tORGMI']
for orgid in orgMIRNAID2MIMAT:

    linestr = str(orgid) + "\t"

    for mimat in orgMIRNAID2MIMAT[orgid]:
        mi = mimat2MI[mimat]
        orgmi = mimat2ORGMI[mimat]
        allCreatedIDs.append( linestr + str(mimat) + "\t" + str(mi) + "\t" + str(orgmi))


printToFile(allCreatedIDs, dataDir + "/miRExplore/orgmir.tsv")