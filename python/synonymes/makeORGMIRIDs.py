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

def getOrgMIRNAID(mirnaid):
    mirna = miRNA(matureID1)

    idPart = mirna.getPart(miRNAPART.ID)

    if idPart == None:
        print("no ID! " + matureID1)
        return None

    idPart = idPart.upper()

    idpartIdx = -1
    if idPart in addedOrgMIRNAIDs:
        idpartIdx = addedOrgMIRNAIDs.index(idPart)
    else:
        idpartIdx = len(addedOrgMIRNAIDs)
        addedOrgMIRNAIDs.append(idPart)

    synid = 'ORGMIR' + str(idpartIdx)

    return synid


addedOrgMIRNAIDs = list()
orgMIRNAID2MIMAT = defaultdict(set)
mimat2MI = {}

for row in mirbase:

    matureAcc1 = None
    matureID1 = None
    matureAcc2 = None
    matureID2 = None


    mi = row['Accession']
    matureID1 = row['Mature1_ID']
    matureAcc1 = row['Mature1_Acc']

    matureID2 = row['Mature2_ID']
    matureAcc2 = row['Mature2_Acc']

    if not (matureID1 == None or matureID1 == 'None'):
        orgID = getOrgMIRNAID(matureID1)
        orgMIRNAID2MIMAT[orgID].add(matureAcc1)
        mimat2MI[matureAcc1] = mi

    if not (matureID2 == None or matureID2 == 'None'):
        orgID = getOrgMIRNAID(matureID2)
        orgMIRNAID2MIMAT[orgID].add(matureAcc2)
        mimat2MI[matureAcc2] = mi




allCreatedIDs = ['ORGMIR\tMIMAT\tMI']
for orgid in orgMIRNAID2MIMAT:

    linestr = str(orgid) + "\t"

    for mimat in orgMIRNAID2MIMAT[orgid]:
        mi = mimat2MI[mimat]
        allCreatedIDs.append( linestr + str(mimat) + "\t" + str(mi))


printToFile(allCreatedIDs, dataDir + "/miRExplore/orgmir.tsv")