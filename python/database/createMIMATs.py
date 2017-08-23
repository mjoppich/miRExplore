from porestat.utils.DataFrame import DataFrame
from database.MIRFamily import MIRFamilyDB
from neo4j.v1 import GraphDatabase, basic_auth

from database.Neo4JInterface import neo4jInterface
from utils.idutils import dataDir

mirbase = DataFrame.parseFromFile(dataDir + "/miRExplore/mirnas_mirbase.csv", bConvertTextToNumber=False)
filename = dataDir + "/miRExplore/miFam.dat"
familyDB = MIRFamilyDB(filename)

print(mirbase.getHeader())

db = neo4jInterface(simulate=True)

def createMatureMIRNAEntry(matureAcc, matureID, preMatureID, familyID):
    db.createNodeIfNotExists( ['MIRNA'], {'id': matureAcc, 'name': matureID, 'family': familyID} )
    db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'fam', ['MIRNA_FAMILY'], {'id': familyID}, ['BELONGS_TO'], None)
    db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'premat', ['MIRNA_PRE'], {'id': preMatureID}, ['MATURE_OF'], None)

mirFamiliesByMI = {}

# CREATE MIRNA FAMILY NODES
for family in familyDB:

    db.createNodeIfNotExists(['MIRNA_FAMILY'], {'id': family.familyID, 'name': family.familyName})

    for (miID, miName) in family.childMIMATs:
        mirFamiliesByMI[miID] = family

for row in mirbase:

    matureAcc1 = None
    matureID1 = None
    matureAcc2 = None
    matureID2 = None

    mirnaAccession = row['Accession']

    # CREATE PRE-MATURE MIRNA NODE
    db.createNodeIfNotExists(['MIRNA_PRE'], {'id': mirnaAccession})


    mirnaFamily = mirFamiliesByMI[mirnaAccession] if mirnaAccession in mirFamiliesByMI else None
    mirnaFamilyID = ""

    if mirnaFamily != None:
        mirnaFamilyID = mirnaFamily.familyID

    matureAcc1 = row['Mature1_Acc']
    matureID1 = row['Mature1_ID']

    matureID2 = row['Mature2_ID']
    matureAcc2 = row['Mature2_Acc']

    if not (matureAcc1 == None or matureAcc1 == 'None'):
        createMatureMIRNAEntry(matureAcc1, matureID1, mirnaAccession, mirnaFamilyID)

    if not (matureAcc2 == None or matureAcc2 == 'None'):
        createMatureMIRNAEntry(matureAcc2, matureID2, mirnaAccession, mirnaFamilyID)

db.close()