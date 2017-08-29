from collections import defaultdict

from porestat.utils.DataFrame import DataFrame
from database.MIRFamily import MIRFamilyDB
from neo4j.v1 import GraphDatabase, basic_auth

from database.Neo4JInterface import neo4jInterface
from database.ORGMIRs import ORGMIRDB
from synonymes.mirnaID import miRNA, miRNAPART
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

    mimatORGMIR = orgmirs.mimat2orgmir.get(matureAcc, None)

    if mimatORGMIR != None:
        db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'org', ['MIRNA_ORG'], {'id': mimatORGMIR}, ['IS_ORG_MIRNA'], None)


# CREATE MIRNA FAMILY NODES
mirFamiliesByMI = {}
for family in familyDB:

    db.createNodeIfNotExists(['MIRNA_FAMILY'], {'id': family.familyID, 'name': family.familyName})

    for (miID, miName) in family.childMIMATs:
        mirFamiliesByMI[miID] = family

# CREATE ORG MIRNAS
orgmirs = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

for orgmir in orgmirs.orgmir2mimat:
    db.createNodeIfNotExists(['MIRNA_ORG'], {'id': orgmir})

    orgFamilies = orgmirs.orgmir2mi.get(orgmir, None)

    if orgFamilies != None:

        for mi in orgFamilies:
            orgFamily = mirFamiliesByMI[mi]

            if orgFamily != None:
                db.createRelationship('org', ['MIRNA_ORG'], {'id': orgmir}, 'mifam', ['MIRNA_FAMILY'], {'id': orgFamily}, ['BELONGS_TO'], None)

for row in mirbase:

    matureAcc1 = None
    matureID1 = None
    matureAcc2 = None
    matureID2 = None

    mirnaAccession = row['Accession']
    mirnaFamily = mirFamiliesByMI[mirnaAccession] if mirnaAccession in mirFamiliesByMI else None
    mirnaFamilyID = ""

    if mirnaFamily != None:
        mirnaFamilyID = mirnaFamily.familyID


    # CREATE PRE-MATURE MIRNA NODE
    db.createNodeIfNotExists(['MIRNA_PRE'], {'id': mirnaAccession})

    db.createRelationship('pre', ['MIRNA_PRE'], {'id': mirnaAccession}, 'org', ['MIRNA_ORG'], {'id': 'ORG' + mirnaAccession}, ['IS_ORG_MIRNA'], None)

    if mirnaFamilyID != None:
        db.createRelationship('pre', ['MIRNA_PRE'], {'id': mirnaAccession}, 'mifam', ['MIRNA_FAMILY'], {'id': mirnaFamilyID}, ['BELONGS_TO'], None)


    matureID1 = row['Mature1_ID']
    matureAcc1 = row['Mature1_Acc']

    matureID2 = row['Mature2_ID']
    matureAcc2 = row['Mature2_Acc']

    if not (matureAcc1 == None or matureAcc1 == 'None'):
        createMatureMIRNAEntry(matureAcc1, matureID1, mirnaAccession, mirnaFamilyID)

    if not (matureAcc2 == None or matureAcc2 == 'None'):
        createMatureMIRNAEntry(matureAcc2, matureID2, mirnaAccession, mirnaFamilyID)

db.close()