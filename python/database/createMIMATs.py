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

db = neo4jInterface(simulate=False, printQueries=False)

db.deleteRelationship('n', None, None, 'm', None, None, ['IS_ORG_MI'], None, 'r')
db.deleteRelationship('n', None, None, 'm', None, None, ['IS_ORG_MIR'], None, 'r')
db.deleteRelationship('n', None, None, 'm', None, None, ['IS_ORG_MIRNA'], None, 'r')

db.deleteRelationship('n', ['MIRNA'], None, 'm', ['MIRNA_PRE'], None, ['MIRNA_MATURE_OF'], None, 'r')

db.deleteRelationship('n', ['MIRNA'], None, 'm', ['MIRNA_FAMILY'], None, ['MIRNA_BELONGS_TO'], None, 'r')
db.deleteRelationship('n', ['MIRNA_PRE'], None, 'm', ['MIRNA_FAMILY'], None, ['MIRNA_BELONGS_TO'], None, 'r')
db.deleteRelationship('n', ['MIRNA_ORGMIR'], None, 'm', ['MIRNA_FAMILY'], None, ['MIRNA_BELONGS_TO'], None, 'r')
db.deleteRelationship('n', ['MIRNA_ORGMI'], None, 'm', ['MIRNA_FAMILY'], None, ['MIRNA_BELONGS_TO'], None, 'r')
db.deleteRelationship('n', ['MIRNA_ORGMI'], None, 'm', ['MIRNA_PRE'], None, ['MIRNA_BELONGS_TO'], None, 'r')

print("Relations deleted")

db.deleteNode(["IS_A_MIRNA"], None)

print("Nodes deleted")

db.createUniquenessConstraint('MIRNA', 'id')
db.createUniquenessConstraint('MIRNA_FAMILY', 'id')
db.createUniquenessConstraint('MIRNA_PRE', 'id')
db.createUniquenessConstraint('MIRNA_ORGMIR', 'id')
db.createUniquenessConstraint('MIRNA_ORGMI', 'id')

print("DB prepared")

def createMatureMIRNAEntry(matureAcc, matureID, preMatureID, familyID):

    db.createNodeIfNotExists( ['MIRNA','MIRNAS'], {'id': matureAcc, 'name': matureID, 'family': familyID}, propsCheck={'id': matureAcc})
    db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'premat', ['MIRNA_PRE'], {'id': preMatureID}, ['MIRNA_MATURE_OF'], None)

    if familyID != None:
        db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'fam', ['MIRNA_FAMILY'], {'id': familyID}, ['MIRNA_BELONGS_TO'], None)

    mimatORGMIR = orgmirDB.mimat2orgmir.get(matureAcc, None)
    if mimatORGMIR != None:
        db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'org', ['MIRNA_ORGMIR'], {'id': mimatORGMIR}, ['IS_ORG_MIR'], None)

    mimatORGMI = orgmirDB.mimat2orgmi.get(matureAcc, None)
    if mimatORGMI != None:
        db.createRelationship('mi', ['MIRNA'], {'id': matureAcc}, 'org', ['MIRNA_ORGMI'], {'id': mimatORGMI}, ['IS_ORG_MI'], None)

print("Starting MIRNA_FAMILY")
# CREATE MIRNA FAMILY NODES
mirFamiliesByMI = {}
for family in familyDB:

    db.createNodeIfNotExists(['MIRNA_FAMILY', 'MIRNAS'], {'id': family.familyID, 'name': family.familyName})

    for (miID, miName) in family.childMIMATs:
        mirFamiliesByMI[miID] = family

print("FINISHED MIRNA_FAMILY")
# CREATE ORG MIRNAS
orgmirDB = ORGMIRDB(dataDir + "/miRExplore/orgmir.tsv")

print("Starting MIRNA_ORG")
for orgmir in orgmirDB.orgmir2mimat:
    db.createNodeIfNotExists(['MIRNA_ORGMIR','MIRNAS'], {'id': orgmir})

    orgFamilies = orgmirDB.orgmir2mi.get(orgmir, None)

    if orgFamilies != None:

        orgmir2families = set()
        for mi in orgFamilies:
            orgFamily = mirFamiliesByMI.get(mi, None)

            if orgFamily != None:
                famID = orgFamily.familyID
                orgmir2families.add(famID)

        for famID in orgmir2families:
            db.createRelationship('org', ['MIRNA_ORGMIR'], {'id': orgmir}, 'mifam', ['MIRNA_FAMILY'], {'id': famID}, ['MIRNA_BELONGS_TO'], None)

print("Finished MIRNA_ORGMIR")

print("Starting MIRNA_ORGMI")
for orgmi in orgmirDB.orgmi2mi:
    db.createNodeIfNotExists(['MIRNA_ORGMI','MIRNAS'], {'id': orgmi})

    orgFamilies = orgmirDB.orgmi2mi.get(orgmi, None)

    if orgFamilies != None:
        orgmi2families = set()

        for mi in orgFamilies:
            orgFamily = mirFamiliesByMI.get(mi, None)

            if orgFamily != None:

                famID = orgFamily.familyID
                orgmi2families.add(famID)

        for famID in orgmi2families:
            db.createRelationship('org', ['MIRNA_ORGMI'], {'id': orgmi}, 'mifam', ['MIRNA_FAMILY'], {'id': famID}, ['MIRNA_BELONGS_TO'], None)

print("Finished MIRNA_ORGMI")

print("Starting MIRNA, MIRNA_PRE")

addedAccessionMIMATs = set()

doneRows = 0
for row in mirbase:
    # CREATE PRE-MATURE MIRNA NODE
    mirnaAccession = row['Accession']
    db.createNodeIfNotExists(['MIRNA_PRE','MIRNAS'], {'id': mirnaAccession})

    mirnaFamily = mirFamiliesByMI[mirnaAccession] if mirnaAccession in mirFamiliesByMI else None
    mirnaFamilyID = None

    if mirnaFamily != None:
        mirnaFamilyID = mirnaFamily.familyID

    #MI 2 ORGMIR
    orgmis = orgmirDB.mi2orgmi.get(mirnaAccession, None)
    if orgmis != None:
        allOrgMIs = set(orgmis)
        for orgmi in allOrgMIs:
            db.createRelationship('pre', ['MIRNA_PRE'], {'id': mirnaAccession}, 'orgmi', ['MIRNA_ORGMI'], {'id': orgmi}, ['IS_ORG_MI'], None)

    orgmirs = orgmirDB.mi2orgmir.get(mirnaAccession, None)
    if orgmirs != None:
        allOrgMIRs = set(orgmirs)
        for orgmir in allOrgMIRs:
            db.createRelationship('pre', ['MIRNA_PRE'], {'id': mirnaAccession}, 'orgmir', ['MIRNA_ORGMI'], {'id': orgmir}, ['IS_ORG_MIR'], None)

    if mirnaFamilyID != None and len(mirnaFamilyID) > 0:
        db.createRelationship('pre', ['MIRNA_PRE'], {'id': mirnaAccession}, 'mifam', ['MIRNA_FAMILY'], {'id': mirnaFamilyID}, ['MIRNA_BELONGS_TO'], None)

    matureID1 = None
    matureID2 = None
    matureAcc1 = None
    matureAcc2 = None

    matureID1 = row.getColumn('Mature1_ID', None)
    matureID2 = row.getColumn('Mature2_ID', None)
    matureAcc1 = row.getColumn('Mature1_Acc', None)
    matureAcc2 = row.getColumn('Mature2_Acc', None)

    if not (matureAcc1 == None or matureAcc1 == 'None' or matureAcc1 in addedAccessionMIMATs):
        addedAccessionMIMATs.add(matureAcc1)
        createMatureMIRNAEntry(matureAcc1, matureID1, mirnaAccession, mirnaFamilyID)

    if not (matureAcc2 == None or matureAcc2 == 'None' or matureAcc2 in addedAccessionMIMATs):
        addedAccessionMIMATs.add(matureAcc2)
        createMatureMIRNAEntry(matureAcc2, matureID2, mirnaAccession, mirnaFamilyID)

    if doneRows % 1000 == 0:
        print("Done MI lines: " + str(doneRows/len(mirbase)))
    doneRows += 1

db.close()