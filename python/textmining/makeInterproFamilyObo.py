import gzip
import os
import shutil

from urllib import request

from collections import defaultdict

from nertoolkit.geneontology.GeneOntology import GeneOntology, GOTerm, GORelationType, GORelation

from utils.idutils import miRExploreDir

"""

Step 1: Download Interpro XML

"""
interproFolder = miRExploreDir + "/interpro/"
interproXMLPath = interproFolder + "/interpro.xml"

if not os.path.isfile(interproXMLPath):

    downloadFile = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz"
    downloadLocation = interproFolder + "/interpro.xml.gz"
    print(downloadFile)

    request.urlretrieve(downloadFile, downloadLocation)


    with gzip.open(downloadLocation, 'rb') as f_in, open(interproXMLPath, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


interproELPath = interproFolder + "/interpro.entry.list"

if not os.path.isfile(interproELPath):

    downloadFile = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list"
    downloadLocation = interproELPath
    print(downloadFile)

    request.urlretrieve(downloadFile, downloadLocation)

"""

STEP 2: build obo

"""

import xml.etree.ElementTree as ET

tree = ET.parse(interproXMLPath)
root = tree.getroot()

class IPRFamilyElement:

    def __init__(self):

        self.id = None
        self.name = None
        self.parents = []
        self.xrefs = []

    def __str__(self):
        return "<IPRFamilyElement id={id} name='{name}' parents='{parents}' xrefs='{xrefs}' >".format( id=self.id, name=self.name, parents=",".join(self.parents), xrefs=",".join(self.xrefs) )

    def to_go_term(self):

        retRes = GOTerm()
        retRes.id = self.id
        retRes.name = self.name
        retRes.xref = self.xrefs

        retRes.is_a = []
        for parent in self.parents:

            oRelation = GORelationType['IS_A']
            rel = GORelation(oRelation, termid=parent, desc=None)
            retRes.is_a.append(rel)

        return retRes


seenElements = {}

for interproID in root.findall('interpro'):

    elemType = interproID.attrib['type']

    if elemType != 'Family':
        continue

    elemID = interproID.attrib['id']

    children = defaultdict(list)
    for x in interproID:
        children[ x.tag ].append(x)

    elemName = [x.text for x in children['name']]
    elemName = "|".join(elemName)

    elemRefs = []
    for parent in children['parent_list']:
        for x in parent:
            elemRefs.append(x.attrib['ipr_ref'])

    elemXRefs = []
    for extDoc in children['external_doc_list']:
        for x in extDoc:
            elemXRefs.append(x.attrib['db'] + ":" + x.attrib['dbkey'])

    iprElem = IPRFamilyElement()
    iprElem.id = elemID
    iprElem.name = elemName
    iprElem.parents = elemRefs
    iprElem.xrefs = elemXRefs

    seenElements[iprElem.id] = iprElem


noParents = 0
for iprID in seenElements:

    iprElem = seenElements[iprID]

    if len(iprElem.parents) == 0:
        noParents += 1

    for id in iprElem.parents:
        if not id in seenElements:
            print(id)

    print(iprElem)

print("Elem Count", len(seenElements))
print("Root Count", noParents)

interProObo = GeneOntology()

for iprID in seenElements:
    iprElem = seenElements[iprID]

    iprGoElem = iprElem.to_go_term()

    interProObo.dTerms[iprID] = iprGoElem

interProObo.linkChildren()
interProObo.saveFile(interproFolder + "/interpro.obo")


