import copy
import os, sys

from synonymes.SynonymFile import Synfile
from textdb.DocOrganismDB import DocOrganismDB
from utils.DataFrame import DataFrame
from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import random
from collections import defaultdict

import os

from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor


class DIANATarbaseEntry(DataBaseRel):

    def __init__(self, intuple, docid, lent, rent, datasource, dataID):
        #(species, cellline, tissue, method, measure, direction)

        self.lent=lent
        self.rent=rent

        self.lontent = lent
        self.rontent = rent

        self.data_source=datasource
        self.data_id = dataID

        self.pubID = None


        self.orgs = intuple[0]

        self.cellline = intuple[1]
        self.tissue = intuple[2]
        self.method = intuple[3]
        self.measure = intuple[4]
        self.direction = intuple[5]


    @property
    def lid(self):
        return self.lent[0]

    @property
    def rid(self):
        return self.rent[0]

    @property
    def ltype(self):
        return self.lent[1]

    @property
    def rtype(self):
        return self.rent[1]

    @property
    def docid(self):
        return self.pubID


    def toJSON(self):

#        self.orgs = intuple[0]

        retJSON =  {
            'tissue': self.tissue,
            'method': self.method,
            'measure': self.measure,
            'direction': self.direction,
            'cellline': self.cellline,
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,
            'lontid': self.lontent[0],
            'rontid': self.rontent[0],
            'data_source': self.data_source,
            'data_id': self.data_id
        }

        if self.orgs != None:
            retJSON['orgs'] = tuple(self.orgs)

        return retJSON

class DIANATarbaseDB(DataBaseDescriptor):


    def __init__(self, ltype, rtype):

        super().__init__()

        self.ltyped = ltype
        self.rtyped = rtype

    @property
    def ltype(self):
        return self.ltyped

    @property
    def rtype(self):
        return self.rtyped

    def types(self):
        return set([self.ltype, self.rtype])

    def getOntologyForType(self, ontType):

        if ontType == self.ltype:
            return self.lontology

        if ontType == self.rtype:
            return self.rontology

        return None

    @classmethod
    def choices(cls, elems, k=1, excluded = []):

        resvec = []

        while len(resvec) < k:

            relem = random.choice(elems)

            if not relem in resvec and not relem in excluded:
                resvec.append(relem)

        return resvec


    def get_lid_rels(self, geneID):

        if not self.l_ont_based:
            return self.ltype2rel.get(geneID, None)
        else:
            allRels = self.ltype2rel.get(geneID, None)

            if allRels == None:
                return None

            rootObo = self.lontology.dTerms.get(geneID, None)

            if rootObo != None:

                allChildren = rootObo.getAllChildren()

                lenBefore = len(allRels)
                for child in allChildren:
                    allRels = allRels.union(self.ltype2rel.get(child.term.id, []))

                print("Added", len(allRels)-lenBefore, "for ontology")


            return self.undoOntID(allRels)

    def get_rid_rels(self, geneID):

        if not self.r_ont_based:
            return self.rtype2rel.get(geneID, None)
        else:
            allRels = self.rtype2rel.get(geneID, None)

            if allRels == None:
                return None

            rootObo = self.rontology.dTerms.get(geneID, None)

            if rootObo != None:

                allChildren = rootObo.getAllChildren()

                lenBefore = len(allRels)
                for child in allChildren:
                    allRels = allRels.union(self.rtype2rel.get(child.term.id, []))

                print("Added", len(allRels) - lenBefore, "for ontology")

            return self.undoOntID(allRels)


    def undoOntID(self, rels):

        for rel in rels:

            lid = rel.lid
            rid = rel.rid

            if self.l_ont_based:
                if lid in self.lontology.dTerms:
                    rel.lent = (self.lontology.dTerms[lid].name, rel.lent[1])

            if self.r_ont_based:
                if rid in self.rontology.dTerms:
                    rel.rent = (self.rontology.dTerms[rid].name, rel.rent[1])

        return rels


    @classmethod
    def loadFromFile(cls, filepath, dbtype='pmid', normGeneSymbols=None):

        syns = Synfile(os.path.dirname(filepath) + "/../textmine/synonyms/celllines.all.syn")
        syngrepFile = os.path.dirname(filepath) + "/celllines.index"

        cellline2obo = defaultdict(set)

        with open(syngrepFile, 'r') as fin:

            for line in fin:
                line = line.strip()
                line = line.split("\t")

                word = line[0].strip()
                syn = line[1].split(":")

                syn = syns.get(int(syn[1]))
                synID = syn.id

                cellline2obo[word].add(synID)


        ret = DIANATarbaseDB("gene", "mirna")

        species2org = {
            'Homo sapiens': "hsa",
            'Mus musculus': "mmu"
        }

        seenRels = set()

        geneSymbolsNormalized = 0

        dianaData = DataFrame.parseFromFile(filepath)

        foundCellInfos = []
        seenMirnas = {}

        for idx,row in enumerate(dianaData):

            """
            geneId  geneName        mirna   species cell_line       tissue  category        method  positive_negative       direct_indirect up_down condition
            0910001A06Rik   0910001A06Rik   mmu-miR-124-3p  Mus musculus    NA      NA      NA      Microarrays     POSITIVE        INDIRECT        UP      NA
            """


            geneID = row['geneName']
            mirna = row['mirna']

            species = row['species']
            cellline = row['cell_line']
            tissue = row['tissue']
            method = row['method']
            measure = row['direct_indirect']
            direction = row['up_down']

            if cellline == "NA":
                cellline = None

            if tissue == "NA":
                tissue = None

            if method == "NA":
                method = None

            if measure == "NA":
                measure = None


            if species not in ['Mus musculus', 'Homo sapiens']:
                continue

            org = None

            if mirna in seenMirnas:
                org, mirna = seenMirnas[mirna]
            else:
                origMirna = mirna
                (org, mirna) = cls.harmonizeMIRNA(mirna)
                seenMirnas[origMirna] = (org, mirna)

            docOrgs = set()

            if org != None:
                docOrgs.add(org)

            if species in species2org:
                docOrgs.add(species2org[species])

            geneID = geneID.upper()
            if geneID in normGeneSymbols:
                geneID = normGeneSymbols[geneID]

            docID = "DIANA:" + str(idx)

            if cellline in cellline2obo:

                for oboID in cellline2obo[cellline]:
                    celllInfo = {
                        'docid': docID,
                        'termid': oboID,
                        'termname': cellline,
                        'evidences': []
                    }

                    foundCellInfos.append(celllInfo)


            entry = DIANATarbaseEntry(
                (tuple(docOrgs), cellline, tissue, method, measure, direction),
                docID,
                (geneID, "gene"),
                (mirna, "mirna"),
                "DIANA",
                idx
            )

            ret.ltype2rel[geneID].add(entry)
            ret.rtype2rel[mirna].add(entry)

            ret.all_ltypes.add(geneID)
            ret.all_rtypes.add(mirna)

        return ret, foundCellInfos



if __name__ == '__main__':


    normGeneSymbols = normalize_gene_names(path="/mnt/c/ownCloud/data/miRExplore/obodir/"+"/hgnc_no_withdrawn.syn")

    ret, celllinfo = DIANATarbaseDB.loadFromFile("/mnt/c/ownCloud/data/miRExplore/diana/hsa_mmu.diana.csv", normGeneSymbols=normGeneSymbols)

    for x in ret.get_rels('gene', 'CXCR4'):
        print(x.toJSON())