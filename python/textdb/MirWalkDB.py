from collections import defaultdict

import os

from synonymes.SynonymFile import Synfile
from textdb.AbstractDBClasses import DataBaseRel, DataBaseDescriptor
from utils.DataFrame import DataFrame


class MirWalkEntry(DataBaseRel):

    def __init__(self, lent, rent, datasource, dataID, bindPos, bindProb, organism, ensID):

        self.lent = lent
        self.rent = rent

        self.data_source=datasource
        self.data_id = dataID

        self.binding_pos = tuple(bindPos)
        self.binding_probability = bindProb
        self.data_type = "predicted"
        self.ensemblID = ensID

        self.orgs = tuple(organism)

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
        return None

    def toJSON(self):

        retJSON= {
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,

            'data_type': self.data_type,
            'data_source': self.data_source,
            'data_id': self.data_id,

            'bind_position': self.binding_pos,
            'bind_probability': self.binding_probability,

            'target_ensembl_id': self.ensemblID

        }

        if self.orgs != None:
            retJSON['orgs'] = tuple(self.orgs)

        return retJSON

class MirWalkDB(DataBaseDescriptor):


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

    @classmethod
    def loadFromFile(cls, filepath,org="mmu", bindSite='3UTR',normGeneSymbols=None):

        ltype = 'gene'
        rtype = 'mirna'

        ret = MirWalkDB(ltype, rtype)
        file_base = os.path.basename(filepath)

        geneSymbolsNormalized=0


        seenMirnas = {}

        with open(filepath) as fin:


            for idx, line in enumerate(fin):

                if idx == 0:
                    continue

                if idx % 1000000 == 0:
                    print("Loaded gene mirnas", idx, filepath)

                line = line.strip()
                line = line.split("\t")
                #mmu-miR-328-3p  ENSMUST00000052172      Cxcr4   1685,1710       0.807692307692


                lid = line[2]

                lid = lid.upper()
                if normGeneSymbols != None and lid in normGeneSymbols:
                    lid = normGeneSymbols[lid]
                    geneSymbolsNormalized = 0

                rid = line[0]

                if rid in seenMirnas:
                    org, rid = seenMirnas[rid]
                else:
                    origRid = rid
                    org, rid = cls.harmonizeMIRNA(rid)
                    seenMirnas[origRid] = (org, rid)

                if not org in ['hsa', 'mmu']:
                    continue

                organism = org
                orgs = set()
                orgs.add(org)
                orgs.add(organism)

                dataID = "mirwalk_"+organism+ "_" + bindSite + "_"+str(idx)
                dataSource = 'mirwalk'

                bindPos = tuple([int(x) for x in line[3].split(",")])
                bindProb = float(line[4])

                gene_ensembl_id = line[1]

                relations = set([
                    MirWalkEntry((lid, ltype), (rid, rtype), dataSource, dataID, bindPos, bindProb, orgs, gene_ensembl_id)
                    ])

                for rel in relations:

                    ret.ltype2rel[lid].add(rel)
                    ret.rtype2rel[rid].add(rel)

                ret.all_ltypes.add(lid)
                ret.all_rtypes.add(rid)

        print("Gene Symbols Normalized", geneSymbolsNormalized)


        return ret


if __name__ == '__main__':


    retDB = MirWalkDB.loadFromFile('/mnt/c/ownCloud/data/miRExplore/mirwalk/mmu_miRWalk_3UTR.txt', org="mmu", bindSite="3UTR")

    allRes = retDB.get_rels("gene", "Cxcr4")

    for x in allRes:
        print(x.toJSON())






