from abc import ABC, abstractmethod
from collections import defaultdict

from synonymes.mirnaID import miRNA
from textdb.SQLiteBase import SQLiteBase


class DataBaseRel(ABC):

    @property
    @abstractmethod
    def docid(self):
        pass

    @property
    @abstractmethod
    def ltype(self):
        pass

    @property
    @abstractmethod
    def rtype(self):
        pass

    @property
    @abstractmethod
    def lid(self):
        pass

    @property
    @abstractmethod
    def rid(self):
        pass

    @property
    def l_id_type(self):
        return (self.lid, self.ltype)

    @property
    def r_id_type(self):
        return (self.rid, self.rtype)

    @abstractmethod
    def toJSON(self):
        pass

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(tuple([(x, self.__dict__[x]) for x in sorted(self.__dict__)]))

    def get_interactor(self, type):

        if self.ltype == type:
            return self.lid
        elif self.rtype == type:
            return self.rid

        return None

class DataBaseDescriptor(ABC, SQLiteBase):

    def __init__(self):

        super(DataBaseDescriptor, self).__init__(infile=None)

        self.ltype2rel = defaultdict(set)
        self.rtype2rel = defaultdict(set)

        self.all_ltypes = set()
        self.all_rtypes = set()

        self.lontology = None
        self.rontology = None

    @property
    def l_ont_based(self):
        return self.lontology != None

    @property
    def r_ont_based(self):
        return self.lontology != None

    def get_evidence_docids(self):

        docIDs = set()

        for lid in self.ltype2rel:
            for ev in self.ltype2rel[lid]:

                if ev.docid != None:
                    docIDs.add(ev.docid)

        return docIDs

    @classmethod
    def harmonizeMIRNA(cls, mirna):
        """

        :param mirna:
        :return: tries to return a normalized name ...

        9761 microRNA
       7958 MicroRNA
       2311 MiRNA
       2191 miRNA
       1844 hsa
       1440 let
        578 miRNAS
        437 MICRORNA
        299 microRNAS
        256 MIRNA
        155 mmu
        125 Micro
        116 micro
        """

        possibleOrgStarts = ['mmu', 'hsa']

        recOrg = None

        for x in possibleOrgStarts:

            if mirna.startswith(x + "-"):
                recOrg = x
                mirna = mirna.replace(x+"-", "", 1)
                break

        possibleStarts = ['microRNA', 'MicroRNA', 'MiRNA', 'miRNA', 'MICRORNA', 'microRNA', 'MIRNA']

        for x in possibleStarts:

            if mirna.startswith(x + "-"):
                mirna = mirna.replace(x, "miR", 1)
                break

        try:
            oMirna = miRNA(mirna)
            mirna = oMirna.normalized_str()

        except:
            pass


        return (recOrg, mirna)


    @property
    @abstractmethod
    def ltype(self):
        pass

    @property
    @abstractmethod
    def rtype(self):
        pass

    def get_rels(self, etype, eid):

        if etype == self.ltype:
            return self.get_lid_rels(eid)
        elif etype == self.rtype:
            return self.get_rid_rels(eid)

        return None


    def get_lid_rels(self, geneID):

        if not self.l_ont_based:
            return self.ltype2rel.get(geneID, None)
        else:
            allRels = self.ltype2rel.get(geneID, None)

            if allRels == None:
                return None

            return self.undoOntID(allRels)

    def get_rid_rels(self, geneID):

        if not self.r_ont_based:
            return self.rtype2rel.get(geneID, None)
        else:
            allRels = self.rtype2rel.get(geneID, None)

            if allRels == None:
                return None

            return self.undoOntID(allRels)


    def undoOntID(self, rels):
        return rels

    @classmethod
    @abstractmethod
    def loadFromFile(cls, filepath, ltype='gene', rtype='mirna'):
        pass