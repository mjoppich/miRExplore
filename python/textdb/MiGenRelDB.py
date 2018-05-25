from collections import defaultdict

import os


class MiRGeneRel:

    def __init__(self, intuple, docid, sameParagraph, sameSentence, lent, rent, datasource, dataID):

        self.assocDir = None
        self.assocDirV = None
        self.assocCat = None
        self.assocFound = None
        self.assocSent = None
        self.assocNeg = None
        self.lPOS = None
        self.rPOS = None
        self.assocPos = None

        self.ltype=lent[1]
        self.rtype=rent[1]

        self.lid = lent[0]
        self.rid = rent[0]

        self.data_source=datasource
        self.data_id = dataID


        if intuple != None:
            self.assocDir = intuple[0]
            self.assocDirV = intuple[1]
            self.assocCat = intuple[2]
            self.assocFound = intuple[3]
            self.assocSent = intuple[4]
            self.assocNeg = intuple[5]
            self.rPOS = intuple[6]
            self.lPOS = intuple[7]
            self.assocPos = intuple[8]

        self.docid = docid
        self.same_sentence = sameSentence
        self.same_paragraph = sameParagraph

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):

        return hash(tuple([(x, self.__dict__[x]) for x in sorted(self.__dict__)]))


    def toJSON(self):

        return {
            'rel_direction': self.assocDir,
            'rel_direction_verb': self.assocDirV,
            'rel_category': self.assocCat,
            'rel_verb': self.assocFound,
            'rel_sentence': self.assocSent,
            'rel_negated': self.assocNeg,
            'lpos': self.lPOS,
            'rpos': self.rPOS,
            'rel_pos': self.assocPos,
            'same_paragraph': self.same_paragraph,
            'same_sentence': self.same_sentence,
            'docid': self.docid,
            'ltype': self.ltype,
            'rtype': self.rtype,
            'rid': self.rid,
            'lid': self.lid,

            'data_source': self.data_source,
            'data_id': self.data_id
        }

class MiGenRelDB:


    def __init__(self, ltype, rtype):

        self.ltype = ltype
        self.rtype = rtype

        self.ltype2rel = defaultdict(set)
        self.rtype2rel = defaultdict(set)

        self.all_ltypes = set()
        self.all_rtypes = set()

    def get_evidence_docids(self):

        docIDs = set()

        for lid in self.ltype2rel:
            for ev in self.ltype2rel[lid]:
                docIDs.add(ev.docid)

        return docIDs

    def get_rels(self, etype, eid):

        if etype == self.ltype:
            return self.get_lid_rels(eid)
        elif etype == self.rtype:
            return self.get_rid_rels(eid)

        return None

    def get_lid_rels(self, geneID):
        return self.ltype2rel.get(geneID, None)

    def get_rid_rels(self, mirnaID):
        return self.rtype2rel.get(mirnaID, None)

    @classmethod
    def loadFromFile(cls, filepath, ltype='gene', rtype='mirna'):


        def harmonizeMIRNA(mirna):
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

            possibleStarts = ['microRNA', 'MicroRNA', 'MiRNA', 'miRNA', 'MICRORNA', 'microRNA', 'MIRNA']


            for x in possibleStarts:

                if mirna.startswith(x +"-"):
                    mirna = mirna.replace(x, "miR", 1)
                    break

            return mirna




        ret = MiGenRelDB(ltype, rtype)

        file_base = os.path.basename(filepath)

        with open(filepath, 'r') as fin:

            #PCBP1	miR-3978	ORGMIR4299	29138420	True	True	[('GM', 'GVM', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (6, 16)), ('GM', 'GMV', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (35, 45))]

            for lineIdx, line in enumerate(fin):

                line = line.strip()

                aline = line.split('\t')

                lid = aline[0]
                rid = aline[1]
                rid = harmonizeMIRNA(rid)

                docid = aline[3]

                sameParagraph = aline[4] == 'True'
                sameSentence = aline[5] == 'True'

                if not sameSentence:
                    continue

                tmRelations = None if aline[6] == 'None' else eval(aline[6])

                if tmRelations != None:
                    allrels = set()

                    for relIdx, rel in enumerate(tmRelations):
                        newrel = MiRGeneRel(rel, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', file_base+"_"+str(lineIdx)+"_"+str(relIdx))

                        allrels.add(newrel)

                    relations = allrels
                else:
                    relations = set([MiRGeneRel(None, docid, sameParagraph, sameSentence, (lid, ltype), (rid, rtype), 'pmid', file_base+"_"+str(lineIdx)+"_0")])

                for rel in relations:

                    ret.ltype2rel[lid].add(rel)
                    ret.rtype2rel[rid].add(rel)

                ret.all_ltypes.add(lid)
                ret.all_rtypes.add(rid)

        return ret

