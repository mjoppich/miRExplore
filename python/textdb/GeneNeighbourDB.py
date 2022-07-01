import json
from collections import defaultdict

import HTSeq

import intervaltree

from textdb.SymbolEnsemblDB import SymbolEnsemblDB


class GeneNeighbourDB:

    def __init__(self, org):

        self.chr2neighbours = defaultdict(lambda: intervaltree.IntervalTree())
        self.id2pos = {}
        self.orgid = org


    def get_neighbours(self, geneid, offset=100000):

        if not geneid in self.id2pos:
            return None


        pos = self.id2pos[geneid]

        tree = self.chr2neighbours[pos[0]]

        elems = tree[pos[1]-offset:pos[2]+offset]

        return [x.data for x in elems]


    def get_neighbour_regions(self, geneids, offset=100000, name=None, sym2ensDB=None):

        allInfo = {}
        for geneid in geneids:
            ret = self.get_neighbours(geneid, offset)

            if ret == None:
                continue

            geneIDInfo = {}

            targetName = geneid

            if name != None:
                targetName = name + " (" + geneid + ")"

            geneIDInfo[geneid] = {
                "chr": self.id2pos[geneid][0],
                "start": self.id2pos[geneid][1],
                "end": self.id2pos[geneid][2],
                "strand": self.id2pos[geneid][3],
                "name": targetName
            }

            geneIDInfo['nb'] = []

            for nb in ret:


                nbinfo = {
                "chr": self.id2pos[nb][0],
                "start": self.id2pos[nb][1],
                "end": self.id2pos[nb][2],
                "strand": self.id2pos[nb][3],
                "name": nb
                }


                if sym2ensDB != None:
                    retval = sym2ensDB.get_symbol_for_ensembl(nb)

                    if retval != None:
                        nbinfo["name"] = retval

                geneIDInfo['nb'].append( nbinfo )

            minStart = geneIDInfo[geneid]['start']

            rangeStart = 0
            rangeEnd = 0

            for elem in geneIDInfo['nb']:
                if elem['start'] < minStart:
                    minStart = elem['start']

            minStart = max([minStart-100, 0])

            geneIDInfo[geneid]['show_start'] = geneIDInfo[geneid]['start'] - minStart
            geneIDInfo[geneid]['show_end'] = geneIDInfo[geneid]['end'] - minStart

            for elem in geneIDInfo['nb']:

                elem['show_start'] = elem['start'] - minStart
                elem['show_end'] = elem['end'] - minStart

                if elem['show_end'] > rangeEnd:
                    rangeEnd = elem['show_end']

            geneIDInfo['range'] = {'start': rangeStart, 'end': rangeEnd}

            allInfo[geneid] = geneIDInfo

        return allInfo


    @classmethod
    def loadFromFile(cls, org, inputgff = "/mnt/c/ownCloud/data/miRExplore/obodir/mm10_primary_assembly_and_lncRNA.gtf"):


        ret = GeneNeighbourDB(org)


        with open(inputgff, 'r') as fin:

            addedGenes = 0


            for line in fin:

                line = line.strip().split()

                #print(line)

                chr = line[0]
                type = line[2]

                if type != 'gene':
                    continue

                attrStr = "\t".join(line[8:])
                allAttr = HTSeq.parse_GFF_attribute_string(attrStr=attrStr)

                geneID = allAttr.get('gene_id', None)

                if geneID != None:

                    origID = geneID

                    if "." in geneID:
                        geneID = geneID[:geneID.index(".")]

                    if geneID in ret.id2pos:
                        print("Duplicate gene ID", geneID, origID)

                    start = int(line[3])
                    end = int(line[4])

                    strand = line[6]

                    ret.chr2neighbours[chr].addi(start, end, geneID)
                    ret.id2pos[geneID] = (chr, start, end, strand)

                    addedGenes += 1

                    if addedGenes % 10000 == 0:
                        print("Added Genes", addedGenes)



            print(addedGenes)

        return ret


if __name__ == '__main__':


    org = "mmu"

    symbol2ensemblDB = SymbolEnsemblDB.loadFromFile("/mnt/c/ownCloud/data/miRExplore/obodir/sym2ens/")
    allensids = symbol2ensemblDB.get_all_genes("CCL2")

    print(allensids)

    if not org in allensids:
        exit(-1)

    orgCollectGenes = allensids[org]

    if orgCollectGenes == None or len(orgCollectGenes) == 0:
        exit(-2)


    rDB = GeneNeighbourDB.loadFromFile("mmu")

    elems = rDB.get_neighbour_regions(orgCollectGenes, name="CCL2", sym2ensDB=symbol2ensemblDB)

    print(json.dumps(elems, indent=2))

    for x in elems:
        print(x, elems[x])
