from collections import defaultdict

import HTSeq

import intervaltree


class GeneNeighbourDB:

    def __init__(self):

        self.chr2neighbours = defaultdict(lambda: intervaltree.IntervalTree())
        self.id2pos = {}


    def get_neighbours(self, geneid, offset=100000):

        if not geneid in self.id2pos:
            return None


        pos = self.id2pos[geneid]

        tree = self.chr2neighbours[pos[0]]

        return tree[pos[1]-offset:pos[2]+offset]


    @classmethod
    def loadFromFile(cls, inputgff = "/mnt/c/ownCloud/data/miRExplore/obodir/mm10_primary_assembly_and_lncRNA.gtf"):


        ret = GeneNeighbourDB()


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

                    start = int(line[3])
                    end = int(line[4])

                    strand = line[6]

                    ret.chr2neighbours[chr].addi(start, end, geneID)
                    ret.id2pos[geneID] = (chr, start, end)

                    addedGenes += 1

                    if addedGenes % 1000 == 0:
                        print("Added Genes", addedGenes)



            print(addedGenes)

        return ret


if __name__ == '__main__':



    rDB = GeneNeighbourDB.loadFromFile()

    elems = rDB.get_neighbours("ENSMUSG00000025907.14")

    print([x.data for x in elems])

