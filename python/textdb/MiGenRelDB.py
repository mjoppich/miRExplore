from collections import defaultdict


class MiGenRel:

    def __init__(self, intuple):

        self.assocDir = intuple[0]
        self.assocDirV = intuple[1]
        self.assocCat = intuple[2]
        self.assocFound = intuple[3]
        self.assocSent = intuple[4]
        self.assocNeg = intuple[5]

        self.mirnaPos = intuple[6]
        self.genePos = intuple[7]
        self.assocPos = intuple[8]


    def toJSON(self):

        return {
            'rel_direction': self.assocDir,
            'rel_direction_verb': self.assocDirV,
            'rel_category': self.assocCat,
            'rel_verb': self.assocFound,
            'rel_sentence': self.assocSent,
            'rel_negated': self.assocNeg,
            'mirna_pos': self.mirnaPos,
            'gene_pos': self.genePos,
            'rel_pos': self.assocPos,
        }

class MiGenRelDB:


    def __init__(self):

        self.mirna2rel = defaultdict(list)
        self.gene2rel = defaultdict(list)


    @classmethod
    def loadFromFile(cls, filepath):


        ret = MiGenRelDB()

        with open(filepath, 'r') as fin:

            #PCBP1	miR-3978	ORGMIR4299	29138420	True	True	[('GM', 'GVM', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (6, 16)), ('GM', 'GMV', 'POS', 'express', '29138420.2.5', False, (26, 34), (0, 5), (35, 45))]

            for line in fin:

                line = line.strip()

                aline = line.split('\t')

                gene = aline[0]
                mirna = aline[1]

                if mirna.startswith('MiRNA'):
                    mirna = mirna.replace('MiRNA', 'miR', 1)

                mirnaID = aline[2]
                docid = aline[3]

                sameParagraph = eval(aline[4])
                sameSentence = eval(aline[5])

                relations = eval(aline[6])

                if relations != None:
                    allrels = []

                    for rel in relations:
                        newrel = MiGenRel(rel)

                        allrels.append(newrel)

                    relations = allrels

                info = {
                    'docid': docid,
                    'gene': gene,
                    'mirna': mirna,
                    'same_paragraph': sameParagraph,
                    'same_sentence': sameSentence,
                    'evidences': relations
                }


                ret.gene2rel[gene].append(info)
                ret.mirna2rel[mirna].append(info)

        return ret

