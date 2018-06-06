from collections import defaultdict


class TestRelLoader:

    def __init__(self):

        self.tester2rels = defaultdict(list)



    @classmethod
    def loadFromFile(cls, file):

        retDB = TestRelLoader()

        with open(file, 'r') as fin:

            for line in fin:

                aline = line.strip().split("\t")

                tester = aline[0]
                lid = aline[1]
                rid = aline[2]
                evSentID = aline[3]

                retDB.tester2rels[tester].append((lid, rid, evSentID))


        return retDB