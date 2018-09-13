from collections import defaultdict


class PMID2PMCDB:


    def __init__(self):

        self.pmid2pmc = defaultdict(set)
        self.pmc2pmid = defaultdict(set)

    def hasID(self, sid):

        if sid in self.pmid2pmc:
            return True

        if str(sid) in self.pmid2pmc:
            return True

        return False

    def getID(self, sid, default=None):

        if sid in self.pmid2pmc:
            return self.pmid2pmc[sid]

        if str(sid) in self.pmid2pmc:
            return self.pmid2pmc[str(sid)]

        return default

    @classmethod
    def loadFromFile(cls, filepath):


        ret = PMID2PMCDB()

        with open(filepath, 'r') as fin:

            # 29113155        PMC22222222
            for line in fin:

                line = line.strip()

                aline = line.split('\t')

                if len(aline) < 2:
                    continue

                pmid = aline[0]
                pmc = aline[1]

                ret.pmid2pmc[pmid].add(pmc)
                ret.pmc2pmid[pmc].add(pmid)


        return ret