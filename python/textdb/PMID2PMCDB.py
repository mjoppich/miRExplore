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


    def getAllPMIDs(self):
        return [x for x in self.pmid2pmc]

    def getAllPMCs(self):
        return [x for x in self.pmc2pmid]

    @classmethod
    def loadFromFile(cls, filepath, PMC2PMID=True):


        ret = PMID2PMCDB()

        with open(filepath, 'r') as fin:

            # 29113155        PMC22222222
            for line in fin:

                aline = line.strip().split('\t')

                assert(len(aline) == 2)

                if PMC2PMID:
                    pmid = aline[1]
                    pmc = aline[0]
                else:
                    pmid = aline[0]
                    pmc = aline[1]

                assert(pmc.startswith("PMC"))

                ret.pmid2pmc[pmid].add(pmc)
                ret.pmc2pmid[pmc].add(pmid)


        return ret