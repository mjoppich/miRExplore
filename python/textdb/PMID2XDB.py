from collections import defaultdict


class PMID2XDB:


    def __init__(self):

        self.pmid2pmc = {}

    def hasPMC(self, pmid):

        if pmid in self.pmid2pmc:
            return True

        if str(pmid) in self.pmid2pmc:
            return True

        return False

    def getPMC(self, pmid, default=None):

        if pmid in self.pmid2pmc:
            return self.pmid2pmc[pmid]

        if str(pmid) in self.pmid2pmc:
            return self.pmid2pmc[str(pmid)]

        return default

    @classmethod
    def loadFromFile(cls, filepath):


        ret = PMID2XDB()

        with open(filepath, 'r') as fin:

            # 29113155        DOID:1389       polyneuropathy  [('29113155.2.1', 0, 14), ('29113155.2.2', 139, 153)]
            for line in fin:

                line = line.strip()

                aline = line.split('\t')

                pmid = aline[0]
                termid = aline[1]
                termname = aline[2]

                loc = eval(aline[3])


                info = {
                    'docid': pmid,
                    'termid': termid,
                    'termname': termname,
                    'evidences': loc
                }


                ret.pmid2info[pmid].append(info)


        return ret

