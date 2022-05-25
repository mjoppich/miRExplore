from collections import defaultdict


class PMID2XDB:


    def __init__(self, assocObo):

        self.docid2info = defaultdict(list)
        self.all_term_names = []

        if assocObo != None:

            for termid in assocObo.dTerms:

                oterm = assocObo.dTerms[termid]
                self.all_term_names.append((oterm.name, oterm.id))


    def add_database(self, obase):

        assert(isinstance(obase, PMID2XDB))

        for docid in obase.docid2info:

            for elem in obase.docid2info[docid]:
                self.docid2info[docid].append(elem)

        self.all_term_names += obase.all_term_names



    def hasDOC(self, docid):

        if docid in self.docid2info:
            return True

        if str(docid) in self.docid2info:
            return True

        return False

    def getDOC(self, docid, default=None):

        if docid in self.docid2info:
            return self.docid2info[docid]

        if str(docid) in self.docid2info:
            return self.docid2info[str(docid)]

        return default

    def getTerms(self):

        return self.all_term_names

    @classmethod
    def loadFromFile(cls, filepath, assocObo, reqDocIDs=None):


        ret = PMID2XDB(assocObo)

        with open(filepath, 'r') as fin:

            # 29113155        DOID:1389       polyneuropathy  [('29113155.2.1', 0, 14), ('29113155.2.2', 139, 153)]
            for line in fin:

                line = line.strip()
                aline = line.split('\t')
                pmid = aline[0]

                if reqDocIDs != None and not pmid in reqDocIDs:
                    continue

                termid = aline[1]
                termname = aline[2]

                loc = eval(aline[3])
                #loc = aline[3]


                info = {
                    'docid': pmid,
                    'termid': termid,
                    'termname': termname,
                    'evidences': loc
                }


                ret.docid2info[pmid].append(info)


        return ret

