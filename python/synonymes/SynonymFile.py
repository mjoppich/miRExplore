import codecs
from collections import OrderedDict, Counter

from synonymes.GeneOntology import GeneOntology, GOTerm, GOSynonyme, GOSynonymeScope
from synonymes.Synonym import Synonym



class Synfile:


    @classmethod
    def syn2goterm(cls, syn, termid):

        ret = GOTerm()

        ret.id = termid
        ret.name = syn.id
        ret.synonym = set()

        for synWord in syn.syns:

            goSyn = GOSynonyme( synWord, GOSynonymeScope.EXACT)

            ret.synonym.add(goSyn)


        return ret



    def __init__(self, sFileLocation):

        self.mSyns = {}
        self.line2syn = {}

        self.synIDs = None
        self.synIDidx = None

        self.location = None

        if sFileLocation != None:
            self._load_file(sFileLocation)

    def _load_file(self, sFileLocation):

        def addSyn(sLine, iLine):

            oSyn = Synonym.parseFromLine(sLine)

            self.mSyns[ oSyn.id ] = oSyn
            self.line2syn[iLine] = oSyn.id

        with codecs.open(sFileLocation, 'r', 'latin1') as infile:
            idx = 0
            for line in infile:
                addSyn(line, idx)
                idx += 1

            self.location = sFileLocation

    def export_flat_obo(self, filepath, termPrefix):

        go = GeneOntology()

        newterms = []

        for synid in self.mSyns:

            syn = self.mSyns[synid]

            got = self.syn2goterm(syn, synid)#termPrefix + str(len(go.dTerms) + len(newterms)))

            if got != None:
                newterms.append(got)

        go.addterms(newterms)

        go.saveFile(filepath)



    def __iter__(self):

        self.synIDs = [x for x in self.mSyns]
        self.synIDidx = 0

        return self

    def __next__(self):

        curIdx = self.synIDidx
        self.synIDidx += 1

        if curIdx < len(self.synIDs):
            return self.mSyns[self.synIDs[curIdx]]

        raise StopIteration()


    def __len__(self):
        return len(self.mSyns)

    def get(self, iSynID):

        return self.mSyns.get(self.line2syn.get(iSynID, None), None)

    def histogramSynonymes(self):

        oSynCounter = Counter()

        for iSynID in self.mSyns:

            oSyn = self.mSyns[iSynID]

            vSyns = set(oSyn.getSynonymes())

            for sSyn in vSyns:

                oSynCounter[sSyn] += 1

        return oSynCounter

    def histogramSynonymesOrdered(self, limit = None):

        synCounter = self.histogramSynonymes()

        if limit != None:
            synOrdered = OrderedDict( synCounter.most_common(limit) )
        else:
            synOrdered = OrderedDict(sorted(synCounter.items(), key=lambda x: x[1]))

        return synOrdered


class AssocSynfile(Synfile):

    def __init__(self, path):
        super().__init__(None)

        self.synid2class = {}
        self.synline2class = {}

        self._load_file(path)

    def _load_file(self, sFileLocation):
        def addSyn(sLine, iLine):

            aline = sLine.strip().split(",")

            regClass = aline[0]
            sLine = aline[1]

            oSyn = Synonym.parseFromLine(sLine)

            self.mSyns[oSyn.id] = oSyn
            self.line2syn[iLine] = oSyn.id

            self.synid2class[oSyn.id] = regClass
            self.synline2class[iLine] = regClass

        with codecs.open(sFileLocation, 'r', 'latin1') as infile:
            idx = 0
            for line in infile:
                addSyn(line, idx)
                idx += 1



if __name__ == '__main__':

    synfile = Synfile("/mnt/c/dev/data/tm_soehnlein/synonyms/messenger.syn")

    synfile.export_flat_obo("/mnt/c/dev/data/tm_soehnlein/obos/messenger.obo", "MSG:")
