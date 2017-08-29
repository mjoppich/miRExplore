import codecs
from collections import OrderedDict, Counter

from synonymes.Synonym import Synonym
class Synfile:

    def __init__(self, sFileLocation):

        self.mSyns = {}
        self.line2syn = {}

        def addSyn(sLine, iLine):

            oSyn = Synonym.parseFromLine(sLine)

            self.mSyns[ oSyn.id ] = oSyn
            self.line2syn[iLine] = oSyn.id

        with codecs.open(sFileLocation, 'r', 'latin1') as infile:
            idx = 0
            for line in infile:
                addSyn(line, idx)
                idx += 1

        self.synIDs = None
        self.synIDidx = None

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