from collections import OrderedDict, Counter

from synonymes.Synonym import Synonym
class Synfile:

    def __init__(self, sFileLocation):

        self.mSyns = {}

        def addSyn(sLine, iLine):

            oSyn = Synonym(sLine)

            self.mSyns[ oSyn.id ] = oSyn

        with open(sFileLocation, 'r') as infile:
            idx = 0
            for line in infile:
                addSyn(line, idx)
                idx += 1

    def __len__(self):
        return len(self.mSyns)

    def get(self, iSynID):

        if iSynID in self.mSyns:

            return self.mSyns[iSynID]

        return None

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