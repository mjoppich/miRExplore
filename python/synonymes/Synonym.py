import csv

from io import StringIO


class Synonym:

    @classmethod
    def parseFromLine(cls, line):
        aLine = line.split(':', 1)

        retSyn = Synonym(aLine[0].strip())

        aSyns = aLine[1].split('|')
        for x in [x.strip() for x in aSyns]:
            retSyn.addSyn(x)

        return retSyn

    def __init__(self, id):

        id = id.strip()
        id = id.replace(':', '_')
        self.id = id
        self.currentIdx = 0
        self.syns = []

    def get(self, synIdx):

        if synIdx >= 0 and synIdx < len(self.syns):
            return self.syns[synIdx]

        return None

    def removeCommonSynonymes(self, commonSyns):

        for x in commonSyns:
            if x in self:
                self.removeSyn(x)

    def getSynonymes(self):

        return self.syns

    def addTextSyns(self, synText):

        if synText == None:
            return

        synText = synText.strip()

        if len(synText) == 0:
            return

        file_like = StringIO(synText)
        csvreader = csv.reader(file_like, delimiter=',', quotechar='"', skipinitialspace=True)

        for line in list(csvreader):
            for x in line:
                self.addSyn(x)

    def addSyn(self, newSyn):

        if newSyn == None:
            return

        newSyn =newSyn.strip()

        if newSyn.startswith("see "):
            newSyn = newSyn[len("see "):]

        if newSyn.endswith("~withdrawn"):
            newSyn = newSyn[:-len("~withdrawn")]

        if len(newSyn) == 0:
            return

        if len(newSyn) < 3:
            return

        if newSyn[0] == newSyn[len(newSyn)-1] and newSyn[0] == '"':
            newSyn = newSyn[1:len(newSyn)-1]

        if not newSyn in self.syns:
            self.syns.append(newSyn)

    def __contains__(self, item):
        return item in self.syns

    def __iter__(self):

        self.currentIdx = 0
        return self

    def __next__(self):

        if self.currentIdx >= len(self.syns):
            raise StopIteration

        self.currentIdx+=1
        return self.syns[self.currentIdx-1]

    def __str__(self):
        return self.id.replace(':', '_') + ":" + "|".join(self.syns)

    def __len__(self):
        return len(self.syns)

    def removeSyn(self, syn):

        if type(syn) == str:

            if syn in self.syns:
                idx = self.syns.index(syn)
                self.removeSyn(idx)

        else:

            isyn = int(syn)

            if isyn >= 0 and isyn < len(self.syns):
                del self.syns[isyn]
