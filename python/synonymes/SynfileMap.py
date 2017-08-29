from synonymes.SynonymFile import Synfile


class SynonymID:

    def __init__(self):

        self.synfile = None
        self.synid = None
        self.synmatch = None

    def __str__(self):
        return "{synfile}:{synid}:{synmatch}".format(synfile=self.synfile, synid=self.synid, synmatch=self.synmatch)

    def __repr__(self):
        return self.__str__()

    @classmethod
    def fromStr(cls, string):

        retObj = SynonymID()

        astring = [x.strip() for x in string.split(':')]

        retObj.synfile = int(astring[0])
        retObj.synid = int(astring[1])
        retObj.synmatch = int(astring[2])

        return retObj

class SynfileMap:

    def __init__(self, file, loadSynFiles=False):
        self.synfiles = {}
        self.synids = {}

        for line in open(file, 'r'):
            aline = line.split(':')

            fileEntry = (aline[0].strip(), int(aline[1].strip()))

            self.synfiles[ fileEntry[1] ] = fileEntry[0]
            self.synids[fileEntry[0]] = fileEntry[1]

        self.loadedSynFiles = {}

        if loadSynFiles:
            self.loadSynFiles()

    def loadSynFiles(self, replacePath = None):

        for x in self.synfiles:
            synPath = self.synfiles[x]

            if replacePath != None:
                synPath = synPath.replace( replacePath[0], replacePath[1])

            self.loadedSynFiles[x] = Synfile(synPath)

    def getSynfiles(self):

        return [self.synfiles[x] for x in self.synfiles]


    def getSynonyme(self, x):

        if not type(x) == SynonymID:
            raise ValueError("x must be SynonymeID")

        synfile = self.loadedSynFiles.get(x.synfile, None)

        if synfile == None:
            return None

        synonyme = synfile.get(x.synid)
        return synonyme

    def getSynFileByID(self, id):

        return self.loadedSynFiles.get(id, None)

    def getSynID(self, x):

        if x in self.synids:
            return self.synids[x]

        return None

    def getSynName(self, x):

        if x in self.synfiles:
            return self.synfiles[x]

        return None