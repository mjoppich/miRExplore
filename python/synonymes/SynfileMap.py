class SynfileMap:

    def __init__(self, file):
        self.synfiles = {}
        self.synids = {}

        for line in open(file, 'r'):
            aline = line.split(':')

            fileEntry = (aline[0].strip(), int(aline[1].strip()))

            self.synfiles[ fileEntry[1] ] = fileEntry[0]
            self.synids[fileEntry[0]] = fileEntry[1]

    def getSynfiles(self):

        return [self.synfiles[x] for x in self.synfiles]


    def getSynID(self, x):

        if x in self.synids:
            return self.synids[x]

        return None

    def getSynName(self, x):

        if x in self.synfiles:
            return self.synfiles[x]

        return None