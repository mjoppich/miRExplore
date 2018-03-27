class SentenceID:

    def __init__(self):

        self.docID = None
        self.parID = None
        self.senID = None

    def __eq__(self, other):

        if not isinstance(other, SentenceID):
            return False

        return other.docID == self.docID and other.parID == self.parID and other.senID == self.senID

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        if self.senID != None:
            return "{docID}.{parID}.{senID}".format(docID=self.docID, parID=self.parID, senID=self.senID)
        else:
            return "{docID}.{parID}".format(docID=self.docID, parID=self.parID)

    @classmethod
    def fromStr(cls, line):

        aline = line.strip().split('.')

        retObj = SentenceID()
        retObj.docID = aline[0]

        if len(aline) > 1:
            retObj.parID = aline[1]

        if len(aline) > 2:
            retObj.senID = aline[2]

        return retObj