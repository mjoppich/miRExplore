from synonymes.SynfileMap import SynonymID
from textmining.SentenceID import SentenceID


class SyngrepHit:

    def __init__(self):

        self.documentID = None
        self.synonymID = None
        self.synonym = None

        self.foundSyn = None
        self.hitSyn = None

        self.perfectHit = False



    @classmethod
    def fromLine(cls, line, synfileMap):


        aline = [x.strip() for x in line.split('\t')]

        if len(aline) < 7:
            return None

        retObj = SyngrepHit()
        retObj.documentID = SentenceID.fromStr(aline[0])

        retObj.synonymID = SynonymID.fromStr(aline[1])

        if synfileMap != None:
            retObj.synonym = synfileMap.getSynonyme(retObj.synonymID)

            if retObj.synonym == None:
                retObj.synonym = synfileMap.getSynonyme(retObj.synonymID)

        retObj.foundSyn = aline[2]
        retObj.hitSyn = aline[5]

        if len(aline) < 7:
            retObj.perfectHit = True
        else:
            retObj.perfectHit = aline[6].upper() == 'TRUE'

        return retObj
