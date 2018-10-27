from textmining.SentenceID import SentenceID


class AllowedSentences:

    def __init__(self, maxSentDist):

        self.maxSentDist = maxSentDist


    def getAllowedSentences(self, evidences, sentDist=None):

        if sentDist == None:
            sentDist = self.maxSentDist

        seenSentences = set()

        for evEv in evidences:

            didSID = evEv[0]

            seenSentences = seenSentences.union(self.getAllowedSentencesForSentID(didSID, sentDist=sentDist))

        return seenSentences

    def getAllowedSentencesForSentID(self, sentID, sentDist=None):

        if sentDist == None:
            sentDist = self.maxSentDist

        allowedSentences = set()

        allowedSentences.add(str(sentID))

        aSent = sentID.split(".")
        aSentNum = int(aSent[-1])

        allowedSentences.add(str(SentenceID.fromArray(aSent)))

        aSent = aSent[0:2]

        allowedSentences.add(str(SentenceID.fromArray([aSent[0], "1", "1"])))

        for i in range(1, sentDist + 1):

            if aSentNum - i < 1 and aSent[1] == "2":
                # TODO know how many sentences this section has ...

                for j in range(1, 100):
                    allowedSentences.add(str(SentenceID.fromArray([aSent[0], "1", str(j)])))

            else:

                allowedSentences.add(str(SentenceID.fromArray(aSent + [aSentNum - i])))


            allowedSentences.add(str(SentenceID.fromArray(aSent + [aSentNum + i])))

        return allowedSentences


    def getDistanceDictBySentID(self, sentID, sentDist=None):

        if sentDist == None:
            sentDist = self.maxSentDist

        allowedSentences = {}

        allowedSentences[sentID] = 0

        aSent = sentID.split(".")
        aSentNum = int(aSent[-1])

        newSentID = str(SentenceID.fromArray(aSent))

        allowedSentences[newSentID] = 0

        aSent = aSent[0:2]

        allowedSentences[str(SentenceID.fromArray([aSent[0], "1", "1"]))] = 0

        for i in range(1, sentDist + 1):

            if aSentNum - i < 1 and aSent[1] == "2":
                # TODO know how many sentences this section has ...

                for j in range(1, 100):
                    allowedSentences[str(SentenceID.fromArray([aSent[0], "1", str(j)]))] = 0

            else:

                allowedSentences[str(SentenceID.fromArray(aSent + [aSentNum - i]))] = -i


            allowedSentences[str(SentenceID.fromArray(aSent + [aSentNum + i]))] = i

        return allowedSentences

if __name__ == '__main__':


    to = AllowedSentences(5)

    print(to.getAllowedSentencesForSentID("59779.2.10"))