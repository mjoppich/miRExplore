import codecs
from collections import defaultdict

import os
from enum import Enum

from textmining.SentenceID import SentenceID

class RegPos(Enum):
    BEFORE=-1
    AFTER=-2
    BEFORE_GM=1
    AFTER_GM=2
    BEFORE_MG=3
    AFTER_MG=4
    BETWEEN=5

class ElemOrder(Enum):
    MG=0
    GM=1

class Sentence:

    def __init__(self, sentid, text):

        self.id = sentid
        self.text = text

    def __str__(self):
        return "{id}\t{txt}".format(id = self.id, txt=self.text)


    def extract_text(self, posM, posG):

        order = None
        text = None

        if posM[1] < posG[0]:
            order = ElemOrder.MG
        elif posG[1] < posM[0]:
            order = ElemOrder.GM
        else:
            order = None

        if order == ElemOrder.MG:
            textBefore = self.text[0:posM[0]]
            text = self.text[ posM[1]:posG[0] ]
            textAfter = self.text[posG[1]:]
        else:
            textBefore = self.text[0:posG[0]]
            text = self.text[ posG[1]: posM[0] ]
            textAfter = self.text[posM[1]:]


        return (textBefore, text, textAfter, order)


class SentenceDB:

    def __init__(self, file):

        self.pubmed2sents = defaultdict(list)
        self.filename = os.path.abspath(file)

        if not os.path.isfile(self.filename):
            raise ValueError("Not a valid filename: " + self.filename)

        self.pubmed2sents = self.loadFile(self.filename)


    def loadFile(self, filename):

        retObj = defaultdict(list)

        with codecs.open(filename, 'r') as infile:

            for line in infile:
                #line = line.decode('latin1')

                line = line.strip()
                if len(line) == 0:
                    continue

                aline = line.split("\t")

                if len(aline) != 2:
                    continue

                sentID = SentenceID.fromStr(aline[0])
                sentText = aline[1]

                retObj[sentID.docID].append(Sentence(sentID, sentText))

        return retObj


    def get_sentences(self, idval, default=None):

        if not idval in self.pubmed2sents:
            return default

        return self.pubmed2sents[idval]

    def get_sentence(self, sentid, default=None):

        allsent = self.get_sentences(sentid.docID)

        if allsent == None:
            return default

        for sent in allsent:

            if sent.id == sentid:
                return sent

        return default
