from collections import defaultdict
import codecs
from textmining.SyngrepHit import SyngrepHit


class SyngrepHitFile:

    def __init__(self, filename, synfileMap = None):

        self.docid2hits = defaultdict(list)

        with codecs.open(filename, 'rb') as infile:

            for line in infile:

                line = line.decode('latin1')

                hit = SyngrepHit.fromLine(line, synfileMap)

                if hit == None:
                    continue

                docID = hit.documentID.docID

                if docID == None:
                    exit(-1)

                self.docid2hits[docID].append(hit)

    def getHitsForDocument(self, docID):
        return self.docid2hits.get(docID, None)

    def __len__(self):
        return len(self.docid2hits)

    def __contains__(self, item):

        return item in self.docid2hits

    def __iter__(self):
        self.allDocIDs = [x for x in self.docid2hits]
        self.currentDocIdx = 0
        return self

    def __next__(self):

        if self.currentDocIdx < len(self.allDocIDs):
            idx = self.currentDocIdx
            self.currentDocIdx += 1
            return self.allDocIDs[idx]

        raise StopIteration()