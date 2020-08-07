import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from textmining.SentenceID import SentenceID
from collections import defaultdict
import datetime

class PubmedDateDB:

    def __init__(self):
        self.docid2date = {}

    def get_document_date(self, docid):

        docDate = self.get_document(docid)

        if docDate == None:
            return None

        ddate = list(docDate)
    
        if ddate[0] == 0:
            return None
            
        if ddate[1] == 0:
            ddate[1] = 1
        if ddate[2] == 0:
            ddate[2] = 1
    
        dtDate = datetime.datetime.strptime("{}-{}-{}".format(*ddate), '%Y-%m-%d')

        return dtDate

    def get_document_timestamp(self, docid):

        dt = self.get_document_date(docid)

        if dt == None:
            return None

        return datetime.datetime.timestamp(dt)

    def get_document(self, docid):

        if isinstance(docid, (str, int)):
            docid = str(docid)

        elif isinstance(docid, (SentenceID,)):
            docid = str(docid.docID)

        return self.docid2date.get(docid, None)


    @classmethod
    def loadFromFile(cls, infile):

        ret = PubmedDateDB()

        print("Loading Dates", file=sys.stderr)

        with open(infile, 'r') as fin:

            for line in fin:

                #14172228        1964    6       0
                line = line.strip().split('\t')

                pmid = line[0]
                year = int(line[1])
                month = int(line[2])
                day = int(line[3])

                ret.docid2date[pmid] = (year, month, day)

        print("Loading Dates Finished", file=sys.stderr)
        return ret


if __name__ == '__main__':

    allDatesFile = "/mnt/d/dev/data/pmid_jun2020/aggregated_pmid/allpmids.date"
    sentDB = PubmedDateDB.loadFromFile(allDatesFile)
    date = sentDB.get_document("459579") # 1979 7 0
    print(date)


