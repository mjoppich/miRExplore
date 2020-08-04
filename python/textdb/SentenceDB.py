import codecs
import glob

import os

from collections import defaultdict, deque

import sys


class SentenceDB:

    def __init__(self):

        self.docid2file= {}

        self.loadedSentFilePath = None
        self.loadedSentFile = None

        self.recent_elements = deque([])

    def add_database(self, odb):

        assert(isinstance(odb, SentenceDB))

        for docid in odb.docid2file:
            self.docid2file[docid] = odb.docid2file[docid]



    def get_sentence(self, docid):

        adocid = docid.split(".")

        filename = self.docid2file.get(adocid[0], None)

        if filename == None:
            return None

        for sent in self.recent_elements:

            if sent[0] == docid:

                self.recent_elements.remove(sent)
                self.recent_elements.append(sent)

                return sent

        if filename != self.loadedSentFilePath:
            loadFile = filename + ".sent"

            sys.stderr.write("Loading file: " + loadFile  + " for docid " + docid + "\n")
            self.loadedSentFile = self.loadFile(loadFile)
            self.loadedSentFilePath = filename

        allDocSents = self.loadedSentFile[adocid[0]]

        for sent in allDocSents:
            if sent[0] == docid:

                self.recent_elements.append( sent )

                if len(self.recent_elements) > 20000:
                    self.recent_elements.popleft() # if too large, pop

                return sent

        return None



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
                    #print("Invalid input line for sent:", line, file=sys.stderr)
                    continue

                sentID = aline[0]
                sentText = aline[1]

                docID = sentID.split(".")[0]

                retObj[docID].append((sentID, sentText))

        return retObj


    @classmethod
    def prepareDB(cls, basepath, outpath):

        docid2filename = defaultdict(set)

        for file in glob.glob(basepath + "/*.sent"):


            filename, fileExt = os.path.splitext(os.path.basename(file))
            with open(file, 'r') as fin:

                print(file)

                for line in fin:

                    aline = line.split("\t")

                    pmid = ".".join(aline[0].split(".")[0:-2])

                    docid2filename[filename].add(pmid)




        with open(outpath, 'w') as fout:

            for file in docid2filename:
                fout.write("{fname}\t{pmids}\n".format(fname=file, pmids=",".join(docid2filename[file])))



    @classmethod
    def loadFromFile(cls, basepath, infile):

        ret = SentenceDB()

        if not os.path.exists(infile):
            print("Building Sentence DB")
            cls.prepareDB(basepath, infile)

            if not os.path.exists(infile):
                print("Error creating sentence DB")
                assert(os.path.exists(infile))


        with open(infile, 'r') as fin:

            for line in fin:

                line = line.strip().split('\t')

                pmids = line[1].split(",")

                for docid in pmids:
                    ret.docid2file[docid] = basepath + "/" + line[0]

        return ret


if __name__ == '__main__':

    sentenceLocation = "/mnt/c/dev/data/pmid/"
    x2sentFile = "/mnt/c/dev/data/pmid2sent"

    #SentenceDB.prepareDB(sentenceLocation, x2sentFile)

    sentDB = SentenceDB.loadFromFile(sentenceLocation, x2sentFile)

    senttxt = sentDB.get_sentence("459579.1.1")
    print(senttxt)


