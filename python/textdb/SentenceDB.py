import codecs
import glob

import os

from collections import defaultdict, deque

import sys


class SentenceDB:

    def __init__(self, basepath):

        self.basepath = basepath
        self.docid2file= {}

        self.loadedSentFilePath = None
        self.loadedSentFile = None

        self.recent_elements = deque([])


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

            sys.stderr.write("Loading file: " + filename  + " for docid " + docid + "\n")
            self.loadedSentFile = self.loadFile(self.basepath + "/" + filename + ".sent")
            self.loadedSentFilePath = filename

        allDocSents = self.loadedSentFile[adocid[0]]

        for sent in allDocSents:
            if sent[0] == docid:

                self.recent_elements.append( sent )

                if len(self.recent_elements) > 10000:
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

                sentID = aline[0]
                sentText = aline[1]

                docID = sentID.split(".")[0]

                retObj[docID].append((sentID, sentText))

        return retObj


    @classmethod
    def prepareDB(cls, basepath, outpath):

        docid2filename = defaultdict(list)

        for file in glob.glob(basepath + "/pubmed18*.title"):


            filename, fileExt = os.path.splitext(os.path.basename(file))
            with open(file, 'r') as fin:

                for line in fin:

                    aline = line.split("\t")

                    pmid = aline[0]

                    docid2filename[filename].append(pmid)




        with open(outpath, 'w') as fout:

            for file in docid2filename:
                fout.write("{fname}\t{pmids}\n".format(fname=file, pmids=",".join(docid2filename[file])))



    @classmethod
    def loadFromFile(cls, basepath, infile):

        ret = SentenceDB(basepath)

        with open(infile, 'r') as fin:

            for line in fin:

                line = line.strip().split('\t')

                pmids = line[1].split(",")

                for docid in pmids:
                    ret.docid2file[docid] = line[0]

        return ret


if __name__ == '__main__':

    #SentenceDB.prepareDB("/home/mjoppich/dev/data/pubmed/", "/home/mjoppich/ownCloud/data/miRExplore/textmine/aggregated_pmid/pmid2sent")

    sentDB = SentenceDB.loadFromFile("/home/mjoppich/dev/data/pubmed/", "/home/mjoppich/ownCloud/data/miRExplore/textmine/aggregated_pmid/pmid2sent")

    senttxt = sentDB.get_sentence("28783539.2.8")

    print(senttxt)


