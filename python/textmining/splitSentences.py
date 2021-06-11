
import argparse
import os,sys
import io

class Sentence:

    def __init__(self, id, text):
        self.sentID = id
        self.text = text


    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "{}\t{}".format(str(self.sentID), self.text)

    @classmethod
    def from_line(cls, line):
        aline = line.split("\t")

        return Sentence(SentenceID.fromStr(aline[0]), aline[1])

class SentenceID:

    def __init__(self):

        self.docID = None
        self.parID = None
        self.senID = None

    def __hash__(self):
        return hash(self.docID) ^ hash(self.parID) ^ hash(self.senID)

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


def load_file(sFileLocation, outpath):

    bname, bext = os.path.splitext(os.path.basename(sFileLocation))

    fOutpath = bname + ".{}.sent"
    

    for encoding in [("utf8", "strict"), ("latin1", "strict"),("utf8", "ignore"), ("latin1", "ignore")]:
        
        curOutFileCount = 0
        curFileDocs = set()

        try:

            outfile = io.open(os.path.join(outpath, fOutpath.format(curOutFileCount)), "w", encoding=encoding[0])

            print("Loading", sFileLocation, "with encoding", encoding)
            with io.open(sFileLocation, 'r', encoding=encoding[0], errors=encoding[1]) as infile:
                idx = 0
                for line in infile:
                    sent = Sentence.from_line(line)

                    curFileDocs.add(sent.sentID.docID)

                    if len(curFileDocs) > args.maxdocs:
                        curFileDocs = set()
                        curOutFileCount += 1
                        outfile = io.open(os.path.join(outpath, fOutpath.format(curOutFileCount)), "w", encoding=encoding[0])

                    print(str(sent), end="", file=outfile)

                    
            return True
        except Exception as e:
            if hasattr(e, 'message'):
                print(e.message)
            else:
                print(e)
            continue



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-s', '--sentfile', type=argparse.FileType('r'), required=True, help='alignment files')
    parser.add_argument('-o', '--outpath', type=str, required=True, help='alignment files')
    parser.add_argument('-n', '--maxdocs', type=int, default=100, help='alignment files')
    args = parser.parse_args()

    if not os.path.isdir(args.outpath):
        raise argparse.ArgumentTypeError("outpath must be valid dir!")


    ret = load_file(args.sentfile.name, args.outpath)

    if not ret:
        raise argparse.ArgumentError("probably encoding error!")