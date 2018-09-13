import sys
from collections import defaultdict

class RelexParser:

    def __init__(self):

        self.sent2res = defaultdict(list)

    @classmethod
    def loadFromFile(cls, filename):

        ret = RelexParser()

        with open(filename) as fin:

            relMode = False
            hitMode = False

            curSent = None
            curHits = {}


            for line in fin:
                line = line.strip()

                if line.startswith(">"):
                    curSent = line[1:]
                    hitMode = False
                    relMode = False
                    continue

                if line.startswith("#HITS"):
                    hitMode = True
                    relMode = False
                    continue

                if line.startswith("#RELATIONS"):
                    hitMode = False
                    relMode = True
                    continue

                if hitMode:
                    aline = line.split("\t")

                    #2	protein	0:1	74	10	mast cells
                    curHits[aline[0]] = (aline[2], aline[3], aline[5])


                if relMode:
                    aline = line.split("\t")

                    #1	2	nsubj	increas	5|5|5|1|3|3
                    hit1 = curHits[aline[0]]
                    hit2 = curHits[aline[1]]

                    neutHit = None
                    cellHit = None

                    if hit1[0] == "neutrophil":
                        neutHit = hit1
                        cellHit = hit2

                    if hit2[0] == "neutrophil":

                        if neutHit != None:
                            print("DOUBLE NEUTROPHIL!", file=sys.stderr)
                            print(curSent, "\t".join(hit1+hit2), file=sys.stderr)
                            continue

                        neutHit = hit2
                        cellHit = hit1

                    if neutHit == None or cellHit == None:
                        print("SOME HIT IS NULL", neutHit, cellHit, file=sys.stderr)
                        print(curSent, "\t".join(hit1+hit2), file=sys.stderr)
                        continue

                    curHit = {
                        "sent": curSent,
                        "neutrophil": neutHit,
                        "cell": cellHit
                    }

                    ret.sent2res[curSent].append(curHit)

        return ret


if __name__ == '__main__':

    retDB = RelexParser.loadFromFile("./relex.out")

    print(len(retDB.sent2res), file=sys.stderr)