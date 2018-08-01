from synonymes.mirnaID import miRNA, miRNAPART


class MI2Mirna:

    def __init__(self):

        self.mi2mirna = {}
        self.mirnaNum2mi = {}


    @classmethod
    def loadFromFile(cls, infile="/mnt/c/ownCloud/data/miRExplore/obodir/mirnas_mirbase.csv"):

        retDB = MI2Mirna()

        with open(infile, 'r') as fin:

            for line in fin:

                line = line.strip()

                if len(line) == 0:
                    continue

                line = line.split("\t")

                miID = line[0]
                mirna = line[1]

                if not (mirna.startswith("mmu") or mirna.startswith("hsa")):
                    continue

                try:


                    miObj = miRNA(mirna)

                    miNum = miObj.getPart(miRNAPART.ID, None)

                    if miNum == None:
                        continue

                    org = mirna[0:3]

                    retDB.mi2mirna[miID] = miNum
                    retDB.mirnaNum2mi[(org, miNum)] = miID


                except:
                    continue

        return retDB


if __name__ == '__main__':

    db = MI2Mirna.loadFromFile()

    print(db.mi2mirna.get('MI0000587', None))
    print(db.mirnaNum2mi.get(("mmu", "103"), None))