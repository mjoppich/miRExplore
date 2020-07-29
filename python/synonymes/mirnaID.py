from collections import OrderedDict, defaultdict
from enum import Enum

import re

from utils.idutils import isNumber, containsDigitAndChars


class miRNAPART(Enum):
    ORGANISM=(0,0)
    MATURE=(1,1)
    ID=(2,2)
    PRECURSOR=(2,3)
    MATURE_SEQS=(3,4)
    ARM=(4,5)

class miRNASynonymeTYPE(Enum):
    MIMAT=0
    MI=1
    FAMILY=2
    MIORG=3
    MIALL=4

class miRNACOMPARISONLEVEL(Enum):
    ORGANISM=[miRNAPART.ORGANISM]
    MATUREID=[miRNAPART.MATURE, miRNAPART.ID]
    PRECURSOR=[miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]


class miRNA:

    @classmethod
    def compositions(cls):

        return {
            miRNASynonymeTYPE.MIALL: [
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.ARM],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.ARM],
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID],
            ],
            miRNASynonymeTYPE.MIMAT: [ [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM] ],
            miRNASynonymeTYPE.MI: [ [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS] ],
            miRNASynonymeTYPE.MIORG: [
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID],

            ],
            miRNASynonymeTYPE.FAMILY: [
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS,miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.MATURE_SEQS, miRNAPART.ARM],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.MATURE_SEQS],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR],
                [miRNAPART.MATURE, miRNAPART.ID, miRNAPART.ARM]
            ],

        }

    def make_strings(self, listOfPartsList):

        allRets = []

        for x in listOfPartsList:
            partRet = self.make_string(x)

            allRets += partRet

        return allRets

    def getPart(self, part, default=None):

        if part in self.parts:
            return self.parts[part]

        return default

    def getStringFromParts(self, partsList, normalized=False):

        collectedParts = {}

        outStr = ""

        for part in partsList:

            if part in self.parts:

                partValue = self.parts[part]

                if normalized and part == miRNAPART.MATURE and self.parts[part].upper() in ['MIR', 'MIRNA', 'MICRORNA']:
                    partValue = 'miR'

                collectedParts[part] = partValue

        hasPrinted = False
        for i in range(0, len(partsList)):

            sep = "-"
            if i == 0 or partsList[i] == miRNAPART.PRECURSOR or not hasPrinted:
                sep = ""

            partVal = collectedParts.get(partsList[i], None)

            if partVal == None:
                continue

            outStr = outStr + sep + collectedParts[partsList[i]]
            hasPrinted = True


        return outStr


    def make_string(self, partsList):


        def makeArmVersions(organismStr, matureStr, idStr, precursorStr, matureSeqStr, armStr, armVersions):
            retVersion = []
            armValues = self.armAlternatives[armStr] if len(armStr) > 0 and armVersions else [armStr]

            delim = '-'

            for selArmStr in armValues:

                joinedFirst = [x for x in [
                    organismStr,
                    matureStr
                ] if len(x) > 0]

                joinedLast = [x for x in [
                    idStr + precursorStr,
                    matureSeqStr
                ] if len(x) > 0]

                joinedFirstStr = delim.join(joinedFirst)
                joinedLastStr = delim.join(joinedLast)

                if selArmStr == '*':
                    joinedLastStr += selArmStr
                elif len(selArmStr) > 0:
                    joinedLastStr += delim + selArmStr

                if len(joinedFirstStr) > 0:
                    retVersion.append(joinedFirstStr + delim + joinedLastStr)
                    retVersion.append(joinedFirstStr + joinedLastStr)
                else:

                    if len(joinedLastStr) > 0:
                        retVersion.append(joinedLastStr)

            return retVersion


        organismStr = self.parts[miRNAPART.ORGANISM] if miRNAPART.ORGANISM in partsList and miRNAPART.ORGANISM in self.parts else ''
        matureStr = self.parts[miRNAPART.MATURE] if miRNAPART.MATURE in partsList and miRNAPART.MATURE in self.parts else ''
        idStr = self.parts[miRNAPART.ID] if miRNAPART.ID in partsList and miRNAPART.ID in self.parts else ''
        precursorStr = self.parts[miRNAPART.PRECURSOR] if miRNAPART.PRECURSOR in partsList and miRNAPART.PRECURSOR in self.parts else ''
        matureSeqStr = self.parts[miRNAPART.MATURE_SEQS] if miRNAPART.MATURE_SEQS in partsList and miRNAPART.MATURE_SEQS in self.parts else ''
        armStr = self.parts[miRNAPART.ARM] if miRNAPART.ARM in partsList and miRNAPART.ARM in self.parts else ''

        matureValues = ['miR','mir', 'MiR', 'Mir']

        createdVersions = []

        addArmVersions = True
        if miRNAPART.ARM not in partsList:
            addArmVersions = False

        if miRNAPART.MATURE not in partsList or matureStr.upper() == 'LET' or not matureStr.upper() in self.validMature:
            matureValues = [matureStr]


        if organismStr == '':
            matureValues += ['microRNA', 'MicroRNA', 'micro-RNA', 'Micro-RNA', 'miRNA', 'MiRNA']

        for selMatureStr in matureValues:
            newVersions = makeArmVersions(organismStr, selMatureStr, idStr, precursorStr, matureSeqStr, armStr, addArmVersions)
            createdVersions += newVersions

            if organismStr == '':
                for suffix in ['-mediated']:
                    altVersions = [x+suffix for x in newVersions]
                    createdVersions += altVersions

        return createdVersions


    def accept(self, other, compLevel=miRNACOMPARISONLEVEL.MATUREID):

        if type(other) == str:
            try:
                other = miRNA(other)
            except:
                return False

        if not isinstance(other, miRNA):
            raise ValueError("Comparison must be either miRNA object or string convertible to miRNA", other)

        for part in compLevel.value:
            if not part in other.part2idx or other.part2idx[part] == None:
                return False

            if not other.getPart(part) == self.getPart(part):
                return False

        return True


    def normalized_str(self):

        thisElems = {}
        for x in self.parts:
            thisElems[x] = self.parts[x]

            if x == miRNAPART.MATURE and thisElems[x].upper() in ['MIR', 'MIRNA']:
                thisElems[x] = 'miR'


        return '-'.join([str(thisElems[x]) for x in thisElems])

    def __str__(self):
        return self.getStringFromParts([miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM])

        #return '-'.join([str(self.parts[x]) for x in self.parts])


    def __repr__(self):
        return self.__str__()

    def __init__(self, mirnaStr):

        self.part2idx = {}
        self.idx2part = {}
        for part in miRNAPART:
            self.part2idx[part] = part.value
            self.idx2part[part.value[1]] = part

        self.mirnaStr = None#mirnaStr
        self.parts = OrderedDict()

        self.armAlternatives = defaultdict(list)

        self.armAlternatives['5p'] = ['5p', 's' , '']
        self.armAlternatives['3p'] = ['3p', 'as', '*']

        self.armAlternatives['s'] = ['5p', 's' , '']
        self.armAlternatives['as'] = ['3p', 'as', '*']

        self.armAlternatives['*'] = ['3p', 'as', '*']

        if mirnaStr != None:
            self.parseMirnaStr(mirnaStr)

    @classmethod
    def parseFromComponents(cls, org=None, mature=None, mirid=None, prec=None, mseq=None, arm=None):

        retMir = miRNA(mirnaStr=None)

        if org != None and len(org) > 0:
            retMir.parts[miRNAPART.ORGANISM] = org

        if mature != None and len(mature) > 0:
            retMir.parts[miRNAPART.MATURE] = mature

        if mirid != None and len(mirid) > 0:
            retMir.parts[miRNAPART.ID] = mirid

        if prec != None and len(prec) > 0:
            retMir.parts[miRNAPART.PRECURSOR] = prec

        if mseq != None and len(mseq) > 0:
            retMir.parts[miRNAPART.MATURE_SEQS] = mseq

        if arm != None and len(arm) > 0:
            retMir.parts[miRNAPART.ARM] = arm

        return retMir
        #miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR, miRNAPART.MATURE_SEQS, miRNAPART.ARM


    def parseMirnaStr(self, mirnaStr):

        if mirnaStr.upper().startswith("MICRO"):

            if mirnaStr.startswith("microRNA-"):
                mirnaStr = mirnaStr.replace("microRNA-", "miR-")
            elif mirnaStr.startswith("microRNA"):
                mirnaStr = mirnaStr.replace("microRNA", "miR-")

            elif mirnaStr.startswith("MicroRNA-"):
                mirnaStr = mirnaStr.replace("MicroRNA-", "miR-")
            elif mirnaStr.startswith("MicroRNA"):
                mirnaStr = mirnaStr.replace("MicroRNA", "miR-")

            elif mirnaStr.startswith("Micro-RNA-"):
                mirnaStr = mirnaStr.replace("Micro-RNA-", "miR-")
            elif mirnaStr.startswith("Micro-RNA"):
                mirnaStr = mirnaStr.replace("Micro-RNA", "miR-")

            elif mirnaStr.startswith("micro-RNA-"):
                mirnaStr = mirnaStr.replace("micro-RNA-", "miR-")
            elif mirnaStr.startswith("micro-RNA"):
                mirnaStr = mirnaStr.replace("micro-RNA", "miR-")

        if mirnaStr.startswith("MiRNA"):
            mirnaStr = mirnaStr.replace("MiRNA", "miR")

        if mirnaStr.startswith("miRNA"):
            mirnaStr = mirnaStr.replace("miRNA", "miR")

        if mirnaStr.upper().startswith("MIR") and len(mirnaStr)>3 and mirnaStr[3].isdigit():
            #print("Going to replace mir ...", mirnaStr)
            mirnaStr = mirnaStr.replace(mirnaStr[0:3], "miR-", 1)

        rePattern = re.search("[miR|let]-([0-9]+-[a-z])([^a-z].+$|$)", mirnaStr)

        if rePattern != None:
            matchStart = rePattern.start(1)
            matchEnd = rePattern.end(1)

            nstr = mirnaStr[0:matchStart]
            nstr += mirnaStr[matchStart:matchEnd].replace("-", "", 1).replace("s", "")
            nstr += mirnaStr[matchEnd:]

            #print("Replacing", mirnaStr, "with", nstr)
            mirnaStr = nstr


        amirna = [x for x in mirnaStr.split("-")]# if len(x) > 0

        if len(amirna) == 1:

            for abbrev in ['MIR', 'MICRORNA', 'MIRNA']:

                abbrevLen = len(abbrev)

                if amirna[0].upper().startswith(abbrev) and len(amirna[0]) > abbrevLen and amirna[0][abbrevLen].isdigit():
                    origval = amirna[0]

                    amirna[0] = abbrev
                    amirna.append(origval[abbrevLen:])

        self.validMature = ['MIR', 'LET', 'MICRORNA']

        if amirna[0].upper() in self.validMature:

            if amirna[0].upper() != 'LET':
                amirna[0] = 'miR'

            self.part2idx.pop(miRNAPART.ORGANISM)

            for key in self.part2idx:
                x = self.part2idx[key]
                self.part2idx[key] = (x[0]-1, x[1]-1)
                self.idx2part[x[1]-1] = key

        mirIdx = self.part2idx[miRNAPART.MATURE][0]
        if len(amirna[mirIdx]) > 3 and amirna[mirIdx].upper()[0:3].upper() in self.validMature:
            amirna.insert(mirIdx+1, amirna[mirIdx][3:])
            amirna[mirIdx] = amirna[mirIdx][:3]

        if not amirna[mirIdx][0:min(3, len(amirna[mirIdx]))].upper() in self.validMature:
            # e.g. in ose-bantam
            self.part2idx.pop(miRNAPART.MATURE)

            for key in self.part2idx:
                x = self.part2idx[key]
                if x[0] > mirIdx:
                    self.part2idx[key] = (x[0] - 1, x[1] - 1)
                    self.idx2part[x[1] - 1] = key
        else:

            if mirIdx+2 < len(amirna) and mirIdx+3 != len(amirna):
                possibleID = containsDigitAndChars(amirna[mirIdx+1])
                possibleMatureSeq = containsDigitAndChars(amirna[mirIdx+2])

                if possibleID and possibleMatureSeq:
                    amirna = amirna[:mirIdx+1] + [amirna[mirIdx+1] + "-" + amirna[mirIdx+2]] + amirna[mirIdx+3:]

        if isNumber(amirna[-1]):
            self.part2idx.pop(miRNAPART.ARM)

        testMatureSeqIdx = self.part2idx[miRNAPART.ID][0]+1
        if testMatureSeqIdx < len(amirna):

            if not isNumber(amirna[testMatureSeqIdx]):

                matureSeqIdx = self.part2idx[miRNAPART.MATURE_SEQS]
                self.part2idx.pop(miRNAPART.MATURE_SEQS)

                for key in self.part2idx:
                    value = self.part2idx[key]

                    if value[0] > matureSeqIdx[0]:
                        self.part2idx[key] = (value[0]-1, value[1]-1)
                        self.idx2part[value[1]-1] = key
        else:

            removeKeys = set()
            for key in self.part2idx:
                idx = self.part2idx[key]

                if idx[0] >= testMatureSeqIdx:
                    removeKeys.add(key)

            for key in removeKeys:
                self.part2idx.pop(key)

        """
        removeKeys = set()
        for key in self.part2idx:
            idx = self.part2idx[key]

            if idx[0] >= len(amirna):
                removeKeys.add(key)

        for key in removeKeys:
            self.part2idx.pop(key)
        """

        self.mirnaStr = mirnaStr
        self.parts = OrderedDict()

        self.armAlternatives = defaultdict(list)

        self.armAlternatives['5p'] = ['5p', 's' , '']
        self.armAlternatives['3p'] = ['3p', 'as', '*']

        self.armAlternatives['s'] = ['5p', 's' , '']
        self.armAlternatives['as'] = ['3p', 'as', '*']

        self.armAlternatives['*'] = ['3p', 'as', '*']

        # fetch organism
        if miRNAPART.ORGANISM in self.part2idx:
            posIdx = self.part2idx[miRNAPART.ORGANISM]
            org = amirna[ posIdx[0] ]
            self.parts[ self.idx2part[posIdx[1]] ] = org

        if miRNAPART.MATURE in self.part2idx:
            posIdx = self.part2idx[miRNAPART.MATURE]
            matureType = amirna[ posIdx[0] ]
            self.parts[ self.idx2part[posIdx[1]] ] = matureType


        if miRNAPART.ID in self.part2idx:
            posIdx = self.part2idx[miRNAPART.ID]
            mirnaid = amirna[ posIdx[0] ]

            if mirnaid[0].isdigit():
                mirnaid = re.sub("[^0-9]", '', mirnaid) # remove any non numbers
            else:
                mirnaid = re.sub("[0-9]", '', mirnaid) # remove any numbers (e.g. ebv-mir-BHRF1-2-5p)

            self.parts[self.idx2part[posIdx[1]]] = mirnaid

        if miRNAPART.PRECURSOR in self.part2idx:
            posIdx = self.part2idx[miRNAPART.PRECURSOR]
            mirnaPrec = amirna[ posIdx[0]]

            if mirnaPrec[0].isdigit():
                mirnaPrec = re.sub("[0-9]", '', mirnaPrec) # remove any non numbers
            else:
                mirnaPrec = re.sub("[^0-9]", '', mirnaPrec)  # remove any numbers (e.g. ebv-mir-BHRF1-2-5p)

            if len(mirnaPrec) > 0:
                self.parts[ self.idx2part[posIdx[1]]] = mirnaPrec

        if miRNAPART.MATURE_SEQS in self.part2idx:
            posIdx = self.part2idx[miRNAPART.MATURE_SEQS]
            matureSeq = amirna[posIdx[0]]
            self.parts[ self.idx2part[posIdx[1]]] = matureSeq

        if miRNAPART.ARM in self.part2idx:
            posIdx = self.part2idx[miRNAPART.ARM]
            armID = amirna[posIdx[0]]
            self.parts[ self.idx2part[posIdx[1]]] = armID

    def printParts(self):
        for part in self.parts:
            print(str(part) + "\t" + str(self.parts[part]))





if __name__ == '__main__':


    def testMIRNA( testName ):


        testMirna = miRNA(testName)
        testMirna.printParts()

        print()
        print()

        testStr = testMirna.getStringFromParts(
                [miRNAPART.ORGANISM, miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR,
                 miRNAPART.MATURE_SEQS,
                 miRNAPART.ARM], normalized=True)

        print(testStr)

        print()
        print()

        allMirnaComps = miRNA.compositions()

        for task in allMirnaComps:

            comps = allMirnaComps[task]
            print(task, len(comps))

            for comp in comps:
                print(str(task) + " " + str(testMirna.make_string(comp)))

        print()
        print()

        print(testMirna.accept("miR-126"))




    testMirna = miRNA("miR-335-g-5p")
    testMirna = miRNA("miR-146-medi")
    testMirna = miRNA("let-7-f-3p")

    testMirna = miRNA("miR-181c-5p")
    testMirna.printParts()

    exit()
    testMIRNA("miR-146-3p")
    testMIRNA("miR-146-a")

    testMIRNA('microRNA-106b')
    #testMIRNA('miRUS')
    testMIRNA('hsa-miR-126a-1-3p')



    #testMIRNA('sv40-miR-S1-5p')
    #testMIRNA('kshv-miR-K12-10a-5p')
    #testMIRNA('ose-bantam-5p')
    #testMIRNA('ebv-mir-BHRF1-2-5p')
    #testMIRNA('ath-MIR399a-5p')