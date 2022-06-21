from collections import OrderedDict, defaultdict
from enum import Enum

import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

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

        #if len(amirna) <= mirIdx:
        #    print("Stupid amirna", amirna)

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

    allTestMirs = ['miR-4251', 'miR-551a', 'miR-7846-3p', 'miR-4632-5p', 'miR-4632-3p', 'miR-4684-5p', 'miR-4684-3p', 'miR-1976', 'miR-4420', 'miR-4254', 'miR-12132', 'miR-3116', 'miR-3117-5p', 'miR-3117-3p', 'miR-7156-5p', 'miR-7156-3p', 'miR-137-5p', 'miR-137-3p', 'miR-2682-5p', 'miR-2682-3p', 'miR-190b-5p', 'miR-190b-3p', 'miR-488-5p', 'miR-488-3p', 'miR-181b-3p', 'miR-181a-5p', 'miR-181b-5p', 'miR-181a-3p', 'miR-6740-5p', 'miR-6740-3p', 'miR-6769b-3p', 'miR-6769b-5p', 'miR-29b-2-5p', 'miR-29c-3p', 'miR-29c-5p', 'miR-29b-3p', 'miR-3122', 'miR-4428', 'miR-559', 'miR-5192', 'miR-3126-3p', 'miR-3126-5p', 'miR-3679-3p', 'miR-3679-5p', 'miR-7157-5p', 'miR-7157-3p', 'miR-6888-5p', 'miR-6888-3p', 'miR-4774-3p', 'miR-4774-5p', 'miR-10b-3p', 'miR-10b-5p', 'miR-561-3p', 'miR-561-5p', 'miR-548f-3p', 'miR-6809-3p', 'miR-6809-5p', 'miR-6811-5p', 'miR-6811-3p', 'miR-885-5p', 'miR-885-3p', 'miR-4270', 'miR-3134', 'miR-128-3p', 'miR-128-2-5p', 'miR-26a-1-3p', 'miR-26a-5p', 'miR-564', 'miR-3921', 'miR-567', 'miR-198', 'miR-5682', 'miR-5002-5p', 'miR-6083', 'miR-5002-3p', 'miR-5092', 'miR-5704', 'miR-6827-5p', 'miR-6827-3p', 'miR-3919', 'miR-28-3p', 'miR-28-5p', 'miR-944', 'miR-4274', 'miR-4798-3p', 'miR-4798-5p', 'miR-574-5p', 'miR-574-3p', 'miR-4802-5p', 'miR-4802-3p', 'miR-4449', 'miR-548ah-5p', 'miR-548ah-3p', 'miR-4450', 'miR-4451', 'miR-3684', 'miR-8066', 'miR-1255a', 'miR-3139', 'miR-4799-3p', 'miR-4799-5p', 'miR-4635', 'miR-581', 'miR-5687', 'miR-4803', 'miR-548p', 'miR-143-5p', 'miR-145-3p', 'miR-145-5p', 'miR-143-3p', 'miR-378a-3p', 'miR-378a-5p', 'miR-218-2-3p', 'miR-585-3p', 'miR-585-5p', 'miR-218-5p', 'miR-378e', 'miR-5003-3p', 'miR-5003-5p', 'miR-340-5p', 'miR-340-3p', 'miR-7853-5p', 'miR-4639-3p', 'miR-4639-5p', 'miR-6891-3p', 'miR-6891-5p', 'miR-6835-3p', 'miR-6835-5p', 'miR-30c-2-3p', 'miR-30a-5p', 'miR-30c-5p', 'miR-30a-3p', 'miR-3145-3p', 'miR-3145-5p', 'miR-4648', 'miR-550a-5p', 'miR-550a-3p', 'miR-1200', 'miR-4649-5p', 'miR-4649-3p', 'miR-3914', 'miR-4652-3p', 'miR-4652-5p', 'miR-6875-3p', 'miR-6875-5p', 'miR-6132', 'miR-592', 'miR-153-5p', 'miR-595', 'miR-153-3p', 'miR-3690', 'miR-6086', 'miR-4768-5p', 'miR-4768-3p', 'miR-3915', 'miR-1468-5p', 'miR-1468-3p', 'miR-4329', 'miR-504-3p', 'miR-504-5p', 'miR-1184', 'miR-6842-5p', 'miR-6842-3p', 'miR-6843-3p', 'miR-4662b', 'miR-7848-3p', 'miR-151a-3p', 'miR-151a-5p', 'miR-6845-5p', 'miR-6845-3p', 'miR-6852-3p', 'miR-6852-5p', 'miR-6853-5p', 'miR-6853-3p', 'miR-3153', 'miR-3910', 'miR-4670-3p', 'miR-4670-5p', 'miR-24-3p', 'miR-27b-5p', 'miR-23b-3p', 'miR-2278', 'miR-27b-3p', 'miR-24-1-5p', 'miR-23b-5p', 'miR-6081', 'miR-6854-5p', 'miR-6854-3p', 'miR-32-5p', 'miR-32-3p', 'miR-4668-3p', 'miR-4668-5p', 'miR-4672', 'miR-6855-5p', 'miR-6855-3p', 'miR-126-3p', 'miR-126-5p', 'miR-6722-3p', 'miR-6722-5p', 'miR-7114-3p', 'miR-7114-5p', 'miR-7847-3p', 'miR-675-3p', 'miR-675-5p', 'miR-483-5p', 'miR-483-3p', 'miR-6124', 'miR-4486', 'miR-4688', 'miR-6745', 'miR-3161', 'miR-6746-5p', 'miR-6746-3p', 'miR-6753-3p', 'miR-6753-5p', 'miR-139-5p', 'miR-139-3p', 'miR-326', 'miR-5579-3p', 'miR-5579-5p', 'miR-708-5p', 'miR-708-3p', 'miR-6716-5p', 'miR-6716-3p', 'miR-6756-5p', 'miR-6756-3p', 'miR-100-3p', 'let-7a-5p', 'miR-125b-5p', 'miR-10526-3p', 'miR-125b-1-3p', 'miR-100-5p', 'let-7a-2-3p', 'miR-4697-5p', 'miR-4697-3p', 'miR-5699-3p', 'miR-5699-5p', 'miR-511-5p', 'miR-511-3p', 'miR-8086', 'miR-938', 'miR-604', 'miR-605-5p', 'miR-605-3p', 'miR-7151-5p', 'miR-7151-3p', 'miR-7152-3p', 'miR-7152-5p', 'miR-4676-3p', 'miR-4676-5p', 'miR-3085-5p', 'miR-3085-3p', 'miR-6507-3p', 'miR-6507-5p', 'miR-936', 'miR-4681', 'miR-614', 'miR-6505-3p', 'miR-6505-5p', 'miR-4701-3p', 'miR-4701-5p', 'miR-9898', 'miR-6757-5p', 'miR-6757-3p', 'miR-548c-5p', 'miR-548c-3p', 'miR-617', 'miR-618', 'miR-3922-5p', 'miR-3922-3p', 'miR-6761-3p', 'miR-6761-5p', 'miR-1178-3p', 'miR-1178-5p', 'miR-2681-3p', 'miR-2681-5p', 'miR-4705', 'miR-8073', 'miR-624-5p', 'miR-624-3p', 'miR-7855-5p', 'miR-4709-3p', 'miR-4709-5p', 'miR-4506', 'miR-2392', 'miR-770-5p', 'miR-12121', 'miR-6765-5p', 'miR-6765-3p', 'miR-4715-5p', 'miR-4715-3p', 'miR-1233-5p', 'miR-1233-3p', 'miR-4510', 'miR-10393-5p', 'miR-10393-3p', 'miR-147b-3p', 'miR-147b-5p', 'miR-7973', 'miR-1266-5p', 'miR-1266-3p', 'miR-190a-3p', 'miR-190a-5p', 'miR-12135', 'miR-11181-5p', 'miR-11181-3p', 'miR-3174', 'miR-4714-5p', 'miR-4714-3p', 'miR-6511b-5p', 'miR-1225-3p', 'miR-6511b-3p', 'miR-1225-5p', 'miR-193b-3p', 'miR-365a-5p', 'miR-193b-5p', 'miR-365a-3p', 'miR-484', 'miR-4518', 'miR-6771-5p', 'miR-6771-3p', 'miR-140-3p', 'miR-140-5p', 'miR-6504-3p', 'miR-7854-3p', 'miR-6504-5p', 'miR-8058', 'miR-3182', 'miR-6774-3p', 'miR-6774-5p', 'miR-6775-5p', 'miR-6775-3p', 'miR-22-3p', 'miR-22-5p', 'miR-1180-5p', 'miR-1180-3p', 'miR-4728-5p', 'miR-4728-3p', 'miR-6783-3p', 'miR-6783-5p', 'miR-6784-3p', 'miR-6784-5p', 'miR-152-3p', 'miR-152-5p', 'miR-10a-5p', 'miR-10a-3p', 'miR-3614-3p', 'miR-3614-5p', 'miR-301a-3p', 'miR-454-5p', 'miR-454-3p', 'miR-301a-5p', 'miR-4729', 'miR-548aa', 'miR-4524a-5p', 'miR-4524a-3p', 'miR-6787-5p', 'miR-6787-3p', 'miR-6718-5p', 'miR-4317', 'miR-5190', 'miR-4526', 'miR-3975', 'miR-924', 'miR-5583-5p', 'miR-5583-3p', 'miR-4743-5p', 'miR-4743-3p', 'miR-4320', 'miR-548av-3p', 'miR-548av-5p', 'miR-6870-5p', 'miR-6870-3p', 'miR-3194-5p', 'miR-3194-3p', 'miR-4758-3p', 'miR-4758-5p', 'miR-6813-5p', 'miR-6813-3p', 'miR-637', 'miR-4747-5p', 'miR-4747-3p', 'miR-6885-5p', 'miR-6885-3p', 'miR-23a-3p', 'miR-24-2-5p', 'miR-27a-3p', 'miR-23a-5p', 'miR-27a-5p', 'miR-6795-5p', 'miR-6795-3p', 'miR-6796-5p', 'miR-6796-3p', 'miR-330-5p', 'miR-330-3p', 'miR-4751', 'miR-8061', 'miR-3198', 'miR-648', 'miR-650', 'miR-548j-5p', 'miR-548j-3p', 'miR-4764-3p', 'miR-4764-5p', 'miR-6819-3p', 'miR-6819-5p', 'miR-1249-3p', 'miR-1249-5p', 'miR-4535', 'miR-3667-3p', 'miR-3667-5p', 'miR-6821-3p', 'miR-6821-5p', 'miR-6508-3p', 'miR-6508-5p', 'miR-6814-5p', 'miR-6814-3p', 'miR-6815-3p', 'miR-6815-5p']
    allTestMirs = ['miR-193b']

    for x in allTestMirs:
        testMIRNA(x)

    exit(0)


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