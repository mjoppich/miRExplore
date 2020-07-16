import csv

from io import StringIO

import re

from utils.idutils import isNumber


class Synonym:

    @classmethod
    def parseFromLine(cls, line):
        aLine = line.split(':', 1)

        retSyn = Synonym(aLine[0].strip())


        if len(aLine) > 1:
            aSyns = aLine[1].split('|')
            for x in [x.strip() for x in aSyns]:
                retSyn.syns.append(x)

        else:
            retSyn.syns.append(aLine[0].strip())

        return retSyn

    def __init__(self, id,idIsSyn=True):

        id = id.strip()
        id = id.replace(':', '_')
        self.id = id
        self.currentIdx = 0

        if idIsSyn:
            self.syns = [self.id]
        else:
            self.syns = []

    def get(self, synIdx):

        if synIdx >= 0 and synIdx < len(self.syns):
            return self.syns[synIdx]

        return None

    def match(self, other):

        if other in self.syns:
            return True

        return False


    def removeCommonSynonymes(self, commonSyns):

        syn2rem = set()
        for x in self:
            if x in commonSyns:
                syn2rem.add(x)

        for x in syn2rem:
            self.removeSyn(x)

    def getSynonymes(self):

        return self.syns

    def addTextSyns(self, synText,captureBrackets=True, addBrackets=True):

        if synText == None or synText == 'None':
            return

        synText = synText.strip()

        if len(synText) == 0:
            return

        vAllWords = self.getAllSplittedSyns(synText, ', ', captureBrackets=captureBrackets, addBrackets=addBrackets)

        for x in vAllWords:
            self.addSyn(x)


    def addSyn(self, newSyn):

        if newSyn == None:
            return

        newSyn = newSyn.strip()

        if len(newSyn) == 0:
            return
        if len(newSyn) < 3:
            return

        if newSyn[0] == newSyn[len(newSyn)-1] and newSyn[0] == '"':
            newSyn = newSyn[1:len(newSyn)-1]

        if len(newSyn) < 3:
            return
        if len(newSyn) == 0:
            return

        if newSyn.startswith('symbol withdrawn, '):
            newSyn = newSyn[len('symbol withdrawn, '):]

        if newSyn.startswith("see "):
            newSyn = newSyn[len("see "):]

        if newSyn.endswith("~withdrawn"):
            newSyn = newSyn[:-len("~withdrawn")]

        if '|' in newSyn:
            newSyn = newSyn.replace('|', '-')

        if newSyn.upper() in ['SYMBOL WITHDRAWN', 'WITHDRAWN', 'PSEUDOGENE', 'RNA', 'ENTRY WITHDRAWN']:
            return

        newSyn = re.sub('\s\s+',' ',newSyn)

        if len(re.sub('[0-9\W]*', '', newSyn)) == 0:
            return

        if newSyn == '1,3':
            print(newSyn)

        if not newSyn in self.syns:
            self.syns.append(newSyn)

    def __contains__(self, item):
        return item in self.syns

    def __iter__(self):

        self.currentIdx = 0
        return self

    def __next__(self):

        if self.currentIdx >= len(self.syns):
            raise StopIteration

        self.currentIdx+=1
        return self.syns[self.currentIdx-1]

    def __str__(self):
        return self.id.replace(':', '_') + ":" + "|".join([x for x in self.syns if len(x) > 0])

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.syns)

    def removeSynUpper(self, dictOfSyns, removeExtraChars=False):

        for excludeName in dictOfSyns:

            listToExclude = dictOfSyns[excludeName]

            syns2remove = set()
            for syn in self.syns:

                usyn = syn
                if not excludeName in ["taxnames", "manual"]:
                    usyn = syn.upper()                   

                if syn in listToExclude or usyn in listToExclude:
                    syns2remove.add(syn)

            for syn in syns2remove:
                self.removeSyn(syn)

    def addAlphaBetaVariants(self):
        self.__searchReplaceVariant('A', 'alpha', [u"\u03B1", '-'+u"\u03B1"])
        self.__searchReplaceVariant('B', 'beta', [u"\u03B2", '-' + u"\u03B2"])

        self.__searchReplaceVariantInner('alpha', 'a', [u"\u03B1"])
        self.__searchReplaceVariantInner('beta','b', [u"\u03B2"])
        self.__searchReplaceVariantInner('gamma','c', [u"\u03B3"])

    def __searchReplaceVariantInner(self, sword, shortw, replaceWith):

        newsyns = set()
        isword = " " + sword + " "

        longHasSWORD = False

        for x in self.syns:
            if isword in x or x.endswith(sword):
                longHasSWORD = True
                break

        if not longHasSWORD:
            return

        for x in self.syns:

            m = re.search(r'{shortw}\d+$'.format(shortw=shortw), x)

            for repWith in replaceWith:

                if isword in x:
                    nsyn = x.replace(sword, repWith)
                    newsyns.add(nsyn)

                    nsyn = x.replace(isword, repWith)
                    newsyns.add(nsyn)

                elif x.endswith(sword):

                    nsyn = x[:-len(sword)] + repWith
                    newsyns.add(nsyn)

                elif m:
                    fnumber = m.group()

                    xo = x[:-len(fnumber)]

                    if len(xo) > 1 and (xo[-1].isdigit() or xo[-1].isalpha()) and xo[-1] not in ["I"]:
                        nnumber = fnumber.replace(shortw, repWith)
                        nsyn = xo + nnumber
                        newsyns.add(nsyn)


        if len(newsyns) > 0:
            print(self.id, newsyns)
            for x in newsyns:
                self.syns.append(x)




    def __searchReplaceVariant(self, endChar, endWord, replaceWith):
        endsWithWord = False

        for x in self.syns:
            if x.endswith(endWord) or x.upper().endswith(endWord.upper()):
                endsWithWord = True
                break

        newsyns = set()
        for syn in self.syns:

            endsWithChar = syn.endswith(endChar)#self.id.endswith(endChar)

            if endsWithChar and endsWithWord:

                for x in replaceWith:
                    aid = syn.rsplit(endChar, 1)
                    aid.append(x)

                    newsyn = "".join(aid)

                    print(self.id + ": " + syn + " ADD SYN " + newsyn)
                    newsyns.add(newsyn)

        for x in newsyns:
            self.syns.append(x)

    def addHyphenVariants(self):

        newSyns = set()
        for x in self.syns:

            if x.startswith("HGNC") or x.startswith("MGI"):
                continue

            m = re.search(r'\d+$', x)

            if m:
                fnumber = m.group()

                xo = x[:-len(fnumber)]
                if not xo[-1].isalpha():
                    continue

                addSyn = xo + "-" + fnumber

                newSyns.add( addSyn )

            m = re.search(r'\d+[AB]$', x)

            if m:
                fnumber = m.group()

                xo = x[:-len(fnumber)]

                if len(xo) == 0:
                    continue

                if not xo[-1].isalpha():
                    continue

                addSyn = xo + "-" + fnumber

                newSyns.add( addSyn )

        for x in newSyns:
            self.syns.append(x)



    def removeNumbers(self):

        i = 0
        while i < len(self.syns):

            if isNumber(self.syns[i]):
                self.removeSyn(i)
            else:
                i+= 1


    def removeSyn(self, syn):

        if type(syn) == str:

            while syn in self.syns:
                idx = self.syns.index(syn)
                self.removeSyn(idx)

        else:

            isyn = int(syn)

            #print("Removing syn", self.syns[isyn])

            if isyn >= 0 and isyn < len(self.syns):
                del self.syns[isyn]

    def splitQuotedDelimited(self, search, delimiter=', ', quotechars=['\"', '\'']):

            allWords = []
            i = 0
            while i <len(search):

                x = search[i]

                if x in quotechars:
                    y = search.find(x, i+1)

                    if y != -1:
                        allWords.append((i,y))
                        i = y

                i += 1


            foundWords = []

            i = -len(delimiter)
            lasti = 0
            lastAdded = 0
            while i < len(search):

                lasti = i+len(delimiter)
                i = search.find(delimiter, lasti)

                if i == -1:
                    i = len(search)
                    break

                if len(allWords) > 0:
                    ignoreDel = False
                    for word in allWords:
                        if word[0] <= i and i <= word[1]:
                            ignoreDel = True
                            break

                    if ignoreDel:
                        continue

                foundWords.append(search[lastAdded:i])
                lastAdded = i+len(delimiter)

            if lastAdded != len(search):
                foundWords.append(search[lastAdded:len(search)])

            return foundWords


    def findNestedMatches(self, search, delimiter=',', quotechars={'(': ')'}):

        nestLevel = 0
        nestStart = 0
        stack = []

        foundNesting = []

        i = 0
        while i < len(search):

            char = search[i]

            if nestLevel == 0:
                nestStart = i

            if char in quotechars:
                nestLevel += 1
                stack.append( quotechars[char] )

            elif len(stack) > 0 and char == stack[-1]:
                stack = stack[0:len(stack)-1]
                nestLevel -= 1

                if nestLevel == 0:

                    nestWord = search[nestStart+1:i]
                    foundNesting.append(nestWord)

            i += 1

        return foundNesting



    def getAllSplittedSyns(self, search, delimiter=', ', quotechars=['\"', '\''], captureBrackets=True, addBrackets=True):

        vAllWords = self.splitQuotedDelimited(search)

        setAllWords = set()

        for word in vAllWords:
            setAllWords.add(word)

            foundBrackets = []

            if captureBrackets:
                foundBracketWords = self.findNestedMatches(word)
                for x in foundBracketWords:

                    foundBrackets.append( "("+x+")" )

                    for y in self.splitQuotedDelimited(x):

                        if y[0]==y[len(y)-1] and len(y) > 0:
                            y = y[1:len(y)-1]

                        if addBrackets:
                            setAllWords.add(y)

            testWord = word
            for bracketWord in foundBrackets:
                testWord = testWord.replace(bracketWord, '').strip()

            if len(testWord) > 0:
                setAllWords.add(testWord)

        return list(setAllWords)

if __name__ == '__main__':

    syn = Synonym('bla')
    syn.getAllSplittedSyns('"Flavin reductase", "biliverdin reductase B (flavin reductase (NADPH))"	SDR43U1	"short chain dehydrogenase/reductase family 43U, member 1", "(flavin reductase (NADPH))"')
    syn.getAllSplittedSyns('"cytochrome P450, subfamily XXIA (steroid 21-hydroxylase, congenital adrenal hyperplasia), polypeptide 2", "cytochrome P450, family 21, subfamily A, polypeptide 2"')
    syn.splitQuotedDelimited('"family with sequence similarity 44, member C", "biorientation of chromosomes in cell division 1 pseudogene"')
    print(syn.splitQuotedDelimited('"non-protein coding RNA 181", "A1BG antisense RNA (non-protein coding)", "A1BG antisense RNA 1 (non-protein coding)"'))
    print(syn.getAllSplittedSyns('"non-protein coding RNA 181", "A1BG antisense RNA (non-protein coding)", "A1BG antisense RNA 1 (non-protein coding)"'))
    syn.splitQuotedDelimited('"C3 and PZP-like, alpha-2-macroglobulin domain containing 9"')

    print(syn.syns)