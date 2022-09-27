import os,sys

"""

PARALLELISM (FORK)

"""

import queue
import threading
import time
import sys
import os
import pickle
import subprocess
import uuid

from pathos import multiprocessing as mp
import glob
import itertools
import uuid

class MyProcessPool(mp.ProcessPool):

    def __init__(self, procs=4):
        super(MyProcessPool, self).__init__(nodes=procs)
        self._id = str(uuid.uuid4())

class MapReduce:

    def __init__(self, procs = 4):
        self.pool = None
        self.nprocs = procs

    def exec(self, oIterable, oFunc, sEnvironment, chunkSize = 1, pReduceFunc = None):

        self.pool = mp.ProcessPool(self.nprocs)
        allResults = []

        resultObj = None

        for x in self.chunkIterable(oIterable, chunkSize):
            allResults.append( self.pool.apipe( oFunc, x, sEnvironment ) )

        self.pool.close()

        while len(allResults) > 0:

            i=0
            while i < len(allResults):

                if allResults[i].ready():

                    result = allResults[i].get()

                    if pReduceFunc != None:

                        resultObj = pReduceFunc(resultObj, result, sEnvironment)

                    else:

                        if resultObj == None:
                            resultObj = []

                        resultObj.append(result)

                    del allResults[i]

                else:
                    i += 1

            time.sleep(0.5)

        self.pool.join()
        self.pool.clear()

        return resultObj



    @classmethod
    def chunkIterable(cls, iterable, size):

        it = iter(iterable)
        item = list(itertools.islice(it, size))

        while item:
            yield item
            item = list(itertools.islice(it, size))


"""

PARALLELISM (THREAD)

"""

class ExtProcThreaded(MapReduce):

    def __init__(self, procs = 4):
        super().__init__(procs=procs)


    def exec(self, oIterable):

        self.pool = mp.ThreadPool(self.nprocs)

        os.makedirs("logs", exist_ok=True)
                
        def call_ext_proc(cmd):
            cmdUUID = uuid.uuid4()
            cmdUUID = str(cmdUUID)[-8:]

            print("Running", cmd, cmdUUID)
            sys.stdout.flush()

            with open("logs/" + cmdUUID + ".out", "w") as outfile, open("logs/" + cmdUUID + ".err", "w") as errfile:
                outfile.write("{}\t{}\n".format(cmdUUID, cmd))
                outfile.flush

                errfile.write("{}\t{}\n".format(cmdUUID, cmd))
                errfile.flush
                subprocess.check_call(cmd, stdout=outfile, stderr=errfile)
                #os.system("sleep 5")
                
            print("Finished", cmdUUID)
            sys.stdout.flush()

            return cmdUUID


        self.pool.map(call_ext_proc, oIterable)
        self.pool.close()
        self.pool.join()
        
        

        return True

"""

MAINTENANCE STUFF


"""

def generate_args_for_subfile(n, newinput, positional=[]):
    """ check an argparse namespace against a module's get_args method.
    Ideally, there would be something built in to argparse, but no such luck.
    This tries to reconstruct the arg list that argparse.parse_args would expect

    Taken from https://stackoverflow.com/questions/36320427/create-arg-string-from-argumentparser-parsed-args-in-python
    """
    arg_list = [[k, v] for k, v in sorted(vars(n).items())]
    argparse_formatted_list = [sys.executable, __file__] #interpreter, current file
    for l in arg_list:
        ####  deal with flag arguments (store true/false)
        if l[1] == True:
            argparse_formatted_list.append("--{}".format(l[0]))
        elif l[1] == False or l[1] is None:
            pass  # dont add this arg
        # add positional argments
        elif l[0] in positional:
            argparse_formatted_list.append(str(l[0]))
        # add the named arguments
        else:

            if l[0] in  ["characters", "c"]:
                continue

            if l[0] == "input":
                l[1] = newinput

            if l[0] == "process_level":
                l[1] = "single"

            if l[0] == "exclude":
                if len(l[1]) == 0:
                    continue

            l[0] = l[0].replace("_", "-")
        
            argparse_formatted_list.append("--{}".format(l[0]))

            if type(l[1]) in [list]:
                for a in l[1]:
                    argparse_formatted_list.append(a)

            else:
                l1str = str(l[1])
                argparse_formatted_list.append(l1str)

    argparse_formatted_list.append("--no-synfile-map")
    
    return argparse_formatted_list

"""

SYNONYM STUFF

"""

class Synonym:

    @classmethod
    def parseFromLine(cls, line, lidx=None):
        aLine = line.split(':', 1)

        retSyn = Synonym(aLine[0].strip())


        if len(aLine) > 1:
            aSyns = aLine[1].split('|')
            for x in [x.strip() for x in aSyns]:
                retSyn.syns.append(x)

        else:
            retSyn.syns.append(aLine[0].strip())

        if lidx != None:
            retSyn.synIdx = lidx

        return retSyn

    def __init__(self, id,idIsSyn=True):

        id = id.strip()
        id = id.replace(':', '_')
        self.id = id
        self.currentIdx = 0
        self.synIdx = -1

        if idIsSyn:
            self.syns = [self.id]
        else:
            self.syns = []

class Synfile:

    def __init__(self, sFileLocation):

        self.mSyns = {}
        self.line2syn = {}
        self.synid2line = {}

        self.synIDs = None
        self.synIDidx = None

        self.location = None

        if sFileLocation != None:
            self._load_file(sFileLocation)

    def _load_file(self, sFileLocation):

        def addSyn(sLine, iLine):

            oSyn = Synonym.parseFromLine(sLine, iLine)

            self.mSyns[ oSyn.id ] = oSyn
            self.line2syn[iLine] = oSyn.id
            self.synid2line[oSyn.id] = iLine

        for encoding in [("utf8", "strict"), ("latin1", "strict"),("utf8", "ignore"), ("latin1", "ignore")]:

            if self.location != None:
                continue

            try:
                print("Loading", sFileLocation, "with encoding", encoding)
                with io.open(sFileLocation, 'r', encoding=encoding[0], errors=encoding[1]) as infile:
                    idx = 0
                    for line in infile:
                        addSyn(line, idx)
                        idx += 1

                    self.location = sFileLocation
            except Exception as e:

                print(e)                
                self.mSyns = {}
                self.line2syn = {}
                self.synid2line = {}

                continue

        if len(self.mSyns) == 0:
            raise ValueError("No Synonymes!")

    def __iter__(self):

        self.synIDs = [x for x in self.mSyns]
        self.synIDidx = 0

        return self

    def __next__(self):

        curIdx = self.synIDidx
        self.synIDidx += 1

        if curIdx < len(self.synIDs):
            return self.mSyns[self.synIDs[curIdx]]

        raise StopIteration()


    def __len__(self):
        return len(self.mSyns)

    def get(self, iSynID):

        return self.mSyns.get(self.line2syn.get(iSynID, None), None)


class StdOut:
    def __enter__(self):
        return sys.stdout

    def __exit__(self, typ, val, trace):
        pass


"""

NER/SYNGREP STUFF

"""


import argparse
import string

import ahocorasick
import io
from collections import defaultdict

if __name__ == '__main__':

#   -np 1 -s .//synonyms/hgnc.syn -i ./pmc/pmc_11.sent -o ./results.pmc.raw//hgnc -nocells -tl 3 -prunelevel none -e excludes/all_excludes.syn

    parser = argparse.ArgumentParser(description='Textmine documents')
    parser.add_argument('-i', '--input', nargs='+', type=str, help='inputfile', required=False)
    parser.add_argument('-f', '--folder', nargs='+', type=str, help='inputfolder', default=None, required=False)

    parser.add_argument('-e', '--exclude', nargs='+', type=str, help='exclude files', required=False, default=[])
    parser.add_argument('-o', '--output', type=str, help='inputfile', required=True)
    parser.add_argument('-s', '--synonyms', nargs='+', type=str, help='syn files', required=True)
    parser.add_argument('-c', '--characters', type=str, required=False, default=' ,.;:-()[]{}=!"§$%&/=?+*\'#-')
    parser.add_argument('-tl', '--trustLength', type=int, required=False, default=4)
    parser.add_argument('-nocells', '--nocells',action="store_true", required=False, default=False)

    parser.add_argument('-nosuffix', '--nosuffix',action="store_true", required=False, default=False)
    parser.add_argument('-olm', '--only-longest-match',action="store_true", required=False, default=False)

    parser.add_argument('-se', '--submatch-exclude',nargs='+', type=argparse.FileType("r"), default=[])
    parser.add_argument('-w', '--test-is-word',action="store_false", required=False, default=True)


    parser.add_argument('-generules', '--generules',action="store_true", required=False, default=False)
    parser.add_argument('-prunelevel', '--prunelevel', type=str, required=False, default="")
    parser.add_argument('-np', '--threads', type=int, required=False, default=8)

    parser.add_argument('-pl', '--process-level', choices=["auto", "single", "threaded", "forked"], required=False, default="auto")
    parser.add_argument('-sm', '--no-synfile-map',action="store_true", required=False, default=False)


    # test run: /usr/bin/python3 /mnt/f/dev/git/miRExplore/python/textmining/textmineDocument.py --input textmining/textmineTestSents.sent --synonyms textmining/textmineTestSyns.syn --output -
    # test run: /usr/bin/python3 /mnt/f/dev/git/miRExplore/python/textmining/textmineDocument.py --input textmineTestSents.sent --synonyms textmineTestSyns.syn --output -
    
    # actual run: /usr/bin/python3 /mnt/f/dev/git/miRExplore/python/textmining/textmineDocument.py --synonyms /mnt/f/dev/data/pmid_jun2020/synonyms/disease.syn --output /tmp/ -tl 5 -prunelevel none -e /mnt/f/dev/data/pmid_jun2020/excludes/all_excludes.syn --input /mnt/f/dev/data/pmid_jun2020/pmc/pmc_118.sent

    args = parser.parse_args()

    #python3 textmineDocument.py --input /mnt/f/dev/data/pmid_jun2020/pmid/pubmed20n0049.sent --synonyms /mnt/f/dev/data/pmid_jun2020/synonyms/hgnc.syn --output - --trustLength 4
    
    if args.input is None and args.folder is None:
        raise ValueError("Either input or folder must be set!")

    if args.folder:
        args.input = []
        for y in args.folder:
            args.input +=  [x for x in glob.glob(y + "/*.sent")]

    if not args.no_synfile_map:
        if args.output != "-":
            assert(os.path.isdir(args.output))
            synfileMapFile = open(os.path.join(args.output, "synfile.map"), 'w')
        else:
            synfileMapFile = sys.stderr
        

    whiteSpaceReplaceChars = args.characters[5:] # leave out blank?   

    synfileMap = {}
    mapSynfile = {}
    for inFile in args.synonyms:
        idx = len(synfileMap)

        inFile = os.path.abspath(inFile)

        synFile = Synfile(inFile)

        synfileMap[idx] = synFile
        mapSynfile[synFile] = idx

        if not args.no_synfile_map:
            print(inFile, idx, sep=": ", file=synfileMapFile)

    if not args.no_synfile_map:
        if args.output != '-':
            synfileMapFile.close()


    def makeUpper(word):
        return word.replace('ß', 'ẞ').upper()

    

    def harmonizeSyn(synword, whiteSpaceReplaceChars):

        newSyn = ""
        repIdx2Orig = {}
        prevCWasWC = False
        for c in synword:

            if c in whiteSpaceReplaceChars:

                if prevCWasWC:
                    continue
                else:
                    newSyn += " "
                    prevCWasWC = True

            else:
                newSyn += c
                prevCWasWC = False

        return newSyn

    A = ahocorasick.Automaton()
    excludeWordAutomaton = ahocorasick.Automaton()

    def addSynwordToAutomaton(A, synWord, data):

        alreadyAdded = A.get(synWord, None)

        if alreadyAdded is None:
            A.add_word( synWord , set([data]) )
        else:
            alreadyAdded.add(data)
            A.add_word( synWord , alreadyAdded )

    excludeWords = set()
    upperExcludeWords = set()

    for exfile in args.submatch_exclude:

        for line in exfile:
            line = line.strip()
            excludeWordAutomaton.add_word(line.upper(), line)


    for exfile in args.exclude:
        with open(exfile, 'r') as fin:

            for line in fin:
                line = line.strip()
                excludeWords.add(line)

    for file_idx in synfileMap:

        synFile = synfileMap[file_idx]

        for synonym in synFile:

            syn_idx = synonym.synIdx

            for synWord in synonym.syns:

                synWord = str(synWord)

                if len(synWord) == 1:
                    # this removes single char synonyms, like a P T ....
                    continue

                # easy cases: word directly
                addSynwordToAutomaton(A, synWord, (file_idx, syn_idx, synWord, synWord, synonym))

                if synWord in excludeWords:
                    continue # we are done here, but added the word for exact matching

                uSyn = makeUpper(synWord)
                if not uSyn in excludeWords:
                    addSynwordToAutomaton(A, uSyn , (file_idx, syn_idx, uSyn, synWord, synonym) )

                # add version without special chars
                fSyn = harmonizeSyn(synWord, whiteSpaceReplaceChars)
                if not fSyn in excludeWords:
                    addSynwordToAutomaton(A, fSyn , (file_idx, syn_idx, fSyn, synWord, synonym) )


                ufSyn = makeUpper(fSyn)
                if not ufSyn in excludeWords:
                    addSynwordToAutomaton(A, ufSyn , (file_idx, syn_idx, ufSyn, synWord, synonym) )


    A.make_automaton()
    excludeWordAutomaton.make_automaton()

    allowedWordChars = set(string.ascii_letters).union({str(x) for x in range(0,10)}).union({x for x in args.characters})

    """
    testSentences = [
        "doc1.2.1\thsa-miR 155 is possibly contained in this sentence of hsa-miR-155, in contrast to HSA----MIR-155.",
        "doc2.2.1\thsa-miR 155 is possibly contained in this sentence, but can it differentiate Ccl2l2 from CCl2?",
        ]

    testResults = [
        ("doc1.2.1", 0,11, 0,0), ("doc1.2.1", 54,11, 0,0), ("doc1.2.1", 82,14, 0,0),
        ("doc2.2.1", 0,11, 0,0), ("doc2.2.1", 77,6, 0,2), ("doc2.2.1", 89,4, 0,1), ("doc2.2.1", 89,4, 0,2),
    ]
    """

    def makeUppercaseSentence(sent, transversion):
        usent = makeUpper(sent)
        if len(usent) != len(sent):
            upIdx2Orig = {}
            usent = ""
            for oidx, c in enumerate(sent):
                uc = c.upper()

                if len(uc) == 1:
                    upIdx2Orig[len(usent)] = transversion[oidx]
                
                elif len(uc) > 1:
                    for i in range(0, len(uc)):
                        upIdx2Orig[len(usent) + i] = transversion[oidx]
                
                usent += uc

            return usent, upIdx2Orig
        else:
            return usent, transversion

    def makePosMapping( sent, whiteSpaceReplaceChars ):

        # removes trailing . from sentence
        sent = sent.rstrip()
        if len(sent) > 0 and sent[-1] in '.:,;':
            sent = sent[:-1]

        posMapped = []
        
        # no changes
        posMapped.append( ( sent, {i:i for i in range(len(sent))} ) )

        #upper case sentence
        nsent, nIdx2Orig = makeUppercaseSentence(sent, {i:i for i in range(len(sent))})
        posMapped.append( (nsent , nIdx2Orig ) )

        newSent = ""
        repIdx2Orig = {}
        prevCWasWC = False
        for oidx, c in enumerate(sent):

            if c in whiteSpaceReplaceChars:

                if prevCWasWC:
                    continue
                else:
                    newSent += " "
                    prevCWasWC = True

            else:
                newSent += c
                prevCWasWC = False
            
            repIdx2Orig[len(newSent)-1] = oidx

        repIdx2Orig[len(newSent)] = len(sent)

        posMapped.append( (newSent, repIdx2Orig) )

        nsent, nIdx2Orig = makeUppercaseSentence(newSent, repIdx2Orig)
        posMapped.append( (nsent , nIdx2Orig ) )

        return posMapped





    def print_results(sentID, foundSyns, foundExcludes, fout):
        for file_idx in foundSyns:

            if args.only_longest_match:

                foundElemsByPos = defaultdict(lambda: defaultdict(set))

                for syn_idx in foundSyns[file_idx]:
                    for (start_index, end_index) in foundSyns[file_idx][syn_idx]:

                        sres = foundSyns[file_idx][syn_idx][(start_index, end_index)]    
                        sres = sorted(sres, key=lambda x: x[5], reverse=False)[0]

                        shouldBeExcluded = False

                        if not foundExcludes is None:
                            shouldBeExcluded = len(foundExcludes[file_idx][syn_idx][start_index]) > 0

                        if len(excludeWordAutomaton) > 0 and shouldBeExcluded:
                            continue
                        
                        (score, text_word, matchWord, synWord, original_value, sentItIdx, prefix, suffix) = sres
                        lres = list(sres)
                        lres.append(file_idx)
                        lres.append(syn_idx)

                        foundElemsByPos[start_index][start_index + len(text_word)].add( tuple(lres) )

                for start_index in foundElemsByPos:

                    max_end_index = max(foundElemsByPos[start_index].keys())
                    
                    for sres in foundElemsByPos[start_index][max_end_index]:
                        (score, text_word, matchWord, synWord, original_value, sentItIdx, prefix, suffix, file_idx, syn_idx) = sres

                        outline = "{sentid}\t{listid}:{synid}\t{matched}\t{start}\t{length}\t{syn}\t{exact}\t{prefix}\t{suffix}\t{sentit}".format(
                            sentid=sentID, listid=file_idx, synid=syn_idx,
                            matched=text_word, start=start_index, length=len(text_word), #end_index-start_index+1
                            syn=synWord, exact=str(text_word==synWord).lower(), prefix=prefix, suffix=suffix, sentit=sentItIdx)

                        print(outline, file=fout)

            else:


                for syn_idx in foundSyns[file_idx]:

                    for (start_index, end_index) in foundSyns[file_idx][syn_idx]:    
                        sres = foundSyns[file_idx][syn_idx][(start_index, end_index)]    
                        sres = sorted(sres, key=lambda x: x[5], reverse=False)[0]
                        (score, text_word, matchWord, synWord, original_value, sentItIdx, prefix, suffix) = sres

                        shouldBeExcluded = False

                        if not foundExcludes is None:
                            shouldBeExcluded = len(foundExcludes[file_idx][syn_idx][start_index]) > 0

                        if len(excludeWordAutomaton) > 0 and shouldBeExcluded:
                            continue

                        # PMC2663906.3.264	0:33862	LRP	61	3	LRP	true	(	−
                        outline = "{sentid}\t{listid}:{synid}\t{matched}\t{start}\t{length}\t{syn}\t{exact}\t{prefix}\t{suffix}\t{sentit}".format(
                            sentid=sentID, listid=file_idx, synid=syn_idx,
                            matched=text_word, start=start_index, length=len(text_word), #end_index-start_index+1
                            syn=synWord, exact=str(text_word==synWord).lower(), prefix=prefix, suffix=suffix, sentit=sentItIdx)

                        print(outline, file=fout)



    def textmineFile(filenames, env):

        for filename in filenames:
            
            outfile = os.path.join(args.output, os.path.splitext(os.path.basename(filename))[0] + ".index")
            print(filename, outfile)

            stdout_copy= StdOut()


            with io.open(filename, 'r') as infile, (io.open(outfile, 'w') if args.output != "-" else stdout_copy) as fout:

                for line in infile:

                    line = line.strip()
                    
                    if len(line) == 0:
                        continue

                    aline = line.split("\t")

                    if len(aline) != 2:
                        continue

                    sentID, sentText = aline[0], aline[1]

                    allTestSentences = makePosMapping(sentText, whiteSpaceReplaceChars)


                    foundSyns = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

                    foundExcludes = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

                    for sentItIdx, sentWIdx in enumerate(allTestSentences):
                        sent, sent2idx = sentWIdx

                        if len(excludeWordAutomaton) > 0:
                            for end_index, matchWord in excludeWordAutomaton.iter(sent):
                                start_index = end_index - len(matchWord) + 1
                                orig_start_index = sent2idx[start_index]
                                orig_end_index = sent2idx[end_index]
                                foundExcludes[file_idx][syn_idx][orig_start_index].add( (matchWord, ) )

                        for end_index, hitKWs in A.iter(sent):
                            
                            for (file_idx, syn_idx, matchWord, synWord, original_value) in hitKWs:

                                start_index = end_index - len(matchWord) + 1
                                test_idx_suffix = end_index+1
                                test_idx_prefix = start_index-1

                                if not start_index in sent2idx or not end_index in sent2idx:
                                    print(sentID)
                                    print(sentText)
                                    print(sent, start_index, end_index)

                                orig_start_index = sent2idx[start_index]
                                orig_end_index = sent2idx[end_index]

                                mod_text_word = sent[start_index:end_index+1]
                                assert(mod_text_word == matchWord)
                                text_word = sentText[orig_start_index:orig_end_index+1]

                                if len(synWord) <= args.trustLength:
                                    if not text_word == synWord:
                                        continue

                                accept = True

                                if args.test_is_word:

                                    if test_idx_suffix < len(sent):
                                        testChar = sent[test_idx_suffix]
                                        accept =  testChar in args.characters
                                    else:
                                        accept = True

                                    if test_idx_prefix >= 0:
                                        testChar = sent[test_idx_prefix]
                                        accept = accept and testChar in args.characters
                                    else:
                                        accept = accept and True

                                prefix = ""
                                suffix = ""



                                addres = (text_word, matchWord, synWord, original_value, sentItIdx)
                                #print(accept, file_idx, syn_idx, (start_index, end_index), addres)
                                if accept:

                                    if test_idx_suffix < len(sent):
                                        suffix = sent[test_idx_suffix: sent.find(" ", test_idx_suffix)]

                                    if args.nosuffix and len(suffix) > 0:
                                        continue

                                    score = 0

                                    if sent == sentText:
                                        score += 1
                                    if matchWord == synWord:
                                        score += 1

                                    addres = (score, text_word, matchWord, synWord, original_value, sentItIdx, prefix, suffix)

                                    foundSyns[file_idx][syn_idx][(orig_start_index, orig_end_index)].add( addres )



                    print_results(sentID, foundSyns, foundExcludes, fout)



    



    print("Starting analysis with TL", args.trustLength)

    processLevel = args.process_level #["auto", "single", "threaded", "forked"]

    if processLevel == "auto":

        fileSizes = []
        for infile in args.input:

            fileSizes.append( os.path.getsize(infile)/(1024*1024) ) # bytes -> kilo-bytes -> mega-bytes

            if len(fileSizes) > 10:
                break

        meanFileSize = sum(fileSizes) / len(fileSizes)
        reason = ""

        if len(args.input) == 1:
            processLevel = "single"
            reason = "Single File"

        else:
            if meanFileSize > 1000: #(~1GB)
                processLevel = "threaded"
                reason = "Mean File Size large"
            else:
                processLevel = "forked"
                reason = "Mean File Size small"

        print("With a mean file size of {}MB in first {} files I choose {} processing for reason: {}".format(meanFileSize, len(fileSizes), processLevel, reason))



    if processLevel == "single":

        print("Single File Processing")

        for inputFile in args.input:
            textmineFile([inputFile], None)

    elif processLevel == "threaded":

        print("Threaded Processing")

        ll = ExtProcThreaded(args.threads)
        allInputs = [generate_args_for_subfile(args, [x]) for x in args.input]      
        ll.exec(allInputs)

    elif processLevel == "forked":
        print("Forked Processing")
        ll = MapReduce(args.threads)
        _ = ll.exec( args.input, textmineFile, None, 1, None)


    print("Done")


#" ".join(['/usr/bin/python3', '/mnt/raidtmp/joppich/pubmed_pmc/pmc/miRExplore/python//textmining/textmineDocument.py', '--exclude', 'excludes/all_excludes.syn', '--input', '../oa_bulk/extracted_sentences/oa_noncomm_xml_PMC002xxxxxx.sent', '--nocells', '--output', './results.pmc.raw//ncit', '--process-level', 'single', '--prunelevel', 'none', '--synonyms', './/synonyms/ncit.syn', '--threads', '14', '--trustLength', '5'])
# /usr/bin/python3 /mnt/raidtmp/joppich/pubmed_pmc/pmc/miRExplore/python//textmining/textmineDocument.py --exclude excludes/all_excludes.syn --input ../oa_bulk/extracted_sentences/oa_noncomm_xml_PMC002xxxxxx.sent --nocells --output ./results.pmc.raw//ncit --process-level single --prunelevel none --synonyms .//synonyms/ncit.syn --threads 14 --trustLength 5