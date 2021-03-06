import glob
import os, sys
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from collections import defaultdict
import shutil

from textmining.SentenceID import SentenceID

if __name__ == '__main__':


    for fileExt in ['.sent', '.author', '.citation', '.title']:

        storagePath = '/mnt/raidtmpbio2/joppich/pmid_jun2020/'

        allfiles = glob.glob(storagePath+'/*'+fileExt)
        allfiles = sorted(allfiles, key=lambda x: os.path.basename(x), reverse=True)

        seenPMIDs = set()
        readInSuffix = "_readin"

        for infileName in allfiles:

            print(infileName)

            foundDocs = defaultdict(list)

            with open(infileName, 'r') as infile:

                for line in infile:
                    aline = line.split('\t')
                    sentID = SentenceID.fromStr(aline[0])

                    try:
                        val = int(sentID.docID)
                    except:
                        val = sentID.docID

                    foundDocs[val].append(line)

            allFilePMIDs = set([x for x in foundDocs])
            pmidIntersect = allFilePMIDs.intersection(seenPMIDs)
            seenPMIDs = seenPMIDs.union(allFilePMIDs)
            print("File PMIDs: {}".format(len(allFilePMIDs)))

            if len(pmidIntersect) > 0:

                print(fileExt + ' : ' + "Removing PMIDs in " + str(infileName) + ": " + str(len(pmidIntersect)))
                docids = sorted([x for x in foundDocs])

                shutil.move(infileName, infileName+readInSuffix)

                with open(infileName, 'w') as outfile:

                    for docID in docids:

                        if docID in pmidIntersect:
                            continue

                        for line in foundDocs[docID]:
                            outfile.write(line)
