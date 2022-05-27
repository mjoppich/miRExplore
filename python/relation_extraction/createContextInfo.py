import argparse
import glob
import sys
import os

from collections import defaultdict

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from synonymes.GeneOntology import GeneOntology

from synonymes.SynfileMap import SynfileMap
from textmining.SentenceDB import SentenceDB
from textmining.SyngrepHitFile import SyngrepHitFile


from utils.parallel import MapReduce


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='aggregate tm results', add_help=False)
    parser.add_argument('-s', '--sentdir', type=str, help='where are the sentences?', required=True)
    parser.add_argument('-r', '--resultdir', type=str, help='where are all the index-files?', required=True)
    parser.add_argument('-d', '--datadir', type=str, help='where is te miRExplore bsae?', required=True)

    parser.add_argument('--obo', type=argparse.FileType('r'), required=True)
    parser.add_argument('--phit', type=int, default=6, required=False)

    parser.add_argument('--threads', type=int, default=8, required=False)

    parser.add_argument('--accept-pmids', type=argparse.FileType('r'), required=False, default=None)
    parser.add_argument('--mine-path', type=str, default="/mnt/e/data/pmid_jun2020/", required=False)

    parser.add_argument('--sentid-no-text', dest='sentid_no_text', action="store_true", required=False, default=False)

    args = parser.parse_args()


    accept_pmids = None

    if args.accept_pmids != None:

        accept_pmids = set()

        for line in args.accept_pmids:

            line = line.strip()

            if len(line) > 0:
                accept_pmids.add(line)


    #resultBase = dataDir + "/miRExplore/textmine/results_pmc/"
    resultBase = args.resultdir

    oboSyns = SynfileMap(resultBase + "/synfile.map")
    oboSyns.loadSynFiles( (args.mine_path, args.datadir) )


    allfiles = glob.glob(resultBase + "/*.index")
    allfileIDs = [os.path.basename(x).replace(".index", "") for x in allfiles]
    allfileIDs = sorted(allfileIDs, reverse=True)

    #allfileIDs = [894]

    celloObo = GeneOntology(args.obo.name)


    def getTerm(synid, obo):

        if synid in obo.dTerms:
            return obo.getID(synid)

        synid = synid.replace('_', ':', 1)

        return obo.getID(synid)

    def analyseFile(splitFileIDs, env):

        fileCoocs = []


        for splitFileID in splitFileIDs:


            indexFile = resultBase + "/"+splitFileID +".index"

            oboHits = SyngrepHitFile(indexFile, oboSyns, sentIDNoText=args.sentid_no_text)
            if len(oboHits) == 0:
                continue

            # not used anyhow ...
            #sentFile = args.sentdir + "/" + splitFileID + ".sent"
            #sentDB = SentenceDB(sentFile)

            sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

            for docID in oboHits:

                if accept_pmids != None:
                    if not docID in accept_pmids:
                        continue

                docHits = oboHits.getHitsForDocument(docID)

                synid2loc = defaultdict(list)

                allSynIDs = set()
                for hit in docHits:

                    if len(hit.foundSyn) < args.phit:
                        if hit.perfectHit != True:
                            continue

                    if "and " in hit.foundSyn:
                        continue

                    if hit.synonym is None:
                        print(hit, file=sys.stderr)
                        sys.stderr.flush()

                    allSynIDs.add(hit.synonym.id)

                    synid2loc[hit.synonym.id].append(
                        (str(hit.documentID), hit.position[0], hit.position[1], hit.hitSyn)
                    )


                removeIDs = set()
                for synID in allSynIDs:

                    gterm = getTerm(synID, celloObo)
                    if gterm == None:
                        sys.stderr.write("Invalid synID: " + synID)
                        removeIDs.add(synID)
                        continue

                    allChildren = gterm.getAllChildren(maxLevel=2)

                    for x in allChildren:
                        removeIDs.add(x)


                allowedIDs = [x for x in allSynIDs if not x in removeIDs]
                #allowedIDs = allSynIDs.remove(removeIDs)

                allterms = []
                for synID in allowedIDs:

                    gterm = getTerm(synID, celloObo)

                    allterms.append(gterm)

                    fileCoocs.append((docID, gterm.id, gterm.name, synid2loc[synID]))



        sys.stderr.write("Found {cnt} elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs)))

        printed = printStuff(None, fileCoocs, None)

        sys.stderr.write("Found {cnt} (printed: {printed}) elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs), printed=printed))


        return None



    threads = args.threads

    if __debug__:
        threads = 1
        sys.stderr.write("Running on threads:" + str(threads) + "\n")

    sys.stderr.write("Debug Mode? " + str(__debug__) + " and threads " + str(threads) + "\n")


    def printStuff(old, fileCoocs, env):

        printed = 0

        allCoocs = set()
        for cooc in fileCoocs:
            allCoocs.add((cooc[0], cooc[1]))


        for cooc in fileCoocs:

            if (cooc[0], cooc[1] + "_excl") in allCoocs:

                sys.stderr.write("Excluded: " + "{pmid}\t{cl}\t{name}\t{pos}\n".format(
                pmid=cooc[0],
                cl=cooc[1],
                name=cooc[2],
                pos=str(cooc[3])
            ) + " for " + cooc[1] + "_excl")

                continue

            print("{pmid}\t{cl}\t{name}\t{pos}\n".format(
                pmid=cooc[0],
                cl=cooc[1],
                name=cooc[2],
                pos=str(cooc[3])
            ), end='', flush=True)

            printed += 1

        return printed




    ll = MapReduce(threads)
    result = ll.exec( allfileIDs, analyseFile, None, 1, None)

