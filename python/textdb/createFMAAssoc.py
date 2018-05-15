import argparse
import glob
import sys
import os

from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from synonymes.SynfileMap import SynfileMap
from textmining.SentenceDB import SentenceDB
from textmining.SyngrepHitFile import SyngrepHitFile


from utils.parallel import MapReduce


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='aggregate tm results', add_help=False)
    parser.add_argument('-s', '--sentdir', type=str, help='where are the sentences?', required=True)
    parser.add_argument('-r', '--resultdir', type=str, help='where are all the index-files?', required=True)
    parser.add_argument('-d', '--datadir', type=str, help='where is te miRExplore bsae?', required=True)

    args = parser.parse_args()

    #resultBase = dataDir + "/miRExplore/textmine/results_pmc/"
    resultBase = args.resultdir


    diseaseSyns = SynfileMap(resultBase + "/model_anatomy/synfile.map")
    diseaseSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', args.datadir) )


    allfiles = glob.glob(resultBase + "/model_anatomy/*.index")
    allfileIDs = [os.path.basename(x).replace(".index", "") for x in allfiles]
    allfileIDs = sorted(allfileIDs, reverse=True)

    #allfileIDs = [894]

    fmaObo = GeneOntology(args.datadir + "miRExplore/foundational_model_anatomy/fma_obo.obo")


    def getTerm(synid, obo):

        if synid in obo.dTerms:
            return obo.getID(synid)

        synid = synid.replace('_', ':', 1)

        return obo.getID(synid)

    def analyseFile(splitFileIDs, env):

        fileCoocs = []


        for splitFileID in splitFileIDs:


            diseaseFile = resultBase + "/model_anatomy/"+splitFileID +".index"

            diseaseHits = SyngrepHitFile(diseaseFile, diseaseSyns)
            if len(diseaseHits) == 0:
                continue

            sentFile = args.sentdir + "/" + splitFileID + ".sent"#"/mnt/c/dev/data/pmc/allsent/"+splitFileID +".sent"
            sentDB = SentenceDB(sentFile)

            sys.stderr.write("Found something in: " + str(splitFileID) + "\n")

            for docID in diseaseHits:

                docHits = diseaseHits.getHitsForDocument(docID)

                synid2loc = defaultdict(list)

                allSynIDs = set()
                for hit in docHits:

                    if "and " in hit.foundSyn:
                        continue

                    allSynIDs.add(hit.synonym.id)

                    synid2loc[hit.synonym.id].append(
                        (str(hit.documentID), hit.position[0], hit.position[1])
                    )


                removeIDs = set()
                for synID in allSynIDs:

                    gterm = getTerm(synID, fmaObo)
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

                    gterm = getTerm(synID, fmaObo)

                    allterms.append(gterm)

                    fileCoocs.append((docID, gterm.id, gterm.name, synid2loc[synID]))



        sys.stderr.write("Found {cnt} elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs)))

        printed = printStuff(None, fileCoocs, None)

        sys.stderr.write("Found {cnt} (printed: {printed}) elems in files {ids}\n".format(cnt=str(len(fileCoocs)), ids=str(splitFileIDs), printed=printed))


        return None



    threads = 6

    if __debug__:
        threads = 1
        sys.stderr.write("Running on threads:" + str(threads) + "\n")

    sys.stderr.write("Debug Mode? " + str(__debug__) + " and threads " + str(threads) + "\n")


    def printStuff(old, fileCoocs, env):

        printed = 0

        for cooc in fileCoocs:

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

