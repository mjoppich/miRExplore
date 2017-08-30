from collections import Counter

from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir, loadExludeWords

resultBase = dataDir + "/miRExplore/textmine/results/"
mirnaSyns = SynfileMap(resultBase + "/mirna/synfile.map")
mirnaSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

hgncSyns = SynfileMap(resultBase + "/hgnc/synfile.map")
hgncSyns.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

hgncFoundSyns = Counter()
mirnaFoundSyns = Counter()
filesDone = 0

excludedSyns = loadExludeWords()

for splitFileID in range(892, 0, -1):

    if filesDone > 100:
        break

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    print(fileID)

    hgncFile = resultBase + "/hgnc/medline17n"+fileID+".index"
    hgncHits = SyngrepHitFile(hgncFile, hgncSyns)

    #mirnaFile = resultBase + "/mirna/medline17n"+fileID+".index"
    #mirnaHits = SyngrepHitFile(mirnaFile, mirnaSyns)

    for doc in hgncHits:

        docHits = hgncHits.getHitsForDocument(doc)

        for hit in docHits:
            hgncFoundSyns[ hit.hitSyn ] += 1


    filesDone += 1

for (syn, cnt) in hgncFoundSyns.most_common(100):
    #if syn in excludedSyns:
    print(str(syn) + " -> " + str(cnt))