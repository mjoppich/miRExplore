from collections import Counter

from synonymes.SynfileMap import SynfileMap
from textmining.SyngrepHitFile import SyngrepHitFile
from utils.idutils import dataDir, loadExludeWords

resultBase = dataDir + "/miRExplore/textmine/results/"
indexFoundSyns = Counter()
excludedSyns = loadExludeWords()

checkResultsFor = 'disease'
analyseFiles = 100
maxFiles = 892

checkSynsMap = SynfileMap(resultBase + "/"+checkResultsFor+"/synfile.map")
checkSynsMap.loadSynFiles( ('/home/users/joppich/ownCloud/data/', dataDir) )

for splitFileID in range(maxFiles, maxFiles-analyseFiles-1, -1):

    fileID = "{:>4}".format(splitFileID).replace(" ", "0")

    print(fileID)

    indexFile = resultBase + "/"+checkResultsFor+"/medline17n"+fileID+".index"
    foundHits = SyngrepHitFile(indexFile, checkSynsMap)

    for doc in foundHits:

        docHits = foundHits.getHitsForDocument(doc)

        for hit in docHits:
            indexFoundSyns[ hit.hitSyn ] += 1

for (syn, cnt) in indexFoundSyns.most_common(100):
    #if syn in excludedSyns:
    print(str(syn) + " -> " + str(cnt))