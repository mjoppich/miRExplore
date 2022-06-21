import re

from textdb.MiGenRelDB import MiGenRelDB
from utils.tmutils import normalize_gene_names

mainPath = "/mnt/d/owncloud/data/miRExplore/"

normGeneSymbols = normalize_gene_names(path=mainPath + "/obodir/" + "/hgnc_no_withdrawn.syn")

mirelPMIDhsa = MiGenRelDB.loadFromFile(mainPath + "/textmine/aggregated_pmid/"+ "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene",
                                       normGeneSymbols=normGeneSymbols, switchLR=True)


print(mirelPMIDhsa.get_rels("mirna", "miR-758"))

exit(0)



def makeListingGroups(baseHits, conjunts):

    resElems = {}

    for baseHit in baseHits:
        spos = baseHit.start()
        epos =baseHit.end()

        curGroup = []

        for conjElem in conjunts:
            if conjElem.start() == epos:
                curGroup.append(conjElem)
                epos = conjElem.end()


        if len(curGroup) > 0:
            resElems[baseHit] = curGroup

    return resElems


text = "we could show that miR-148/126/340 is important. On the other hand, miR-148/-122/-190 do not seem to play a role"

baseHits = [x for x in re.finditer("miR-[0-9]*", text)]
fconj = [x for x in re.finditer("\/[\-]*([0-9]+)", text)]

print(baseHits)
print(fconj)

resElems = makeListingGroups(baseHits, fconj)

for x in resElems:
    print(x, resElems[x])