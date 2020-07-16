
import os, sys
from collections import defaultdict, Counter

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

from synonymes.mirnaID import miRNA, miRNAPART
from textmining.MirGeneRelCheck import SentenceRelationChecker

from utils.tmutils import normalize_gene_names


scaiBase = "/mnt/d/owncloud/data/miRExplore/scai_corpus/"


if sys.argv[1].upper() == "TRAIN":
    scaiFile = "miRNA_train_fixed.xml"
elif sys.argv[1].upper() == "TEST":
    scaiFile = "miRNA_test_fixed.xml"
else:
    exit(-1)

sentFile = open(sys.argv[2], 'w')
synFile = open(sys.argv[3], 'w')

print(sentFile.name)
print(synFile.name)

normGeneSymbols = normalize_gene_names(path=scaiBase + "/../obodir/" + "/hgnc_no_withdrawn.syn")

relexAccepted = []

with open("relexfiles/scai_"+sys.argv[1].lower()+"_relex.out") as fin:

    wasRelation = False
    curSentID = None
    for line in fin:

        if line.startswith(">"):
            curSentID = line.strip()[1:]

        if line.startswith("#RELATIONS:"):
            wasRelation = True
            continue

        if wasRelation:
            relexAccepted.append(curSentID)

        wasRelation = False

    print("Relex Hits", len(relexAccepted))


        

    """
    >365.0.0
    #HITS:hitno:	type	idlist	start	length	idlist	hittext
    1	protein	365.0.0_Specific_miRNAs	8	7	miR-34a
    2	protein	365.0.0_Genes/Proteins	109	4	bcl2
    #RELATIONS:hitno_from	hitno_to	rule	keywords	characterization
    1	2	nsubj	inhibit	5|5|1|1|3|1
    """


from lxml import etree

correctIdentified = 0
incorrectIdentified = 0
totalChecks = 0
incorrectClass = Counter()

relationNum = 0
elemCaseCounter = Counter()


with open(os.path.join(scaiBase, scaiFile), 'r') as fin:
    tree = etree.parse(fin)
    root = tree.getroot()
    scaiPairs = []

    for elem in root.findall(".//document"):

        pmid = elem.attrib['origId']

        for sentElem in elem:

            allEntities = sentElem.findall(".//entity")
            allPairs = sentElem.findall(".//pair")

            sentText = sentElem.attrib["text"]

            entId2elem = {}

            for entity in allEntities:
                entId = entity.attrib['id']
                entText = entity.attrib['text']
                entType = entity.attrib['type']
                entOffset = tuple([int(x) for x in entity.attrib['charOffset'].split("-")])

                if entType in ["Specific_miRNAs", "Genes/Proteins"]:

                    if "Genes" in entType:
                        if entText in normGeneSymbols:
                            entText = normGeneSymbols[entText]
                        elif entText.upper() in normGeneSymbols:
                            gene = normGeneSymbols[entText.upper()]
                    else:
                        try:
                            entText = miRNA(entText).getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])
                        except:
                            pass

                    entTuple = (entText, entType, (entOffset[0], entOffset[1]+1))
                    entId2elem[entId] = entTuple


                    sentEntText = sentText[entTuple[2][0]:entTuple[2][1]]

            for pair in allPairs:

                validInteraction = pair.attrib['interaction'].lower() == "true"
                pairE1 = pair.attrib['e1']
                pairE2 = pair.attrib['e2']

                #if pairInt == 'true':
                if pairE1 in entId2elem and pairE2 in entId2elem:

                    totalChecks += 1

                    e1 = entId2elem[pairE1]
                    e2 = entId2elem[pairE2]

                    if not e1[1] in ["Specific_miRNAs"]:

                        tmp=e1
                        e1=e2
                        e2=tmp


                    relationID = pair.attrib["id"]
                    relationID = "{}.0.0".format(relationNum)

                    relationNum += 1
                    print(relationID, sentText.strip().rstrip(".").replace("/", ","), sep="\t", file=sentFile)

                    #8652807.2.10    protein 0       6:16202 55-60   tumor
                    e1Word = sentText[e1[2][0]:e1[2][1]]
                    e2Word = sentText[e2[2][0]:e2[2][1]]

                    print(relationID, len(relationID) + e1[2][0]+1, e1[2][1]-e1[2][0], e1Word, "protein", relationID + "_" + e1[1], sep="\t", file=synFile)
                    print(relationID, len(relationID) + e2[2][0]+1, e2[2][1]-e2[2][0], e2Word, "protein", relationID + "_" + e2[1], sep="\t", file=synFile)


                    validInteraction = pair.attrib['interaction'].lower() == "true"
                    acceptInteraction = relationID in relexAccepted

                    elemCase = (acceptInteraction, validInteraction)
                    elemCaseCounter[elemCase] += 1


    def printStats(outfile):
        print("Total:     ", totalChecks, file=outfile)
        print(file=outfile)
        print(file=outfile)
        print(file=outfile)
        print(file=outfile)
        print("T,T", elemCaseCounter[(True, True)], file=outfile)
        print("T,F", elemCaseCounter[(True, False)], file=outfile)
        print("F,T", elemCaseCounter[(False, True)], file=outfile)
        print("F,F", elemCaseCounter[(False, False)], file=outfile)

        precision = elemCaseCounter[(True, True)] / (elemCaseCounter[(True, True)]+elemCaseCounter[(True, False)])
        recall = elemCaseCounter[(True, True)] / (elemCaseCounter[(True, True)]+elemCaseCounter[(False, True)])

        f1 = 2* precision * recall / (precision+recall)

        specificity = elemCaseCounter[(False, False)] / (elemCaseCounter[(True, False)] + elemCaseCounter[(False, False)])

        print("precision", precision, file=outfile)
        print("recall", recall, file=outfile)
        print("specificity", specificity, file=outfile)
        print("f1", f1, file=outfile)

    printStats(sys.stdout)


    sentFile.close()
    synFile.close()
                    