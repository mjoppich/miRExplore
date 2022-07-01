
import os, sys
from collections import defaultdict

scaiBase = "/mnt/d/owncloud/data/miRExplore/scai_corpus/"
trainFile = "miRNA-Train-Corpus.xml"
testFile = "miRNA-Test-Corpus.xml"

from lxml import etree

sentNum = 11000

for scaiFiles in [(trainFile, "/corpus_train/ORIG_TRAIN_FILE.TXT"), (testFile, "/corpus_test/ORIG_TEST_FILE.TXT")]:

    scaiFile = scaiFiles[0]
    with open(sys.argv[1] + scaiFiles[1], 'w') as fout:


        with open(os.path.join(scaiBase, scaiFile), 'r') as fin:
            tree = etree.parse(fin)
            root = tree.getroot()
            scaiPairs = []

            for elem in root.findall(".//document"):

                pmid = elem.attrib['origId']

                for sentElem in elem:

                    sentence = sentElem.attrib["text"]

                    allEntities = sentElem.findall(".//entity")
                    allPairs = sentElem.findall(".//pair")

                    entId2elem = {}

                    for entity in allEntities:
                        entId = entity.attrib['id']
                        entText = entity.attrib['text']
                        entType = entity.attrib['type']
                        entOffset = [int(x) for x in entity.attrib['charOffset'].split("-")]

                        if entType in ["Specific_miRNAs", "Genes/Proteins"]:

                            """
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
                                    
                            """

                            entname = "e1" if entType in ["Specific_miRNAs"] else "e2"

                            entId2elem[entId] = (entText, entType, entOffset[0], entOffset[1], entname)


                    for pair in allPairs:

                        pairInt = pair.attrib['interaction'].lower()
                        pairE1 = pair.attrib['e1']
                        pairE2 = pair.attrib['e2']

                        if pairE1 in entId2elem and pairE2 in entId2elem:
                            e1 = entId2elem[pairE1]
                            e2 = entId2elem[pairE2]

                            if not e1[2] < e2[2]:
                                tmpE = e1
                                e1 = e2
                                e2 = tmpE

                            #print(sentence)

                            ent1 = sentence[e1[2]:e1[3]+1]
                            ent2 = sentence[e2[2]:e2[3]+1]

                            assert(e1[4] != e2[4])
                            assert(e1[2] < e2[2])

                            sentNum += 1

                            """
                            if e1[4] == "e1":
                                ent1 = "REGULATOR"
                            else:
                                ent1 = "REGULATEE"

                            if e2[4] == "e1":
                                ent2 = "REGULATOR"
                            else:
                                ent2 = "REGULATEE"
                            """

                            esentence = "".join([sentence[:e1[2]], "<{}>".format(e1[4]), ent1,"</{}>".format(e1[4]), sentence[e1[3]+1:e2[2]], "<{}>".format(e2[4]), ent2,"</{}>".format(e2[4]), sentence[e2[3]+1:]])

                            assert(all(["<e1>" in esentence, "<e2>" in esentence, "</e1>" in esentence, "</e2>" in esentence]))

                            #print(ent1, ent2)
                            print(sentNum, "\"" + esentence + "\"", sep="\t" , file=fout)

                            if pairInt == "true":
                                print("interaction({},{})".format(e1[4], e2[4]), file=fout)
                            else:
                                print("no_interaction", file=fout)

                            if e1[4] == "e1":
                                print("Comment: {}:::{}".format(ent1, ent2), file=fout)
                            else:
                                print("Comment: {}:::{}".format(ent2, ent1), file=fout)
                            print(file=fout)


                        """
                        
                        1	"The system as described above has its greatest application in an arrayed <e1>configuration</e1> of antenna <e2>elements</e2>."
                        Component-Whole(e2,e1)
                        Comment: Not a collection: there is structure here, organisation.
        
                        """