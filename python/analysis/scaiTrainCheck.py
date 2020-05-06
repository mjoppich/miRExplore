
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

normGeneSymbols = normalize_gene_names(path=scaiBase + "/../obodir/" + "/hgnc_no_withdrawn.syn")


from lxml import etree
import spacy

nlp = spacy.load('/mnt/d/spacy/models/en_core_sci_lg-0.2.4/en_core_sci_lg/en_core_sci_lg-0.2.4/')


relChecker = SentenceRelationChecker(nlp)

correctIdentified = 0
incorrectIdentified = 0
totalChecks = 0
incorrectClass = Counter()

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


                    relRes = relChecker.check_sentence(sentText
                                                        , {"entity_type": "mirna", "entity_type_token": "e1", "entity_location": e1[2]}
                                                        , {"entity_type": "gene", "entity_type_token": "e2", "entity_location": e2[2]}
                                                        , fix_special_chars=False
                                                        )


                    fullsentence = relRes['full_sentence']
                    acceptInteraction = relRes['accept_relation']

                    if not acceptInteraction == validInteraction:

                        incorrectClass[(validInteraction, acceptInteraction)]+=1
                        incorrectIdentified += 1
                        print("Incorrect Sentence Start")
                        print("SCAI:", validInteraction, "MIREXPLORE:", acceptInteraction)
                        print(sentText)
                        print(e1)
                        print(e2)
                        print(fullsentence)

                    
                        relRes2 = relChecker.check_sentence(sentText
                                                           , {"entity_type": "mirna", "entity_type_token": "e1", "entity_location": e1[2]}
                                                           , {"entity_type": "gene", "entity_type_token": "e2", "entity_location": e2[2]}
                                                           , fix_special_chars=False
                                                           , verbose=True)

                        relRes2["full_sentence"] = ""
                        #relRes2["entity1"] = None
                        #relRes2["entity2"] = None
                        print(relRes2)
                        if relRes["accept_relation"] != relRes2["accept_relation"]:
                            print("DiffRes")
                            print(relRes)
                        print("Incorrect Sentence End")
                        print()
                        print()
                    else:
                        correctIdentified += 1

            """
            <entity id="miRNA-corp.d1.s0.e0" text="down-regulation" type="Relation_Trigger" charOffset="23-37"/>
            <entity id="miRNA-corp.d1.s0.e1" text="miRNA-128" type="Specific_miRNAs" charOffset="42-50"/>
            <entity id="miRNA-corp.d1.s0.e2" text="glioma" type="Diseases" charOffset="70-75"/>
            <entity id="miRNA-corp.d1.s0.e3" text="GBM" type="Diseases" charOffset="81-83"/>
            <entity id="miRNA-corp.d1.s0.e4" text="up-regulating" type="Relation_Trigger" charOffset="111-123"/>
            <entity id="miRNA-corp.d1.s0.e5" text="ARP5" type="Genes/Proteins" charOffset="125-128"/>
            <entity id="miRNA-corp.d1.s0.e6" text="ANGPTL6" type="Genes/Proteins" charOffset="131-137"/>
            <entity id="miRNA-corp.d1.s0.e7" text="Bmi-1" type="Genes/Proteins" charOffset="141-145"/>
            <entity id="miRNA-corp.d1.s0.e8" text="GBM" type="Diseases" charOffset="210-212"/>
            <pair id="miRNA-corp.d1.s0.p0" type="Specific_miRNAs-Diseases" interaction="True" e2="miRNA-corp.d1.s0.e2" e1="miRNA-corp.d1.s0.e1"/>
            <pair id="miRNA-corp.d1.s0.p1" type="Specific_miRNAs-Diseases" interaction="True" e2="miRNA-corp.d1.s0.e3" e1="miRNA-corp.d1.s0.e1"/>
            <pair id="miRNA-corp.d1.s0.p2" type="Specific_miRNAs-Genes/Proteins" interaction="True" e2="miRNA-corp.d1.s0.e5" e1="miRNA-corp.d1.s0.e1"/>
            <pair id="miRNA-corp.d1.s0.p3" type="Specific_miRNAs-Genes/Proteins" interaction="True" e2="miRNA-corp.d1.s0.e6" e1="miRNA-corp.d1.s0.e1"/>
            <pair id="miRNA-corp.d1.s0.p4" type="Specific_miRNAs-Genes/Proteins" interaction="True" e2="miRNA-corp.d1.s0.e7" e1="miRNA-corp.d1.s0.e1"/>
            <pair id="miRNA-corp.d1.s0.p5" type="Specific_miRNAs-Diseases" interaction="False" e2="miRNA-corp.d1.s0.e8" e1="miRNA-corp.d1.s0.e1"/>
            """

    def printStats(outfile):
        print("Total:     ", totalChecks, file=outfile)
        print("Correct:   ", correctIdentified, correctIdentified/totalChecks, file=outfile)
        print("Incorrect: ", incorrectIdentified, incorrectIdentified/totalChecks, file=outfile)
        print("classes", incorrectClass, file=outfile)

    printStats(sys.stdout)
    printStats(sys.stderr)