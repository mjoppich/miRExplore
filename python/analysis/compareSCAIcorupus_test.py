
import os
from collections import defaultdict

from synonymes.mirnaID import miRNA, miRNAPART

from utils.tmutils import normalize_gene_names


scaiBase = "/mnt/d/owncloud/data/miRExplore/scai_corpus/"
scaiFile = "miRNA-Train-Corpus.xml"

normGeneSymbols = normalize_gene_names(path=scaiBase + "/../obodir/" + "/hgnc_no_withdrawn.syn")


from lxml import etree

def processFile(fin, org):
    foundRels = defaultdict(list)
    for line in fin:

        line = line.strip().split("\t")

        # miR-29-b        Mir-29b MIRNA   GRN     PGRN    GENE    20479936        True    True    [('12', '1V2', 'NEG', 'downregulat', '20479936.2.4', False, (0, 7), (74, 78), (8, 21), 'all_rels', 1, 1, 2, 0)]

        mirna = line[0]
        gene = line[3]
        docid = line[6]

        evs = eval(line[9])

        if gene in normGeneSymbols:
            gene = normGeneSymbols[gene]
        elif gene.upper() in normGeneSymbols:
            gene = normGeneSymbols[gene.upper()]

        for ev in evs:

            relAcc = ev[-4] != 0 or ev[-3] != 0 or ev[-2] != 0

            if relAcc:

                try:
                    rel = {
                        'mirna': miRNA(mirna).getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                        'gene': gene,
                        'docid': docid,
                        'org': org
                    }

                    foundRels[docid].append(rel)

                except:
                    continue

    return foundRels


with open(os.path.join(scaiBase, scaiFile), 'r') as fin:
    tree = etree.parse(fin)
    root = tree.getroot()
    scaiPairs = []

    for elem in root.findall(".//document"):

        pmid = elem.attrib['origId']

        for sentElem in elem:

            allEntities = sentElem.findall(".//entity")
            allPairs = sentElem.findall(".//pair")

            entId2elem = {}

            for entity in allEntities:
                entId = entity.attrib['id']
                entText = entity.attrib['text']
                entType = entity.attrib['type']
                entOffset = entity.attrib['charOffset']

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

                    entId2elem[entId] = (entText, entType, entOffset)


            for pair in allPairs:

                pairInt = pair.attrib['interaction'].lower()
                pairE1 = pair.attrib['e1']
                pairE2 = pair.attrib['e2']

                if pairInt == 'true':
                    if pairE1 in entId2elem and pairE2 in entId2elem:

                        e1 = entId2elem[pairE1]
                        e2 = entId2elem[pairE2]

                        scaiPairs.append(
                            (pmid, e1[0], e2[0])
                        )


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


    docs = set()
    scaiRels = defaultdict(set)
    for x in scaiPairs:
        docs.add(x[0])
        scaiRels[x[0]].add(x)

    with open(scaiBase + "/scaidocs_train", 'w') as fout:
        for x in docs:
            print(x,file=fout )

    foundRels = defaultdict(list)

    with open(scaiBase + "/mirexplore.train.hsa.pmid", "r") as fin:

        fr = processFile(fin, "hsa")

        for x in fr:
            for y in fr[x]:
                foundRels[x].append(y)

    with open(scaiBase + "/mirexplore.train.mmu.pmid", "r") as fin:

        fr = processFile(fin, "mmu")

        for x in fr:
            for y in fr[x]:
                foundRels[x].append(y)


    #invalid docs
    invalidDocs = {
        "18262516": "neither human nor mouse (c. elegans)",
        "20607356": "no hgnc gene names"
    }

    for x in scaiRels:

        if x in invalidDocs:
            continue

        if not x in foundRels:
            print("scai has more", x)

    relevantDocs = set()

    for docid in scaiRels:

        if docid in invalidDocs:
            continue

        scaiEvs = scaiRels[docid]

        mirEvs = foundRels[docid]

        #print(scaiEvs)
        #print(mirEvs)

        scaiInts = set([(x[1], x[2]) for x in scaiEvs])
        mirInts = set([(x['mirna'], x['gene']) for x in mirEvs])

        #print(scaiInts)
        #print(mirInts)

        resInts = (scaiInts.difference(mirInts), mirInts.difference(scaiInts), mirInts.intersection(scaiInts))

        #print(docid, resInts[0], resInts[1], resInts[2], sep="\t")
        print(docid, scaiInts, mirInts, sep="\t")
        relevantDocs.add(docid)

    print("\n".join(relevantDocs))