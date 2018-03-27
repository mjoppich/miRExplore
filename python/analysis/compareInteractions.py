import pickle
import tempfile

from collections import defaultdict

import os

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from analysis.TargetScanDB import TargetScanDB
from analysis.miRWalk3DB import miRWalk3DB
from analysis.miRecordDB import miRecordDB
from database.Neo4JInterface import neo4jInterface
import networkx as nx

from synonymes.mirnaID import miRNA, miRNAPART
from utils.cytoscape_grapher import CytoscapeGrapher


if __name__ == '__main__':

    selChemokines = ['CXCL9', 'CXCL10', 'CCL22', 'CCL4', 'CXCL13', 'CXCR2', 'CXCL7', 'CXCL5', 'CXCR4', 'CXCL12', 'CCR7']
    selChemokines = ['CXCR2','CCL9', 'CXCL5', 'CXCL1', 'CXCL13', 'CXCL7', 'CXCL12', 'CXCL14', 'CCR7', 'CXCL9', 'CCL2', 'CCR7', 'CCL7', 'CCL3', 'CXCL10', 'CCL22', 'CCL22', 'CCL4', 'CXCR4', 'CX3CL1']

    #test = InteractionRetriever(chemokines=selChemokines)
    #test = InteractionRetriever(chemokines=['CXCR2','CCL9', 'CXCL5', 'CXCL1', 'CXCL13', 'CXCL7', 'CXCL12', 'CXCL14', 'CCR7', 'CXCL9', 'CCL2', 'CCR7', 'CCL7', 'CCL3', 'CXCL10', 'CCL22', 'CCL22', 'CCL4', 'CXCR4', 'CX3CL1'])

    interactions = {
        'CCL9': ['miR-30d-3p', 'miR-3473c', 'let-7g-5p'],
        'CXCL5': ['miR-204-5p', 'let-7g-5p', 'miR-362-3p', 'miR-155-5p'],
        'CXCL1': ['miR-194-2-3p', 'miR-128-3p', 'miR-194-5p', 'miR-199b-5p', 'miR-467g', 'miR-122-5p'],
        'CXCL13': ['miR-122-5p'],
        'CXCL14': ['miR-301b-3p'],
        'CXCR2': ['let-7g-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7a-5p', 'let-7i-5p', 'miR-98-5p'],
        'CXCL7': ['let-7g-5p'],
        'CCL2': ['let-7a-5p', 'let-7b-5p', 'let-7f-5p', 'let-7c-5p', 'let-7g-5p', 'let-7i-5p', 'miR-181a-5p'],
        'CXCL9': ['miR-1935'],
        'CCL3': ['miR-30a-5p','miR-30b-5p','miR-30c-5p','miR-30d-5p','miR-30e-5p'],
        'CCL7': ['miR-181a-5p', 'miR-322-5p', 'miR-29a-5p', 'miR-29b-1-5p'],
        'CCL22': [  'miR-34a-5p'],
        'CXCL10': ['miR-503-3p', 'miR-186-5p'],
        'CCR5': ['miR-186-5p', 'miR-669j', 'miR-21-5p', 'miR-146a-5p', 'miR-150-5p', 'miR-146b-5p', 'miR-669k-3p', 'miR-142-3p', 'miR-34a-5p'],
        'CCL4': ['miR-27b-3p', 'miR-27a-3p', 'miR-21-3p', 'miR-467f'],
        'CX3CL1': ['miR-15a-5p', 'miR-322-5p', 'miR-706', 'miR-762', 'miR-665-3p', 'miR-758-3p', 'miR-381-3p'],
        'CXCR4': ['miR-381-3p', 'miR-21-3p', 'miR-467a-5p', 'miR-467h', 'miR-218-5p', 'miR-1a-3p', 'miR-181d-5p', 'miR-206-3p', 'miR-181b-5p', 'miR-9-5p', 'miR-132-3p', 'miR-25-3p', 'miR-467d-5p', 'miR-669k-3p', 'miR-146b-5p', 'miR-467b-5p', 'miR-467e-5p', 'miR-467f', 'miR-146a-5p'],
        'CCR7': ['let-7g-5p', 'miR-23b-3p', 'miR-669p-5p', 'miR-23a-5p', 'let-7e-5p', 'miR-669l-5p', 'miR-15a-5p', 'miR-467e-5p', 'miR-21-5p', 'miR-16-5p', 'let-7d-5p', 'miR-669n', 'miR-98-5p', 'let-7b-5p', 'let-7a-5p', 'let-7i-5p', 'let-7c-5p', 'miR-15b-5p', 'miR-467h'],
        'CXCL12': [
            'miR-532-5p', 'miR-130b-3p', 'miR-222-3p', 'miR144-3p', 'miR-542-3p', 'miR-149-5p', 'miR-330-3p', 'miR-532-3p', 'miR-3470b', 'miR-125b-5p', 'miR-221-3p', 'miR-19b-3p', 'miR-301b-3p',
            'miR-34b-5p', 'miR-125a-3p', 'miR-126-3p', 'miR-16-1-3p', 'miR-882', 'miR-497-5p', 'miR-26a-5p', 'miR-124-3p', 'miR-26b-5p', 'miR-5620-3p', 'mIR-19a-3p', 'miR-130a-3p', 'miR-690',
            'miR-185-5p', 'miR-31-5p', 'miR-340-5p', 'miR-1843-5p', 'miR-466f-3p', 'miR-301a-3p', 'miR-101a-3p', 'miR-210-3p', 'miR-107-3p', 'miR-706', 'miR-23b-3p', 'miR-146a-5p', 'miR-467f',
            'miR-322-5p', 'miR-15a-5p', 'miR-29b-1-5p', 'let-7e-5p', 'miR-23a-3p', 'miR-338-3p', 'miR-103-3p', 'miR-362-3p', 'let-7g-5p', 'miR-155-5p', 'miR-140-5p', 'miR-122-5p', 'miR-22-3p', 'miR-3470a', 'let-7d-5p'
        ]

    }

    for gene in interactions:

        mirlist = interactions[gene]
        mirids = [miRNA(x) for x in mirlist]
        interactions[gene] = mirids

    graph = nx.Graph()

    foundInteractions = defaultdict(set)


    graphConnections = defaultdict(list)
    with open('/tmp/mirtex/mirel', 'r') as fin:

        for line in fin:

            if 'something' in line:
                continue

            aline = line.split('\t')

            gene = aline[0]

            if not gene in interactions:
                continue

            mirna = aline[1]
            evidences = eval(aline[3])

            for x in evidences:
                graphConnections[(gene, mirna)].append(('PUBMED', x[0], x[1]))

    mirecords = miRecordDB.from_xslx()
    for elem in mirecords.elems:
        graphConnections[(elem[0].upper(), elem[1])].append(('MIRECORD', elem[2]))

    targetscandb = TargetScanDB.from_tsv()
    targetscandb.make_dictionary()



    mirwalk = miRWalk3DB()#.from_xslx(targetGenes=[x for x in interactions], minScore=0.95)
    loadedElems = 0
    for elem in mirwalk.elems:
        if elem[0].upper() in interactions:

            tsres = targetscandb.gene2mirnas[elem[0].upper()]

            accept = False
            if tsres == None or len(tsres) == 0:
                accept=True

            else:
                mirwalk_mirna = miRNA(elem[1].replace('hsa-', '').replace('mmu-', ''))

                for tselem in tsres:
                    if mirwalk_mirna.accept(tselem[1]):
                        accept=True
                        break

            if accept:
                graphConnections[(elem[0].upper(), elem[1])].append(('MIRWALK', str(elem[2])))
                loadedElems += 1



    print("Loaded miRWalk 3:", len(mirwalk.elems), loadedElems)

    pickleFile = '/home/mjoppich/chemokines.upd.graph.pickle'
    if os.path.isfile(pickleFile):

        with open(pickleFile, 'rb') as infile:
            dbconns = pickle.load(infile)

            for elem in dbconns:
                gene = elem[0].upper()
                mirna = elem[1]

                evs = list(elem[2])

                for ev in evs:
                    if ev[0] == 'MIRTARBASE':
                        graphConnections[(gene, mirna)].append(ev)

    gene2keys = defaultdict(list)

    for x in graphConnections:
        #print(x)

        foundInteractions[x[0]].add(x[1].replace('hsa-', '').replace('mmu-', ''))

        gene = x[0]
        mirna = x[1].replace('hsa-', '').replace('mmu-', '')

        gene2keys[gene].append(x)

        graph.add_edge(x[0], x[1], attr_dict={'color': "#FF0000"})

    dirTemp = tempfile.mkdtemp()
    print("Graph Data located in " + dirTemp)

    print()
    print()
    print()

    def findEdgeInfo( geneID, mirnaID, mustBeType = None):

        foundResults = defaultdict(set)

        allkeys = gene2keys[geneID]

        for x in allkeys:
            if x[0] == geneID and mirnaID.accept(x[1]):

                for ev in graphConnections[x]:
                    if ev[0] in mustBeType:
                        foundResults[ev[0]].add(ev[1])

        return  foundResults

    #for x in foundInteractions:
    #        print(x, foundInteractions[x])

    print()
    print()
    print()

    missingInteractions = defaultdict(set)
    additionalInteractions = defaultdict(set)

    foundAcceptedInteractions = defaultdict(set)


    for x in interactions:

        defInts = interactions[x]

        if x in foundInteractions:
            fInts = foundInteractions[x]
        else:
            fInts = set()

        for mirna in defInts:

            mirnaFound = False
            for foundMirna in fInts:
                if mirna.accept(foundMirna):
                    mirnaFound = True
                    break

            if mirnaFound == False:
                missingInteractions[x].add(mirna)
            else:
                foundAcceptedInteractions[x].add(mirna)


        for mirna in fInts:

            mirnaFound = False

            for defMirna in defInts:
                if defMirna.accept(mirna):
                    mirnaFound = True
                    break

            if mirnaFound == False:
                additionalInteractions[x].add(miRNA(mirna))


    missingDF = DataFrame()
    missingDF.addColumns(['chemokine', 'miRNA Group', 'miRNA', 'Original Network', 'PubMed', 'MIRECORD', 'MIRTARBASE', 'MIRWALK'])


    linkedDF = DataFrame()
    linkedDF.addColumns(['chemokine', 'miRNA Group', 'miRNA', 'Original Network', 'PubMed', 'MIRECORD','MIRTARBASE', 'MIRWALK'])

    totalMissing = 0
    print("Missing miRNAs")
    for x in missingInteractions:
        print(x, len(missingInteractions[x]), len(interactions[x]), missingInteractions[x])

        totalMissing += len(missingInteractions[x])

        selInts = missingInteractions[x]

        geneCap = x.upper()
        geneLow = x.lower()

        if not geneLow[0].isdigit():
            geneLow = geneLow.capitalize()

        dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
            geneCap=geneCap, geneLow=geneLow)

        for mirna in selInts:

            addRow = {
                'chemokine': x,
                'miRNA Group': mirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                'miRNA': mirna,
                'Original Network': "False</br>"
                                    "Missing</br>"
                                    "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term="+x+"+"+mirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID])+"\">Search PUBMED</a>"
                                    "</br><a href=\""+dianaLink+"\">Search DIANA</a>"
                ,
                'PubMed': "",
                'MIRECORD': '',
                'MIRTARBASE': "",
                'MIRWALK': ""
            }

            row = DataRow.fromDict(addRow)
            missingDF.addRow(row)
            linkedDF.addRow(row)

    print("Total Missing miRNAs", totalMissing)


    def makeLinks( mylist, lnkFnc):

        return [lnkFnc(x) for x in mylist]

    print("Accepted miRNAs")
    for x in foundAcceptedInteractions:

        geneCap = x.upper()
        geneLow = x.lower()

        if not geneLow[0].isdigit():
            geneLow = geneLow.capitalize()

        dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
            geneCap=geneCap, geneLow=geneLow)

        for interact in foundAcceptedInteractions[x]:
            #print(x, interact, findEdgeInfo(x, interact))

            allEdgeInfos = findEdgeInfo(x, interact, ['PUBMED', 'MIRECORD', 'MIRTARBASE', 'MIRWALK'])

            pubmedEvs = sorted([x for x in allEdgeInfos['PUBMED']])
            mirecordEvs = sorted([x for x in allEdgeInfos['MIRECORD']])
            mirtarbaseEvs = sorted([x for x in allEdgeInfos['MIRTARBASE']])
            mirwalkEvs = set([x for x in allEdgeInfos['MIRWALK']])

            addRow = {
                'chemokine': x,
                'miRNA Group': interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                'miRNA': interact,
                'Original Network': True,
                'PubMed': ", ".join(pubmedEvs),
                'MIRECORD': ", ".join(mirecordEvs),

                'MIRTARBASE': ", ".join(mirtarbaseEvs),
                'MIRWALK': ", ".join(mirwalkEvs)
            }

            row = DataRow.fromDict(addRow)
            missingDF.addRow(row)

            mirtarbaseEvs = makeLinks(mirtarbaseEvs,
                    lambda
                        lid: "<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid=" + lid + "\">" + lid + "</a>")

            if len(mirtarbaseEvs) > 0:
                mirtarbaseEvs = ["<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/search.php?opt=search_box&kw=" + x + "&sort=id\">" + x + "</a>"] + mirtarbaseEvs

            addRow = {
                'chemokine': x,
                'miRNA Group': interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                'miRNA': interact,
                'Original Network': "True</br>"
                                    "Accepted</br>"
                                    "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term="+x+"+"+interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID])+"\">Search PUBMED</a>"
                                    "</br><a href=\""+dianaLink+"\">Search DIANA</a>",
                'PubMed': ", ".join(makeLinks(pubmedEvs, lambda lid: "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/"+lid+"\">"+lid+"</a>")),
                'MIRECORD': ", ".join(makeLinks(mirecordEvs, lambda
                    lid: "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/" + lid + "\">" + lid + "</a>")),

                'MIRTARBASE': ", ".join(mirtarbaseEvs),
                'MIRWALK': ", ".join(mirwalkEvs)
            }

            row = DataRow.fromDict(addRow)
            linkedDF.addRow(row)


    totalAdditional = 0
    print("Additional miRNAs")
    for x in additionalInteractions:

        print(x, len(additionalInteractions[x]), len(interactions[x]), additionalInteractions[x])
        totalAdditional += len(additionalInteractions[x])

        geneCap = x.upper()
        geneLow = x.lower()

        if not geneLow[0].isdigit():
            geneLow = geneLow.capitalize()

        dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
            geneCap=geneCap, geneLow=geneLow)

        for interact in additionalInteractions[x]:
            #print(x, interact, findEdgeInfo(x, interact))

            allEdgeInfos = findEdgeInfo(x, interact, ['PUBMED', 'MIRECORD', 'MIRTARBASE', 'MIRWALK'])

            pubmedEvs = sorted([x for x in allEdgeInfos['PUBMED']])
            mirecordEvs = sorted([x for x in allEdgeInfos['MIRECORD']])
            mirtarbaseEvs = sorted([x for x in allEdgeInfos['MIRTARBASE']])
            mirwalkEvs = set([x for x in allEdgeInfos['MIRWALK']])

            addRow = {
                'chemokine': x,
                'miRNA Group': interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                'miRNA': interact,
                'Original Network': False,
                'PubMed': ", ".join(pubmedEvs),
                'MIRECORD': ", ".join(mirecordEvs),
                'MIRTARBASE': ", ".join(mirtarbaseEvs),
                'MIRWALK': ", ".join(mirwalkEvs)

            }

            row = DataRow.fromDict(addRow)
            missingDF.addRow(row)

            mirtarbaseEvs = makeLinks(mirtarbaseEvs,
                    lambda
                        lid: "<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid=" + lid + "\">" + lid + "</a>")

            if len(mirtarbaseEvs) > 0:
                mirtarbaseEvs = ["<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/search.php?opt=search_box&kw=" + x + "&sort=id\">" + x + "</a>"] + mirtarbaseEvs


            addRow = {
                'chemokine': x,
                'miRNA Group': interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                'miRNA': interact,
                'Original Network': "False</br>"
                                    "Additional</br>"
                                    "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term="+x+"+"+interact.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID])+"\">Search PUBMED</a></br>"
                                    "<a href=\""+dianaLink+"\">Search DIANA</a>",
                'PubMed': ", ".join(makeLinks(pubmedEvs, lambda
                    lid: "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/" + lid + "\">" + lid + "</a>")),
                'MIRECORD': ", ".join(makeLinks(mirecordEvs, lambda
                    lid: "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/" + lid + "\">" + lid + "</a>")),

                'MIRTARBASE': ", ".join(mirtarbaseEvs),
                'MIRWALK': ", ".join(mirwalkEvs)

            }

            row = DataRow.fromDict(addRow)
            linkedDF.addRow(row)

    print("Total Additional miRNAs", totalAdditional)


    missingDF.export("/mnt/c/ownCloud/data/miRExplore/overview_weber_tmn.xlsx", ExportTYPE.XLSX)
    linkedDF.export("/mnt/c/ownCloud/data/miRExplore/overview_weber_tmn.html", ExportTYPE.HTML)

    #CytoscapeGrapher.showGraph(graph, location=dirTemp, name='chemokines', title='Chemokines', nodeLabel=lambda x: x, edgeLabel=lambda x: '')
