import pickle
import tempfile

from collections import defaultdict, Counter
import networkx as nx

from networkx.drawing.nx_agraph import graphviz_layout, pygraphviz_layout
import os, sys

from synonymes.GeneOntology import GeneOntology
from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname("/mnt/d/dev/git/poreSTAT/")))

from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE

from synonymes.mirnaID import miRNA, miRNAPART, miRNACOMPARISONLEVEL
from textdb.makeNetworkView import DataBasePlotter
from utils.cytoscape_grapher import CytoscapeGrapher

import matplotlib.pyplot as plt
from natsort import natsorted
if __name__ == '__main__':

    cellObo = GeneOntology("/mnt/d/owncloud/data/miRExplore/obodir/meta_cells.obo")

    cellTypeName2Terms = {
        "EC": ["META:52"],
        "MC": ["META:148", "META:99"],
        "FC": ["CL:0000891"],
        "SMC": ["META:83"],
    }

    cellType2AccTerms = {}
    for cellT in cellTypeName2Terms:

        cellType2AccTerms[cellT] = set()

        for et in cellTypeName2Terms[cellT]:

            oboT = cellObo.getID(et)

            if oboT != None:
                cellType2AccTerms[cellT].add(et)
                for x in oboT.getAllChildren():
                    cellType2AccTerms[cellT].add(x.termid)
                    print(cellT, x.term.name)

            else:
                print("No such obo term:", et)

    for ct in cellType2AccTerms:
        print(ct, len(cellType2AccTerms[ct]))


    networks = {}

    # endothelial cell activation

    targetMirsECA = [
        'miR-21',
        'miR-92a',
        'miR-217',
        'miR-663',
        'miR-712',
        'miR-7g',
        'let-7g',
        'miR-10a',
        'miR-17-3p',
        'miR-31',
        'miR-124a',
        'miR-125',
        'miR-126',
        'miR-126-5p',
        'miR-143',
        'miR-145',
        'miR-146',
        'miR-155',
        'miR-181b',
        'miR-221',
        'miR-222']

    networks['targetMirsECA'] = targetMirsECA

    # monocyte

    targetMirsMonocyte = [
        'miR-222',
        'miR-323',
        'miR-503',
        'miR-125b',
        'miR-155',
        'miR-342-5p',
        'miR-17',
        'miR-20a',
        'miR-106a',
        'miR-9',
        'miR-21',
        'miR-124',
        'miR-125a-5p',
        'miR-146a',
        'miR-146b',
        'miR-147',
        'miR-223']

    networks['targetMirsMonocyte'] = targetMirsMonocyte

    # foam cell formation
    targetMirsFCF = [
        'miR-9',
        'miR-125a-5p',
        'miR-146a-5p',
        'miR-155'
    ]

    networks['targetMirsFCF'] = targetMirsFCF


    # Angiogenesis
    targetMirsAngio = [
        'let-7f',
        'miR-7f',
        'miR-23',
        'miR-24',
        'miR-27',
        'miR-126',
        'miR-130a',
        'miR-132',
        'miR-150',
        'miR-210',
        'miR-218',
        'miR-378',
        'miR-15b',
        'miR-16',
        'miR-20a',
        'miR-21',
        'miR-26a',
        'miR-17',
        'miR-92',
        'miR-100',
        'miR-200',
        'miR-221',
        'miR-222',
        'miR-223']

    networks['targetMirsAngio'] = targetMirsAngio


    # Vascular remodeling
    targetMirsVasRemod = [
        'miR-21',
        'miR-155',
        'miR-222',
        'miR-126',
        'miR-143',
        'miR-145']

    networks['targetMirsVasRemod'] = targetMirsVasRemod


    # T - cell differentiation and activation

    targetMirsTCell = [
        'miR-17',
        'miR-92',
        'miR-146a',
        'miR-155',
        'miR-182',
        'miR-326',
        'miR-125b',
        'miR-181a']

    networks['targetMirsTCell'] = targetMirsTCell


    # Cholestrol efflux

    targetMirsCholEfflux = [
        'miR-10b',
        'miR-26',
        'miR-27',
        'miR-33a',
        'miR-106b',
        'miR-144',
        'miR-145',
        'miR-155',
        'miR-302a',
        'miR-758',
        'miR-223',
        'miR-378']


    networks['targetMirsCholEfflux'] = targetMirsCholEfflux


    # SMC proliferation / migration
    targetMirsSMCProlif = [
        'miR-24',
        'miR-26a',
        'miR-31',
        'miR-146a',
        'miR-155',
        'miR-208',
        'miR-221',
        'miR-222',
        'miR-7d',
        'let-7d',
        'miR-1',
        'miR-10a',
        'miR-21',
        'miR-29',
        'miR-100',
        'miR-132',
        'miR-133',
        'miR-143',
        'miR-145',
        'miR-195',
        'miR-204',
        'miR-424',
        'miR-638',
        'miR-663']

    networks['targetMirsSMCProlif'] = targetMirsSMCProlif

    summaryDF = DataFrame()
    summaryDF.addColumns(["Network", "Accepted miRNAs", 'Additional miRNAs', "Missing miRNAs"])

    networkGraphs = {}
    makeStory = [

    ]

    allNetworks = [x for x in networks]
    print(allNetworks)

    #exit()

    ignoreNetworks = []

    networkRestrictions = {
        'targetMirsECA': {
            "cells": [
                {"group": "cells", "name": "endothelial cell", "termid": "META:52"}
            ]
            #, "go": [{"group": "go", "name": "", "termid": "GO:0006915"},{"group": "go", "name": "", "termid": "GO:0001775"},{"group": "go", "name": "", "termid": "GO:0006954"}]
        },
        'targetMirsMonocyte': {
            "cells": [
                {"group": "cells", "name": "monocyte", "termid": "META:148"},
                {"group": "cells", "name": "macrophage", "termid": "META:99"}
            ]
            #, "go": [{"group": "go", "name": "", "termid": "GO:0030224"}, {"group": "go", "name": "", "termid": "GO:0042116"}]
        },
        'targetMirsFCF': {
            "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"}]
            #, "go": [{"group": "go", "name": "", "termid": "GO:0090077"}]
        },
        'targetMirsAngio': {
            #"cells": [{"group": "cells", "name": "blood vessel", "termid": "UBERON:0001981"}, {"group": "cells", "name": "blood vessel elastic tissue", "termid": "UBERON:0003614"} , {"group": "cells", "name": "arterial blood vessel", "termid": "UBERON:0003509"}],
            "go": [{"group": "go", "name": "angiogenesis", "termid": "GO:0001525"}]
        },
        'targetMirsVasRemod': {
            #"disease": [], #{"group": "disease", "name": "vascular disease", "termid": "DOID:178"}
            "go": [
                    {"group": "go", "name": "tissue remodeling", "termid": "GO:0048771"},
                    {"group": "go", "name": "regulation of tissue remodeling", "termid": "GO:0034103"},
                {"group": "go", "name": "regulation of blood vessel remodeling", "termid": "GO:0060312"}
                   ]
        }
        ,'targetMirsTCell': {
            "cells": [{"group": "cells", "name": "T cell", "termid": "META:44"}],
            #,"go": [{"group": "go", "name": "", "termid": "GO:0030217"}]
        },
        'targetMirsCholEfflux': {
            "cells": [{"group": "cells", "name": "foam cell", "termid": "CL:0000891"}]
            #,"go": [{"group": "go", "name": "", "termid": "GO:0033344"}]
        },
        'targetMirsSMCProlif': {
            "cells": [{"group": "cells", "name": "smooth muscle cell", "termid": "META:83"}]
            #,"go": [{"group": "go", "name": "", "termid": "GO:0048659"},{"group": "go", "name": "", "termid": "GO:0014909"}]
        }

    }

    networkToTitle = {
        "targetMirsECA": "Endothelial cell activation\\\\and inflammation",
        "targetMirsMonocyte": "Monocyte differentiation\\\\Macrophage activation",
        "targetMirsFCF": "Foam cell formation",
        "targetMirsAngio": "Angiogenesis",
        "targetMirsVasRemod": "Vascular remodeling",
        "targetMirsTCell": "T-cell differentiation\\\\and activation",
        "targetMirsCholEfflux": "Cholesterol efflux",
        "targetMirsSMCProlif": "SMC proliferation\\\\SMC migration"

    }


    restrictDF = DataFrame()
    restrictDF.addColumns(["Network", "Cells", "Disease", "Other"], "")

    for x in networkRestrictions:

        nrestricts = defaultdict(list)

        for rt in networkRestrictions[x]:
            nrestricts[rt] = networkRestrictions[x][rt]

        nrestricts['disease'] += [{'group':'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}]
        restricts = nrestricts


        networkDRdict = defaultdict(str)
        networkDRdict["Network"] =  networkToTitle[x]

        diseaseElems = []
        cellElems = []
        otherElems = []

        for restrictType in restricts:

            if restrictType == "sentences":
                continue

            if restrictType in ["disease"]:
                for elem in restricts[restrictType]:
                    diseaseElems.append( elem['name'] + " ("+elem['termid']+")")

            elif restrictType in ["cells"]:
                for elem in restricts[restrictType]:
                    cellElems.append( elem['name'] + " ("+elem['termid']+")")

            else:
                for elem in restricts[restrictType]:
                    otherElems.append( elem['name'] + " ("+elem['termid']+")")


        networkDRdict['Cells'] =  "\makecell[l]{" +"\\\\".join(sorted(cellElems)) + "}"
        networkDRdict['Disease'] = "\makecell[l]{" + "\\\\".join(sorted(diseaseElems)) + "}"
        networkDRdict['Other'] = "\makecell[l]{" + "\\\\".join(sorted(otherElems)) + "}"

        dr = DataRow.fromDict(networkDRdict)
        restrictDF.addRow(dr)

    print(restrictDF._makeLatex())

    #exit()



    allMissing = {}
    figidx = 0

    mirna2cellOut = open("/mnt/d/yanc_network/important_process.txt", 'w')

    for network in networks:
        figidx+= 1

        networkGraph = nx.Graph()

        if network in ignoreNetworks:
            continue


        interactions = defaultdict(set)

        acceptedInteractions = defaultdict(set)
        typeByGene = defaultdict(lambda: Counter())
        elemsByGene = defaultdict(lambda: defaultdict(set))

        allMirna = set(networks[network])

        miStr2mirna = {}
        allTargetMirna = []
        mirnaObj2str = {}

        mirna2evs = defaultdict(set)

        newAllMirna = set()

        for x in allMirna:
            try:

                oMirna = miRNA(x)

                allTargetMirna.append( oMirna )
                miStr = oMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

                miStr2mirna[miStr] = oMirna
                mirnaObj2str[oMirna] = miStr

                newAllMirna.add( miStr )

            except:
                pass

        allMirna = newAllMirna
        #allMirna = set([str(x) for x in allTargetMirna])
        requestData = None


        if network in networkRestrictions:
            requestData = networkRestrictions[network]
        else:
            requestData = {
                            'sentences': "false",
                            }

        requestData['sentences'] = "false"


        requestData["mirna"]= list(allMirna)
        print(allMirna)

        if not 'disease' in requestData and not network in ['targetMirsVasRemod']:
            requestData['disease'] = [{'group': 'disease', 'termid': 'DOID:1936', 'name': 'atherosclerosis'}]#[{'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},{'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}]

        #requestData['disease'] += [
        #            {'group': 'disease', 'termid': 'DOID:1287', 'name': 'cardiovascular system disease'},
        #            {'group': 'disease', 'termid': 'DOID:2349', 'name': 'arteriosclerosis'}
        #        ]


        print(requestData)

        graph, nodeCounter, edge2datasourceCount, jsonRes = DataBasePlotter.fetchGenes(requestData)

        print(len(jsonRes['rels']))

        htmlDF = DataFrame()
        htmlDF.addColumns(['gene rel', 'gene', 'miRNA Group', 'miRNA', 'Original Network', 'PubMed', 'MIRECORD', 'MIRTARBASE', 'DIANA', 'Disease', 'Cells', 'GO'])

        for rel in jsonRes['rels']:

            orderedEdge = [None, None]

            if rel['ltype'] == "gene":
                orderedEdge[0] = rel['lid']
            elif rel['ltype'] == "mirna":
                orderedEdge[1] = rel['lid']

            if rel['rtype'] == "gene":
                orderedEdge[0] = rel['rid']
            elif rel['rtype'] == "mirna":
                orderedEdge[1] = rel['rid']

            orderedEdges = set()

            if orderedEdge[1].startswith("microRNAS"):
                continue

            wasAccepted = False
            """
            for tMirna in allTargetMirna:

                if tMirna.accept(orderedEdge[1]):

                    wasAccepted = True

                    orderedEdges.add(
                        (orderedEdge[0], str(tMirna))
                    )
            """

            if not wasAccepted:
                orderedEdges.add(tuple(orderedEdge))


            for oEdge in orderedEdges:

                origEdge = tuple(oEdge)
                edgeStatus = None
                oEdge = list(oEdge)

                wasFound = False
                for miObj in mirnaObj2str:

                    if miObj.accept(oEdge[1], compLevel=miRNACOMPARISONLEVEL.PRECURSOR):
                        oEdge[1] = mirnaObj2str[miObj]
                        wasFound = True
                        break

                if not wasFound:

                    try:
                        miObj = miRNA(oEdge[1])
                        miStr = miObj.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR])

                        miRNA(miStr)

                        oEdge[1] = miStr

                    except:
                        print("Could not read/load", oEdge, miStr)
                        continue

                allGeneMirna = interactions[oEdge[0]]
                miAccepted = False

                allAcceptedStr = set()


                for strMirna in allGeneMirna:
                    miObj = miStr2mirna.get(strMirna, None)

                    if miObj == None:
                        continue

                    miObjAccepts = miObj.accept(oEdge[1], compLevel=miRNACOMPARISONLEVEL.PRECURSOR)
                    miAccepted = miAccepted or miObjAccepts


                    if miObjAccepts:
                        acceptedInteractions[oEdge[0]].add(strMirna)
                        edgeStatus = "accepted"

                        allAcceptedStr.add(strMirna)

                        #print(oEdge[0], oEdge[1], strMirna)

                if not miAccepted:
                    edgeStatus = "additional"

                networkGraph.add_edge(oEdge[0], oEdge[1], color= 'g' if edgeStatus == "accepted" else "b")

                typeByGene[oEdge[0]][edgeStatus] += 1
                elemsByGene[oEdge[0]][edgeStatus].add(oEdge[1])

                objMirna = miRNA(oEdge[1])

                pmidEvs = set()
                mirtarbaseEvs = set()
                mirecordsEvs = set()
                dianaEvs = set()

                docDiseases = set()
                docCells = set()
                docGOs = set()

                for ev in rel['evidences']:

                    docid = ev['docid']

                    mirna2evs[oEdge[1]].add(docid)

                    disEvs = jsonRes['pmidinfo'].get('disease', {}).get(docid, {})

                    for disEv in disEvs:
                        did = disEv['termid']
                        dname = disEv['termname']

                        docDiseases.add((dname, did, docid))

                    cellEvs = jsonRes['pmidinfo'].get('cells', {}).get(docid, {})

                    for cellEv in cellEvs:
                        did = cellEv['termid']
                        dname = cellEv['termname']

                        for ct in cellType2AccTerms:

                            ctTerms = cellType2AccTerms[ct]

                            if did in ctTerms:
                                print(network, oEdge[0], oEdge[1], ct, docid, sep="\t", file=mirna2cellOut)

                        docCells.add((dname, did, docid))

                    goEvs = jsonRes['pmidinfo'].get('go', {}).get(docid, {})

                    for goEv in goEvs:
                        did = goEv['termid']
                        dname = goEv['termname']

                        docGOs.add((dname, did, docid))

                    if ev['data_source'] == "DIANA":
                        dianaEvs.add( (ev['method'], ev['direction']) )

                    elif ev['data_source'] == "miRTarBase":
                        mirtarbaseEvs.add(
                            (ev['data_id'], ",".join(ev['exp_support']), ev['functional_type'], ev['docid'])
                        )
                    elif ev['data_source'] == "pmid":
                        pmidEvs.add((ev['docid'],))
                    elif ev['data_source'] == "mirecords":
                        mirecordsEvs.add((ev['docid']))
                    else:
                        print("Unhandled data source", ev['data_source'])


                dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
                    geneCap=oEdge[0].upper(), geneLow=oEdge[1].upper())


                pmidStr = "<br/>".join(
                    [
                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}\" target=\"_blank\">{pmid}</a>".format(pmid=elem[0]) for elem in pmidEvs
                    ]
                )

                mirtarbaseStr = "<br/>".join(
                    [
                        "<a href=\"http://mirtarbase.mbc.nctu.edu.tw/php/detail.php?mirtid={mtbid}\">{mtbid}</a>".format(mtbid=elem[0]) for elem in mirtarbaseEvs
                    ]
                )

                mirecordStr = "<br/>".join(
                    [
                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/{pmid}\" target=\"_blank\">{pmid}</a>".format(pmid=elem[0]) for elem in mirecordsEvs
                    ]
                )

                dianaStr = "<br/>".join(
                    [
                        "{method} {direction}".format(method=elem[0], direction=elem[1]) for elem in dianaEvs
                    ]
                )

                goStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for
                        elem
                        in docGOs
                    ]
                )

                cellStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for
                        elem
                        in docCells
                    ]
                )

                diseaseStr = "<br/>".join(
                    [
                        "{method} ({direction}, {docid})".format(method=elem[0], direction=elem[1], docid=elem[2]) for
                        elem
                        in docDiseases
                    ]
                )

                addRow = {
                    'gene rel': oEdge[0] + "<br/>" + oEdge[1],
                    'gene': oEdge[0],
                    'miRNA Group': objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                    'miRNA': "<br/>".join(allAcceptedStr),
                    'Original Network': "{edgestate}</br>".format(edgestate=edgeStatus) +
                                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term={miRes}+{miShort}\">Search PUBMED</a>".format(miRes=oEdge[1], miShort=objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID]))+
                                        "</br><a href=\"{dianaLink}\">Search DIANA</a>".format(dianaLink=dianaLink)
                    ,
                    'PubMed': pmidStr,
                    'MIRECORD': mirecordStr,
                    'MIRTARBASE': mirtarbaseStr,
                    'DIANA': dianaStr,
                    'Disease': diseaseStr,
                    'Cells': cellStr,
                    'GO': goStr
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)

        for gene in interactions:

            for mirna in interactions[gene]:

                edgeWasFound = mirna in acceptedInteractions[gene]

                if edgeWasFound:
                    continue

                edgeStatus = "missing"

                networkGraph.add_edge(gene, mirna, color='r')


                typeByGene[gene][edgeStatus] += 1
                elemsByGene[gene][edgeStatus].add(mirna)

                objMirna = miRNA(mirna)
                dianaLink = "http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=&genes%5B%5D={geneCap}&genes%5B%5D={geneLow}&sources%5B%5D=1&sources%5B%5D=7&sources%5B%5D=9&publication_year=&prediction_score=&sort_field=&sort_type=&query=1".format(
                    geneCap=gene.upper(), geneLow=gene.upper())

                addRow = {
                    'gene rel': gene,
                    'gene': gene,
                    'miRNA Group': objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID, miRNAPART.PRECURSOR]),
                    'miRNA': mirna,
                    'Original Network': "{edgestate}</br>".format(edgestate=edgeStatus) +
                                        "<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/?term={gene} {miRes}+{miShort}\">Search PUBMED</a>".format(gene=gene,miRes=mirna, miShort=objMirna.getStringFromParts([miRNAPART.MATURE, miRNAPART.ID]))+
                                        "</br><a href=\"{dianaLink}\">Search DIANA</a>".format(dianaLink=dianaLink)
                    ,
                    'PubMed': "",
                    'MIRECORD': "",
                    'MIRTARBASE': "",
                    'DIANA': ""
                }


                row = DataRow.fromDict(addRow)
                htmlDF.addRow(row)



        elemsByMirna = defaultdict(set)

        for gene in elemsByGene:

            for expMirna in allMirna:

                for category in elemsByGene[gene]:

                    if category == 'missing':
                        continue

                    for foundMirna in elemsByGene[gene][category]:

                        elemsByMirna[foundMirna].add(gene)


        foundMirnas = set([x for x in elemsByMirna])

        minEvs = 1

        while True:

            addMirnas = [x for x in foundMirnas.difference(allMirna) if len(mirna2evs[x]) >= minEvs]

            if len(addMirnas) > 50:
                minEvs += 1

            else:
                break

        print(network)
        print("Found Mirnas", len(foundMirnas), list(foundMirnas))
        print("Expected Mirnas", len(allMirna), list(allMirna))
        print("Intersected Mirnas", len(foundMirnas.intersection(allMirna)), list(foundMirnas.intersection(allMirna)))
        print("Missing Mirnas", len(allMirna.difference(foundMirnas)), allMirna.difference(foundMirnas))
        print("Additional Mirnas", len(foundMirnas.difference(allMirna)), foundMirnas.difference(allMirna))
        print("Additional Mirnas filtered", len(addMirnas))
        print("Filter level", minEvs)

        allMissing[network] = allMirna.difference(foundMirnas)

        rowDict = {}
        rowDict['Network'] = "\makecell[l]{"+networkToTitle[network] + "\\\\(min evidences: "+str(minEvs) + ", additionals: "+str(len(foundMirnas.difference(allMirna)))+")" + "}"
        rowDict['Accepted miRNAs'] = "\makecell[l]{" +"\\\\".join(natsorted(foundMirnas.intersection(allMirna), key=lambda x: x.split("-")[1])) + "}"
        rowDict['Additional miRNAs'] = "\makecell[l]{" + "\\\\".join(natsorted(addMirnas, key=lambda x: x.split("-")[1])) + "}"
        rowDict['Missing miRNAs'] = "\makecell[l]{"  + "\\\\".join(natsorted(allMirna.difference(foundMirnas), key=lambda x: x.split("-")[1])) + "}"

        newRow = DataRow.fromDict(rowDict)

        #["Network", "Accepted miRNAs", "Missing miRNAs"]
        summaryDF.addRow( newRow )

        if False:

            print(network)
            for gene in sorted([x for x in typeByGene]):
                print(gene, typeByGene[gene], elemsByGene[gene]['missing'])

            print()
            print()

            print(network)
            for gene in sorted([x for x in typeByGene]):
                print("Gene:", gene, "Status: ", ", ".join([": ".join([x, str(typeByGene[gene][x])]) for x in typeByGene[gene]]), "Missing miRNAs: "+",".join(elemsByGene[gene]['missing']))

            print()
            print()
            print()
            print()


        networkGraphs[network] = networkGraph

        htmlDF.export("/mnt/d/yanc_network/" + network.replace(" ", "_") + ".html", ExportTYPE.HTML)
        htmlDF.export("/mnt/d/yanc_network/" + network.replace(" ", "_") + ".tsv", ExportTYPE.TSV)



    figidx = 0
    for stages in makeStory:

        mergedGraph = networkGraphs[stages[0]]

        for i in range(1, len(stages)):
            mergedGraph = nx.compose(mergedGraph, networkGraphs[stages[i]])


        hasLargeStage = any(['large' in stage for stage in stages])

        pos = nx.spring_layout(mergedGraph)

        for stage in stages:


            networkGraph = networkGraphs[stage]

            edges = networkGraph.edges()
            colors = [networkGraph[u][v]['color'] for u, v in edges]

            d = nx.degree(networkGraph)
            nodes = networkGraph.nodes()
            nodeColors = []
            nodeSizes = []

            allSizes = [x[1] for x in d]
            minSize = min(allSizes)
            maxSize = max(allSizes)
            diffSize = maxSize-minSize

            fontSize = 16
            minNodeSize = 1200
            figSize = (20, 14)
            edgeWidth = 3

            if hasLargeStage:
                fontSize = 8
                minNodeSize = 100
                figSize = (20,30)
                edgeWidth = 0.75


            plt.figure(figidx, figsize=figSize)
            figidx += 1

            maxNodeSize = 3000
            diffNodeSize = maxNodeSize-minNodeSize
            nodeList = []






            for x in nodes:
                if any([x.lower().startswith(y) for y in ['mir', 'let']]):
                    nodeColors.append('blue')
                else:
                    nodeColors.append('green')

                nodeDegree = d[x]

                nodeDegree -= minSize
                nodeDegree = nodeDegree / diffSize

                nodeSize = minNodeSize + diffNodeSize * nodeDegree

                nodeSizes.append( nodeSize )
                nodeList.append(x)

            nx.draw(networkGraph, pos, font_size=fontSize, with_labels=False, node_color=nodeColors, edges=edges, edge_color=colors, nodelist=nodeList, node_size=nodeSizes, width=edgeWidth, font_weight='bold', dpi=1000)


            for p in pos:  # raise text positions

                clist = list(pos[p])

                if p in nodeList:
                    if nodeSizes[nodeList.index(p)] < 1000:
                        clist[1] = clist[1] + 0.005
                    else:
                        clist[1] = clist[1] + 0.02
                pos[p] = tuple(clist)

            nx.draw_networkx_labels(networkGraph, pos, font_weight='bold', font_size=fontSize)

            plt.suptitle(stage)

            plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".png")
            plt.savefig("/mnt/d/yanc_network/" + stage.replace(" ", "_") + ".pdf")

    #plt.show()

    print(summaryDF._makeLatex())


    for x in allMissing:
        for mirna in allMissing[x]:
            print(x, mirna)

        print()
        print()


