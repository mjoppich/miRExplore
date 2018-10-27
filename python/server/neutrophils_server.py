import argparse
import datetime
import pickle
import re
import shlex
import regex
import sys
import os

from neutrophils.allowedSentences import AllowedSentences
from synonymes.SynonymFile import Synfile
from textmining.SentenceID import SentenceID
from utils.HashDict import HashDict
from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")

sys.stdout = sys.stderr

from natsort import natsorted


import time

from synonymes.GeneOntology import GeneOntology
from textdb.AbstractDBClasses import DataBaseDescriptor
from textdb.MiGenRelDB import MiGenRelDB
from textdb.MirTarBaseDB import MirTarBaseDB
from textdb.MirandaRelDB import MirandaRelDB
from textdb.PMID2PMCDB import PMID2PMCDB
from textdb.PMID2XDB import PMID2XDB
from textdb.SentenceDB import SentenceDB
from textdb.feedback_db import feedbackDB

from textdb.TestRelLoader import TestRelLoader
from analysis.miRecordDB import miRecordDB


from io import StringIO

from flask import Flask, jsonify, request, redirect, url_for, send_from_directory
import json
import pprint
from collections import defaultdict, Counter

from flask_cors import CORS

dataurl = str(os.path.dirname(os.path.realpath(__file__))) + "/../../" + 'frontend_neutrophils/src/static/'

app = Flask(__name__, static_folder=dataurl, static_url_path='/static')

app.debug = True
CORS(app)

app.config['DEBUG'] = False
app.config['UPLOAD_FOLDER'] = ""

allOrgInfos = [{
        "name": "Homo sapiens",
        "group": "organism",
        "termid": "Homo sapiens",
        "syns": ['human', 'hsa']
    },
    {
        "name": "Mus musculus",
        "group": "organism",
        "termid": "Mus musculus",
        "syns": ['mouse', 'mmu']
    }

]


# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']



@app.route('/')
def root():

    retFile = 'index.html'

    return app.send_static_file(retFile)

@app.route('/test', methods=['GET', 'POST'])
def test():

    return "<html><body>neutrophils Server v0.01</body></html>", 200, None


@app.route('/help', methods=['GET', 'POST'])
def help():
    res = "<html><body><ul>"

    for x in [rule.rule for rule in app.url_map.iter_rules() if rule.endpoint !='static']:
        res += "<li>"+str(x)+"</li>"

    res +="</body></html>"

    return res, 200, None



@app.route('/oboinfo', methods=['POST'])
def oboinfo():
    interactReq = request.get_json(force=True, silent=True)

    global categoriesObo
    global messengersObo
    global cellsObo

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'oboname' in interactReq:
        return app.make_response((jsonify( {'error': 'must include oboname'} ), 400, None))

    print(interactReq)

    oboname = interactReq['oboname']


    selObo = {
        'cells': cellsObo,
        'messengers': messengersObo,
        'categories': categoriesObo
    }

    if not oboname in selObo:
        return app.make_response((jsonify( {'error': 'oboname must be in ' + ", ".join([x for x in selObo])} ), 400, None))

    obo = selObo[oboname]
    retObj = []


    for termid in sorted([x for x in obo.dTerms]):

        node = obo.dTerms[termid]

        nodeEntry = {}

        nodeEntry['termid'] = node.id
        nodeEntry['termname'] = node.name

        allChildren = []


        if node.children != None:
            for elem in node.children:
                allChildren.append({
                    'id': elem.term.id,
                    'name': elem.term.name
                })

        nodeEntry['children'] = sorted(allChildren, key=lambda x: x['name'])

        allSyns = []


        if node.synonym != None:
            for syn in node.synonym:

                if syn == None:
                    continue
                allSyns.append(syn.syn)

        nodeEntry['synonymes'] = sorted(allSyns)

        retObj.append(nodeEntry)


    return app.make_response((jsonify(retObj), 200, None))



@app.route('/sankey', methods=['POST'])
def sankey_network():
    interactReq = request.get_json(force=True, silent=True)

    global categoriesObo
    global messengersObo

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'elements' in interactReq:
        return app.make_response((jsonify( {'error': 'must include CELLS'} ), 400, None))

    print(interactReq)

    obolevel = interactReq.get('obolevel', 3)
    messenger_obolevel = interactReq.get('messenger_obolevel', 3)

    selected = interactReq.get('elements')
    selectedCells = [x['termid'] for x in selected if x['group'] == 'CELLS']


    graphObj = make_plot_info(interactReq.get('elements', None), interactReq.get('categories', None),
                              interactReq.get('messengers', None), organisms=interactReq.get('organisms', None),
                              majorParents=True, onlyCells=False, obolevel=obolevel, messenger_obolevel=messenger_obolevel, selectedElems=selectedCells )

    return app.make_response((jsonify(graphObj), 200, None))


@app.route('/interaction_network', methods=['POST'])
def interaction_network():
    interactReq = request.get_json(force=True, silent=True)

    global categoriesObo
    global messengersObo

    if interactReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    if not 'elements' in interactReq:
        return app.make_response((jsonify({'error': 'must include CELLS'}), 400, None))

    print(interactReq)
    obolevel = interactReq.get('obolevel', 3)
    messenger_obolevel = interactReq.get('messenger_obolevel', 3)


    graphObj = make_plot_info(interactReq.get('elements', None), interactReq.get('categories', None),
                              interactReq.get('messengers', None), organisms=interactReq.get('organisms', None),
                              majorParents=True, onlyCells=True, obolevel=obolevel, messenger_obolevel=messenger_obolevel)


    for edge in graphObj['links']:
        edge['group1'] = edge['messengers']
        edge['group2'] = edge['categories']


    return app.make_response((jsonify(graphObj), 200, None))



def getOboID2Major(selectedTermIDs, hitsPerOboID, obo, primaryRoots, remainingRoots, oboLevel = 3, minOccurrences=10):


    selectedIDs = {}
    for termID in selectedTermIDs:

        if termID in obo.dTerms:

            selectedIDs[termID] = termID

            ecTerm = obo.dTerms[termID]
            ac = ecTerm.getAllChildren()

            for child in ac:
                selectedIDs[child.term.id] = termID


    def prepareObo2Major( base, overwrite=True ):

        if overwrite or not base.id in oboid2majorid:
            oboid2majorid[base.id] = base.id
            oboid2rootdist[base.id] = -1000

        baseChildren = base.getAllChildren(maxLevel=oboLevel)
        baseChildrenIDs = [x.term.id for x in baseChildren]


        for child in baseChildren:
            if not child.term.id in oboid2rootdist or oboid2rootdist[child.term.id] > 0:
                if overwrite or not child.term.id in oboid2majorid:

                    if child.term.id in selectedIDs:

                        oboid2majorid[child.term.id] = selectedIDs[child.term.id]
                        oboid2rootdist[child.term.id] = -1000

                    else:


                        if hitsPerOboID[child.term.id] > minOccurrences:
                            oboid2majorid[child.term.id] = child.term.id
                        else:

                            if overwrite and child.term.id.startswith("CL"):
                                pass#print(child.term.name, child.term.id, hitsPerOboID[child.term.id])

                            #find best parent to assign this child to

                            level2parent = child.term.getAllParentyByLevel()

                            parentSet = False

                            for i in range(min(level2parent), max(level2parent)):
                                if i in level2parent:

                                    for parent in level2parent[i]:
                                        if hitsPerOboID[parent.id] > minOccurrences:
                                            oboid2majorid[child.term.id] = parent.id
                                            parentSet = True
                                            break

                                if parentSet:
                                    break

                            if parentSet == False:
                                oboid2majorid[child.term.id] = base.id
                                print("Base Ontology Term for", child.term.id, child.term.name)

                        if child.term.id == 'CL:1000497':
                            print("Too few hits for", child.term.name, oboid2majorid[child.term.id], hitsPerOboID[child.term.id], base.id)


                        oboid2rootdist[child.term.id] = -1000

        allchildren = base.getChildrenAtLevel(oboLevel, withLevel=True)


        for child, l in allchildren:

            newc = child.getAllChildren(withLevel=True)

            if not child.id in oboid2majorid:
                oboid2majorid[child.id] = child.id
                oboid2rootdist[child.id] = -1000

            for nc, cl in newc:

                if not nc.term.id in oboid2rootdist or oboid2rootdist[nc.term.id] > cl:

                    if overwrite or not nc.term.id in oboid2majorid:

                        if nc.term.id in selectedIDs:

                            oboid2majorid[nc.term.id] = selectedIDs[nc.term.id]
                            oboid2rootdist[nc.term.id] = -1000

                        else:

                            oboid2majorid[nc.term.id] = oboid2majorid[child.id]

                            if nc.term.id == 'CL:1000497':
                                print("Too few hits for", child.name, nc.term.name, oboid2majorid[child.id], oboid2majorid[nc.term.id], base.id)

                            oboid2rootdist[nc.term.id] = cl



    oboid2majorid = {}
    oboid2rootdist = {}

    print([x for x in primaryRoots])

    for root in primaryRoots:
        prepareObo2Major(root)


    if remainingRoots != None:
        for root in remainingRoots:
            prepareObo2Major(root, overwrite=False)

    return oboid2majorid




def make_plot_info(cells, categories, messengers, organisms, majorParents=True, onlyCells=False, obolevel=2, messenger_obolevel=1, selectedElems=None):

    global messengerID2Name
    global ccPMID

    elemJSON = makeInteractionObject(cells, categories, messengers, organisms, fetchSentences=False)

    rels = elemJSON['rels']
    addinfo = elemJSON['pmidinfo']

    docInfos = defaultdict(lambda: defaultdict(list))

    for x in addinfo:

        alladdinfo = addinfo[x]
        for pmid in alladdinfo:
            for ev in alladdinfo[pmid]:
                termid = (ev['termid'], ev['evidences'])

                docInfos[x][pmid].append(termid)




    """
    
    count evidence-hits per obo id
    
    """

    hitsPerOboID = defaultdict(lambda: Counter())
    messengerHitsPerOboID = Counter()



    lOntology = ccPMID.lontology
    rOntology = ccPMID.rontology

    evCount = 10

    if majorParents:
        evCount = 0
        for x in rels:

            for evidence in x['evidences']:

                evCount += 1

                src = evidence.get('lontid', evidence['lid'])
                tgt = evidence.get('rontid', evidence['rid'])

                hitsPerOboID[ccPMID.ltype][src] += 1

                if src in lOntology.dTerms:
                    oboNode = lOntology.dTerms[src]

                    for child in oboNode.getAllParents():
                        hitsPerOboID[ccPMID.ltype][child.id] += 1

                hitsPerOboID[ccPMID.rtype][tgt] += 1

                if tgt in rOntology.dTerms:
                    oboNode = rOntology.dTerms[tgt]

                    for child in oboNode.getAllParents():
                        hitsPerOboID[ccPMID.rtype][child.id] += 1

                evDocId = evidence['docid']
                messengerIDs = docInfos['messengers'].get(evDocId, [])

                for messengerID,nev in messengerIDs:

                    if messengerID in messengersObo.dTerms:
                        oboNode = messengersObo.dTerms[messengerID]

                        messengerHitsPerOboID[oboNode.id] += 1

                        for child in oboNode.getAllParents():

                            if child == None:
                                print("null child", oboNode.id)
                                continue

                            messengerHitsPerOboID[child.id] += 1


    """
    
    In order to avoid too many single nodes, assign each child to a major parent ....
    
    """
    ontType2OboMajor = {}

    print(ccPMID.types())

    for ontType in set(ccPMID.types()):

        ont = ccPMID.getOntologyForType(ontType)

        allRoots = ont.getRoots()

        rootNodes = None
        if ontType == "CELLS":
            rootNodes = [ont.dTerms.get('CL:0000255', None)]


        if rootNodes == None or rootNodes == [None]:
            rootNodes = ont.getRoots()

        selectedCellElems = [x['termid'] for x in cells if x['group'] == ontType]
        oboid2majorid = getOboID2Major(selectedCellElems, hitsPerOboID[ontType], ont, rootNodes, allRoots, minOccurrences=-1)

        ontType2OboMajor[ontType] = oboid2majorid



    if messengers != None:
        selectedMessengerElems = [x['termid'] for x in messengers if x['group'] == 'messengers']
    else:
        selectedMessengerElems = []

    messengerMinOccs = int(evCount * 0.1)
    print("Messenger Minimum Occurrences", messengerMinOccs)
    messengerOboid2Major = getOboID2Major(selectedMessengerElems, messengerHitsPerOboID, messengersObo, messengersObo.getRoots(), None, messenger_obolevel, minOccurrences=-1)


    oldname2newname = {}

    allNodes = set()


    for x in rels:

        if majorParents:

            for evidence in x['evidences']:

                src = evidence.get('lontid', evidence['lid'])
                tgt = evidence.get('rontid', evidence['rid'])

                srctype = evidence.get('ltype', "CELLS")
                tgtype = evidence.get('rtype', "CELLS")

                newsrc = ontType2OboMajor[srctype].get(src, src)
                newtgt = ontType2OboMajor[tgtype].get(tgt, tgt)


                if newsrc != None:
                    srcTerm = ccPMID.lontology.dTerms[newsrc]
                    srcName = srcTerm.name

                    allNodes.add(srcName)
                else:
                    print("No oboid2majorid for", src)
                    srcName = x['lid']

                if newtgt != None:
                    tgtTerm = ccPMID.rontology.dTerms[newtgt]
                    tgtName = tgtTerm.name

                    allNodes.add(tgtName)
                else:
                    print("No oboid2majorid for", tgt)
                    tgtName = x['rid']


                oldname2newname[src] = srcName
                oldname2newname[tgt] = tgtName

        else:
            allNodes.add(x['lid'])
            allNodes.add(x['rid'])


    for termid in messengersObo.dTerms:
        obonode = messengersObo.dTerms[termid]

        if obonode.name in allNodes:
            print("Duplicate messenger term", obonode.name)

        allNodes.add(obonode.name)

    for termid in categoriesObo.dTerms:
        obonode = categoriesObo.dTerms[termid]

        if obonode.name in allNodes:
            print("Duplicate category term", obonode.name)

        allNodes.add(obonode.name)

    allNodes.add('mUnknown')
    allNodes.add('cUnknown')


    allNodes = list(sorted(allNodes))

    """
    now we need to add edges :)
    """

    allEdges = Counter()
    messengerEdgesCounter = Counter()
    categoryEdgesCounter = Counter()

    cellSrc = set()
    cellTgt = set()


    def check_circ(lid, rid):

        edgeCellCell = (allNodes.index(lid), allNodes.index(rid))

        ccSrc = edgeCellCell[1] in cellSrc
        ccTgt = edgeCellCell[0] in cellTgt

        if ccSrc or ccTgt:
            edgeCellCell = tuple(reversed(edgeCellCell))
        else:
            return edgeCellCell

        ccSrc = edgeCellCell[1] in cellSrc
        ccTgt = edgeCellCell[0] in cellTgt

        if ccSrc or ccTgt:
            return None

        return edgeCellCell

    name2obo = defaultdict(set)

    """
    
    BACKUP of original code: cell cell messenger effect
    
    for x in rels:

        for evidence in x['evidences']:

            lid = evidence.get('lontid', evidence['lid'])
            rid = evidence.get('rontid', evidence['rid'])

            if majorParents:

                lid = oldname2newname.get(lid, evidence['lid'])
                rid = oldname2newname.get(rid, evidence['rid'])

            name2obo[lid].add(evidence.get('lontid', evidence['lid']))
            name2obo[rid].add(evidence.get('rontid', evidence['rid']))


            if (lid == rid):
                rid = rid + "_c"

                if not rid in allNodes:
                    allNodes.append(rid)

            edgeCellCell = check_circ(lid, rid)

            if edgeCellCell == None:


                addL = lid + "_c" in allNodes
                addR = rid + "_c" in allNodes

                if addL and not addR:
                    lid = lid + "_c"
                elif not addL and addR:
                    rid = rid + "_c"
                elif addL and addR:
                    print("lid and rid modified in allnodes?", lid, rid, allNodes)
                elif not addL and not addR:
                    lid = lid + "_c"
                    allNodes.append(lid)

                edgeCellCell = check_circ(lid, rid)

                if edgeCellCell == None:

                    print("Circular", lid, rid, allNodes)
                    continue

            cellSrc.add(edgeCellCell[0])
            cellTgt.add(edgeCellCell[1])


            evDocId = evidence['docid']

            messengerEdges = docInfos['messengers'].get(evDocId, ['mUnknown'])
            categoryEdges = docInfos['categories'].get(evDocId, ['cUnknown'])

            assert( len(messengerEdges) > 0 )

            messengerEdgesCounter[edgeCellCell] += len(messengerEdges)
            categoryEdgesCounter[edgeCellCell] += len(categoryEdges)

            if onlyCells:
                allEdges[edgeCellCell] += 1

            for messengerID in messengerEdges:


                if messengerID == 'mUnknown':
                    messengerName = messengerID
                else:

                    messengerID = messengerOboid2Major[messengerID]
                    messengerName = messengerID2Name[messengerID]

                mNodeIdx = allNodes.index(messengerName)

                if not onlyCells:
                    allEdges[(edgeCellCell[1], mNodeIdx)] += len(categoryEdges)
                    allEdges[edgeCellCell] += len(categoryEdges)

                for categoryID in categoryEdges:

                    cNodeIdx = allNodes.index(categoryID)

                    if not onlyCells:
                        allEdges[(mNodeIdx, cNodeIdx)] += 1
    """

    """
    
    Generating Sankey Chart: CELL MESSAGE CELL EFFECT
    
    """
    allInteracts = 0

    for x in rels:

        for evidence in x['evidences']:

            lid = evidence.get('lontid', evidence['lid'])
            rid = evidence.get('rontid', evidence['rid'])

            sentid = evidence['rel_sentence']

            if majorParents:

                lid = oldname2newname.get(lid, evidence['lid'])
                rid = oldname2newname.get(rid, evidence['rid'])

            if rid == 'root':
                print("ROOT", evidence)

            name2obo[lid].add(evidence.get('lontid', evidence['lid']))
            name2obo[rid].add(evidence.get('rontid', evidence['rid']))


            if (lid == rid):
                rid = rid + "_c"

                if not rid in allNodes:
                    allNodes.append(rid)

            edgeCellCell = check_circ(lid, rid)

            if edgeCellCell == None:


                addL = lid + "_c" in allNodes
                addR = rid + "_c" in allNodes

                if addL and not addR:
                    lid = lid + "_c"
                elif not addL and addR:
                    rid = rid + "_c"
                elif addL and addR:
                    print("lid and rid modified in allnodes?", lid, rid, allNodes)
                elif not addL and not addR:
                    lid = lid + "_c"
                    allNodes.append(lid)

                edgeCellCell = check_circ(lid, rid)

                if edgeCellCell == None:

                    print("Circular", lid, rid, allNodes)
                    continue

            cellSrc.add(edgeCellCell[0])
            cellTgt.add(edgeCellCell[1])


            evDocId = evidence['docid']

            sentence = evidence.get("sentence", "")

            aSent = sentid.split(".")
            aSentNum = int(aSent[-1])
            aSent = aSent[0:2]
            allowedSentIDs = set()
            allowedSentIDs.add(str(SentenceID.fromArray(aSent + [aSentNum - 1])))
            allowedSentIDs.add(str(SentenceID.fromArray(aSent + [aSentNum + 1])))

            effect = docInfos['categories'].get(evDocId, [None])
            message = docInfos['messengers'].get(evDocId, [None])

            foundEffects = set()
            if effect != None:

                for termID, evidences in effect:
                    for x in evidences:
                        if x[0] == sentid:
                            foundEffects.add(termID)

                if len(foundEffects) == 0:
                    for termID, evidences in effect:
                        for x in evidences:
                            if x[0] in allowedSentIDs:
                                foundEffects.add(termID)

            foundMessages = set()
            if message != None:

                for termID, evidences in message:
                    for x in evidences:
                        if x[0] == sentid:
                            foundMessages.add(termID)

                if len(foundMessages) == 0:
                    for termID, evidences in message:
                        for x in evidences:
                            if x[0] in allowedSentIDs:
                                foundMessages.add(termID)



            messengerEdgesCounter[edgeCellCell] += len(foundMessages)
            categoryEdgesCounter[edgeCellCell] += len(foundEffects)

            for messengerID in foundMessages:


                if messengerID == 'mUnknown':
                    messengerName = messengerID
                else:

                    origMessID = messengerID
                    messengerID = messengerOboid2Major[messengerID]
                    messengerName = messengerID2Name[messengerID]

                    if "root" in [origMessID, messengerID, messengerName]:
                        print("root sighted in", [origMessID, messengerID, messengerName])

                mNodeIdx = allNodes.index(messengerName)


                for categoryID in foundEffects:

                    cNodeIdx = allNodes.index(categoryID)

                    if not onlyCells:


                        allEdges[(edgeCellCell[0], mNodeIdx)] += 1
                        allEdges[(mNodeIdx, edgeCellCell[1])] += 1
                        allEdges[(edgeCellCell[1], cNodeIdx)] += 1

                        allInteracts +=1

    print("Max interacts", allInteracts)

    """
    
    Fetch actually used edges
    
    """

    usedNodeIndices = {}
    usedNodes = []
    for src, tgt in allEdges:

        print(src, tgt, allNodes[src], allNodes[tgt], allEdges[(src, tgt)])

        nodename = allNodes[src]
        if not nodename in usedNodes:

            usedNodes.append(nodename)
            usedNodeIndices[src] = usedNodes.index(nodename)

        nodename = allNodes[tgt]
        if not nodename in usedNodes:
            usedNodes.append(nodename)
            usedNodeIndices[tgt] = usedNodes.index(nodename)

    """

    reindex edges

    """

    usedEdges = {}
    usedMessengerEdgesCounter = Counter()
    usedCategoryEdgesCounter = Counter()

    for src, tgt in allEdges:

        newsrc = usedNodeIndices[src]
        newtgt = usedNodeIndices[tgt]

        assert(newsrc != newtgt)

        usedMessengerEdgesCounter[(newsrc, newtgt)] = messengerEdgesCounter[(src, tgt)]
        usedCategoryEdgesCounter[(newsrc, newtgt)] = categoryEdgesCounter[(src, tgt)]

        usedEdges[(newsrc, newtgt)] = allEdges[(src, tgt)]



    graphObj = {
        'nodes': [{'name': x, 'obo': list(name2obo[x])} for x in usedNodes],
        'links': [{'source': x[0],
                   'target': x[1],
                   'srcname': usedNodes[x[0]],
                   'tgtname': usedNodes[x[1]],
                   'value': usedEdges[x],
                   'messengers': usedMessengerEdgesCounter[x],
                   'categories': usedCategoryEdgesCounter[x]
                   }
                  for x in usedEdges]
    }


    debug=False

    if debug:
        print("nodes")

        for x in graphObj['nodes']:
            print(x)


        print("links")
        for x in graphObj['links']:
            print(x)

    return graphObj


@app.route("/query", methods=['POST'])
def queryDB():
    interactReq = request.get_json(force=True, silent=True)

    global categoriesObo
    global messengersObo

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'elements' in interactReq:
        return app.make_response((jsonify( {'error': 'must include CELLS'} ), 400, None))

    print(interactReq)

    fetchSentences=True
    if 'sentences' in interactReq:
        if interactReq['sentences'].upper() == 'FALSE':
            fetchSentences = False

    maxSentDist=1
    if 'sentence_distance' in interactReq:
        maxSentDist = int(interactReq['sentence_distance'])

    retInfo = makeInteractionObject(interactReq.get('elements', None), interactReq.get('categories', None),
                              interactReq.get('messengers', None), None, fetchSentences=fetchSentences, maxSentDist=maxSentDist)

    return app.make_response((jsonify(retInfo), 200, None))



@app.route('/find_interactions', methods=['GET', 'POST'])
def findInteractions():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'elements' in interactReq:
        return app.make_response((jsonify( {'error': 'must include CELLS'} ), 400, None))

    print(interactReq)


    fetchSentences=True
    if 'sentences' in interactReq:
        if interactReq['sentences'].upper() == 'FALSE':
            fetchSentences = False

    maxSentDist=1
    if 'sentence_distance' in interactReq:
        maxSentDist = int(interactReq['sentence_distance'])


    return returnInteractions(interactReq.get('elements', None),
                              interactReq.get('categories', None),
                              interactReq.get('messengers', None),
                              organisms=interactReq.get('organisms', None),
                              fetchSentences=fetchSentences,
                              maxSentDist=maxSentDist
                              )

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def makeInteractionObject(cells=None, categories=None, messengers=None, organisms=None, fetchSentences=True, maxSentDist=1):

    global relDBs
    global sentDB
    global pmid2Categories
    global pmid2Messenger
    global categoriesObo
    global messengersObo

    cells = cells if cells!=None else []
    categories = categories if categories!=None else []
    messengers = messengers if messengers!=None else []

    foundRels = defaultdict(list)

    allDocIDs = set()
    allRels = []

    if True:
        #proceed as usual

        allRelsByType = defaultdict(list)

        for relDB in relDBs:

            for celltype in cells:

                dbrels = relDB.get_rels(celltype['group'], celltype['termid'])

                if dbrels == None:
                    continue

                adbrels = [x for x in dbrels if x.assocSent != None]

                allRelsByType[celltype['group']] += adbrels


        for etype in allRelsByType:

            allRels += allRelsByType[etype]


    allRels = natsorted(allRels, key=lambda x: x.docid)

    allRelsEx = None
    if len(allRels) > 0:
        allRelsEx=allRels[0]

    print("Loading", len(allRels), "relations like", allRelsEx)

    loadedSents = 0

    if organisms != None:

        organisms = [x['termid'] for x in organisms]

        neworgs = set()

        if 'Homo sapiens' in organisms:
            neworgs.add("hsa")

        if 'Mus musculus' in organisms:
            neworgs.add("mmu")

        organisms = neworgs

        print("Must include any organism", organisms)

    allowedIDs = {}
    if categories != None:

        allowedTermIDs = []
        for delem in categories:
            elemTerm = categoriesObo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        allowedIDs['categories'] = set(allowedTermIDs)

    if messengers != None:

        allowedTermIDs = []
        for delem in messengers:
            elemTerm = messengersObo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        allowedIDs['messengers'] = set(allowedTermIDs)


    """
    
    FETCH ALL POSSIBLY RELEVANT DOCIDs
    
    """

    allDocIDs = set()
    for rel in allRels:

        if organisms != None:

            if rel.orgs == None:
                continue

            if not any([x in rel.orgs for x in organisms]):
                continue

        if rel.docid != -1:
            allDocIDs.add(rel.docid)

    """

    CHECK ALL POSSIBLY RELEVANT DOCIDs for RELEVANCE

    """



    addInfo = {}
    if pmid2Categories:

        addInfo['categories'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2Categories.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['categories'] = allDocInfos

    if pmid2Messenger:

        addInfo['messengers'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2Messenger.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['messengers'] = allDocInfos


    """
    
    FOUND CATEGORIES FOR ALL DOCIDs
    
    
    GOING TO LOAD SENTENCES
        
    """

    seenTuples = set()

    oneDebug = 0
    takenPMIDs = set()

    missingGroupsCounter = Counter()
    missingSentenceCounter = Counter()
    missingTrust = 0
    incorrectVerbDir = 0
    seenEvidences = 0
    duplicateDescriptor=0
    reviewPMIDCount = 0
    noOrganismCount = 0
    incorrectOrganismCount = 0

    allowedSentsGen = AllowedSentences(maxSentDist)

    print("Loading sentences")
    for rel in allRels:

        seenEvidences += 1

        if rel.docid in reviewPMIDs:
            reviewPMIDCount += 1
            continue


        if organisms != None:

            if rel.orgs == None:
                noOrganismCount += 1
                continue

            if not any([x in rel.orgs for x in organisms]):
                incorrectOrganismCount += 1
                continue

        if rel.docid != -1:
            allDocIDs.add(rel.docid)

        evJSON = rel.toJSON()

        """
        
        THIS PART WILL REJECT EVIDENCES IF INCORRECT SELECTOR
        
        """

        if categories != None and len(allowedIDs['categories']) > 0:

            docID = evJSON['docid']

            docIDDiseaseInfo = addInfo['categories'].get(docID, [])

            if len(docIDDiseaseInfo) == 0:
                #print("docIDDisease")
                continue
            else:

                acceptEv = any([x['termid'] in allowedIDs['categories'] for x in docIDDiseaseInfo])

                if not acceptEv:
                    #print("docIDDisease ev")
                    continue

        if messengers != None and len(allowedIDs['messengers']) > 0:

            docID = evJSON['docid']

            docIDGOInfo = addInfo['messengers'].get(docID, [])

            if len(docIDGOInfo) == 0:
                #print("docIDmessenger")
                continue
            else:

                acceptEv = any([x['termid'] in allowedIDs['messengers'] for x in docIDGOInfo])
                #print("docIDmessenger ev")

                if not acceptEv:
                    continue


        """
        
        THIS SECTION IS ONLY PASSED BY EVIDENCES WHICH ARE SELECTED
        
        """



        directionalAction = True
        if directionalAction == True:

            if not evJSON['rel_direction'] == "12":
                incorrectVerbDir += 1
                continue


        if True:

            docID = evJSON['docid']
            rejectDoc = False


            missingGroups = []

            for infoGroup in ['messengers', 'categories']:

                docIDInfo = addInfo[infoGroup].get(docID, [])

                if len(docIDInfo) == 0:
                    rejectDoc = True
                    missingGroups.append(infoGroup)

            missingGroupsCounter[tuple(sorted(missingGroups))] += 1

            if rejectDoc:
                if evJSON['docid'] == "22391089":
                    print("Rejected infoGroup", docID, infoGroup)
                continue

            sameSentence = True


            if sameSentence:

                evSentID = evJSON['rel_sentence']
                rejectDoc = False
                missingSentenceGroup = set()

                for infoGroup in ['messengers', 'categories']:
                    docIDInfo = addInfo[infoGroup].get(docID, [])

                    #docIDInfo != 0
                    seenSentences = set()

                    for docEvs in docIDInfo:

                        newSents = allowedSentsGen.getAllowedSentences(docEvs['evidences'])

                        if newSents != None:
                            seenSentences = seenSentences.union(newSents)


                    if evJSON['docid'] == "313430":
                        print("DocID info check", docID, docIDInfo, infoGroup, seenSentences)


                    if not evSentID in seenSentences:
                        rejectDoc = True

                        missingSentenceGroup.add(infoGroup)

                        if oneDebug == 0:

                            print(evSentID, seenSentences, docIDInfo)
                            oneDebug += 1

                missingSentenceCounter[tuple(sorted(missingSentenceGroup))] += 1

                if len(missingSentenceGroup) > 0:
                    if 'messengers' in missingSentenceGroup or 'categories' in missingSentenceGroup:
                        continue


                evTrust = evJSON['trust']

                trustDoc = evTrust['conj'] == 0

                if evTrust['stack']:
                    trustDoc = True

                if evTrust['verb']:
                    trustDoc = True

                if evTrust['relex']:
                    trustDoc = True

                if not trustDoc:
                    rejectDoc = True

                if rejectDoc:

                    missingTrust += 1

                    if evJSON['docid'] == "22391089":
                        print("Rejected trust")
                    continue


        tupleDescriptor = (evJSON['docid'], evJSON['rel_sentence'], evJSON['lid'], evJSON['rid'])

        if tupleDescriptor in seenTuples:
            duplicateDescriptor += 1
            continue


        if fetchSentences:
            evSent = evJSON.get('rel_sentence', None)

            if evSent != None:
                evSentTxt = sentDB.get_sentence(evSent)

                if evSentTxt != None:

                    loadedSents += 1
                    evJSON['sentence'] = evSentTxt[1]

        relID = (rel.l_id_type, rel.r_id_type)

        docID = evJSON['docid']
        takenPMIDs.add(docID)

        foundRels[relID].append(evJSON)

        seenTuples.add(tupleDescriptor)

    print("Loaded Sentences", loadedSents)
    print("Presenting results from", len(takenPMIDs), "PMIDs")


    takenEvs = 0
    allrels = []
    for lent, rent in foundRels:

        lrEvs = foundRels[(lent, rent)]

        takenEvs += len(lrEvs)

        allrels.append({'lid':lent[0],
                        'rid': rent[0],
                        'ltype': lent[1],
                        'rtype': rent[1],
                        'evidences': lrEvs
                        })

    returnObj = {
        'rels': allrels,
        'pmidinfo': addInfo
    }

    print("")
    print("")
    print("")
    print("")
    print("")
    print("")
    print("")
    print("Max Sentence Distance", maxSentDist)
    print("Seen Evidences", seenEvidences)
    print("incorrect Verb Direction", incorrectVerbDir)
    print("missing Trust", missingTrust)
    print("missing groups", missingGroupsCounter)
    print("duplicate descriptor", duplicateDescriptor)
    print("review Doc", reviewPMIDCount)
    print("no org", noOrganismCount)
    print("wrong org", incorrectOrganismCount)
    print("Missing group for sentence", missingSentenceCounter)
    print("taken rels", len(allrels))
    print("taken evidences", takenEvs)
    print("")
    print("")
    print("")
    print("")
    print("")
    print("")
    print("")
    print("")

    reviewPMIDCount = 0
    noOrganismCount = 0
    incorrectOrganismCount = 0

    return returnObj


def returnInteractions(cells=None, categories=None, messengers=None, organisms=None, fetchSentences=True, maxSentDist=1):

    returnObj = makeInteractionObject(cells, categories, messengers, organisms, fetchSentences=fetchSentences, maxSentDist=maxSentDist)

    return app.make_response((jsonify(returnObj), 200, None))

@app.route("/feedback_info", methods=['GET', 'POST'])
def getFeedbackInfo():

    sendJSON = request.get_json(force=True, silent=True)
    dataID = sendJSON['data_id']

    print(sendJSON)

    dbRes = mirFeedback.get_feedback(dataID)

    print(dbRes)

    if dbRes == None or len(dbRes) == 0:
        dbRes = []


    feedbackPos = 0
    feedbackNeg = 0

    for fback in dbRes:
        if fback[4] == True:
            feedbackPos += 1
        else:
            feedbackNeg -= 1

    feedbackCum = feedbackPos+feedbackNeg

    return app.make_response((jsonify({
        'data_id': dataID,
        'feedbacks': len(dbRes),
        'feedback_pos': feedbackPos,
        'feedback_neg': feedbackNeg,
        'feedback_cum': feedbackCum
    }), 200, None))


@app.route('/organisms', methods=['GET', 'POST'])
def getOrganisms():
    #return app.make_response((jsonify(['Homo sapiens', 'Mus musculus']), 200, None))

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) <= 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for orgInfo in allOrgInfos:

        termName = orgInfo['name']
        termGroup = orgInfo['group']
        termID = orgInfo['termid']
        syns = orgInfo['syns']

        for word in [termName]+syns:
            if reMatch.match(word):
                jsonResult.append(
                    {
                        'name': word,
                        'termid': termID,
                        'group': 'organism'
                    }
                )

        if len(jsonResult) > 100:
            break

    return app.make_response((jsonify(jsonResult), 200, None))



@app.route('/autocomplete', methods=['GET', 'POST'])
def findID():

    jsonResult = set()

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    print("Autocomplete")
    print(searchWord)

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify( [] ), 200, None))

    reMatch = regex.compile(searchWord+'{e<=3}', regex.IGNORECASE)

    jsonResultByType = defaultdict(set)

    for relDB in relDBs:

        ltype = relDB.ltype
        rtype = relDB.rtype

        if relDB.l_ont_based:

            lont = relDB.lontology

            for termid in lont.dTerms:

                goterm = lont.dTerms[termid]

                if len(jsonResultByType[ltype]) <= 100 and reMatch.match(termid):
                    jsonResultByType[ltype].add(HashDict({ 'name': goterm.name, 'termid': goterm.id, 'group': ltype }))

                if len(jsonResultByType[ltype]) <= 100 and reMatch.match(goterm.name):
                    jsonResultByType[ltype].add(HashDict({ 'name': goterm.name, 'termid': goterm.id, 'group': ltype }))

        else:
            for entName in relDB.all_ltypes:
                if len(jsonResultByType[ltype]) <= 100 and reMatch.match(entName):
                    jsonResultByType[rtype].add(HashDict({'name': entName, 'termid': entName, 'group': ltype}))

        if relDB.r_ont_based:

            ront = relDB.rontology

            for termid in ront.dTerms:

                goterm = ront.dTerms[termid]

                if len(jsonResultByType[ltype]) <= 100 and reMatch.match(termid):
                    jsonResultByType[ltype].add(HashDict({'name': goterm.name, 'termid': goterm.id, 'group': rtype}))

                if len(jsonResultByType[ltype]) <= 100 and reMatch.match(goterm.name):
                    jsonResultByType[ltype].add(HashDict({'name': goterm.name, 'termid': goterm.id, 'group': rtype}))


        else:
            for entName in relDB.all_rtypes:
                if len(jsonResultByType[rtype]) <= 100 and reMatch.match(entName):
                    jsonResultByType[rtype].add(HashDict({'name': entName, 'termid': entName, 'group': rtype}))

    jsonResult = []
    for type in jsonResultByType:

        jsonResult += list(
            jsonResultByType[type]
        )


    return app.make_response((jsonify( jsonResult ), 200, None))


@app.route('/category_ac', methods=['GET', 'POST'])
def disease_autocomplete():

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in pmid2Categories.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'category'
                }
            )

        if len(jsonResult) > 100:
            break

    return app.make_response((jsonify(jsonResult), 200, None))

@app.route('/messengers_ac', methods=['GET', 'POST'])
def go_autocomplete():

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in pmid2Messenger.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'messengers'
                }
            )

        if len(jsonResult) > 100:
            break

    return app.make_response((jsonify(jsonResult), 200, None))

@app.route('/cells_ac', methods=['GET', 'POST'])
def cells_autocomplete():
    global ccPMID

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in ccPMID.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'CELLS'
                }
            )

        if len(jsonResult) > 100:
            break

    return app.make_response((jsonify(jsonResult), 200, None))

@app.route("/relation_feedback", methods=["GET", "POST"])
def rel_feedback():

    global mirFeedback

    feedbackInfo = request.get_json(force=True, silent=True)

    print(feedbackInfo)

    docID = feedbackInfo['docid']
    relSentID = feedbackInfo.get('rel_sentence', None)

    dataSource = feedbackInfo['data_source']
    dataID = feedbackInfo['data_id']

    lid = feedbackInfo['lid']
    rid = feedbackInfo['rid']

    lontid = feedbackInfo.get('lontid', None)
    rontid = feedbackInfo.get('rontid', None)

    ltypePOS = feedbackInfo.get('lpos', None)
    rtypePOS = feedbackInfo.get('rpos', None)

    ltype = feedbackInfo['ltype']
    rtype = feedbackInfo['rtype']

    approve = feedbackInfo['approve']

    searchTermGenes = tuple(feedbackInfo['search_terms'].get('gene', []))

    print(dataID, approve)

    mirFeedback.add_feedback( (dataSource, dataID, docID, relSentID, approve, lid, rid, ltype, rtype, ltypePOS, rtypePOS, searchTermGenes, lontid, rontid) )

    return app.make_response((jsonify({}), 200, None))



relDBs = None
sentDB = None
pmid2Categories = None
pmid2Messenger = None
categoriesObo = None
messengersObo = None
pmid2pmcDB = None
mirFeedback = None
ccPMID = None
cellsObo = None
neutrophilsObo = None
messengerID2Name = None
reviewPMIDs = None


def start_app_from_args(args):

    global relDBs
    global sentDB
    global pmid2Categories
    global pmid2Messenger
    global categoriesObo
    global messengersObo
    global ccPMID
    global cellsObo
    global mirFeedback
    global pmid2pmcDB
    global messengerID2Name
    global reviewPMIDs


    with open(args.obodir + "/reviewPMID") as fin:
        reviewPMIDs = set()

        for line in fin:
            line = line.strip()

            aline = line.split("\t")

            if aline[1] == "REVIEW":
                reviewPMIDs.add(aline[0])



    pmidBase = args.textmine + '/aggregated/'
    pmcBase = args.textmine + '/pmc_aggregation/'

    normGeneSymbols = normalize_gene_names(path=args.obodir + "/hgnc_no_withdrawn.syn")

    print("Loading Interactions")

    # allInteractions = defaultdict(list)
    allDBS = None

    print(datetime.datetime.now(), "Loading PMID2PMC")
    pmid2pmcDB = PMID2PMCDB.loadFromFile(pmidBase + '/pmid2pmc')
    print(datetime.datetime.now(), "Loading mirel")

    testRels = None#TestRelLoader.loadFromFile(pmidBase + "/test_rels_4")
    neutrophilsObo = GeneOntology(args.obodir + "/neutrophils.obo")
    cellsObo = GeneOntology(args.obodir + "/target_cells.obo")

    ccPMC = MiGenRelDB.loadFromFile(pmcBase + "/neutrophils_cell.pmid", ltype="NEUTROPHIL", rtype="CELLS",
                                     normGeneSymbols=normGeneSymbols, lontology=neutrophilsObo, rontology=cellsObo,
                                     lReplaceSc=False, rReplaceSc=False, switchLR=False)

    allPMCDocIDs = ccPMC.get_evidence_docids()

    rejectPMIDs = set()
    for pmc in allPMCDocIDs:

        if pmc in pmid2pmcDB.pmc2pmid:

            for pmid in pmid2pmcDB.pmc2pmid[pmc]:
                rejectPMIDs.add(pmid)


    print("Rejectful PMIDs", len(rejectPMIDs))

    """
    Here we should ensure to ignore any PMID which has already been listed in PMC ...
    
    """

    print("REMOVE PMCs from PMID !!!!", file=sys.stderr)

    ccPMID = MiGenRelDB.loadFromFile(pmidBase + "/neutrophils_cell.pmid", ltype="NEUTROPHIL", rtype="CELLS",
                                     normGeneSymbols=normGeneSymbols, lontology = neutrophilsObo, rontology = cellsObo,
                                     lReplaceSc=False, rReplaceSc=False, switchLR=False, ignoreDocIDs=rejectPMIDs)


    relDBs = [ccPMID, ccPMC]
    relDBs = [x for x in relDBs if x != None]

    print(datetime.datetime.now(), "Finished mirel")

    mirFeedback = feedbackDB(args.feedback)

    print(datetime.datetime.now(), "Loading sents")
    sentDB = SentenceDB.loadFromFile(args.sentdir, "/mnt/c/dev/data/pmid2sent")
    sentDBPMC = SentenceDB.loadFromFile(args.pmcsentdir, "/mnt/c/dev/data/pmc2sent")

    sentDB.add_database(sentDBPMC)

    print(datetime.datetime.now(), "Finished sents")

    requiredDOCIDs = set()
    for rdb in relDBs:

        assert (isinstance(rdb, DataBaseDescriptor))

        for rpmid in rdb.get_evidence_docids():
            requiredDOCIDs.add(rpmid)

    pickleFile = args.textmine + "/dbs.pickle"

    if os.path.isfile( pickleFile ):
        print(datetime.datetime.now(), "Loading pickle", pickleFile)
        with open(pickleFile, 'rb') as fin:
            allDBS = pickle.load(fin)

        pmid2Categories = allDBS[0]
        pmid2Messenger = allDBS[1]

        print(datetime.datetime.now(), "Loading pickle ended", pickleFile)

    print(datetime.datetime.now(), "Loading ontologies")

    categoriesObo = GeneOntology(args.obodir + "/effects.obo")
    messengersObo = GeneOntology(args.obodir + "/messages.obo")

    messengerID2Name = {}

    for termID in messengersObo.dTerms:
        term= messengersObo.dTerms[termID]
        messengerID2Name[term.id] = term.name


    print(datetime.datetime.now(), "Loading ontologies finished")


    if allDBS == None:

        print(datetime.datetime.now(), "Loading GO")
        pmid2Categories = PMID2XDB.loadFromFile(pmidBase + "/effects.pmid", categoriesObo, requiredDOCIDs)
        pmc2Categories = PMID2XDB.loadFromFile(pmcBase + "/effects.pmid", categoriesObo, requiredDOCIDs)


        print(datetime.datetime.now(), "Loading Disease")
        pmid2Messenger = PMID2XDB.loadFromFile(pmidBase + "/messages.pmid", messengersObo, requiredDOCIDs)
        pmc2Messenger = PMID2XDB.loadFromFile(pmcBase + "/messages.pmid", messengersObo, requiredDOCIDs)

        pmid2Categories.add_database(pmc2Categories)
        pmid2Messenger.add_database(pmc2Messenger)


        allDBS = (pmid2Categories, pmid2Messenger)

        print(datetime.datetime.now(), "Writing Pickle", pickleFile)

        with open(pickleFile, 'wb') as fout:
            pickle.dump(allDBS, fout)

    print(datetime.datetime.now(), "Loading finished", pickleFile)


def getCLParser():
    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str,
                        help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-sp', '--pmcsentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=False)
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    return parser

if __name__ == '__main__':

    parser = getCLParser()

    args = parser.parse_args()

    for x in args.__dict__:
        print(x, args.__dict__[x])

    start_app_from_args(args)

    print("Starting Flask on port", args.port)
    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)


def gunicorn_start(datadir="/home/j/joppich/tm_soehnlein/", sentdir="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/pmid_sent", feedbackFile=None):

    parser = getCLParser()

    argstr = "--textmine {datadir} --obodir {datadir}/obodir --sentdir {sentdir} --feedback {feedbackFile}".format(datadir=datadir, sentdir=sentdir, feedbackFile=feedbackFile)

    print("Starting app with")
    print(argstr)

    args = parser.parse_args(shlex.split(argstr))

    start_app_from_args(args)

    return app
