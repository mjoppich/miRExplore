import argparse
import datetime
import pickle
import re
import shlex
import regex
import sys
import os

from synonymes.SynonymFile import Synfile
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


    graphObj = make_plot_info(interactReq.get('elements', None), interactReq.get('categories', None),
                              interactReq.get('messengers', None), organisms=interactReq.get('organisms', None),
                              majorParents=True, onlyCells=False, obolevel=obolevel )

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


    graphObj = make_plot_info(interactReq.get('elements', None), interactReq.get('categories', None),
                              interactReq.get('messengers', None), organisms=interactReq.get('organisms', None),
                              majorParents=True, onlyCells=True, obolevel=obolevel)

    return app.make_response((jsonify(graphObj), 200, None))



def make_plot_info(cells, categories, messengers, organisms, majorParents=True, onlyCells=False, obolevel=2):

    if obolevel == None:
        print("Error: obolevel was none", obolevel)
        obolevel = 4

    elemJSON = makeInteractionObject(cells, categories, messengers, organisms, fetchSentences=False)

    rels = elemJSON['rels']
    addinfo = elemJSON['pmidinfo']

    docInfos = defaultdict(lambda: defaultdict(set))

    for x in addinfo:

        alladdinfo = addinfo[x]
        for pmid in alladdinfo:
            for ev in alladdinfo[pmid]:
                termid = ev['termid']

                docInfos[x][pmid].add(termid)


    """
    
    In order to avoid too many single nodes, assign each child to a major parent ....
    
    """

    allRoots = cellsObo.getRoots()
    cellRoot = cellsObo.dTerms.get('CL:0000000', None)
    if cellRoot != None:
        #allRoots = []
        allRoots.append( cellRoot )

    oboid2majorid = {}
    oboid2rootdist = {}

    for root in allRoots:

        oboid2majorid[root.id] = root.id
        oboid2rootdist[root.id] = -1000

        allchildren = root.getChildrenAtLevel(obolevel, withLevel=True)

        print("Root node", root.id, root.name, [(x.id, x.name, l) for (x,l) in allchildren])

        for child,l in allchildren:

            newc = child.getAllChildren(withLevel=True)
            oboid2majorid[child.id] = child.id
            oboid2rootdist[child.id] = -1000

            for nc,cl in newc:
                if 'CL:0000893' in nc.term.id or 'CL:0000893' in child.id:
                    print( nc.term.id, nc.term.name, cl, "=>", child.name, child.id, l)


                if not nc.term.id in oboid2rootdist or oboid2rootdist[nc.term.id] > cl:
                    oboid2majorid[nc.term.id] = child.id
                    oboid2rootdist[nc.term.id] = cl

                    if 'CL:0000893' in nc.term.id or 'CL:0000893' in child.id:
                        print( nc.term.id, nc.term.name, cl, "==>", child.name, child.id, l)



    oldname2newname = {}

    allNodes = set()

    for x in rels:

        if majorParents:

            for evidence in x['evidences']:

                src = evidence.get('lontid', evidence['lid'])
                tgt = evidence.get('rontid', evidence['rid'])

                newsrc = oboid2majorid.get(src, src)
                newtgt = oboid2majorid.get(tgt, tgt)

                if newsrc != None:
                    srcTerm = cellsObo.dTerms[newsrc]
                    srcName = srcTerm.name

                    allNodes.add(srcName)
                else:
                    print("No oboid2majorid for", src)
                    srcName = x['lid']

                if newtgt != None:
                    tgtTerm = cellsObo.dTerms[newtgt]
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

    print("allNodes")
    for x in allNodes:
        print("obo node", x)

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


    for x in rels:

        for evidence in x['evidences']:

            lid = evidence.get('lontid', evidence['lid'])
            rid = evidence.get('rontid', evidence['rid'])

            if majorParents:

                lid = oldname2newname.get(lid, evidence['lid'])
                rid = oldname2newname.get(rid, evidence['rid'])

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

                mNodeIdx = allNodes.index(messengerID)

                if not onlyCells:
                    allEdges[(edgeCellCell[1], mNodeIdx)] += len(categoryEdges)
                    allEdges[edgeCellCell] += len(categoryEdges)

                for categoryID in categoryEdges:

                    cNodeIdx = allNodes.index(categoryID)

                    if not onlyCells:
                        allEdges[(mNodeIdx, cNodeIdx)] += 1


    """
    
    Fetch actually used edges
    
    """

    usedNodeIndices = {}
    usedNodes = []
    for src, tgt in allEdges:

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
        'nodes': [{'name': x} for x in usedNodes],
        'links': [{'source': x[0],
                   'target': x[1],
                   'value': usedEdges[x],
                   'messengers': usedMessengerEdgesCounter[x],
                   'categories': usedCategoryEdgesCounter[x]
                   }
                  for x in usedEdges]
    }

    print("nodes")

    for x in graphObj['nodes']:
        print(x)


    print("links")
    for x in graphObj['links']:
        print(x)

    return graphObj

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

    return returnInteractions(interactReq.get('elements', None), interactReq.get('categories', None), interactReq.get('messengers', None), organisms=interactReq.get('organisms', None), fetchSentences=fetchSentences)

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def makeInteractionObject(cells=None, categories=None, messengers=None, organisms=None, fetchSentences=True):

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




    print("Loading sentences")
    for rel in allRels:

        if organisms != None:

            if rel.orgs == None:
                continue

            if not any([x in rel.orgs for x in organisms]):
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
                continue
            else:

                acceptEv = any([x['termid'] in allowedIDs['categories'] for x in docIDDiseaseInfo])

                if not acceptEv:
                    continue

        if messengers != None and len(allowedIDs['messengers']) > 0:

            docID = evJSON['docid']

            docIDGOInfo = addInfo['messengers'].get(docID, [])

            if len(docIDGOInfo) == 0:
                continue
            else:

                acceptEv = any([x['termid'] in allowedIDs['messengers'] for x in docIDGOInfo])

                if not acceptEv:
                    continue


        """
        
        THIS SECTION IS ONLY PASSED BY EVIDENCES WHICH ARE SELECTED
        
        """



        if fetchSentences:
            evSent = evJSON.get('rel_sentence', None)

            if evSent != None:
                evSentTxt = sentDB.get_sentence(evSent)

                if evSentTxt != None:

                    loadedSents += 1
                    evJSON['sentence'] = evSentTxt[1]

        relID = (rel.l_id_type, rel.r_id_type)

        foundRels[relID].append(evJSON)

    print("Loaded Sentences", loadedSents)


    allrels = []
    for lent, rent in foundRels:

        lrEvs = foundRels[(lent, rent)]

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

    return returnObj


def returnInteractions(cells=None, categories=None, messengers=None, organisms=None, fetchSentences=True):


    returnObj = makeInteractionObject(cells, categories, messengers, organisms, fetchSentences=fetchSentences)

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

    reMatch = regex.compile(searchWord+'{e<=3}')

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

            ront = relDB.lontology

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


    pmidBase = args.textmine + '/aggregated/'

    normGeneSymbols = normalize_gene_names(path=args.obodir + "/hgnc_no_withdrawn.syn")

    print("Loading Interactions")

    # allInteractions = defaultdict(list)
    allDBS = None

    print(datetime.datetime.now(), "Loading PMID2PMC")
    pmid2pmcDB = PMID2PMCDB.loadFromFile(pmidBase + '/pmid2pmc')
    print(datetime.datetime.now(), "Loading mirel")

    testRels = None#TestRelLoader.loadFromFile(pmidBase + "/test_rels_4")

    cellsObo = GeneOntology(args.obodir + "/cl.obo")
    ccPMID = MiGenRelDB.loadFromFile(pmidBase + "/cell_cell.pmid", ltype="CELLS", rtype="CELLS", normGeneSymbols=normGeneSymbols, lontology = cellsObo, rontology = cellsObo)


    relDBs = [ccPMID]
    relDBs = [x for x in relDBs if x != None]

    print(datetime.datetime.now(), "Finished mirel")

    mirFeedback = feedbackDB(args.feedback)

    print(datetime.datetime.now(), "Loading sents")
    sentDB = SentenceDB.loadFromFile(args.sentdir, pmidBase+"/pmid2sent")
    print(datetime.datetime.now(), "Finished sents")

    requiredPMIDs = set()
    for rdb in relDBs:

        assert (isinstance(rdb, DataBaseDescriptor))

        for rpmid in rdb.get_evidence_docids():
            requiredPMIDs.add(rpmid)

    if os.path.isfile(pmidBase + "/dbs.pickle"):
        print(datetime.datetime.now(), "Loading pickle")
        with open(pmidBase + "/dbs.pickle", 'rb') as fin:
            allDBS = pickle.load(fin)

        pmid2Categories = allDBS[0]
        pmid2Messenger = allDBS[1]

        print(datetime.datetime.now(), "Loading pickle ended")

    print(datetime.datetime.now(), "Loading ontologies")

    categoriesObo = GeneOntology(args.obodir + "/categorization.obo")
    messengersObo = GeneOntology(args.obodir + "/messenger.obo")

    print(datetime.datetime.now(), "Loading ontologies finished")


    if allDBS == None:

        print(datetime.datetime.now(), "Loading GO")
        pmid2Categories = PMID2XDB.loadFromFile(pmidBase + "/categories.pmid", categoriesObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading Disease")
        pmid2Messenger = PMID2XDB.loadFromFile(pmidBase + "/messenger.pmid", messengersObo, requiredPMIDs)


        allDBS = (pmid2Categories, pmid2Messenger)

        print(datetime.datetime.now(), "Writing Pickle")

        with open(pmidBase + "/dbs.pickle", 'wb') as fout:
            pickle.dump(allDBS, fout)

    print(datetime.datetime.now(), "Loading finished")


def getCLParser():
    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str,
                        help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=True)
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


def gunicorn_start():

    parser = getCLParser()

    argstr = "--textmine /home/j/joppich/tm_soehnlein/ --obodir /home/j/joppich/tm_soehnlein/obos --sentdir /home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/pmid_sent --feedback /home/j/joppich/tm_soehnlein/feedback --port 65522"

    args = parser.parse_args(shlex.split(argstr))

    start_app_from_args(args)

    return app
