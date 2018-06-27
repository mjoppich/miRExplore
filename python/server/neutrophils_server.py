import argparse
import datetime
import pickle
import re
import shlex
import regex
import sys
import os


sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from synonymes.SynonymFile import Synfile
from utils.HashDict import HashDict
from utils.tmutils import normalize_gene_names

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
from collections import defaultdict

from flask_cors import CORS


app = Flask(__name__)
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


@app.route('/find_interactions', methods=['GET', 'POST'])
def findInteractions():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'elements' in interactReq:
        return app.make_response((jsonify( {'error': 'must include CELLS'} ), 400, None))

    print(interactReq)

    return returnInteractions(interactReq.get('elements', None), interactReq.get('categories', None), interactReq.get('messengers', None), organisms=interactReq.get('organisms', None))

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def returnInteractions(cells=None, categories=None, messengers=None, organisms=None):


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

    for rel in allRels:

        if organisms != None:

            if rel.orgs == None:
                continue

            if not any([x in rel.orgs for x in organisms]):
                continue

        if rel.docid != -1:
            allDocIDs.add(rel.docid)

        evJSON = rel.toJSON()

        evSent = evJSON.get('rel_sentence', None)

        if evSent != None:
            evSentTxt = sentDB.get_sentence(evSent)

            if evSentTxt != None:

                loadedSents += 1
                evJSON['sentence'] = evSentTxt[1]

        relID = (rel.l_id_type, rel.r_id_type)

        foundRels[relID].append(evJSON)

    print("Loaded Sentences", loadedSents)

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


    allrels = []
    for lent, rent in foundRels:

        lrEvs = foundRels[(lent, rent)]

        okEvs = []

        for jsonEV in lrEvs:

            if categories != None and len(allowedIDs['categories']) > 0:

                docID = jsonEV['docid']

                docIDDiseaseInfo = addInfo['categories'].get(docID, [])

                if len(docIDDiseaseInfo) == 0:
                    continue
                else:

                    acceptEv = any([x['termid'] in allowedIDs['categories'] for x in docIDDiseaseInfo])

                    if not acceptEv:
                        continue

            if messengers != None and len(allowedIDs['messengers']) > 0:

                docID = jsonEV['docid']

                docIDGOInfo = addInfo['messengers'].get(docID, [])

                if len(docIDGOInfo) == 0:
                    continue
                else:

                    acceptEv = any([x['termid'] in allowedIDs['messengers'] for x in docIDGOInfo])

                    if not acceptEv:
                        continue

            okEvs.append(jsonEV)


        if len(okEvs) > 0:



            allrels.append({'lid':lent[0],
                            'rid': rent[0],
                            'ltype': lent[1],
                            'rtype': rent[1],
                            'evidences': okEvs
                            })

    returnObj = {
        'rels': allrels,
        'pmidinfo': addInfo
    }

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


    if reMatch.match("Tester"):

        for testerID in testRels.tester2rels:
            if not testerID.startswith("Tester"):
                continue
            jsonResult.append({'name': testerID, 'group': 'gene'})

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

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c.docid) for c in re.split('(\d+)', text) ]


def start_app_from_args(args):

    global relDBs
    global categoriesObo
    global messengersObo
    global pmid2Categories
    global pmid2Messenger
    global mirFeedback
    global sentDB


    pmidBase = args.textmine + '/aggregated/'

    normGeneSymbols = normalize_gene_names(args.obodir + "/hgnc_no_withdrawn.syn")

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
        pmid2Categories = None
        pmid2Messenger = None

        print(datetime.datetime.now(), "Loading GO")
        pmid2Categories = PMID2XDB.loadFromFile(pmidBase + "/categories.pmid", categoriesObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading Disease")
        pmid2Messenger = PMID2XDB.loadFromFile(pmidBase + "/messenger.pmid", messengersObo, requiredPMIDs)


        allDBS = (pmid2Categories, pmid2Messenger)

        print(datetime.datetime.now(), "Writing Pickle")

        with open(pmidBase + "/dbs.pickle", 'wb') as fout:
            pickle.dump(allDBS, fout)

    print(datetime.datetime.now(), "Loading finished")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str,
                        help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=True)
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    args = parser.parse_args()

    for x in args.__dict__:
        print(x, args.__dict__[x])

    start_app_from_args(args)

    print("Starting Flask on port", args.port)
    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)


def gunicorn_start():
    argstr = "--textmine /home/j/joppich/tm_soehnlein/ --obodir /home/j/joppich/tm_soehnlein/obos --sentdir /home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/pmid_sent/ --feedback /home/j/joppich/tm_soehnlein/feedback --port 65522"

    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str,
                         help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=True)
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    args = parser.parse_args(shlex.split(argstr))

    start_app_from_args(args)

    return app
