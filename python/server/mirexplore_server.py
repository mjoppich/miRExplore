import argparse
import datetime
import pickle
import regex
import sys
import os

from synonymes.SynonymFile import Synfile
from textdb.SymbolEnsemblDB import SymbolEnsemblDB
from utils.tmutils import normalize_gene_names

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


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


dataurl = str(os.path.dirname(os.path.realpath(__file__))) + "/../../" + 'frontend/src/static/'
app = Flask(__name__, static_folder=dataurl, static_url_path='/static')
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

    return "<html><body>miRExplore Server v0.01</body></html>", 200, None


@app.route('/help', methods=['GET', 'POST'])
def help():
    res = "<html><body><ul>"

    for x in [rule.rule for rule in app.url_map.iter_rules() if rule.endpoint !='static']:
        res += "<li>"+str(x)+"</li>"

    res +="</body></html>"

    return res, 200, None



@app.route('/get_features/<geneID>')
def getGeneFeatures(geneID):

    if geneID == None:
        return app.make_response((jsonify( {'error': 'invalid geneid'} ), 400, None))

    return getGeneMirnaFeatures(geneID, None)

@app.route('/get_features/<geneID>/<mirnaID>')
def getGeneMirnaFeatures(geneID, mirnaID):

    if geneID == None:
        return app.make_response((jsonify( {'error': 'invalid geneid'} ), 400, None))

    return app.make_response((jsonify({'gene_id': geneID, 'mirna_id': mirnaID }), 200, None))




@app.route('/find_interactions', methods=['GET', 'POST'])
def findInteractions():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'gene' in interactReq and not 'mirna' in interactReq:
        return app.make_response((jsonify( {'error': 'must include gene or mirna'} ), 400, None))

    print(interactReq)

    return returnInteractions(interactReq.get('gene', None), interactReq.get('mirna', None), interactReq.get('lncrna', None), organisms=interactReq.get('organisms', None), diseaseRestrictions=interactReq.get('disease', None), goRestrictions=interactReq.get('go', None), cellRestrictions=interactReq.get('cells', None))

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def returnInteractions(genes=None, mirnas=None, lncrnas=None, organisms=None, diseaseRestrictions=None, goRestrictions=None, cellRestrictions=None):

    genes = genes if genes!=None else []
    mirnas = mirnas if mirnas!=None else []
    lncrnas = lncrnas if lncrnas!=None else []

    foundRels = defaultdict(list)

    allDocIDs = set()
    allRels = []

    if any([x.startswith("Tester") for x in genes]):

        testerID = [x for x in genes if x.startswith("Tester")][0]

        print("Retrieving testing Interactions for tester", testerID)

        testerRels = testRels.tester2rels[testerID]

        if "all" in testRels.tester2rels:
            testerRels += testRels.tester2rels["all"]

        print(testerID, len(testerRels), testerRels)

        targetsByLid = defaultdict(set)

        for (lid, rid, sid) in testerRels:
            targetsByLid[lid].add( (rid, sid) )


        if mirelPMID != None:

            for lid in targetsByLid:

                dbrels = mirelPMID.get_rels('gene', lid)

                if dbrels == None:
                    continue

                takenRels = []
                allrids = targetsByLid[lid]

                for rel in dbrels:
                    if (rel.rid, rel.assocSent) in allrids:
                        takenRels.append(rel)

                allRels += takenRels

    else:
        #proceed as usual

        allRelsByType = defaultdict(list)

        for relDB in relDBs:

            for gene in genes:

                dbrels = relDB.get_rels('gene', gene)

                if dbrels == None:
                    continue

                allRelsByType['gene'] += dbrels

            for mirna in mirnas:

                dbrels = relDB.get_rels('mirna', mirna)

                if dbrels == None:
                    continue

                allRelsByType['mirna'] += dbrels

            for lncrna in lncrnas:

                dbrels = relDB.get_rels('lncrna', lncrna)

                if dbrels == None:
                    continue

                allRelsByType['lncrna'] += dbrels


        for etype in allRelsByType:

            allRels += allRelsByType[etype]


    allRels = sorted(allRels, key=lambda x: x.docid)
    print("Loading", len(allRels), "relations")

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

        foundRels[(rel.l_id_type, rel.r_id_type)].append(evJSON)

    print("Loaded Sentences", loadedSents)

    addInfo = {}

    if pmid2disease:

        addInfo['disease'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2disease.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['disease'] = allDocInfos

    if pmid2go:

        addInfo['go'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2go.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['go'] = allDocInfos

    if pmid2cell:

        addInfo['cells'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2cell.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['cells'] = allDocInfos

    allowedIDs = defaultdict(list)

    if diseaseRestrictions != None:

        allowedTermIDs = []
        for delem in diseaseRestrictions:
            elemTerm = diseaseObo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        allowedIDs['disease'] = set(allowedTermIDs)

    if goRestrictions != None:

        allowedTermIDs = []
        for delem in goRestrictions:
            elemTerm = goObo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        allowedIDs['go'] = set(allowedTermIDs)

    if cellRestrictions != None:

        allowedTermIDs = []
        for delem in cellRestrictions:
            elemTerm = cellObo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        allowedIDs['cells'] = set(allowedTermIDs)


    allrels = []
    for lent, rent in foundRels:

        lrEvs = foundRels[(lent, rent)]

        okEvs = []

        for jsonEV in lrEvs:

            if diseaseRestrictions != None and len(allowedIDs['disease']) > 0:

                docID = jsonEV['docid']

                docIDDiseaseInfo = addInfo['disease'].get(docID, [])

                if len(docIDDiseaseInfo) == 0:
                    continue
                else:

                    acceptEv = any([x['termid'] in allowedIDs['disease'] for x in docIDDiseaseInfo])

                    if not acceptEv:
                        continue

            if goRestrictions != None and len(allowedIDs['go']) > 0:

                docID = jsonEV['docid']

                docIDGOInfo = addInfo['go'].get(docID, [])

                if len(docIDGOInfo) == 0:
                    continue
                else:

                    acceptEv = any([x['termid'] in allowedIDs['go'] for x in docIDGOInfo])

                    if not acceptEv:
                        continue

            if cellRestrictions != None and len(allowedIDs['cells']) > 0:

                docID = jsonEV['docid']

                docIDCellInfo = addInfo['cells'].get(docID, [])

                if len(docIDCellInfo) == 0:
                    continue
                else:

                    acceptEv = any([x['termid'] in allowedIDs['cells'] for x in docIDCellInfo])

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

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify( [] ), 200, None))

    reMatch = regex.compile(searchWord+'{e<=3}')

    jsonResultByType = defaultdict(set)

    for relDB in relDBs:

        ltype = relDB.ltype
        rtype = relDB.rtype

        for entName in relDB.all_ltypes:
            if len(jsonResultByType[ltype]) <= 100 and reMatch.match(entName):
                jsonResultByType[ltype].add(entName)

        for entName in relDB.all_rtypes:
            if len(jsonResultByType[rtype]) <= 100 and reMatch.match(entName):
                jsonResultByType[rtype].add(entName)

    jsonResult = []
    for type in jsonResultByType:

        jsonResult += list(
            [{'name': interact, 'group': type} for interact in jsonResultByType[type]]
        )


    if reMatch.match("Tester"):

        for testerID in testRels.tester2rels:
            if not testerID.startswith("Tester"):
                continue
            jsonResult.append({'name': testerID, 'group': 'gene'})

    return app.make_response((jsonify( jsonResult ), 200, None))


@app.route('/disease_ac', methods=['GET', 'POST'])
def disease_autocomplete():

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in pmid2disease.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'disease'
                }
            )

        if len(jsonResult) > 100:
            break

    return app.make_response((jsonify(jsonResult), 200, None))

@app.route('/go_ac', methods=['GET', 'POST'])
def go_autocomplete():

    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in pmid2go.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'go'
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

    for (termName, termID) in pmid2cell.getTerms():

        if reMatch.match(termName):
            jsonResult.append(
                {
                    'name': termName,
                    'termid': termID,
                    'group': 'cells'
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

    ltypePOS = feedbackInfo.get('lpos', None)
    rtypePOS = feedbackInfo.get('rpos', None)

    ltype = feedbackInfo['ltype']
    rtype = feedbackInfo['rtype']

    approve = feedbackInfo['approve']

    searchTermGenes = tuple(feedbackInfo['search_terms'].get('gene', []))

    print(dataID, approve)

    mirFeedback.add_feedback( (dataSource, dataID, docID, relSentID, approve, lid, rid, ltype, rtype, ltypePOS, rtypePOS, searchTermGenes) )

    return app.make_response((jsonify({}), 200, None))



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str, help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=True)
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    args = parser.parse_args()

    for x in args.__dict__:
        print(x, args.__dict__[x])

    pmidBase = args.textmine + '/aggregated_pmid/'
    pmcBase = args.textmine + '/aggregated_pmc/'

    normGeneSymbols = normalize_gene_names()

    print("Loading Interactions")

    # allInteractions = defaultdict(list)

    symbol2ensemblDB = SymbolEnsemblDB.loadFromFile(args.obodir + "/sym2ens/")

    mirandaDB_mm10 = MirandaRelDB.loadFromFile(filepath=args.obodir + "/mm10_interactionsAllGenes.txt", symbol2ens=symbol2ensemblDB)
    #mirandaDB_hg38 = MirandaRelDB.loadFromFile(filepath=args.obodir + "/hg38_interactionsAllGenes.txt")

    recordsDB = miRecordDB.loadFromFile(filelocation=args.obodir + "/mirecords_v4.xlsx", normGeneSymbols=normGeneSymbols)
    mirtarbaseDB = MirTarBaseDB.loadFromFile(filepath=args.obodir + "/miRTarBase.csv", normGeneSymbols=normGeneSymbols)

    allDBS = None

    print(datetime.datetime.now(), "Loading PMID2PMC")
    pmid2pmcDB = PMID2PMCDB.loadFromFile(args.textmine + '/pmid2pmc')
    print(datetime.datetime.now(), "Loading mirel")

    testRels = None#TestRelLoader.loadFromFile(pmidBase + "/test_rels_4")

    mirelPMID = MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.hsa.pmid", ltype="gene", rtype="mirna", normGeneSymbols=normGeneSymbols)
    lncMirPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/lncrna_mirna.cur.pmid", ltype="lncrna", rtype="mirna")
    geneLncPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/gene_lncrna.cur.pmid", ltype="gene", rtype="lncrna")

    relDBs = [recordsDB, mirtarbaseDB, mirelPMID, lncMirPMID, geneLncPMID, mirandaDB_mm10]

    relDBs = [x for x in relDBs if x != None]

    print(datetime.datetime.now(), "Finished mirel")

    mirFeedback = feedbackDB(args.feedback)

    print(datetime.datetime.now(), "Loading sents")
    sentDB = SentenceDB.loadFromFile(args.sentdir,
                                    args.textmine+ "/aggregated_pmid/pmid2sent")
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

        pmid2go = allDBS[0]
        pmid2disease = allDBS[1]
        pmid2fma = allDBS[2]
        pmid2cell = allDBS[3]

        print(datetime.datetime.now(), "Loading pickle ended")

    print(datetime.datetime.now(), "Loading ontologies")

    diseaseObo = GeneOntology(args.obodir + "/doid.obo")
    goObo = GeneOntology(args.obodir + "/go.obo")
    cellObo = GeneOntology(args.obodir + "/meta_cells.obo")

    print(datetime.datetime.now(), "Loading ontologies finished")


    if allDBS == None:
        pmid2go = None
        pmid2disease = None
        pmid2fma = None
        pmid2cell = None

        print(datetime.datetime.now(), "Loading GO")
        pmid2go = PMID2XDB.loadFromFile(pmidBase + "/go.pmid", goObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading Disease")
        pmid2disease = PMID2XDB.loadFromFile(pmidBase + "/disease.pmid", diseaseObo, requiredPMIDs)
        # print(datetime.datetime.now(), "Loading FMA")
        # pmid2fma = PMID2XDB.loadFromFile(pmidBase + "/model_anatomy.pmid")
        print(datetime.datetime.now(), "Loading cellline")
        pmid2cell = PMID2XDB.loadFromFile(pmidBase + "/celllines.pmid", cellObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading mirna")

        allDBS = (pmid2go, pmid2disease, pmid2fma, pmid2cell)

        print(datetime.datetime.now(), "Writing Pickle")

        with open(pmidBase + "/dbs.pickle", 'wb') as fout:
            pickle.dump(allDBS, fout)

    print(datetime.datetime.now(), "Loading finished")


    print("Starting Flask on port", args.port)

    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)
