import argparse
import datetime
import pickle
import regex
import sys
import os

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


from io import StringIO

from flask import Flask, jsonify, request, redirect, url_for, send_from_directory
import json
import pprint
from collections import defaultdict

from flask_cors import CORS

from analysis.miRecordDB import miRecordDB

app = Flask(__name__)
CORS(app)

app.config['DEBUG'] = False
app.config['UPLOAD_FOLDER'] = ""


# For a given file, return whether it's an allowed type or not
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

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


@app.route('/get_interactions/gene/<geneID>')
def getHomCluster(geneID):

    if geneID == None:
        return app.make_response((jsonify( {'error': 'invalid homID'} ), 400, None))

    return returnInteractions([geneID])


@app.route('/find_interactions', methods=['GET', 'POST'])
def findInteractions():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'gene' in interactReq and not 'mirna' in interactReq:
        return app.make_response((jsonify( {'error': 'must include gene or mirna'} ), 400, None))

    print(interactReq)

    return returnInteractions(interactReq.get('gene', None), interactReq.get('mirna', None), diseaseRestrictions=interactReq.get('disease', None))

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def returnInteractions(genes=None, mirnas=None, lncrnas=None, diseaseRestrictions=None, goRestrictions=None, cellRestrictions=None):

    genes = genes if genes!=None else []
    mirnas = mirnas if mirnas!=None else []
    lncrnas = lncrnas if lncrnas!=None else []

    foundRels = defaultdict(list)

    allDocIDs = set()


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


    allRels = []
    for etype in allRelsByType:

        allRels += allRelsByType[etype]


    allRels = sorted(allRels, key=lambda x: x.docid)

    for rel in allRels:

        if rel.docid != -1:
            allDocIDs.add(rel.docid)

        evJSON = rel.toJSON()
        evSent = evJSON.get('rel_sentence', None)

        if evSent != None:
            evSentTxt = sentDB.get_sentence(evSent)

            if evSentTxt != None:
                evJSON['sentence'] = evSentTxt[1]

        foundRels[(rel.lid, rel.rid)].append(evJSON)

    addInfo = {}

    if pmid2disease:

        addInfo['disease'] = None

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmid2disease.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        addInfo['disease'] = allDocInfos


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


    allrels = []
    for lid, rid in foundRels:

        lrEvs = foundRels[(lid, rid)]

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

            okEvs.append(jsonEV)


        if len(okEvs) > 0:

            allrels.append({'lid':lid,
                            'rid': rid,
                            'evidences': okEvs
                            })

    returnObj = {
        'rels': allrels,
        'pmidinfo': addInfo
    }

    return app.make_response((jsonify(returnObj), 200, None))


@app.route('/organisms')
def getOrganisms():
    return app.make_response((jsonify(['human', 'mouse', 'Homo sapiens', 'Mus musculus']), 200, None))

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

    mirFeedback.add_feedback( (dataSource, dataID, docID, relSentID, approve, lid, rid, ltype, rtype, ltypePOS, rtypePOS) )

    return app.make_response((jsonify({}), 200, None))



if __name__ == '__main__':


    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str, help='Base for Textmining. Includes aggregated_ and results folder', required=True)
    parser.add_argument('-o', '--obodir', type=str, help='Path to all obo-files/existing databases', required=True)
    parser.add_argument('-s', '--sentdir', type=str, help='Path to sentences', required=True)
    parser.add_argument('-f', '--feedback', type=str, help="Path for feedback stuff", required=True)
    parser.add_argument('-p', '--port', type=int, help="port to run on", required=False, default=5000)

    args = parser.parse_args()

    pmidBase = args.textmine + '/aggregated_pmid/'
    pmcBase = args.textmine + '/aggregated_pmc/'

    print("Loading Interactions")

    # allInteractions = defaultdict(list)

    mirandaDB = MirandaRelDB.loadFromFile(filepath=args.obodir + "/miranda_test.tsv")


    recordsDB = None#miRecordDB.loadFromFile(filelocation=args.obodir + "/mirecords_v4.xlsx")
    mirtarbaseDB = None#MirTarBaseDB.loadFromFile(filepath=args.obodir + "/miRTarBase.csv")

    # for elem in mirecords.elems:
    #    allInteractions[(elem[0].upper(), elem[1])].append(('MIRECORD', elem[2]))

    allDBS = None

    print(datetime.datetime.now(), "Loading PMID2PMC")
    pmid2pmcDB = PMID2PMCDB.loadFromFile(args.textmine + '/pmid2pmc')
    print(datetime.datetime.now(), "Loading mirel")

    mirelPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.spacy.pmid", ltype="gene", rtype="mirna")
    lncMirPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/lncrna_mirna.cur.pmid", ltype="lncrna", rtype="mirna")
    geneLncPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/gene_lncrna.cur.pmid", ltype="gene", rtype="lncrna")

    relDBs = [recordsDB, mirtarbaseDB, mirelPMID, lncMirPMID, geneLncPMID, mirandaDB]

    relDBs = [x for x in relDBs if x != None]

    print(datetime.datetime.now(), "Finished mirel")

    mirFeedback = feedbackDB(args.feedback)

    print(datetime.datetime.now(), "Loading sents")
    sentDB = SentenceDB.loadFromFile(args.sentdir,
                                    args.textmine+ "aggregated_pmid/pmid2sent")
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

    diseaseObo = GeneOntology(args.obodir + "/doid.obo")

    if allDBS == None:
        pmid2go = None
        pmid2disease = None
        pmid2fma = None
        pmid2cell = None

        # print(datetime.datetime.now(), "Loading GO")
        # pmid2go = PMID2XDB.loadFromFile(pmidBase + "/go.pmid")
        print(datetime.datetime.now(), "Loading Disease")
        pmid2disease = PMID2XDB.loadFromFile(pmidBase + "/disease.pmid", diseaseObo, requiredPMIDs)
        # print(datetime.datetime.now(), "Loading FMA")
        # pmid2fma = PMID2XDB.loadFromFile(pmidBase + "/model_anatomy.pmid")
        # print(datetime.datetime.now(), "Loading cellline")
        # pmid2cell = PMID2XDB.loadFromFile(pmidBase + "/cellline.pmid")
        print(datetime.datetime.now(), "Loading mirna")

        allDBS = (pmid2go, pmid2disease, pmid2fma, pmid2cell)

        print(datetime.datetime.now(), "Writing Pickle")

        with open(pmidBase + "/dbs.pickle", 'wb') as fout:
            pickle.dump(allDBS, fout)

    print(datetime.datetime.now(), "Loading finished")


    print("Starting Flask on port", args.port)

    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)
