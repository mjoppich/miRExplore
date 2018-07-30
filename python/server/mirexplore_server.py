import argparse
import datetime
import pickle
import regex
import sys
import os
import shlex


from synonymes.SynonymFile import Synfile
from synonymes.mirnaID import miRNA, miRNAPART
from textdb.SymbolEnsemblDB import SymbolEnsemblDB
from textdb.featureviewer import FeatureViewer
from textdb.rfamDB import RFamDB
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
from collections import defaultdict, Counter

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



@app.route('/hfi_general', methods=['POST'])
def hfi_general():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify({'error': 'invalid json'}), 400, None))

    entType = interactReq.get('type', "gene")

    if not entType.upper() in ['GENE', 'MIRNA']:
        return app.make_response((jsonify({'error': 'invalid enttity type ' + entType}), 400, None))

    mirnas = []
    genes = []

    if entType.upper() == 'GENE':
        entName = interactReq.get('name', None)
        if entName != None:
            genes.append(entName)

    elif entType.upper() == 'MIRNA':

        entName = interactReq.get('name', None)
        if entName != None:
            mirnas.append(entName)

    if all([len(x) == 0 for x in [mirnas, genes]]):
        return app.make_response((jsonify({'error': 'no entity names given'}), 400, None))

    relObj = returnInteractions(genes, mirnas, None, loadSentences=False)
    rels = relObj['rels']

    seenInteractors = defaultdict(set)
    seenPMIDs = set()
    seenEvidenceTypes = Counter()

    for rel in rels:
        """
                    allrels.append({'lid':lent[0],
                            'rid': rent[0],
                            'ltype': lent[1],
                            'rtype': rent[1],
                            'evidences': okEvs
                            })
        """

        if entType.upper() == 'GENE':

            for ev in rel['evidences']:

                docid = ev.get('docid', None)

                if docid != None:
                    seenPMIDs.add(docid)

                data_source = ev.get('data_source', None)

                seenEvidenceTypes[data_source] += 1


            if rel['ltype'].upper() == entType.upper():
                seenInteractors[rel['rtype']].add( rel['rid'] )

            else:
                seenInteractors[rel['ltype']].add( rel['lid'] )



    seenMirnaIDs = set()
    if entType.upper() == 'GENE':

        for x in seenInteractors['mirna']:

            txtMirna = x

            try:
                mirnaID = miRNA(txtMirna)

                mirID = mirnaID.getPart(miRNAPART.ID, None)

                mirID = int(mirID)

                seenMirnaIDs.add(mirID)

            except:
                seenMirnaIDs.add(x)




    answer = {
        'search': genes,
        'interactor_count': len(seenInteractors),
        'interactor_types': [x for x in seenInteractors],
        'data_source_count': len(seenEvidenceTypes),
        'data_sources': [x for x in seenEvidenceTypes],
        'evidence_count': sum([seenEvidenceTypes[x] for x in seenEvidenceTypes]),
        'interactors': list(seenMirnaIDs)
    }

    return app.make_response((jsonify(answer), 200, None))







@app.route('/get_features/<geneID>')
def getGeneFeatures(geneID):

    if geneID == None:
        return app.make_response((jsonify( {'error': 'invalid geneid'} ), 400, None))

    return getGeneMirnaFeatures(geneID, None)

@app.route('/get_features/<geneID>/<mirnaID>')
def getGeneMirnaFeatures(geneID, mirnaID):

    global featureViewer
    global symbol2ensemblDB


    if geneID == None:
        return app.make_response((jsonify( {'error': 'invalid geneid'} ), 400, None))


    resObj = {
        'gene_id': geneID
    }


    if mirnaID != None:
        resObj['mirna_id'] = mirnaID

    ret = []

    collectGenes = symbol2ensemblDB.get_all_genes(geneID)

    for org in collectGenes:
        geneEnsIDs = collectGenes[org]

        if geneEnsIDs == None:
            continue


        for geneEnsID in geneEnsIDs:

            lret = featureViewer.getFeatures(geneEnsID, mirnaID)

            if len(lret) > 0:
                ret += list(lret)

    if ret != None:
        resObj['features'] = ret


    return app.make_response((jsonify(resObj), 200, None))




@app.route('/find_interactions', methods=['GET', 'POST'])
def findInteractions():
    interactReq = request.get_json(force=True, silent=True)

    if interactReq == None:
        return app.make_response((jsonify( {'error': 'invalid json'} ), 400, None))

    if not 'gene' in interactReq and not 'mirna' in interactReq:
        return app.make_response((jsonify( {'error': 'must include gene or mirna'} ), 400, None))

    print(interactReq)


    gene = interactReq.get('gene', None)
    mirna = interactReq.get('mirna', None)
    lncrna = interactReq.get('lncrna', None)
    organisms = interactReq.get('organisms', None)
    diseases = interactReq.get('disease', None)
    goes = interactReq.get('go', None)
    cells = interactReq.get('cells', None)
    ncits = interactReq.get('ncits', None)

    loadSentences = interactReq.get('sentences', "true").upper() == "TRUE"

    retObj = returnInteractions(gene, mirna, lncrna, organisms=organisms, diseaseRestrictions=diseases, goRestrictions=goes, cellRestrictions=cells, ncitRestrictions=ncits, loadSentences=loadSentences)

    return app.make_response((jsonify(retObj), 200, None))


@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def returnInteractions(genes=None, mirnas=None, lncrnas=None, organisms=None, diseaseRestrictions=None, goRestrictions=None, cellRestrictions=None, ncitRestrictions=None, loadSentences=True):

    global mirFeedback
    global mirandaDB_mm10
    global relDBs
    global diseaseObo
    global goObo
    global cellObo
    global ncitObo
    global pmid2go
    global pmid2disease
    global pmid2fma
    global pmid2cell
    global sentDB

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

        if evSent != None and loadSentences:
            evSentTxt = sentDB.get_sentence(evSent)

            if evSentTxt != None:

                loadedSents += 1
                evJSON['sentence'] = evSentTxt[1]

        foundRels[(rel.l_id_type, rel.r_id_type)].append(evJSON)

    print("Loaded Sentences", loadedSents)

    addInfo = {}

    def makeAddInfo( pmidFile ):

        allDocInfos = {}

        for docid in allDocIDs:
            docInfo = pmidFile.getDOC(docid)

            if docInfo != None:
                allDocInfos[docid] = docInfo

        return allDocInfos

    if pmid2disease:
        addInfo['disease'] = makeAddInfo(pmid2disease)


    if pmid2go:
        addInfo['go'] = makeAddInfo(pmid2go)


    if pmid2cell:
        addInfo['cells'] = makeAddInfo(pmid2cell)

    if pmid2ncit:
        addInfo['ncits'] = makeAddInfo(pmid2ncit)

    #if pmid2fma:
    #    addInfo['fma'] = makeAddInfo(pmid2fma)


    def getAllowedTermIDs(restrictions, obo):

        allowedTermIDs = []
        for delem in restrictions:
            elemTerm = obo.getID(delem['termid'])

            if elemTerm == None:
                continue

            elemTerms = [x.term.id for x in elemTerm.getAllChildren()] + [elemTerm.id]
            allowedTermIDs += elemTerms

        return set(allowedTermIDs)



    allowedIDs = defaultdict(list)

    if diseaseRestrictions != None:
        allowedIDs['disease'] = getAllowedTermIDs(diseaseRestrictions, diseaseObo)

    if goRestrictions != None:
        allowedIDs['go'] = getAllowedTermIDs(goRestrictions, goObo)

    if cellRestrictions != None:
        allowedIDs['cells'] = getAllowedTermIDs(cellRestrictions, cellObo)

    if ncitRestrictions != None:
        allowedIDs['ncits'] = getAllowedTermIDs(ncitRestrictions, ncitObo)



    allrels = []
    for lent, rent in foundRels:

        lrEvs = foundRels[(lent, rent)]

        okEvs = []

        for jsonEV in lrEvs:


            def checkEvID(infoID):

                if len(addInfo[infoID]) == 0:
                    return True

                docID = jsonEV['docid']

                infoIDInfos = addInfo[infoID].get(docID, [])

                if len(infoIDInfos) == 0:
                    return False
                else:

                    acceptEv = any([x['termid'] in allowedIDs[infoID] for x in infoIDInfos])

                    if not acceptEv:
                        return False

                return True



            if diseaseRestrictions != None and not checkEvID('disease'):
                continue

            if goRestrictions != None and not checkEvID('go'):
                continue

            if cellRestrictions != None and not checkEvID('cells'):
                continue

            if ncitRestrictions != None and not checkEvID('ncits'):
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

    return returnObj


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

    global relDBs
    global testRels

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

    global pmid2disease

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


@app.route('/ncit_ac', methods=['GET', 'POST'])
def ncit_autocomplete():

    global pmid2ncit


    searchWords = request.get_json(force=True, silent=True)
    searchWord = searchWords['search']

    if searchWord == None or len(searchWord) < 2:
        return app.make_response((jsonify([]), 200, None))

    reMatch = regex.compile(searchWord + '{e<=3}')

    jsonResult = list()

    for (termName, termID) in pmid2ncit.getTerms():

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

@app.route('/go_ac', methods=['GET', 'POST'])
def go_autocomplete():

    global pmid2go


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

    global pmid2cell

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

    global mirFeedback

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


mirFeedback = None
mirandaDB_mm10 = None
relDBs = None
diseaseObo = None
goObo = None
cellObo = None
ncitObo = None
pmid2go = None
pmid2disease = None
pmid2fma = None
pmid2cell = None
pmid2ncit = None
testRels = None
mirelPMID = None
sentDB=None
featureViewer = None
symbol2ensemblDB = None


def start_app_from_args(args):


    global mirFeedback
    global mirandaDB_mm10
    global relDBs
    global diseaseObo
    global goObo
    global cellObo
    global pmid2go
    global pmid2disease
    global pmid2fma
    global pmid2cell
    global testRels
    global mirelPMID
    global sentDB
    global featureViewer
    global symbol2ensemblDB
    global pmid2ncit
    global ncitObo

    pmidBase = args.textmine + '/aggregated_pmid/'
    pmcBase = args.textmine + '/aggregated_pmc/'

    normGeneSymbols = normalize_gene_names(path=args.obodir + "/hgnc_no_withdrawn.syn")

    print("Loading Interactions")

    # allInteractions = defaultdict(list)

    symbol2ensemblDB = SymbolEnsemblDB.loadFromFile(args.obodir + "/sym2ens/")

    mirandaDB_mm10 = MirandaRelDB.loadFromFile(filepath=args.obodir + "/mm10_interactionsAllGenes.txt", symbol2ens=symbol2ensemblDB)
    # mirandaDB_hg38 = MirandaRelDB.loadFromFile(filepath=args.obodir + "/hg38_interactionsAllGenes.txt")

    recordsDB = miRecordDB.loadFromFile(filelocation=args.obodir + "/mirecords_v4.xlsx", normGeneSymbols=normGeneSymbols)
    mirtarbaseDB = MirTarBaseDB.loadFromFile(filepath=args.obodir + "/miRTarBase.csv", normGeneSymbols=normGeneSymbols)

    allDBS = None

    print(datetime.datetime.now(), "Loading PMID2PMC")
    #pmid2pmcDB = PMID2PMCDB.loadFromFile(args.textmine + '/pmid2pmc')
    print(datetime.datetime.now(), "Loading mirel")

    testRels = None  # TestRelLoader.loadFromFile(pmidBase + "/test_rels_4")

    mirelPMID = MiGenRelDB.loadFromFile(pmidBase + "/mirna_gene.hsa.pmid", ltype="mirna", rtype="gene", normGeneSymbols=normGeneSymbols, switchLR=True)
    lncMirPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/lncrna_mirna.pmid", ltype="lncrna", rtype="mirna")
    geneLncPMID = None#MiGenRelDB.loadFromFile(pmidBase + "/gene_lncrna.pmid", ltype="gene", rtype="lncrna")

    relDBs = [recordsDB, mirtarbaseDB, mirelPMID, lncMirPMID, geneLncPMID, mirandaDB_mm10]

    relDBs = [x for x in relDBs if x != None]

    print(datetime.datetime.now(), "Finished mirel")

    mirFeedback = feedbackDB(args.feedback)

    print(datetime.datetime.now(), "Loading sents")
    sentDB = SentenceDB.loadFromFile(args.sentdir, pmidBase + "/pmid2sent")
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
        pmid2ncit = allDBS[4]

        print(datetime.datetime.now(), "Loading pickle ended")

    print(datetime.datetime.now(), "Loading ontologies")

    diseaseObo = GeneOntology(args.obodir + "/doid.obo")
    goObo = GeneOntology(args.obodir + "/go.obo")
    cellObo = GeneOntology(args.obodir + "/meta_cells.obo")
    ncitObo = GeneOntology(args.obodir + "/ncit.obo")
    fmaObo = GeneOntology(args.obodir + "/fma_obo.obo")


    print(datetime.datetime.now(), "Loading ontologies finished")

    if allDBS == None:
        pmid2go = None
        pmid2disease = None
        pmid2fma = None
        pmid2cell = None
        pmid2ncit = None

        print(datetime.datetime.now(), "Loading GO")
        pmid2go = PMID2XDB.loadFromFile(pmidBase + "/go.pmid", goObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading Disease")
        pmid2disease = PMID2XDB.loadFromFile(pmidBase + "/disease.pmid", diseaseObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading FMA")
        pmid2fma = PMID2XDB.loadFromFile(pmidBase + "/model_anatomy.pmid", fmaObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading cellline")
        pmid2cell = PMID2XDB.loadFromFile(pmidBase + "/celllines.pmid", cellObo, requiredPMIDs)
        print(datetime.datetime.now(), "Loading ncit")
        pmid2ncit = PMID2XDB.loadFromFile(pmidBase + "/ncit.pmid", ncitObo, requiredPMIDs)

        allDBS = (pmid2go, pmid2disease, pmid2fma, pmid2cell, pmid2ncit)

        print(datetime.datetime.now(), "Writing Pickle")

        with open(pmidBase + "/dbs.pickle", 'wb') as fout:
            pickle.dump(allDBS, fout)

        print(datetime.datetime.now(), "Finished Writing Pickle")

    print(datetime.datetime.now(), "Loading Features")
    rfDB = RFamDB.loadFromFile(pmidBase + '/rfam.regions.mirexplore')
    featureViewer = FeatureViewer('mmu', args.obodir, rfamDB=rfDB)

    print(datetime.datetime.now(), "Loading finished")
def getCLParser():
    parser = argparse.ArgumentParser(description='Start miRExplore Data Server', add_help=False)
    parser.add_argument('-t', '--textmine', type=str, help='Base for Textmining. Includes aggregated_ and results folder', required=True)
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

    print("Starting Flask on port", args.port)

    start_app_from_args(args)

    print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint != 'static'])
    app.run(threaded=True, host="0.0.0.0", port=args.port)

def gunicorn_start( datadir="/home/proj/biosoft/ws/projekte/dataintegration/yancDB",
                    sentdir="/home/proj/biocluster/praktikum/neap_ss18/neapss18_noncoding/pmid_sent",
                    feedbackFile=None):

    parser = getCLParser()

    argstr = "--textmine {datadir} --obodir {datadir}/obodir --sentdir {sentdir} --feedback {feedbackFile}".format(datadir=datadir, sentdir=sentdir, feedbackFile=feedbackFile)

    print("Starting app with")
    print(argstr)

    args = parser.parse_args(shlex.split(argstr))

    start_app_from_args(args)

    return app