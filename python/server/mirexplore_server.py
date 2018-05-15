import re
import sys
import os

from textdb.PMID2PMCDB import PMID2PMCDB
from textdb.PMID2XDB import PMID2XDB

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../../")

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

print("Loading Interactions")

allInteractions = defaultdict(list)

mirecords = miRecordDB.from_xslx()
for elem in mirecords.elems:
    allInteractions[(elem[0].upper(), elem[1])].append(('MIRECORD', elem[2]))


pmidBase = '/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmid/'
pmcBase = '/mnt/c/ownCloud/data/miRExplore/textmine/aggregated_pmc/'


pmid2pmcDB = PMID2PMCDB.loadFromFile('/mnt/c/ownCloud/data/miRExplore/pmid2pmc')

pmid2go = PMID2XDB.loadFromFile(pmidBase + "/go.pmid")
pmid2disease = PMID2XDB.loadFromFile(pmidBase + "/disease.pmid")
pmid2fma = PMID2XDB.loadFromFile(pmidBase + "/model_anatomy.pmid")
pmid2cell = PMID2XDB.loadFromFile(pmidBase + "/cellline.pmid")



mirelPMID = MiGenRelDB.loadFromFile()


print("Loading finished")



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

    return returnInteractions(interactReq.get('gene', None), interactReq.get('mirna', None))

@app.route('/status')
def getStatus():
    return app.make_response((jsonify({'error': 'must include homid'}), 400, None))


def returnInteractions(genes=None, mirnas=None):

    jsonResult = defaultdict(list)
    jsonResultByDatabase=defaultdict(list)

    if genes != None and mirnas == None:

        for interaction in allInteractions:

            if interaction[0] in genes:

                for elem in allInteractions[interaction]:


                    thisElem = {
                        'gene': interaction[0],
                        'mirna': interaction[1],
                        'database': elem[0],
                        'database_id': elem[1]
                    }

                    jsonResult[interaction[0]].append( thisElem )
                    jsonResultByDatabase[thisElem['database']].append(thisElem)


    returnObj = {
        'databaseview': jsonResultByDatabase,
        'elemview': jsonResult
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

    reMatch = re.compile(searchWord)

    jsonResultGene = set()
    jsonResultMIRNA = set()

    for interact in allInteractions:

        if reMatch.match(interact[0]):
            jsonResultGene.add(interact[0])

        if reMatch.match(interact[1]):
            jsonResultMIRNA.add(interact[1])

    jsonResult = list([{'name': interact, 'group': 'gene'} for interact in jsonResultGene])
    jsonResult += list([{'name': interact, 'group': 'mirna'} for interact in jsonResultMIRNA])


    return app.make_response((jsonify( jsonResult ), 200, None))



if __name__ == '__main__':

   print([rule.rule for rule in app.url_map.iter_rules() if rule.endpoint !='static'])

   app.run(threaded=True)
