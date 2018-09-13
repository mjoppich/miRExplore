import sys, os
sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


from collections import Counter, defaultdict

import requests
import json

from synonymes.mirnaID import miRNA, miRNAPART


class DBAccess:

    @classmethod
    def fetchGenes(cls):
        serverAddress = "https://turingwww.bio.ifi.lmu.de"
        serverPort = None
        serverPath = "yancDB"

        serverAddress = "http://localhost"
        serverPort = "5000"
        serverPath = "/"

        def makeServerAddress(address, port, path):

            ret = address

            if port != None:
                ret += ":" + str(port)

            if path != None:
                ret += "/" + path + "/"

            return ret

        r = requests.post(makeServerAddress(serverAddress, serverPort, serverPath) + "/genes")

        jsonRes = r.json()

        return jsonRes['gene']

    @classmethod
    def fetchGeneInfos(cls, requestDict, acceptRel):
        serverAddress = "https://turingwww.bio.ifi.lmu.de"
        serverPort = None
        serverPath = "yancDB"

        serverAddress = "http://localhost"
        serverPort = "5000"
        serverPath = "/"

        def makeServerAddress(address, port, path):

            ret = address

            if port != None:
                ret += ":" + str(port)

            if path != None:
                ret += "/" + path + "/"

            return ret

        requestDict['sentences'] = "false"

        r = requests.post(makeServerAddress(serverAddress, serverPort, serverPath) + "/find_interactions",
                          data=json.dumps(requestDict))


        jsonRes = r.json()



        filteredRels = []

        for rel in jsonRes['rels']:

            if acceptRel(rel):
                filteredRels.append(rel)



        return filteredRels



if __name__ == '__main__':

    allGenes = DBAccess.fetchGenes()


    def acceptRel(rel):


        evByType = Counter()

        for ev in rel['evidences']:

            ds = ev['data_source']
            evByType[ds] += 1

        if 'mirwalk' in evByType and len(evByType) >= 2:
            return True

        return False



    print("Running queries for", len(allGenes), "genes", file=sys.stderr)

    for i in range(0, len(allGenes), 1000):

        chunkGenes = set()

        for j in range(0, 1000):

            if i+j >= len(allGenes):
                continue

            chunkGenes.add(allGenes[i+j])

        print("Running query for", i, "to", i+len(chunkGenes), file=sys.stderr)

        requestDict = {
            "gene": list(chunkGenes),
            "organisms": [{'termid': "Mus musculus"}]
        }

        rels = DBAccess.fetchGeneInfos(requestDict, acceptRel)

        for rel in rels:

            pmids = set()
            mirTarBases = set()
            mirecordsCount = 0
            dianaCount = 0
            orgs = set()

            mwBS = []
            mwBP = []
            mwEns = []

            for ev in rel['evidences']:

                if 'orgs' in ev:

                    for x in ev['orgs']:
                        orgs.add(x)

                if ev['data_source'] == "DIANA":
                    dianaCount += 1
                elif ev['data_source'] == "miRTarBase":
                    mirTarBases.add(ev['data_id'])
                elif ev['data_source'] == "pmid":
                    pmids.add(ev['docid'])
                elif ev['data_source'] == "mirecords":
                    mirecordsCount += 1
                elif ev['data_source'] == "mirwalk":

                    mwBS.append( [str(x) for x in ev['bind_position']] )
                    mwBP.append( str(ev['bind_probability']) )
                    mwEns.append( ev['target_ensembl_id'])

                else:
                    print("Unhandled data source", ev['data_source'])

            dr = {

                'gene': rel['lid'],
                'mirna': rel['rid'],
                'pmid_count': len(pmids),
                'pmids': ";".join(pmids),
                'mirtarbase_count': len(mirTarBases),
                'mirecords': mirecordsCount,
                'diana_count': dianaCount,
                'orgs': ";".join(orgs),
                'bsite': ";".join( ["-".join(x) for x in mwBS]),
                'bprob': ";".join([str(x) for x in mwBP]),
                'ensembl': ";".join(mwEns)
            }

            print(dr['gene'], dr['mirna'], dr['pmid_count'], dr['mirtarbase_count'], dr['mirecords'], dr['diana_count'], dr['orgs'], dr['bsite'], dr['bprob'], dr['ensembl'], dr['pmids'], sep="\t")
