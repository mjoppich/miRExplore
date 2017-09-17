from Bio import Entrez

class CoCitationStore:

    def __init__(self):
        Entrez.email = "joppich@bio.ifi.lmu.de"


    def _extractLinkSetData( self, resultData ):

        retData = dict()

        for i in range(0, len(resultData)):
            linkData = resultData[i]
            queryID = linkData['IdList'][0]

            linkSetDB = linkData["LinkSetDb"]

            if linkSetDB == None or len(linkSetDB) == 0:
                retData[queryID] = []
            else:
                linkIDs = [link["Id"] for link in linkData["LinkSetDb"][0]["Link"]]
                retData[queryID] = linkIDs

        return retData

    def _makeChunks(self, someiterable, n=1000):

        somelist = list(someiterable)

        allRetLists = []

        for i in range(0, len(somelist), n):
            end = min([i+n, len(somelist)])
            allRetLists.append( somelist[i:end] )

        if not len(somelist) == sum([len(x) for x in allRetLists]):
            print("Error creating chunks")
            exit(-1)

        return allRetLists

    def _queryEntrez(self, linkname, pmids):

        chunks = self._makeChunks(pmids)
        finalReturn = {}

        print("Received " + str(len(chunks)) + " chunks")
        for chunk in chunks:
            cites = Entrez.read(Entrez.elink(dbfrom="pubmed", LinkName=linkname, from_uid=chunk))
            citesIDs = self._extractLinkSetData(cites)

            for x in citesIDs:
                finalReturn[x] = citesIDs[x]

        return finalReturn

    def getCites(self, pmids):

        if type(pmids) == str:
            pmids = [pmids]

        citesIDs =self._queryEntrez('pubmed_pubmed_refs', pmids)
        return citesIDs

    def getCitedBy(self, pmids):

        if type(pmids) == str:
            pmids = [pmids]

        citesIDs =self._queryEntrez('pubmed_pubmed_citedin', pmids)
        return citesIDs


if __name__ == '__main__':

    store = CoCitationStore()

    pmid = ['25958325', '27984116']

    print(store.getCites(pmid))
    print(store.getCitedBy(pmid))