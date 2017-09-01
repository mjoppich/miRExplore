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

    def getCites(self, pmids):

        if type(pmids) == str:
            pmids = [pmids]

        cites = Entrez.read(Entrez.elink(dbfrom="pubmed", LinkName="pubmed_pubmed_refs", from_uid=pmids))
        citesIDs = self._extractLinkSetData(cites)

        return citesIDs

    def getCitedBy(self, pmids):

        if type(pmids) == str:
            pmids = [pmids]

        cites = Entrez.read(Entrez.elink(dbfrom="pubmed", LinkName="pubmed_pubmed_citedin", from_uid=pmids))
        citesIDs = self._extractLinkSetData(cites)

        return citesIDs


if __name__ == '__main__':

    store = CoCitationStore()

    pmid = ['25958325', '27984116']

    print(store.getCites(pmid))
    print(store.getCitedBy(pmid))