import argparse
import json
import sys, os
from lxml import etree
from xml.dom import minidom

sys.path.insert(0, str(os.path.dirname(os.path.realpath(__file__))) + "/../")


import requests

from textmining.SentenceID import SentenceID

parser = argparse.ArgumentParser(description='db query', add_help=False)
parser.add_argument('-o', '--output', type=argparse.FileType("w"), help='outfile', default=None, required=False)

args = parser.parse_args()


query = {"elements": [{"group": "NEUTROPHIL", "name": "neutrophils", "termid": "neutrophils"}], "sentences": "true", "obolevel": 3, "messenger_obolevel": 1}
r = requests.post("http://localhost:65522/query", data=json.dumps(query))

res = json.loads(r.content.decode())

sep = "\t"

allPMIDs = set()
allPMC = set()

print("received", len(res['rels']), "relations", file=sys.stderr)

for rel in res['rels']:

    for ev in rel['evidences']:

        docid = ev['docid']

        if docid == None:
            continue

        if docid.startswith("PMC"):
            allPMC.add(docid.replace('PMC', ""))
        else:
            allPMIDs.add(docid)

print("PMIDs", len(allPMIDs))
print("PMCs", len(allPMC))

from Bio import Entrez

Entrez.email = "joppich@bio.ifi.lmu.de"

pmid2info = {}

for dbName, allIDs in [('pmc', allPMC), ('pubmed', allPMIDs)]:

    epostResult = Entrez.read(Entrez.epost(dbName, id=",".join(allIDs)))
    webEnv = epostResult['WebEnv']
    queryKey = epostResult['QueryKey']

    handle = Entrez.efetch(db=dbName, webenv=webEnv, query_key=queryKey, retmode='XML')

    allAuthorArticles = []

    if dbName == 'pmc':
        root = etree.parse(handle)

        for article in root.findall('.//article/'):

            articleIDs = [x.text for x in article.findall(".//article-id[@pub-id-type='pmc']")]

            if len(articleIDs) == 0:
                continue

            link = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC" + articleIDs[0] if len(articleIDs) > 0 else ""


            articleTitles = [x for x in article.findall('.//title-group')]
            articleTitle = ""

            for titleGroup in articleTitles:

                tgText = "".join([x.strip() for x in titleGroup.itertext()])
                articleTitle += tgText


            articleJournals = article.findall('.//journal-title')
            articleJournal = "" if len(articleJournals) == 0 else articleJournals[0].text

            articleYears = article.findall('.//pub-date/year')
            articleYear = -1 if len(articleYears) == 0 else articleYears[0].text


            pmid2info[articleIDs[0]] = (link, articleTitle, articleJournal, articleYear)


    else:

        record = Entrez.read(handle)

        for article in record['PubmedArticle']:

            pubmedID = article['PubmedData']['ArticleIdList'][0] if len(
                article['PubmedData']['ArticleIdList']) > 0 else "-1"
            pubID = int(pubmedID)

            artInfo = article['MedlineCitation']['Article']
            articleTitle = artInfo['ArticleTitle']
            articleJournal = artInfo['Journal']['Title'] if 'Journal' in artInfo else ''

            artDate = artInfo["ArticleDate"]

            articleYear = -1

            if len(artDate) > 0:
                articleYear = artDate[0]["Year"]

            if articleYear == -1 and 'Journal' in artInfo:

                articleJournalInfo = artInfo["Journal"]
                if "JournalIssue" in articleJournalInfo:
                    if 'PubDate' in articleJournalInfo['JournalIssue']:
                        if 'Year' in articleJournalInfo['JournalIssue']['PubDate']:
                            articleYear = articleJournalInfo['JournalIssue']['PubDate']['Year']


            infotuple = ("https://www.ncbi.nlm.nih.gov/pubmed/"+pubmedID, articleTitle, articleJournal, articleYear)

            pmid2info[pubmedID] = infotuple

elemcount = 0

from openpyxl import Workbook
wb = Workbook()
ws = wb.active

header = ["PMID", "Sent ID", "Sentence", "Verb Structure", "Left ID", "Left Ontology ID", "Right ID", "Right Ontology ID", "Link", "Title", "Journal", "Year", "Message ID", "Message Ontology ID", "Effect ID", "Effect Ontology ID", "Verb", "Stack", "Relex", "Conj"]
print(header, sep=sep)
ws.append(header)

relEvCount = 0
evCount = 0
evPMID = set()
evSent = set()

for rel in res['rels']:

    for ev in rel['evidences']:

        docid = ev['docid']

        sentid = ev['rel_sentence']
        sentence = ev.get("sentence", "")

        aSent = sentid.split(".")
        aSentNum = int(aSent[-1])

        aSent = aSent[0:2]

        allowedSentIDs = set()
        allowedSentIDs.add(str(SentenceID.fromArray(aSent + [aSentNum - 1])))
        allowedSentIDs.add(str(SentenceID.fromArray(aSent + [aSentNum + 1])))


        verbdir = ev['rel_direction_verb']

        lid = ev['lid']
        loid = ev['lontid']

        rid = ev['rid']
        roid = ev['rontid']

        trusts = (ev['trust']['verb'], ev['trust']['stack'], ev['trust']['relex'], ev['trust']['conj'])

        #if len([x for x in trusts[0:3] if x > 0]) == 0:
        #    continue


        effect = res['pmidinfo']['categories'].get(docid, [None])
        message = res['pmidinfo']['messengers'].get(docid, [None])

        foundEffects = set()
        if effect != None:

            for termEffect in effect:
                for x in termEffect['evidences']:
                    if x[0] == sentid:
                        foundEffects.add((termEffect['termid'], termEffect['termname']))

            if len(foundEffects) == 0:
                for termEffect in effect:
                    for x in termEffect['evidences']:
                        if x[0] in allowedSentIDs:
                            foundEffects.add((termEffect['termid'], termEffect['termname']))

        foundMessages = set()
        if message != None:

            for termMessage in message:

                for x in termMessage['evidences']:
                    if x[0] == sentid:
                        foundMessages.add((termMessage['termid'], termMessage['termname']))


            if len(foundMessages) == 0:
                for termMessage in message:
                    for x in termMessage['evidences']:
                        if x[0] in allowedSentIDs:
                            foundMessages.add((termMessage['termid'], termMessage['termname']))

        trustStr = sep.join([str(x) for x in trusts])

        infoTuple = pmid2info.get(docid, ("", "","",-2))

        infoTuple = [str(x) for x in infoTuple]

        elemcount += len(foundEffects) * len(foundMessages)

        if len(foundEffects) == 0 or len(foundMessages) == 0:
            print("INC EFF MESS", len(foundEffects), len(foundMessages))

        relEvCount += 1

        for messageElem in foundMessages:
            for effectElem in foundEffects:


                allData = [docid, sentid, sentence, verbdir, lid, loid, rid, roid]+ list(infoTuple) + [messageElem[0], messageElem[1], effectElem[0], effectElem[1]] + list(trusts)
                ws.append(allData)

                ws.cell(row=ws._current_row, column=9).hyperlink = infoTuple[0]
                ws.cell(row=ws._current_row, column=9).value = infoTuple[0]
                ws.cell(row=ws._current_row, column=9).style = "Hyperlink"

                print(allData, sep=sep)

                evCount += 1
                evPMID.add(docid)
                evSent.add(sentid)


print(relEvCount)
print(evCount)
print(len(evPMID))
print(len(evSent))

# Save the file

if args.output != None:

    wb.save(args.output.name)

#print(elemcount, file=sys.stderr)