from Bio import Entrez
from porestat.utils.DataFrame import DataFrame

from utils.idutils import miRExploreDir


dbData = DataFrame.parseFromFile(miRExploreDir + "/miR2Disease/AllEntries.txt", ['mirna', 'disease', 'effect', 'measurement', 'year', 'title'], bConvertTextToNumber=False)

pmidTitleIdx = dbData.getColumnIndex('title')

allTitles = []

for row in dbData:

    title = row['title']

    if title == None:
        continue

    title = title.strip()

    if not title[-1] == '.':
        title += "."

    allTitles.append(title)

allTitles = list(set(allTitles))

print(len(allTitles))

titlesSearch = []
for title in allTitles:

    titlesSearch.append( title )


Entrez.email = 'joppich@bio.ifi.lmu.de'

stepSize=1

fetchIDs = set()

for i in range(0, len(titlesSearch), stepSize):

    maxI = min(len(titlesSearch), i+stepSize)

    print("Fetching", i, maxI)

    chunk = titlesSearch[i:maxI]

    if len(chunk) > 1:
        chunk = [ "(" + x + "[Title])" for x in chunk ]

    for title in chunk:
        print(title)

    searchTerm=" OR ".join(chunk)

    print(searchTerm)

    handle = Entrez.esearch(db='pubmed', term=searchTerm, field='title')
    record = Entrez.read(handle)

    elementCount = int(record['Count'])

    if elementCount > 0:
        for x in record['IdList']:
            fetchIDs.add(x)

print("Fetched IDs", len(fetchIDs))

epostResult = Entrez.read(Entrez.epost('pubmed', id=",".join(fetchIDs)))
webEnv = epostResult['WebEnv']
queryKey = epostResult['QueryKey']

handle = Entrez.efetch(db="pubmed", webenv=webEnv, query_key=queryKey, retmode='XML')
record = Entrez.read(handle)

artTitle2ID = {}

for article in record['PubmedArticle']:
    pubmedID = article['PubmedData']['ArticleIdList'][0] if len(article['PubmedData']['ArticleIdList']) > 0 else "-1"
    pubID = int(pubmedID)

    artInfo = article['MedlineCitation']['Article']
    articleTitle = artInfo['ArticleTitle']

    if not articleTitle in allTitles:
        print("Not in all titles:", pubID, articleTitle)
    else:
        artTitle2ID[articleTitle] = pubID


print("Downloaded Article Info", len(artTitle2ID))

for title in artTitle2ID:
    print(title, artTitle2ID[title])

for title in allTitles:
    if not title in artTitle2ID:
        print("not found", title)


pmidIdx = dbData.addColumn('pmid', -1)

def addPMIDfunc(olddata):

    if olddata[pmidTitleIdx] == None:
        return olddata

    artTitle = olddata[pmidTitleIdx]
    getPMID = artTitle2ID.get(artTitle, None)

    if getPMID == None:
        return olddata

    olddata[pmidIdx] = getPMID

    return olddata


dbData.applyToRow(addPMIDfunc)

dbData.export(miRExploreDir + "/miR2Disease/mirna_disease.tsv")