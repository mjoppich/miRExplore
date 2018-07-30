import matplotlib.pyplot as plt
from itertools import combinations
import pyvenny
from Bio import Entrez
from jinja2 import Environment, PackageLoader
from porestat.utils.DataFrame import DataFrame, DataRow, ExportTYPE


def fifteen_counts(a, b, c, d):
    union = 1

    # get sizes for regions delimited by set areas
    # _ is for - as in substract, i is for intersection, u is for union
    return { 'A_BuCuD': "%0.2f" % (len(a - set.union(b,c,d))/union),
             'B_AuCuD': "%0.2f" % (len(b - set.union(a,c,d))/union),
             'C_AuBuD': "%0.2f" % (len(c - set.union(a,b,d))/union),
             'D_AuBuC': "%0.2f" % (len(d - set.union(a,b,c))/union),

             'AiB_CuD': "%0.2f" % (len(set.intersection(a,b) - set.union(c,d))/union),
             'AiC_BuD': "%0.2f" % (len(set.intersection(a,c) - set.union(b,d))/union),
             'AiD_BuC': "%0.2f" % (len(set.intersection(a,d) - set.union(b,c))/union),
             'BiC_AuD': "%0.2f" % (len(set.intersection(b,c) - set.union(a,d))/union),
             'BiD_AuC': "%0.2f" % (len(set.intersection(b,d) - set.union(a,c))/union),
             'CiD_AuB': "%0.2f" % (len(set.intersection(c,d) - set.union(a,b))/union),

             'AiBiC_D': "%0.2f" % (len(set.intersection(a,b,c) - d)/union),
             'BiCiD_A': "%0.2f" % (len(set.intersection(b,c,d) - a)/union),
             'AiCiD_B': "%0.2f" % (len(set.intersection(a,c,d) - b)/union),
             'AiBiD_C': "%0.2f" % (len(set.intersection(a,b,d) - c)/union),

             'AiBiCiD': "%0.2f" % (len(set.intersection(a,b,c,d))/union) }

def render_four_set_venn(a, b, c, d, titles={'A': 'A',
                                             'B': 'B',
                                             'C': 'C',
                                             'D': 'D'}):
    sections = fifteen_counts(a, b, c, d)

    env = Environment(loader=PackageLoader('pyvenny', 'templates'))
    template = env.get_template('four_set_venn.svg')

    mergedItems = {}

    items = [x for x in sections.items()]
    items = items + [x for x in  titles.items()]

    return template.render( dict(items) )

allelems = {}

with open("/tmp/tm_soehnlein", 'r') as fin:
    for line in fin:

        line = line.strip().split("\t")

        descr = line[0]
        pmids = eval(line[1])

        allelems[descr] = set(pmids)

        pmids = None

for x in allelems:
    print(x, len(allelems[x]))

print(len(allelems['TISSUES'].intersection(allelems['DOID'])))


ntd = set.intersection(allelems['NEUTROPHIL'], allelems['TISSUES'], allelems['DOID'])
ntg = set.intersection(allelems['NEUTROPHIL'], allelems['TISSUES'], allelems['INFLAMM'])


def getPMIDTitles(pmids):
    Entrez.email = "joppich@bio.ifi.lmu.de"

    pmid2title = {}

    if len(pmids) == 0:
        return pmid2title

    epostResult = Entrez.read(Entrez.epost('pubmed', id=",".join([str(x) for x in pmids])))
    webEnv = epostResult['WebEnv']
    queryKey = epostResult['QueryKey']

    handle = Entrez.efetch(db="pubmed", webenv=webEnv, query_key=queryKey, retmode='XML')
    record = Entrez.read(handle)

    pmid2title = {}

    for article in record['PubmedArticle']:
        pubmedID = article['PubmedData']['ArticleIdList'][0] if len(
            article['PubmedData']['ArticleIdList']) > 0 else "-1"
        pubID = int(pubmedID)

        artInfo = article['MedlineCitation']['Article']
        articleTitle = artInfo['ArticleTitle']
        articleJournal = artInfo['Journal']['Title'] if 'Journal' in artInfo else ''

        pmid2title[pubID] = articleTitle

    return pmid2title


res = DataFrame()
res.addColumns(["SET", "PMID_ID", "PMID_TITLE", 'Common'])

print(ntd)
print("NTD", len(ntd))

pmidt = getPMIDTitles(ntd)
for x in sorted([x for x in pmidt]):

    dataDict = {
        'SET': 'NTinfect',
        'PMID_ID': "<a href='https://www.ncbi.nlm.nih.gov/pubmed/"+str(x)+"' target='_blank'>"+str(x)+"</a>",
        'PMID_TITLE': pmidt[x],
        'Common': x in ntg
    }

    res.addRow(DataRow.fromDict(dataDict))

    print(x, pmidt[x])

print("\n\n\n\n\n")

print(ntg)
print("NTG", len(ntg))
pmidt = getPMIDTitles(ntg)
for x in sorted([x for x in pmidt]):
    dataDict = {
        'SET': 'NTinflamm',
        'PMID_ID': "<a href='https://www.ncbi.nlm.nih.gov/pubmed/" + str(x) + "' target='_blank'>" + str(x) + "</a>",
        'PMID_TITLE': pmidt[x],
        'Common': x in ntd

    }

    res.addRow(DataRow.fromDict(dataDict))

    print(x, pmidt[x])

res.export("/home/mjoppich/win/Desktop/soehnlein_overview.html", ExportTYPE.HTML)

with open('/home/mjoppich/win/Desktop/soehnlein.svg', 'w') as f:
    f.write(render_four_set_venn(
        allelems['NEUTROPHIL'],
        allelems['TISSUES'],
        allelems['DOID'],
        allelems['INFLAMM'],

        titles={'A': 'NEUTROPHIL',
                'B': 'TISSUES',
                'C': 'Infection',
                'D': 'Inflammation'}
    )
    )