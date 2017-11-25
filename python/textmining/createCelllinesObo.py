from collections import defaultdict
from nertoolkit.geneontology.GeneOntology import GeneOntology, GOTerm, GOSynonymeScope, GOSynonyme

from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords
from utils.idutils import dataDir, loadExludeWords, printToFile, speciesName2TaxID
import math
celloObo = GeneOntology(dataDir + "miRExplore/cellosaurus/cellosaurus.obo")
tax2cells = defaultdict(set)

allowedTaxIDs = set([str(speciesName2TaxID[x]) for x in speciesName2TaxID])

xref2cello = defaultdict(set)

considerXRefs = ['CL', 'CLO', 'DOID', 'UBERON']

for cellID in celloObo.dTerms:

    oboNode = celloObo.dTerms[cellID]
    oboXRefs = oboNode.xref

    taxID = {'all'}
    if oboXRefs != None:
        for xref in oboXRefs:
            if xref.startswith('NCBI_TaxID'):

                newTaxID = xref.split(' ')[0].split(':')[1]

                if newTaxID in allowedTaxIDs:
                    taxID.add(newTaxID)

            else:

                aref = xref.split('!')[0].split(':')

                if aref[0] in considerXRefs:

                    if '_' in aref[1]:
                        refid = aref[1].replace('_', ':')
                    else:
                        refid = aref[0] + ":" + aref[1]


                    xref2cello[refid].add( oboNode.id )



for oboID in xref2cello:

    if len(xref2cello[oboID]) > 1:
        print("BLA" + oboID)
    else:
        pass
        #print(oboID)

celllObo = GeneOntology(dataDir + "miRExplore/cellline_ontology/clo.obo")
newid = 0

metaObo = GeneOntology.mergeOntologies([celloObo, celllObo])
metaObo.saveFile("/tmp/test.obo")


originalTermIDs = [x for x in metaObo.dTerms]

print("Number of merges: " , len(xref2cello))

digits = int(math.log10(len(xref2cello)))+1
formatStr = '{:0'+str(digits)+'}'

for termID in originalTermIDs:

    if not termID in metaObo.dTerms:
        continue

    if termID in xref2cello:
        allXrefs = xref2cello[ termID ]

        allCelllIDs = set(allXrefs)
        listIDs = list(allCelllIDs) + [termID]

        idstr = 'META:' + formatStr.format(newid)
        newid += 1

        if newid % 1000 == 0:
            print("Merging: " + str(allCelllIDs) + " " + termID + " into " + idstr)
            print("Elements in obo: ", len(metaObo))

        metaObo.mergeTerms( listIDs, idstr, linkChildren=False )

        #replace any occurrence of termID and celloID in xref2cello with idstr
        for elem in xref2cello:
            xrefs = set(xref2cello[elem])

            intersection = allCelllIDs.intersection(xrefs)
            if len(intersection) == 0:
                continue

            xrefs = xrefs.difference(intersection)
            xrefs.add(idstr)

            xref2cello[elem].clear()

            for x in xrefs:
                xref2cello[elem].add(x)

            xref2cello[elem] = xrefs

metaObo.linkChildren()

metaObo.saveFile(dataDir + "/miRExplore/meta_cells.obo")