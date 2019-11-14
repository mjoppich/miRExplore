from collections import Counter
import sys

from utils.idutils import loadExludeWords

sys.path.insert(0, "/mnt/d/dev/git/nameConvert/")


from nameconvert.stores import GeneIdentity
from nameconvert.stores import UniprotStore
from nameconvert.utils.DataFrame import DataFrame
from nameconvert.utils.Files import printToFile, fileExists
from synonymes.Synonym import Synonym



if __name__ == '__main__':


    MGIdata = DataFrame.parseFromFile("/mnt/d/owncloud/data/miRExplore/MRK_Sequence.rpt")

    mgiID = MGIdata.getColumnIndex("MGI Marker Accession ID")
    mgiSym = MGIdata.getColumnIndex("Marker Symbol")
    #mgiSyn = MGIdata.getColumnIndex("Marker Synonyms (pipe-separated)")
    mgiUniprot = MGIdata.getColumnIndex("UniProt IDs")

    MGIinfo = {}
    foundUniprotIDs = set()
    uniprotID2MGI = {}

    locID2sym = {}

    for elem in MGIdata.vElements:

        locID = str(elem[mgiID])
        locSym = str(elem[mgiSym])
        locUniprot = elem[mgiUniprot]

        locID2sym[locID.replace(':', '_', 1)] = locSym.upper()

        if locUniprot == None:
            locUniprot = set()
        else:
            locUniprot = set(locUniprot.split("|"))


        MGIinfo[ locID ] = { 'sym': set([locSym]) }

        for uniprotID in locUniprot:

            if not uniprotID in uniprotID2MGI:
                uniprotID2MGI[uniprotID] = set()

            uniprotID2MGI[ uniprotID ].add( locID )

            foundUniprotIDs.add(uniprotID)


    print(len(foundUniprotIDs))



    uniprotConvFile = "/mnt/d/owncloud/data/miRExplore/uniprotConv_mgi.tsv"

    if fileExists(uniprotConvFile):
        convData = DataFrame.parseFromFile(uniprotConvFile )

        allHeaders = convData.dHeader
        newHeaders = {}

        for x in allHeaders:

            enumKey = x

            if x.startswith('GeneIdentity.'):
                enumKey = x[len('GeneIdentity.'):]

            fieldName = GeneIdentity[enumKey]

            newHeaders[fieldName] = allHeaders[x]

            convData.dHeader = newHeaders

    else:

        odb = UniprotStore()
        convData = odb.fetch(GeneIdentity.UNIPROT, list(foundUniprotIDs),
                             [GeneIdentity.ALIAS, GeneIdentity.PROTEIN_NAMES])

        printToFile([str(convData)], uniprotConvFile)

    uniprotIdx = convData.getColumnIndex( GeneIdentity.UNIPROT )
    geneNamesIdx = convData.getColumnIndex( GeneIdentity.ALIAS )
    protNamesIdx = convData.getColumnIndex( GeneIdentity.PROTEIN_NAMES )

    for result in convData.vElements:

        uniprotID = result[uniprotIdx]

        mgiIDs = uniprotID2MGI[uniprotID]

        if mgiIDs == None or len(mgiIDs) == 0:
            continue

        for mgiID in mgiIDs:

            if not mgiID in MGIinfo:
                continue

            data = MGIinfo[mgiID]

            if result[geneNamesIdx] != None:

                geneNames = result[geneNamesIdx].split(" ")

                for geneName in geneNames:
                    data['sym'].add(geneName)

            if result[protNamesIdx] != None:

                geneNames = result[protNamesIdx].split(") (")

                if len(geneNames) > 1:
                    geneNames[len(geneNames) - 1] = geneNames[len(geneNames) - 1].replace(')', '')

                initGeneName = geneNames[0].split(" (")
                geneNames = [x for x in initGeneName] + geneNames[1:]


                for geneName in geneNames:
                    data['sym'].add(geneName)

            MGIinfo[mgiID] = data


    vAllSyns = []
    for mgiID in MGIinfo:

        if None in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove(None)

        if "None" in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove("None")


        allSyms = []

        for x in MGIinfo[mgiID]['sym']:

            if x.startswith("EC "):
                continue

            allSyms.append(x)

        synline = mgiID.replace(':', '_') + ":" + "|".join( allSyms)

        synonyme = Synonym.parseFromLine(synline)

        vAllSyns.append(synonyme)

    synCounter = Counter()
    for synonyme in vAllSyns:

        for syn in synonyme:
            synCounter[syn] += 1

    setCommonWords = set()

    for synWordCount in synCounter.most_common(10):
        setCommonWords.add(synWordCount[0])
        print(synWordCount[0] + " " + str(synWordCount[1]))

    exWords = loadExludeWords(common=True, generic=True, disease=False, taxnames=False, cell_co=False)

    for cat in exWords:
        for w in exWords[cat]:
            setCommonWords.add(w)

    vPrintSyns = []

    print("family" in setCommonWords)
    print("DNA" in setCommonWords)

    for synonyme in vAllSyns:

        synonyme.removeCommonSynonymes(setCommonWords)

        #if len(synonyme.syns) == 1:
        #    continue

        #if len(synonyme) > 0:
        #    newID = locID2sym[synonyme.id]
        #    synonyme.id = newID

        vPrintSyns.append(synonyme)

    printToFile( vPrintSyns, "/mnt/d/owncloud/data/miRExplore/textmine/synonyms/mgi.syn" )
