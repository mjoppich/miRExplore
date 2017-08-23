from collections import Counter

from nameconvert.stores import GeneIdentity
from nameconvert.stores import UniprotStore
from nameconvert.utils.DataFrame import DataFrame
from nameconvert.utils.Files import printToFile, fileExists
from scripts.hgncSynCount import Synonyme


if __name__ == '__main__':


    MGIdata = DataFrame.parseFromFile("../MRK_Sequence.rpt")

    mgiID = MGIdata.getColumnIndex("MGI Marker Accession ID")
    mgiSym = MGIdata.getColumnIndex("Marker Symbol")
    #mgiSyn = MGIdata.getColumnIndex("Marker Synonyms (pipe-separated)")
    mgiUniprot = MGIdata.getColumnIndex("UniProt IDs")

    MGIinfo = {}
    foundUniprotIDs = set()
    uniprotID2MGI = {}

    for elem in MGIdata.vElements:

        locID = str(elem[mgiID])
        locSym = str(elem[mgiSym])
        locUniprot = elem[mgiUniprot]

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





    if fileExists("../uniprotConv.tsv"):
        convData = DataFrame.parseFromFile("../uniprotConv.tsv" )

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

        printToFile([str(convData)], "../uniprotConv.tsv")

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

        synline = mgiID.replace(':', '_') + ":" + "|".join( MGIinfo[mgiID]['sym'] )

        synonyme = Synonyme.parseFromLine(synline)

        vAllSyns.append(synonyme)

    synCounter = Counter()
    for synonyme in vAllSyns:

        for syn in synonyme:
            synCounter[syn] += 1

    setCommonWords = set()

    for synWordCount in synCounter.most_common(66):
        setCommonWords.add(synWordCount[0])
        print(synWordCount[0] + " " + str(synWordCount[1]))

    vPrintSyns = []

    for synonyme in vAllSyns:

        synonyme.removeCommonSynonymes(setCommonWords)

        if len(synonyme) > 0:
            vPrintSyns.append(synonyme)

    printToFile( vPrintSyns, "../mgi.syn" )
