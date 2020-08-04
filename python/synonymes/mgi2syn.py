from collections import Counter
import sys, os

sys.path.append(os.path.dirname(__file__) + "/../")

from utils.idutils import loadExludeWords

sys.path.insert(0, "/mnt/d/dev/git/nameConvert/")


from utils.DataFrame import DataFrame

from utils.idutils import printToFile
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords



if __name__ == '__main__':


    MGIdata = DataFrame.parseFromFile("/mnt/d/owncloud/data/miRExplore/MRK_Sequence.rpt")

    mgiID = MGIdata.getColumnIndex("MGI Marker Accession ID")
    mgiSym = MGIdata.getColumnIndex("Marker Symbol")
    #mgiSyn = MGIdata.getColumnIndex("Marker Synonyms (pipe-separated)")
    mgiUniprot = MGIdata.getColumnIndex("UniProt IDs")

    mgiGeneType = 20

    MGIinfo = {}
    foundUniprotIDs = set()
    uniprotID2MGI = {}

    locID2sym = {}

    mirIDs = set()

    for elem in MGIdata.data:

        locID = str(elem[mgiID])
        locSym = str(elem[mgiSym])
        locUniprot = elem[mgiUniprot]

        if len(elem) >= 21:
            if elem[20] in ["miRNA gene"]:
                print("miRNA gene", locID, locSym, locUniprot, file=sys.stderr)
                mirIDs.add(locID)

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

    if os.path.exists(uniprotConvFile):
        convData = DataFrame.parseFromFile(uniprotConvFile )

    else:

        odb = UniprotStore()
        convData = odb.fetch(GeneIdentity.UNIPROT, list(foundUniprotIDs),
                             [GeneIdentity.ALIAS, GeneIdentity.PROTEIN_NAMES])

        printToFile([str(convData)], uniprotConvFile)

    uniprotIdx = convData.getColumnIndex( "GeneIdentity.UNIPROT" )
    geneNamesIdx = convData.getColumnIndex( "GeneIdentity.ALIAS" )
    protNamesIdx = convData.getColumnIndex( "GeneIdentity.PROTEIN_NAMES" )

    for result in convData.data:

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

        if mgiID in mirIDs:
            print("Skipping mir ID", mgiID)
            continue

        if None in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove(None)

        if "None" in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove("None")


        allSyms = []

        for x in MGIinfo[mgiID]['sym']:

            if x.startswith("EC "):
                continue

            allSyms.append(x)

        printID = mgiID.replace(':', '_',1)
        if printID in locID2sym and len(locID2sym[printID]) > 0:
            allSyms.append(mgiID)
            allSyms.append(mgiID.replace(':', '_'))
            printID = locID2sym[printID]

        synline = printID + ":" + "|".join( allSyms)
        synonyme = Synonym.parseFromLine(synline)

        vAllSyns.append(synonyme)


    #exWords = loadExludeWords(common=True, generic=True, disease=False, taxnames=False, cell_co=False)

    exWords = loadExludeWords(cell_co=False, common=False, generic=True, syngrep=False)
    vPrintSyns = handleCommonExcludeWords(vAllSyns, exWords, mostCommonCount=500, maxCommonCount=7, addAlphaBeta=True, addHyphenGene=True, removeSyn=lambda synonym: synonym.id.startswith('MIR') and not synonym.id.endswith('HG'))

    printToFile(vPrintSyns, "/mnt/d/dev/data/pmid_jun2020/synonyms/mgi.syn", codec="utf8")

    """
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
    """
