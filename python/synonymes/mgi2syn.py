from collections import Counter
import sys, os

sys.path.append(os.path.dirname(__file__) + "/../")

from utils.idutils import loadExludeWords
from utils.DataFrame import DataFrame

from utils.idutils import printToFile
from synonymes.Synonym import Synonym
from synonymes.SynonymUtils import handleCommonExcludeWords



if __name__ == '__main__':
    import argparse


    parser = argparse.ArgumentParser(description='Convert Medline XML to miRExplore base files')
    parser.add_argument('-r', '--rpt', type=argparse.FileType("r"), required=True, help="input genes file")
    parser.add_argument('-s', '--syn', type=argparse.FileType("w"), required=True, help="output synonym file")
    args = parser.parse_args()


    MGIdata = DataFrame.parseFromFile(args.rpt.name)

    mgiID = MGIdata.getColumnIndex("MGI Marker Accession ID")
    mgiSym = MGIdata.getColumnIndex("Marker Symbol")
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



    uniprotConvFile = args.rpt.name + "_uniprot_conversion.tsv"

    if not os.path.exists(uniprotConvFile):
        
        import biomart
        import pandas as pd
        import io
        server = biomart.BiomartServer( "ensembl.org/biomart" )
        
        mds = server.datasets["mmusculus_gene_ensembl"]
        
        response = mds.search({"attributes":["mgi_description", "mgi_symbol", "uniprotswissprot"]}, header=1)
        rawData = pd.read_csv(io.StringIO("\n".join([x.decode() for x in response.raw.readlines()])), sep="\t")
        rawData.columns = ["GeneIdentity.ALIAS", "GeneIdentity.PROTEIN_NAMES", "GeneIdentity.UNIPROT" ]
        rawData.to_csv(uniprotConvFile, sep="\t", header=True, index=False)

    convData = DataFrame.parseFromFile( uniprotConvFile )

    uniprotIdx = convData.getColumnIndex( "GeneIdentity.UNIPROT" )
    geneNamesIdx = convData.getColumnIndex( "GeneIdentity.ALIAS" )
    protNamesIdx = convData.getColumnIndex( "GeneIdentity.PROTEIN_NAMES" )

    for result in convData.data:

        uniprotID = result[uniprotIdx]
        
        if uniprotID is None:
            continue

        mgiIDs = uniprotID2MGI.get(uniprotID, None)

        if mgiIDs == None or len(mgiIDs) == 0:
            continue

        for mgiID in mgiIDs:

            if not mgiID in MGIinfo:
                continue

            data = MGIinfo[mgiID]

            if result[geneNamesIdx] != None:

                geneNames = result[geneNamesIdx]#.split(" ")
                data['sym'].add(geneNames)

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
            print("Skipping MGI ID", mgiID)
            continue

        if None in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove(None)

        if "None" in MGIinfo[mgiID]['sym']:
            MGIinfo[mgiID]['sym'].remove("None")


        allSyms = []

        for x in MGIinfo[mgiID]['sym']:

            if x.startswith("EC "):
                continue

            if x.upper() in ["TH1"]:
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


    for syn in vAllSyns:

        removeSyns = []
        for synword in syn.syns:

            if len(synword) == 1:
                removeSyns.append(synword)

        if len(removeSyns) > 0:
            print(syn.id, removeSyns)

            syn.removeSyn(removeSyns)

    #exWords = loadExludeWords(common=True, generic=True, disease=False, taxnames=False, cell_co=False)

    exWords = loadExludeWords(cell_co=False, common=False, generic=True, syngrep=False)
    vPrintSyns = handleCommonExcludeWords(vAllSyns, exWords, mostCommonCount=500, maxCommonCount=7, addAlphaBeta=True, addHyphenGene=True, removeSyn=lambda synonym: synonym.id.startswith('MIR') and not synonym.id.endswith('HG'))

    printToFile(vPrintSyns, args.syn.name, codec="utf8")

