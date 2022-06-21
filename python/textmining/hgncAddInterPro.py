import gzip
import os

from collections import defaultdict
from porestat.utils.DataFrame import DataFrame

from utils.idutils import miRExploreDir

hgncData = DataFrame.parseFromFile(miRExploreDir + "/hgnc.tsv")


allUniprotIDs = set()
for row in hgncData:

    uniprotVals = row['UniProt ID(supplied by UniProt)']

    if uniprotVals == None:
        continue

    uniprotVals = uniprotVals.strip()
    uniprotIDs = uniprotVals.split(', ')

    for x in uniprotIDs:
        allUniprotIDs.add(x)

print(len(allUniprotIDs))
allUniprotIDs = sorted(allUniprotIDs)

uniprot2ipr = defaultdict(set)
neededUniprotIDs = miRExploreDir + "/interpro/relevant.uniprot.list"

# zgrep -f relevant.uniprot.list > relevant.uniprot.ipr.list
iprList = miRExploreDir + '/interpro/relevant.uniprot.ipr.list'
iprArchive = miRExploreDir + '/interpro/protein2ipr.dat.gz'

if not os.path.isfile(neededUniprotIDs) or not os.path.isfile(iprList):

    with open(neededUniprotIDs, 'w') as outfile:
        for uniprotID in allUniprotIDs:
            outfile.write("^" + uniprotID + '\t\n')

    exit(0)

    command = " ".join(["zgrep", "-f", neededUniprotIDs, iprArchive, ">", iprList])

    print("Executing command:", command)
    os.subprocess.call(command, shell=True)

with open(iprList, 'r') as f_in:

    cnt = 0
    cntRows = 0
    for line in f_in:
        if line[-1] == '\n':
            line = line[0:len(line)-1]

        aline = line.split('\t')

        uniprotID = aline[0]

        cntRows += 1
        if cntRows % 100000 == 0:
            print("Searched Rows", cntRows)


        if not uniprotID in allUniprotIDs:
            continue

        iprID = aline[1]

        uniprot2ipr[uniprotID].add(iprID)

        cnt += 1

        if cnt % 1000 == 0:
            print("Found uniprot ids", cnt, "of", len(allUniprotIDs))

for x in uniprot2ipr:
    print(x, uniprot2ipr[x])


uniprotColIdx = hgncData.getColumnIndex('UniProt ID(supplied by UniProt)')
iprColIdx = hgncData.addColumn('ipr', None)

def addIPRs( oldData ):


    uniprotVals = oldData[uniprotColIdx]

    if uniprotVals == None:
        return oldData

    uniprotVals = uniprotVals.strip()
    uniprotIDs = uniprotVals.split(', ')

    iprs = set()
    for x in uniprotIDs:
        if x in uniprot2ipr:
            iprs = iprs.union(uniprot2ipr[x])

    iprValues = ", ".join(iprs)

    oldData[iprColIdx] = iprValues
    return oldData

hgncData.applyToRow( addIPRs )
hgncData.export(miRExploreDir + "/interpro/hgnc.interpro.tsv")